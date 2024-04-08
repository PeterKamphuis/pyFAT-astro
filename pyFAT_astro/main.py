# -*- coding: future_fstrings -*-

# This is the python version of FAT
import numpy as np
import os
import psutil
import pyFAT_astro
import pyFAT_astro.Support.support_functions as sf
import pyFAT_astro.Support.read_functions as rf
import sys
import traceback
import warnings
import threading


from datetime import datetime
from multiprocessing import Pool,get_context,Lock,Manager
from omegaconf import OmegaConf
from pyFAT_astro.FAT_Galaxy_Loop import FAT_Galaxy_Loop,MP_initialize_sofia,\
                                        MP_Fitting_Loop
from pyFAT_astro.config.defaults import defaults
from pyFAT_astro.Support.fat_errors import ProgramError,BadCatalogueError
from pyFAT_astro.Support.write_functions import reorder_output_catalogue
from pyFAT_astro.Support.log_functions import full_system_tracking

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))
try:
    from importlib.resources import files as import_pack_files
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    # For Py<3.9 files is not available
    from importlib_resources import files as import_pack_files


# String syntax ''' '''for multiline strings. " " for string without break and ' ' for indexing dictionaries

#from memory_profiler import profile
#@profile
def main(argv):
    try:
        #Get default settings
        print(f"This is version {pyFAT_astro.__version__} of pyFAT.")
        if pyFAT_astro.__branch__:
            print(f"This is a github distribution and we are on the branch {pyFAT_astro.__branch__}.")

        if '-v' in argv or '--version' in argv:
            #print(f"This is version {pyFAT_astro.__version__} of the program.")
            #if pyFAT_astro.__branch__:
            #    print(f"This is a github distribution and we are on the branch {pyFAT_astro.__branch__}.")
            sys.exit()


        help_message = '''
        Use pyFAT in this way for batch fitting:

            pyFAT configuration_file=FAT_Input.yml

        where configuration_file specifies a yaml file with specific settings
        such as the catalog.

        For fitting use pyFAT in this way:

            pyFAT cube_name=Input_Cube.fits

        Where Input_Cube.fits is the observation to be fitted. In this mode
        configuration_file can still be used to specify fit settings but
        catalogue and location setting will be ignored.

            pyFAT -h

        prints this message.

            pyFAT print_examples=True

        prints a yaml file (FAT_defaults.yml)  with the default values for all
        possible fitting parameters and an example input catalogue (FAT_Example_Catalogue.txt).
        The files are printed in the current working directory. In the yaml
        file values designated ??? indicated values without defaults.

        All config parameters can be set directly from the command e.g:

            pyFAT file_name=Input_Cube.fits fitting.ring_size=1.5 'fitting.fixed_parameters=[INCL,SDIS]'

        You can test your installation with:

            pyFAT installation_check=True

        '''

        if '-h' in argv or '--help' in argv:
            print(help_message)
            sys.exit()

        cfg = OmegaConf.structured(defaults)
        if cfg.ncpu == psutil.cpu_count():
            cfg.ncpu -= 1

        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        inputconf = OmegaConf.from_cli(argv)
        cfg_input = OmegaConf.merge(cfg,inputconf)
        if cfg_input.print_examples:
            no_cube = OmegaConf.masked_copy(cfg, ['ncpu','input','output',\
                'fitting','advanced'])
            with open('FAT_defaults.yml','w') as default_write:
                default_write.write(OmegaConf.to_yaml(no_cube))
            my_resources = import_pack_files('pyFAT_astro.config')
            data = (my_resources / 'FAT_Input_Catalogue.txt').read_bytes()
            with open('FAT_Example_Catalogue.txt','w+b') as default_write:
                default_write.write(data)

            print(f'''We have printed the file FAT_defaults.yml FAT_Input_Catalogue.txt in {os.getcwd()}.
''')
            sys.exit()

        if cfg_input.configuration_file:
            succes = False
            while not succes:
                try:
                    yaml_config = OmegaConf.load(cfg_input.configuration_file)
            #merge yml file with defaults
                    cfg = OmegaConf.merge(cfg,yaml_config)
                    succes = True
                except FileNotFoundError:
                    cfg_input.configuration_file = input(f'''
    You have provided a config file ({cfg_input.configuration_file}) but it can't be found.
    If you want to provide a config file please give the correct name.
    Else press CTRL-C to abort.
    configuration_file = ''')



        cfg = OmegaConf.merge(cfg,inputconf)

        if not any([cfg.cube_name, cfg.configuration_file, cfg.installation_check\
                    ,cfg.print_examples,cfg.input.catalogue]):
            print(help_message)
            sys.exit()
        # if we set more cpus than available we limit to the available cpus
        try:
            if cfg.ncpu > len(psutil.Process().cpu_affinity()):
                cfg.ncpu  = len(psutil.Process().cpu_affinity())
        except AttributeError:
            if cfg.ncpu > psutil.cpu_count():
                cfg.ncpu  = psutil.cpu_count()

        #Let's write and input example to the main directory
        if cfg.output.debug:
            with open(f'{cfg.input.main_directory}/FAT_Inputs-Run_{datetime.now().strftime("%d-%m-%Y")}.yml','w') as default_write:
                default_write.write(OmegaConf.to_yaml(cfg))


        #Transform all to a Configuration dictionary
        Original_Configuration = sf.setup_configuration(cfg)

        if cfg.output.debug:
            warnings.showwarning = warn_with_traceback

        #First we check for sofia and TiRiFiC
        Original_Configuration['SOFIA2'] = sf.find_program(Original_Configuration['SOFIA2'], "SoFiA 2")
        Original_Configuration['TIRIFIC'] = sf.find_program(Original_Configuration['TIRIFIC'], "TiRiFiC")

        if cfg.cube_name:
            Full_Catalogue = sf.Proper_Dictionary({})
            Full_Catalogue['ENTRIES'] = ['ENTRIES','ID','DISTANCE','DIRECTORYNAME','CUBENAME']
            Full_Catalogue['ID'] = [f"{os.path.splitext(cfg.cube_name.split('/')[-1])[0]}"]
            Full_Catalogue['DISTANCE'] = [-1.]
            Full_Catalogue['DIRECTORYNAME'] = ['./']
            Full_Catalogue['CUBENAME'] = [f"{os.path.splitext(cfg.cube_name.split('/')[-1])[0]}"]
        elif 'sofia_catalogue' in Original_Configuration['FITTING_STAGES']:
            Full_Catalogue = rf.sofia_input_catalogue(Original_Configuration)
        else:
            Full_Catalogue = rf.catalogue(Original_Configuration['CATALOGUE'],split_char= cfg.advanced.catalogue_split_character)
        # Get the longest directory name to format the output directory properlyFit_Tirific_OSC
        for directory in Full_Catalogue['DIRECTORYNAME']:
            if directory == './':
                directory = Original_Configuration['MAIN_DIRECTORY'].split('/')[-2]
            if len(directory) > Original_Configuration['MAXIMUM_DIRECTORY_LENGTH']:
                Original_Configuration['MAXIMUM_DIRECTORY_LENGTH'] = len(directory)

        # Create a file to write the results to if if required
        if Original_Configuration['OUTPUT_CATALOGUE']:
            if not os.path.exists(Original_Configuration['OUTPUT_CATALOGUE']) or Original_Configuration['NEW_OUTPUT']:
                if os.path.exists(Original_Configuration['OUTPUT_CATALOGUE']) and Original_Configuration['NEW_OUTPUT']:
                    os.rename(Original_Configuration['OUTPUT_CATALOGUE'],f"{os.path.splitext(Original_Configuration['OUTPUT_CATALOGUE'])[0]}_Prev.txt")
                with open(Original_Configuration['OUTPUT_CATALOGUE'],'w') as output_catalogue:
                    comment = 'Comments on Fit Result'
                    AC1 = 'OS'
                    output_catalogue.write(f"{'Directory Name':<{Original_Configuration['MAXIMUM_DIRECTORY_LENGTH']}s} {AC1:>6s} {comment}\n")

        if Original_Configuration['TIMING']:

            with open(Original_Configuration['MAIN_DIRECTORY']+'Timing_Result.txt','w') as timing_result:
                timing_result.write("Timing results for every section of the fit process for all galaxies.  \n")
            # If we do this we should have 1 cpu to keep going
            Original_Configuration['NCPU'] -= 1
            system_monitor = full_system_tracking(Original_Configuration)
            fst = threading.Thread(target=system_monitor.start_monitoring)
            fst.start()

            print(f"We are using {Original_Configuration['NCPU']} cpus for fitting and 1 for timing.")
        else:
            print(f"We are using {Original_Configuration['NCPU']} cpus.")
        #if start_galaxy not negative then it is catalogue ID
        if Original_Configuration['CATALOGUE_START_ID'] in ['-1','-1.']:
            Original_Configuration['CATALOGUE_START_ID'] = int(0)
        else:
            Original_Configuration['CATALOGUE_START_ID'] = int(np.where(Original_Configuration['CATALOGUE_START_ID'] == np.array(Full_Catalogue['ID'],dtype=str))[0][0])
        # If the end galaxy is -1 fit the whole catalogue
        if Original_Configuration['CATALOGUE_END_ID'] in ['-1','-1.']:
            Original_Configuration['CATALOGUE_END_ID'] = int(len(Full_Catalogue['ID']))
            if Original_Configuration['CATALOGUE_END_ID'] == 0:
                Original_Configuration['CATALOGUE_END_ID'] = 1
        else:
            Original_Configuration['CATALOGUE_END_ID'] = int(np.where(Original_Configuration['CATALOGUE_END_ID'] == np.array(Full_Catalogue['ID'],dtype=str))[0][0])+1
        # start the main fitting loop
        if float(Original_Configuration['CATALOGUE_START_ID']) > float(Original_Configuration['CATALOGUE_END_ID']):
            raise BadCatalogueError(f''' Your starting galaxy (Line nr = {Original_Configuration['CATALOGUE_START_ID']}) is listed after your ending galaxy (Line nr = {Original_Configuration['CATALOGUE_END_ID']}), maybe you have double catalogue ids?''')
            sys.exit(1)

        if Original_Configuration['MULTIPROCESSING']:
            Original_Configuration['VERBOSE_SCREEN'] = False
            #output_catalogue = copy.deepcopy(Original_Configuration['OUTPUT_CATALOGUE'])
            #Original_Configuration['OUTPUT_CATALOGUE'] = None
            no_processes,sofia_processes = sf.calculate_number_processes(Original_Configuration)
            Configs_and_Locks = []

            with Manager() as loop_manager:
                timing_lock = loop_manager.Lock()
                catalogue_lock = loop_manager.Lock()
                #In case of multiprocessing we want to make sure to start with
                #The big galaxies
                #Setup an array of configs with locks
                for current_galaxy_index in range(Original_Configuration['CATALOGUE_START_ID'], Original_Configuration['CATALOGUE_END_ID']):
                    Configs_and_Locks.append([sf.set_individual_configuration(current_galaxy_index,Full_Catalogue,Original_Configuration),timing_lock,catalogue_lock])

                #Get all intitial setups
                with get_context("spawn").Pool(processes=sofia_processes) as pool:
                    print(f'Starting size estimates with {sofia_processes} processes')
                    initial_setups = pool.starmap(MP_initialize_sofia, Configs_and_Locks)

                initial_setups = [x for x in initial_setups if x['Succes']]
                sizes = np.array([np.mean(x['Size']) for x in initial_setups]\
                    ,dtype=float)
                if len(sizes) > 0.:
                    sorted_ind = np.flip(sizes.argsort())
                    sorted_initial_setups = [[initial_setups[x],timing_lock,catalogue_lock] \
                        for x in sorted_ind]
                    initial_setups =[]
                    with get_context("spawn").Pool(processes=no_processes) as pool:
                        print(f'Starting fitting with {no_processes} processes')
                        finals = pool.starmap(MP_Fitting_Loop, sorted_initial_setups)

                else:
                    print(f'All galaxies can not be fitted')

            #For clarity we reorder the output results to match the input
            reorder_output_catalogue(Original_Configuration,Full_Catalogue)
            #Stitch all temporary outpu catalogues back together
            #with open(output_catalogue,'a') as catalogue:
            #    for x in results:
            #        catalogue.writelines(x)

        else:
            Original_Configuration['PER_GALAXY_NCPU'] = sf.set_limits(Original_Configuration['NCPU'],1,20)
            for current_galaxy_index in range(Original_Configuration['CATALOGUE_START_ID'], Original_Configuration['CATALOGUE_END_ID']):
                Configuration = sf.set_individual_configuration(current_galaxy_index,Full_Catalogue,Original_Configuration)
                catalogue_line = FAT_Galaxy_Loop(Configuration)
        
        if Original_Configuration['TIMING']:
            system_monitor.stop_monitoring()
            fst.join()

    except SystemExit:
        try:
            system_monitor.stop_monitoring()
            fst.join()
        except:
            pass
        pass
    except KeyboardInterrupt:
        traceback.print_exception(*sys.exc_info())
        try:
            system_monitor.stop_monitoring()
            fst.join()
        except:
            pass
        pass
    except:
        try:
            system_monitor.stop_monitoring()
            fst.join()
        except:
            pass
        raise ProgramError(f'''Something went wrong in the main. This should not happen. Please list an issue on github.''')

main.__doc__ = '''
 NAME:
     main
 PURPOSE:
      Fit Tilted Ring Models with Tirific in a fully automated manner
 CATEGORY:
      Main for fitting galaxies. Tirific still requires interactive fitting this code attempts
      to remedy that

 CALLING SEQUENCE:
     see pyFAT -h

 INPUTS:
    see pyFAT -h

 OUTPUTS:
     See Readme or just run the code

 EXAMPLE:
     pyFAT  configuration_file=/home/your_computer/FAT_dir/FAT_INPUT.yml'
'''
