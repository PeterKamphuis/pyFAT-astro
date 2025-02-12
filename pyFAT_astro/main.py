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

from pyFAT_astro.FAT_Galaxy_Loop import FAT_Galaxy_Loop,MP_initialize_sofia,\
                                        MP_Fitting_Loop
from pyFAT_astro.config.defaults import process_input
from pyFAT_astro.Support.fat_errors import ProgramError,BadCatalogueError,InputError
from pyFAT_astro.Support.write_functions import reorder_output_catalogue
from pyFAT_astro.Support.log_functions import full_system_tracking

def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))



# String syntax ''' '''for multiline strings. " " for string without break and ' ' for indexing dictionaries

#from memory_profiler import profile
#@profile
def main_trace():
    from viztracer import VizTracer
    with VizTracer(output_file="FAT_Run_Viztracer.json",min_duration=1000) as tracer:
        main()

def main():
    argv = sys.argv[1:]
    try:
        #Get default settings
        print(f"This is version {pyFAT_astro.__version__} of pyFAT.")
        if pyFAT_astro.__branch__:
            print(f"This is a github distribution and we are on the branch {pyFAT_astro.__branch__}.")
        cfg =process_input(argv)
        
        # if we set more cpus than available we limit to the available cpus
        try:
            if cfg.ncpu > len(psutil.Process().cpu_affinity()):
                cfg.ncpu  = len(psutil.Process().cpu_affinity())
        except AttributeError:
            if cfg.ncpu > psutil.cpu_count():
                cfg.ncpu  = psutil.cpu_count()

    

        #Transform all to a Configuration dictionary
        Original_Configuration = sf.setup_configuration(cfg)

        if cfg.output.debug:
            warnings.showwarning = warn_with_traceback
        if Original_Configuration['RP_SECTION'] == 'START':
            #First we check for sofia and TiRiFiC
            Original_Configuration['SOFIA2'] = \
                sf.find_program(Original_Configuration['SOFIA2'], "SoFiA 2")
            Original_Configuration['TIRIFIC'] =\
                sf.find_program(Original_Configuration['TIRIFIC'], "TiRiFiC")
            if not cfg.cube_name is None:
                Full_Catalogue = sf.Proper_Dictionary({})
                Full_Catalogue['ENTRIES'] =\
                    ['ENTRIES','ID','DISTANCE','DIRECTORYNAME','CUBENAME']
                Full_Catalogue['ID'] = [f"{cfg.cube_name.split('/')[-1]}"]
                Full_Catalogue['DISTANCE'] = [-1.]
                Full_Catalogue['DIRECTORYNAME'] = ['./']
                Full_Catalogue['CUBENAME'] = [cfg.cube_name.split('/')[-1]]
            elif 'sofia_catalogue' in Original_Configuration['FITTING_STAGES']:
                Full_Catalogue = rf.sofia_input_catalogue(Original_Configuration)
            elif not Original_Configuration['CATALOGUE'] is None:
                Full_Catalogue = rf.catalogue(Original_Configuration['CATALOGUE']\
                    ,split_char= cfg.advanced.catalogue_split_character)
            else:
                raise InputError(f'''                             
{'':8s}We could not find any of the following input: 
{'':8s}cube_name= (We found {cfg.cube_name})
{'':8s}input.catalogue= (We found { Original_Configuration['CATALOGUE']})
{'':8s}and 'sofia_catalogue' was not in the fitting stages (We found { Original_Configuration['FITTING_STAGES']})
{'':8s}Please add any of these when calling pyFAT or at them in the correct manner to the yml configuration file.
''')                        
        
            # Get the longest directory name to format the output directory properlyFit_Tirific_OSC
            for directory in Full_Catalogue['DIRECTORYNAME']:
                if directory == './':
                    directory = Original_Configuration['MAIN_DIRECTORY'].split('/')[-2]
                if len(directory) > Original_Configuration['MAXIMUM_DIRECTORY_LENGTH']:
                    Original_Configuration['MAXIMUM_DIRECTORY_LENGTH'] = len(directory)

            # Create a file to write the results to if if required
            if Original_Configuration['OUTPUT_CATALOGUE']:
                if not os.path.exists(Original_Configuration['OUTPUT_CATALOGUE']) or Original_Configuration['NEW']:
                    if os.path.exists(Original_Configuration['OUTPUT_CATALOGUE']) and Original_Configuration['NEW']:
                        os.rename(Original_Configuration['OUTPUT_CATALOGUE'],f"{os.path.splitext(Original_Configuration['OUTPUT_CATALOGUE'])[0]}_Prev.txt")
                    with open(Original_Configuration['OUTPUT_CATALOGUE'],'w') as output_catalogue:
                        comment = 'Comments on Fit Result'
                        AC1 = 'OS'
                        output_catalogue.write(f"{'Directory Name':<{Original_Configuration['MAXIMUM_DIRECTORY_LENGTH']}s} {AC1:>6s} {comment}\n")

          

              
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
       
        if Original_Configuration['TIMING']:
            if Original_Configuration['RP_SECTION'] == 'START':
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
           
        
        #if we do recovery multiprocessing is turned off but the points should still be written       
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
            # as in recovery we might be switching from MP to single reset PER_GALAXY_NCPU
            Original_Configuration['PER_GALAXY_NCPU'] = sf.set_limits(Original_Configuration['NCPU'],1,20)
            for current_galaxy_index in range(Original_Configuration['CATALOGUE_START_ID'], Original_Configuration['CATALOGUE_END_ID']):
                if Original_Configuration['RP_SECTION'] == 'START':
                    Configuration = sf.set_individual_configuration(current_galaxy_index,Full_Catalogue,Original_Configuration)
                else:
                    Configuration = Original_Configuration
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
    except InputError:
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
