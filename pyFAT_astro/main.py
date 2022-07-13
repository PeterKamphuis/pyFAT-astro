# -*- coding: future_fstrings -*-

# This is the python version of FAT

import numpy as np
import os
import pyFAT_astro.Support.support_functions as sf
import pyFAT_astro.Support.read_functions as rf
import sys
import traceback
import warnings

from datetime import datetime
from multiprocessing import Pool,get_context
from omegaconf import OmegaConf
from pyFAT_astro.FAT_Galaxy_Loops import FAT_Galaxy_Loops
from pyFAT_astro.config.defaults import defaults
from pyFAT_astro.Support.fat_errors import ProgramError

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

def main(argv):



    try:
        #Get default settings
        if '-v' in argv or '--version' in argv:
            print(f"This is version {pyFAT_astro.__version__} of the program.")
            sys.exit()

        help_message = '''
        Use pyFAT in this way for batch fitting:

            pyFAT configuration_file=FAT_Input.yml

        where configuration_file specifies a yaml file with specific settings
        such as the catalog.

        For fitting a single galaxy use pyFAT in this way:

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

        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        inputconf = OmegaConf.from_cli(argv)
        cfg_input = OmegaConf.merge(cfg,inputconf)
        if cfg_input.print_examples:
            no_cube = OmegaConf.masked_copy(cfg, ['input','output','fitting'])
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
        #Let's write and input example to the logs directory
        if cfg.output.debug:
            with open(f'{cfg.input.main_directory}/FAT_Inputs-{datetime.now().strftime("%d-%m-%Y")}.yml','w') as default_write:
                default_write.write(OmegaConf.to_yaml(cfg))

        #Transform all to a Configuration dictionary
        Original_Configuration = sf.setup_configuration(cfg)

        if Original_Configuration['DEBUG']:
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
                timing_result.write("This file contains the system start and end times for the fitting of each galaxy")

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
            raise CatalogError(f''' Your starting galaxy (Line nr = {Original_Configuration['CATALOGUE_START_ID']}) is listed after your ending galaxy (Line nr = {Original_Configuration['CATALOGUE_END_ID']}), maybe you have double catalogue ids?''')
            sys.exit(1)
        if Original_Configuration['MULTIPROCESSING']:
            Original_Configurations,no_processes = sf.calculate_number_processes(Original_Configuration)
            to_maps = [(x,Full_Catalogue ) for x in Original_Configurations]
            with get_context("spawn").Pool(processes=no_processes) as pool:
                results = pool.starmap(FAT_Galaxy_Loops, to_maps)
            #Stitch all temporary outpu catalogues back together
            with open(Original_Configuration['OUTPUT_CATALOGUE'],'a') as catalogue:
                for x in Original_Configurations:
                    with open(x['OUTPUT_CATALOGUE']) as tmp:
                        lines = tmp.readlines()
                    catalogue.writelines(lines[1:])
                    #clean up
                    os.remove(x['OUTPUT_CATALOGUE'])
        else:
            Original_Configuration['PER_GALAXY_NCPU'] = sf.set_limits(Original_Configuration['NCPU'],1,20)
            FAT_Galaxy_Loops(Original_Configuration,Full_Catalogue)
    except SystemExit:
        pass
    except:
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
