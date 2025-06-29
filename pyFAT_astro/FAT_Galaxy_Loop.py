# -*- coding: future_fstrings -*-
# This module contains the main loop over the fitting per galaxy

import copy
import os
import pyFAT_astro
import warnings

import psutil as psu
import pyFAT_astro.Support.clean_functions as cf
import pyFAT_astro.Support.fits_functions as ff
import pyFAT_astro.Support.read_functions as rf
import pyFAT_astro.Support.run_functions as runf
import pyFAT_astro.Support.support_functions as sf
import pyFAT_astro.Support.write_functions as wf

from datetime import datetime
from omegaconf import OmegaConf
from pyFAT_astro.Support.fat_errors import BadCatalogueError
from pyFAT_astro.Support.log_functions import print_log,enter_recovery_point\
    ,update_statistics

def FAT_Galaxy_Loop(Configuration):

    try:
        registered_exception = None
        if Configuration['RP_SECTION'] == 'START':
            current_run = 'Not Initialized'
            Fits_Files = initialize_loop(Configuration)
        else:
            if Configuration['TIMING']:
                Configuration['FAT_PSUPROCESS'] = psu.Process(
                    Configuration['FAT_PID'])
            Fits_Files = copy.deepcopy(Configuration['Fits_Files'])
            cf.cleanup(Configuration,Fits_Files, recovery = True)

        if Configuration['RP_SECTION'] == 'AFTER-INITIALIZATION':
            loop_sofia(Configuration,Fits_Files)
        
        # If you add any make sure that the fitstage  starts with 'Fit_' if Configuration['RP_SECTION'] == 'AFTER-INITIALIZATION':
        if (not Configuration['USED_FITTING'] is None) and \
            (Configuration['RP_SECTION'] == 'AFTER-SOFIA'):
            # Process the found source in sofia to set up the proper fitting and make sure source can be fitted

            Initial_Parameters = runf.check_source(
                    Configuration, Fits_Files)

            #sf.sofia_output_exists(Configuration,Fits_Files)
            print_log(f'''FAT_GALAXY_LOOPS: The source is well defined and we will now setup the initial tirific file
''', Configuration)
        else:
            try:
                Initial_Parameters = Configuration['Initial_Parameters']
            except:
                Initial_Parameters = None
                pass
            if (not Configuration['USED_FITTING'] is None) and \
                Configuration['OPTIMIZED']:

                 #We want to check if the cube has a decent number of pixels per beam.
              
                if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['OPTIMIZED_CUBE']}"):
                    ff.optimized_cube(Configuration, Fits_Files)
              
      

        if 'fit_tirific_osc' in Configuration['FITTING_STAGES']:
            if  Configuration['RP_SECTION'] in ['ITERATION',\
                'AFTER-PREPARATIONS','BEFORE-SMOOTHING']:
                current_run,Tirific_Template = runf.fitting_osc(
                    Configuration, Fits_Files,  Initial_Parameters)
            else:
                current_run = 'Not set'
                Tirific_Template = copy.deepcopy(Configuration['Tirific_Template'])
        elif 'fit_make_your_own' in Configuration['FITTING_STAGES']:
            print_log(f'''FAT_GALAXY_LOOPS: If you add any fitting routine make sure that the fit stage  starts with Fit_
''',Configuration)
            sf.create_directory(Configuration['USED_FITTING'],Configuration['FITTING_DIR'])
            Configuration['FINAL_COMMENT'] = 'This example does not work'
            catalogue_line = cf.finish_galaxy(Configuration)
            return catalogue_line
        else:
            Configuration['FINAL_COMMENT'] = 'You have chosen not to do any fitting'
            catalogue_line = cf.finish_galaxy(Configuration)
            return catalogue_line

        #if all the fitting has gone properly we create nice errors
        if  Configuration['RP_SECTION'] in ['AFTER-SMOOTHING']:
            if Configuration['QUANTITY'] != 5:
                if 'tirshaker' in Configuration['FITTING_STAGES']:
                    runf.tirshaker_call(
                        Configuration,Fits_Files)

                Configuration['FINAL_COMMENT'] = 'The fit has converged succesfully'
        

        catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run,
                         Fits_Files=Fits_Files,exiting=registered_exception)
        if Configuration['QUANTITY'] != 5:
            DHI = rf.get_DHI(
                Configuration, Model=Configuration['USED_FITTING'])
            Totflux = rf.get_totflux(
                Configuration, f"/Finalmodel/Finalmodel_mom0.fits")
            wf.basicinfo(Configuration, template=Tirific_Template, Tot_Flux=Totflux, DHI=[
                         DHI, Configuration['BEAM'][0]*Configuration['RING_SIZE']])
            if Configuration['INSTALLATION_CHECK']:
                cf.installation_check(
                    Configuration)
        print_log(f'''FAT_GALAXY_LOOP: This galaxy has succesfully ran through all parts of the fitting
''', Configuration)
    except Exception as e:
        if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
            Configuration['QUANTITY'] = 5
        else:
            Configuration['QUANTITY'] = 'error'
        registered_exception = e
        Configuration['FINAL_COMMENT'] = e
        Configuration['QUANTITY'] = 'error'
        try:
            tmp=Fits_Files
        except:
            Fits_Files= None
        
        catalogue_line = cf.finish_galaxy(Configuration,
                          current_run=current_run,Fits_Files=Fits_Files,
                          exiting=registered_exception)
       
    return catalogue_line


def initialize_loop(Configuration):
    Configuration['INITIALIZATION_TIME'][0] = datetime.now()
    # First check the starttime
    if Configuration['FITTING_DIR'][-2:] == '//':
        Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'

    if not Configuration['SOFIA_DIR']:
        Configuration['SOFIA_DIR'] = Configuration['FITTING_DIR']
    #ini_mode_factor =25
    # We initially set the variations to fixed for all parameters
    #let's see what happens if we immediately

    #Make a dictionary for the fitsfiles we use
    Fits_Files = {
        'ORIGINAL_CUBE': Configuration['INPUT_CUBE']}

    # If we have a fitting log we start writing
    log_statement = f'''This file is a log of the fitting process run at {Configuration ['FULL_TIME'][0]}.
{"":8s}This is version {pyFAT_astro.__version__} of the program.
'''
    if pyFAT_astro.__branch__:
        log_statement =f'''{log_statement} This is a github distribution and we are on the branch {pyFAT_astro.__branch__}.
'''


    # Adapt configuration to hold some specifics to this galaxy

    Configuration['LOG_DIRECTORY'] = f'''{Configuration['FITTING_DIR']}{Configuration['LOG_DIRECTORY']}/'''
    sf.create_directory(
        Configuration['LOG_DIRECTORY'], Configuration['FITTING_DIR'])
    Configuration['OUTPUTLOG'] = f"{Configuration['LOG_DIRECTORY']}{Configuration['LOG_FILE']}"

    #If it exists move the previous Log
    if os.path.exists(Configuration['OUTPUTLOG']):
        os.rename(
            Configuration['OUTPUTLOG'], f"{os.path.splitext(Configuration['OUTPUTLOG'])[0]}_Prev.txt")

    with open(Configuration['OUTPUTLOG'], 'w') as log:
        log.write(log_statement)
    #Save the original input file
    inputyml =  f'{Configuration["LOG_DIRECTORY"]}/FAT_Inputs-Run_{datetime.now().strftime("%H:%M:%S_%d-%m-%Y")}.yml'
    with open(inputyml,'w') as default_write:
        default_write.write(OmegaConf.to_yaml( Configuration['ORIGINAL_CONFIGURATION']))


    Fits_Files['FITTING_CUBE'] = f"{Configuration['BASE_NAME']}.fits"
    Fits_Files['OPTIMIZED_CUBE'] = f"{Configuration['BASE_NAME']}_opt.fits"
    Fits_Files['MOMENT0'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom0.fits"
    Fits_Files['MOMENT1'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom1.fits"
    Fits_Files['MOMENT2'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom2.fits"
    Fits_Files['MASK'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mask.fits"
    Fits_Files['CHANNEL_MAP'] = f"Sofia_Output/{Configuration['BASE_NAME']}_chan.fits"
    Fits_Files['TIR_RUN_CUBE'] = copy.deepcopy(Fits_Files['FITTING_CUBE'])
    if 'create_fat_cube' in Configuration['FITTING_STAGES']:
        if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['ORIGINAL_CUBE']}"):
            raise BadCatalogueError(
                f'''We can not find the file {Fits_Files['ORIGINAL_CUBE']}. This is likely to be due to a typo in your catalog {Configuration['CATALOGUE']}.''')
    else:
        if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}"):
            raise BadCatalogueError(
                f'''We can not find the file {Fits_Files['FITTING_CUBE']}. This is likely to be due to a typo in your catalog {Configuration['CATALOGUE']}.''')

    # run cleanup
    cf.cleanup(Configuration, Fits_Files)


    if Configuration['DEBUG']:
        from numpy import __version__ as npversion
        from scipy import __version__ as spversion
        from astropy import __version__ as apversion
        from make_moments import __version__ as mmversion
        from TRM_errors import __version__ as teversion
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            from matplotlib import __version__ as mpversion
        #from subprocess import __version__ as subversion

        print_log(f'''FAT_GALAXY_LOOPS: We are using the following versions
{'':8s}NumPy {npversion}
{'':8s}SciPy {spversion}
{'':8s}AstroPy {apversion}
{'':8s}Matplotlib {mpversion}
{'':8s}make_moments {mmversion}
{'':8s}TRM_errors {teversion}

''', Configuration)

    log_statement = f'''We are starting the catalogue entry {Configuration['ID']} in the directory {Configuration['SUB_DIR']}.\n'''
    print_log(
        log_statement, Configuration)

    if Configuration['TIMING']:
        Configuration['FAT_PSUPROCESS'] = psu.Process(
            Configuration['FAT_PID'])
        update_statistics(
            Configuration, message="Creating a CPU RAM Log for analysis.")


    # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
    if not os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"):
        ff.create_fat_cube(Configuration, Fits_Files = Fits_Files)
    # Get a bunch of info from the cube
    rf.read_cube(
        Configuration, Fits_Files['FITTING_CUBE'])
    Configuration['INITIALIZATION_TIME'][1] = datetime.now()
    enter_recovery_point(Configuration,Fits_Files=Fits_Files,\
                        message='Creating the entry point at the end of the initialization.',\
                        point_ID='AFTER-INITIALIZATION')

    return Fits_Files

def loop_sofia(Configuration,Fits_Files):
     
    Configuration['SOFIA_TIME'][0] = datetime.now()
    #If we have Sofia Preprocessed Output request make sure it all exists
    if 'external_sofia' in Configuration['FITTING_STAGES']:
        sf.copy_homemade_sofia(
            Configuration)
    elif 'run_sofia' in Configuration['FITTING_STAGES']:
        runf.sofia(Configuration, Fits_Files)
    else:
        sf.sofia_output_exists(
            Configuration, Fits_Files)
    Configuration['SOFIA_TIME'][1] = datetime.now()
    enter_recovery_point(Configuration,Fits_Files=Fits_Files,\
                        message='Creating the entry point after running sofia.',\
                        point_ID='AFTER-SOFIA')


def MP_initialize_sofia(Configuration,timing_lock,catalogue_lock):
    registered_exception = None
    current_run = 'Not Initialized'
    succes = False
    catalogue_line = 'unset'
    Initial_Parameters = []
    Fits_Files = {}
    try:
        Fits_Files = initialize_loop(Configuration)
        if Configuration['RP_SECTION'] == 'AFTER-INITIALIZATION':
            loop_sofia(Configuration,Fits_Files)
        else:
            Initial_Parameters = copy.deepcopy(Configuration['Initial_Parameters'])

        if (not Configuration['USED_FITTING'] is None) and \
            (Configuration['RP_SECTION'] == 'AFTER-SOFIA'):
            # Process the found source in sofia to set up the proper fitting and make sure source can be fitted

            Initial_Parameters = runf.check_source(Configuration,\
                Fits_Files)
            succes = True
            #sf.sofia_output_exists(Configuration,Fits_Files)
            print_log(f'''FAT_GALAXY_LOOPS: The source is well defined and we will now setup the initial tirific file
''', Configuration)
        else:
            try:
                Initial_Parameters = Configuration['Initial_Parameters']
            except:
                pass
            if Configuration['USED_FITTING'] is None:
                catalogue_line = cf.finish_galaxy(Configuration,
                            current_run=current_run,Fits_Files=Fits_Files,
                            timing_lock=timing_lock, catalogue_lock = catalogue_lock,
                            exiting=registered_exception)
    except Exception as e:
        print(f"What is going on? {e}, {e.__class__.__name__} ")
        succes =False
        if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
            Configuration['QUANTITY'] = 5
        else:
            Configuration['QUANTITY'] = 'error'
        registered_exception = e
        Configuration['FINAL_COMMENT'] = e
        catalogue_line = cf.finish_galaxy(Configuration,
                          current_run=current_run,Fits_Files=Fits_Files,
                          timing_lock=timing_lock, catalogue_lock = catalogue_lock,
                          exiting=registered_exception)

   
    sofia_output = {'Succes': succes, 'Configuration': Configuration,
                      'catalogue_line': catalogue_line, 'Fits_Files':Fits_Files,
                      'Size': Configuration['SIZE_IN_BEAMS'],
                      'Initial_Parameters': Initial_Parameters}
    
    update_statistics(
            Configuration, message="Pause at")  
    Configuration['FAT_PSUPROCESS'] = None
    return sofia_output


def MP_Fitting_Loop(input,timing_lock,catalogue_lock):
    try:
       
        Configuration = input['Configuration']
        Configuration['FAT_PSUPROCESS'] = psu.Process(
            Configuration['FAT_PID'])
        update_statistics(
            Configuration, message="Start fitting loop")  
        Initial_Parameters = input['Initial_Parameters']
        Fits_Files = input['Fits_Files']
        registered_exception = None
        current_run = 'Not Initialized'
       
            # If you add any make sure that the fitstage  starts with 'Fit_'


        if 'fit_tirific_osc' in Configuration['FITTING_STAGES']:   
            current_run,Tirific_Template = runf.fitting_osc(
                Configuration, Fits_Files, Initial_Parameters)
        elif 'fit_make_your_own' in Configuration['FITTING_STAGES']:
            print_log(f'''FAT_GALAXY_LOOPS: If you add any fitting routine make sure that the fit stage  starts with Fit_
''',Configuration)
            sf.create_directory(Configuration['USED_FITTING'],Configuration['FITTING_DIR'])
            Configuration['FINAL_COMMENT'] = 'This example does not work'
            catalogue_line = cf.finish_galaxy(Configuration)
            return catalogue_line
        else:
            Configuration['FINAL_COMMENT'] = 'You have chosen not to do any fitting'
            catalogue_line = cf.finish_galaxy(Configuration)
            return catalogue_line

        #if all the fitting has gone properly we create nice errors

        if Configuration['QUANTITY'] != 5:
            if 'tirshaker' in Configuration['FITTING_STAGES']:
                runf.tirshaker_call(
                    Configuration,Fits_Files)

            Configuration['FINAL_COMMENT'] = 'The fit has converged succesfully'

    #Only
        catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run,
                         Fits_Files=Fits_Files,timing_lock=timing_lock,
                         catalogue_lock = catalogue_lock,
                         exiting=registered_exception)
        if Configuration['QUANTITY'] != 5:
            DHI = rf.get_DHI(
                Configuration, Model=Configuration['USED_FITTING'])
            Totflux = rf.get_totflux(
                Configuration, f"/Finalmodel/Finalmodel_mom0.fits")
            wf.basicinfo(Configuration, template=Tirific_Template, Tot_Flux=Totflux,\
                DHI=[DHI, Configuration['BEAM'][0]*Configuration['RING_SIZE']])
            if Configuration['INSTALLATION_CHECK']:
                cf.installation_check(
                    Configuration)
        print_log(f'''FAT_GALAXY_LOOP: This galaxy has succesfully ran through all parts of the fitting
''', Configuration)
    except Exception as e:
        if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
            Configuration['QUANTITY'] = 5
        else:
            Configuration['QUANTITY'] = 'error'
        registered_exception = e
        Configuration['FINAL_COMMENT'] = e
        Configuration['QUANTITY'] = 'error'
        catalogue_line = cf.finish_galaxy(Configuration,
                          current_run=current_run,Fits_Files=Fits_Files,
                          timing_lock=timing_lock, catalogue_lock = catalogue_lock,
                          exiting=registered_exception)
    return catalogue_line
