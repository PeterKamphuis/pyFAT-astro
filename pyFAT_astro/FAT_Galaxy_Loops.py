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
from pyFAT_astro.Support.fat_errors import BadCatalogueError

def FAT_Galaxy_Loops(id,Proc_Configuration, Full_Catalogue,progress):
    try:
        # we can make a single line when id is set
        if id:
            Proc_Configuration['CATALOGUE_START_ID'] = id
            Proc_Configuration['CATALOGUE_END_ID'] = id+1
            Proc_Configuration['OUTPUT_CATALOGUE'] = None




        for current_galaxy_index in range(Proc_Configuration['CATALOGUE_START_ID'], Proc_Configuration['CATALOGUE_END_ID']):
            registered_exception = None
            current_run = 'Not Initialized'
            Configuration = copy.deepcopy(Proc_Configuration)
            Configuration['START_TIME'] = datetime.now()
            # First check the starttime
            Configuration['ID'] = Full_Catalogue['ID'][current_galaxy_index]
            if Full_Catalogue['DISTANCE'][current_galaxy_index] != -1.:
                Configuration['DISTANCE'] = Full_Catalogue['DISTANCE'][current_galaxy_index]
            Configuration['SUB_DIR'] = Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]
            Configuration['BASE_NAME'] = Full_Catalogue['CUBENAME'][current_galaxy_index]+'_FAT'
            if not Configuration['SOFIA_BASENAME']:
                if 'BASENAME' in Full_Catalogue['ENTRIES']:
                    Configuration['SOFIA_BASENAME'] = Full_Catalogue['BASENAME'][current_galaxy_index]
                else:
                    Configuration['SOFIA_BASENAME'] = Configuration['BASE_NAME']
            #Add our fitting directory to the Configuration
            #Maindir always ends in slash already
            if Full_Catalogue['DIRECTORYNAME'][current_galaxy_index] == './':
                Configuration['FITTING_DIR'] = f"{Configuration['MAIN_DIRECTORY']}"
            else:
                Configuration['FITTING_DIR'] = f"{Configuration['MAIN_DIRECTORY']}{Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]}/"
            if Configuration['FITTING_DIR'][-2:] == '//':
                Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'

            if not Configuration['SOFIA_DIR']:
                Configuration['SOFIA_DIR'] = Configuration['FITTING_DIR']
            #ini_mode_factor =25
            # We initially set the variations to fixed for all parameters
            #let's see what happens if we immediately

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {
                'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}

            # If we have a fitting log we start writing
            log_statement = f'''This file is a log of the fitting process run at {Configuration ['START_TIME']}.
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

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {
                'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}

            #if we skip the create_fat _cube stage peaople could give the fat cube itself

            if Fits_Files['ORIGINAL_CUBE'][-9:] == '_FAT.fits':
                if 'create_fat_cube' in Configuration['FITTING_STAGES']:
                    sf.print_log(f'''FAT_GALAXY_LOOPS: Your input cube ends in _FAT.fits indicating it is a FAT processed cube.
Therefore we remove the Create_FAT_Cube stages from the loop.
''', Configuration)
                    Configuration['FITTING_STAGES'].remove('create_fat_cube')
                fat_ext = ''
            else:
                fat_ext = '_FAT'
            Fits_Files['FITTING_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}{fat_ext}.fits"
            Fits_Files['OPTIMIZED_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}{fat_ext}_opt.fits"
            Fits_Files['MOMENT0'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom0.fits"
            Fits_Files['MOMENT1'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom1.fits"
            Fits_Files['MOMENT2'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mom2.fits"
            Fits_Files['MASK'] = f"Sofia_Output/{Configuration['BASE_NAME']}_mask.fits"
            Fits_Files['CHANNEL_MAP'] = f"Sofia_Output/{Configuration['BASE_NAME']}_chan.fits"
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

            # then we want to read the template
            Tirific_Template = sf.tirific_template()
            if Configuration['DEBUG']:
                from numpy import __version__ as npversion
                from scipy import __version__ as spversion
                from astropy import __version__ as apversion
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    from matplotlib import __version__ as mpversion
                #from subprocess import __version__ as subversion

                sf.print_log(f'''FAT_GALAXY_LOOPS: We are using the following versions
{'':8s}NumPy {npversion}
{'':8s}SciPy {spversion}
{'':8s}AstroPy {apversion}
{'':8s}Matplotlib {mpversion}
''', Configuration)

            log_statement = f'''We are starting the catalogue entry {Configuration['ID']} in the directory {Configuration['SUB_DIR']}.\n'''
            sf.print_log(
                log_statement, Configuration)

            if Configuration['TIMING']:
                Configuration['FAT_PSUPROCESS'] = psu.Process(
                    Configuration['FAT_PID'])
                sf.update_statistic(
                    Configuration, message="Creating a CPU RAM Log for analysis.")


            # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
            if not os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"):
                try:
                    ff.create_fat_cube(Configuration, Fits_Files)
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
                        Configuration['OUTPUT_QUANTITY'] = 5
                    else:
                        Configuration['OUTPUT_QUANTITY'] = 'error'
                    catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run,
                                         exiting=e)
                    continue

            # Get a bunch of info from the cube
            rf.read_cube(
                Configuration, Fits_Files['FITTING_CUBE'])

            #If we have Sofia Preprocessed Output request make sure it all exists
            if Configuration['DEBUG']:
                sf.write_config(
                    f'{Configuration["LOG_DIRECTORY"]}CFG_Before_Sofia.txt', Configuration)

            if 'external_sofia' in Configuration['FITTING_STAGES']:
                sf.copy_homemade_sofia(
                    Configuration)
            elif 'run_sofia' in Configuration['FITTING_STAGES']:
                # Run sofia2
                try:
                    runf.sofia(Configuration, Fits_Files)
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
                        Configuration['OUTPUT_QUANTITY'] = 5
                    else:
                        Configuration['OUTPUT_QUANTITY'] = 'error'
                    catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run, exiting=e)
                    continue
            else:
                sf.sofia_output_exists(
                    Configuration, Fits_Files)
                # We assume sofia is ran and created the proper files
            try:

                # If you add any make sure that the fitstage  starts with 'Fit_'
                if Configuration['USED_FITTING']:
                    # Process the found source in sofia to set up the proper fitting and make sure source can be fitted
                    Initial_Parameters = runf.check_source(
                        Configuration, Fits_Files)
                    #sf.sofia_output_exists(Configuration,Fits_Files)
                    sf.print_log(f'''FAT_GALAXY_LOOPS: The source is well defined and we will now setup the initial tirific file
''', Configuration)
                    #Add your personal fitting types here
                    if Configuration['DEBUG']:
                        sf.write_config(
                            f'{Configuration["LOG_DIRECTORY"]}CFG_Before_Fitting.txt', Configuration)

                if 'fit_tirific_osc' in Configuration['FITTING_STAGES']:
                    current_run = runf.fitting_osc(
                        Configuration, Fits_Files, Tirific_Template, Initial_Parameters)
                elif 'fit_make_your_own' in Configuration['FITTING_STAGES']:
                    sf.print_log(f'''FAT_GALAXY_LOOPS: If you add any fitting routine make sure that the fit stage  starts with Fit_
''',Configuration)
                    sf.create_directory(Configuration['USED_FITTING'],Configuration['FITTING_DIR'])
                    Configuration['FINAL_COMMENT'] = 'This example does not work'
                    catalogue_line = cf.finish_galaxy(Configuration)
                    continue
                else:
                    Configuration['FINAL_COMMENT'] = 'You have chosen not to do any fitting'
                    catalogue_line = cf.finish_galaxy(Configuration)
                    continue

                #if all the fitting has gone properly we create nice errors

                if Configuration['OUTPUT_QUANTITY'] != 5:
                    if 'tirshaker' in Configuration['FITTING_STAGES']:
                        runf.tirshaker_call(
                            Configuration)

                    Configuration['FINAL_COMMENT'] = 'The fit has converged succesfully'

            except Exception as e:
                registered_exception = e
                Configuration['FINAL_COMMENT'] = e
                if e.__class__.__name__ in Configuration['STOP_INDIVIDUAL_ERRORS']:
                    Configuration['OUTPUT_QUANTITY'] = 5
                else:
                    Configuration['OUTPUT_QUANTITY'] = 'error'
            #Only
            catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run,
                             Fits_Files=Fits_Files, exiting=registered_exception)
            if Configuration['OUTPUT_QUANTITY'] != 5:
                DHI = rf.get_DHI(
                    Configuration, Model=Configuration['USED_FITTING'])
                Totflux = rf.get_totflux(
                    Configuration, f"/Finalmodel/Finalmodel_mom0.fits")
                wf.basicinfo(Configuration, template=Tirific_Template, Tot_Flux=Totflux, DHI=[
                             DHI, Configuration['BEAM'][0]*Configuration['RING_SIZE']])
                if Configuration['INSTALLATION_CHECK']:
                    cf.installation_check(
                        Configuration)
    except Exception as e:
        registered_exception = e
        Configuration['FINAL_COMMENT'] = e
        Configuration['OUTPUT_QUANTITY'] = 'error'
        catalogue_line = cf.finish_galaxy(Configuration, current_run=current_run,Fits_Files=Fits_Files,
                          exiting=registered_exception)
    if progress:
        print(f'We are at {progress}% of processing the catalogue')
    return catalogue_line
