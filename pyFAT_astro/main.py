# -*- coding: future_fstrings -*-

# This is the python version of FAT
import sys
import os
import copy
import numpy as np
from omegaconf import OmegaConf,MissingMandatoryValue
import traceback
import warnings
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res
from datetime import datetime
import pyFAT_astro
import pyFAT_astro.Support.read_functions as rf
import pyFAT_astro.Support.support_functions as sf
# Functions that run external programs such as tirific and sofia
import pyFAT_astro.Support.run_functions as runf
# function that keep things orderly and nicely
import pyFAT_astro.Support.clean_functions as cf
# Functions that modify or produce fat fits file
import pyFAT_astro.Support.fits_functions as ff
#functions that write files
import pyFAT_astro.Support.write_functions as wf
#from pyFAT.Support.constants import initialize
from pyFAT_astro.Support.modify_template import write_new_to_template
from pyFAT_astro.config.defaults import defaults
class MissingProgramError(Exception):
    pass
class CatalogError(Exception):
    pass
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))


# String syntax ''' '''for multiline strings. " " for string without break and ' ' for indexing dictionaries

def main(argv):

    warnings.showwarning = warn_with_traceback
    try:




        #Get the directory we are running from, This is for the Installation Check
        start_dir = os.getcwd()
        #Get default settings
        if '-v' in argv or '--version' in argv:
            print(f"This is version {pyFAT_astro.__version__} of the program.")
            sys.exit()

        if '-h' in argv or '--help' in argv:
            print('''
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

        ''')
            sys.exit()

        cfg = OmegaConf.structured(defaults)

        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        inputconf = OmegaConf.from_cli(argv)
        cfg_input = OmegaConf.merge(cfg,inputconf)
        if cfg_input.print_examples:
            no_cube = OmegaConf.masked_copy(cfg, ['input','output','fitting'])
            with open('FAT_defaults.yml','w') as default_write:
                default_write.write(OmegaConf.to_yaml(no_cube))
            my_resources = import_res.files('pyFAT_astro.config')
            data = (my_resources / 'FAT_Input_Catalogue.txt').read_bytes()
            with open('FAT_Example_Catalogue.txt','w+b') as default_write:
                default_write.write(data)

            print(f'''We have printed the file FAT_defaults.yml FAT_Input_Catalogue.txt in {os.getcwd()}.
Exiting moments.''')
            sys.exit()

        if cfg_input.configuration_file:
            succes = False
            while not succes:
                try:
                    yaml_config = OmegaConf.load(cfg_input.configuration_file)
            #merge yml file with defaults
                    cfg = OmegaConf.merge(cfg,yaml_config)
                    succes = True
                except:
                    inputconf.configuration_file = input(f'''
                            You have provided a config file ({inputconf.configuration_file}) but it can't be found.
                            If you want to provide a config file please give the correct name.
                            Else press CTRL-C to abort.
configuration_file = ''')

        cfg = OmegaConf.merge(cfg,inputconf)
        #Add none user mutable input
        #OmegaConf.update(cfg, 'start_directory', f'{os.getcwd()}', force_add=True)

        Original_Configuration = sf.setup_configuration(cfg)

        #First we check for sofia and TiRiFiC
        Original_Configuration['SOFIA2'] = sf.find_program(Original_Configuration['SOFIA2'], "SoFiA 2")
        Original_Configuration['TIRIFIC'] = sf.find_program(Original_Configuration['TIRIFIC'], "TiRiFiC")

        if cfg.cube_name:
            Full_Catalogue = sf.Proper_Dictionary({})
            Full_Catalogue['ENTRIES'] = ['ENTRIES','NUMBER','DISTANCE','DIRECTORYNAME','CUBENAME']
            Full_Catalogue['NUMBER'] = ['0']
            Full_Catalogue['DISTANCE'] = [-2.]
            Full_Catalogue['DIRECTORYNAME'] = ['./']
            Full_Catalogue['CUBENAME'] = [f"{os.path.splitext(cfg.cube_name.split('/')[-1])[0]}"]
        elif 'sofia_catalogue' in Original_Configuration['FITTING_STAGES']:
            Full_Catalogue = rf.sofia_input_catalogue(Original_Configuration)
        else:
            Full_Catalogue = rf.catalogue(Original_Configuration['CATALOGUE'])
        stop_individual_errors = ['SmallSourceError','BadSourceError','SofiaFaintError','BadHeaderError','BadCubeError','BadMaskError','BadCatalogueError']
        # Get the longest directory name to format the output directory properlyFit_Tirific_OSC
        dirname = 'Directory Name'
        maximum_directory_length = len(dirname)
        for directory in Full_Catalogue['DIRECTORYNAME']:
            if directory == './':
                directory = Original_Configuration['MAIN_DIRECTORY'].split('/')[-2]
            if len(directory) > maximum_directory_length:
                maximum_directory_length = len(directory)

        # Create a file to write the results to if if required
        if Original_Configuration['OUTPUT_CATALOGUE']:
            if not os.path.exists(Original_Configuration['OUTPUT_CATALOGUE']) or Original_Configuration['NEW_OUTPUT']:
                with open(Original_Configuration['OUTPUT_CATALOGUE'],'w') as output_catalogue:
                    comment = 'Comments on Fit Result'
                    AC1 = 'OS'
                    output_catalogue.write(f"{dirname:<{maximum_directory_length}s} {AC1:>6s} {comment}\n")

        if Original_Configuration['TIMING']:
            timing_result = open(Original_Configuration['MAIN_DIRECTORY']+'Timing_Result.txt','w')
            timing_result.write("This file contains the system start and end times for the fitting of each galaxy")
            timing_result.close()
        #if start_galaxy not negative then it is catalogue ID
        if -1 != Original_Configuration['CATALOGUE_START_ID']:
            Original_Configuration['CATALOGUE_START_ID'] = np.where(Original_Configuration['CATALOGUE_START_ID'] == Full_Catalogue['NUMBER'])[0][0]
        else:
            Original_Configuration['CATALOGUE_START_ID'] = 0
        # If the end galaxy is -1 fit the whole catalogue
        if Original_Configuration['CATALOGUE_END_ID'] == -1:
            Original_Configuration['CATALOGUE_END_ID'] = len(Full_Catalogue['NUMBER'])
            if Original_Configuration['CATALOGUE_END_ID'] == 0:
                Original_Configuration['CATALOGUE_END_ID'] = 1
        else:
            Original_Configuration['CATALOGUE_END_ID'] = np.where(Original_Configuration['CATALOGUE_END_ID'] == Full_Catalogue['NUMBER'])[0][0]
        # start the main fitting loop

        if float(Original_Configuration['CATALOGUE_START_ID']) > float(Original_Configuration['CATALOGUE_END_ID']):
            raise CatalogError(f''' Your starting galaxy (Line nr = {Original_Configuration['CATALOGUE_START_ID']}) is listed after your ending galaxy (Line nr = {Original_Configuration['CATALOGUE_END_ID']}), maybe you have double catalogue ids?''')
            sys.exit(1)


        for current_galaxy_index in range(Original_Configuration['CATALOGUE_START_ID'],Original_Configuration['CATALOGUE_END_ID']):
            registered_exception = None
            current_run = 'Not Initialized'
            Configuration = copy.deepcopy(Original_Configuration)
            Configuration['START_TIME'] = datetime.now()
            # First check the starttime
            Configuration['ID_NR'] = Full_Catalogue['NUMBER'][current_galaxy_index]
            if  Full_Catalogue['DISTANCE'][current_galaxy_index] != 2.:
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
            ini_mode_factor =25
            # We initially set the variations to fixed for all parameters
            #let's see what happens if we immediately

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}

            # If we have a fitting log we start writing
            log_statement = f'''This file is a log of the fitting process run at {Configuration ['START_TIME']}.
{"":8s}This is version {pyFAT_astro.__version__} of the program.
'''


            # Adapt configuration to hold some specifics to this galaxy

            Configuration['LOG_DIRECTORY']=f'''{Configuration['FITTING_DIR']}{Configuration['LOG_DIRECTORY']}/'''
            sf.create_directory(Configuration['LOG_DIRECTORY'],Configuration['FITTING_DIR'])
            Configuration['OUTPUTLOG'] = f"{Configuration['LOG_DIRECTORY']}{Configuration['LOG_FILE']}"

                #If it exists move the previous Log
            if os.path.exists(Configuration['OUTPUTLOG']):
                os.rename(Configuration['OUTPUTLOG'],f"{Configuration['LOG_DIRECTORY']}Previous_Log.txt")

            with open(Configuration['OUTPUTLOG'],'w') as log:
                log.write(log_statement)

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}
            #if we skip the create_fat _cube stage peaople could give the fat cube itself

            if Fits_Files['ORIGINAL_CUBE'][-9:] == '_FAT.fits':
                sf.print_log(f''' Your input cube ends in _FAT.fits indicating it is a FAT processed cube.
''',Configuration['OUTPUTLOG'],screen=True,debug=Configuration['DEBUG'])
                if 'create_fat_cube' in Configuration['FITTING_STAGES']:
                    sf.print_log(f''' We are remove the Create_FAT_Cube stages from the loop.
''',Configuration['OUTPUTLOG'],screen=True,debug=Configuration['DEBUG'])
                    Configuration['FITTING_STAGES'].remove('create_fat_cube')
                fat_ext = ''
            else:
                fat_ext = '_FAT'
            Fits_Files['FITTING_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}{fat_ext}.fits"
            Fits_Files['OPTIMIZED_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}{fat_ext}_opt.fits"
            Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
            Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
            Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
            Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
            Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"
            if 'create_fat_cube' in Configuration['FITTING_STAGES']:
                if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['ORIGINAL_CUBE']}"):
                    raise CatalogError(f'''We can not find the file {Fits_Files['ORIGINAL_CUBE']}. This is likely to be due to a typo in your catalog.''' )
            else:
                if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}"):
                    raise CatalogError(f'''We can not find the file {Fits_Files['FITTING_CUBE']}. This is likely to be due to a typo in your catalog.''' )

            # run cleanup
            cf.cleanup(Configuration,Fits_Files,debug=Configuration['DEBUG'])

            # then we want to read the template
            Tirific_Template = rf.tirific_template()
            if Configuration['DEBUG']:
                from numpy import __version__ as npversion
                from scipy import __version__ as spversion
                from astropy import __version__ as apversion
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    from matplotlib import __version__ as mpversion
                #from subprocess import __version__ as subversion

                sf.print_log(f'''MAIN: We are using the following versions
{'':8s}NumPy {npversion}
{'':8s}SciPy {spversion}
{'':8s}AstroPy {apversion}
{'':8s}Matplotlib {mpversion}
''',Configuration['OUTPUTLOG'])



            log_statement = f'''We are in loop {current_galaxy_index}. This is catalogue number {Configuration['ID_NR']} and the directory {Configuration['SUB_DIR']}.\n'''
            sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)



            if Configuration['TIMING']:
                with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'w') as file:
                    file.write("Creating a CPU RAM Log for analysis. \n")
            # Check if the input cube exists




            # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
            if not os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"):
                try:
                    ff.create_fat_cube(Configuration, Fits_Files)
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['OUTPUT_QUANTITY'] = 5
                    else:
                        Configuration['OUTPUT_QUANTITY'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'],exiting=e)
                    continue

            # Get a bunch of info from the cube
            rf.read_cube(Configuration,Fits_Files['FITTING_CUBE'],debug =Configuration['DEBUG'] )


            #If we have Sofia Preprocessed Output request make sure it all exists
            if Configuration['DEBUG']:
                wf.write_config(f'{Configuration["LOG_DIRECTORY"]}CFG_Before_Sofia.txt',Configuration,debug = True)
            if 'existing_sofia' in  Configuration['FITTING_STAGES']:
                sf.copy_homemade_sofia(Configuration,debug=Configuration['DEBUG'])
            elif 'run_sofia' in Configuration['FITTING_STAGES']:
                # Run sofia2
                try:
                    runf.sofia(Configuration, Fits_Files,debug=Configuration['DEBUG'])
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['OUTPUT_QUANTITY'] = 5
                    else:
                        Configuration['OUTPUT_QUANTITY'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'],exiting=e)
                    continue
            else:
                sf.sofia_output_exists(Configuration,Fits_Files, debug = Configuration['DEBUG'])
                    # We assume sofia is ran and created the proper files
            try:

                # Process the found source in sofia to set up the proper fitting and make sure source can be fitted
                Initial_Parameters = runf.check_source(Configuration, Fits_Files,debug=Configuration['DEBUG'])
                #sf.sofia_output_exists(Configuration,Fits_Files)

                sf.print_log(f'''The source is well defined and we will now setup the initial tirific file
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
                #Add your personal fitting types here
                wf.write_config(f'{Configuration["LOG_DIRECTORY"]}CFG_Before_Fitting.txt',Configuration,debug = True)
                # If you add any make sure that the fitstage  starts with 'Fit_'
                if 'fit_tirific_osc' in Configuration['FITTING_STAGES']:
                    current_run = runf.fitting_osc(Configuration,Fits_Files,Tirific_Template,Initial_Parameters)
                elif 'fit_make_your_own' in Configuration['FITTING_STAGES']:
                    print_log(f'If you add any fiiting routine make sure that the fit stage  starts with Fit_')
                    Configuration['FINAL_COMMENT'] = 'This example does not work'
                    cf.finish_galaxy(Configuration,maximum_directory_length,debug=Configuration['DEBUG'])
                    continue
                else:
                    Configuration['FINAL_COMMENT'] = 'You have chosen not to do any fitting'
                    cf.finish_galaxy(Configuration,maximum_directory_length,debug=Configuration['DEBUG'])
                    continue
                #cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
                #continue
                Configuration['FINAL_COMMENT'] = 'The fit has converged succesfully'


            except Exception as e:
                registered_exception = e
                Configuration['FINAL_COMMENT'] = e
                if e.__class__.__name__ in stop_individual_errors:
                    Configuration['OUTPUT_QUANTITY'] = 5
                else:
                    Configuration['OUTPUT_QUANTITY'] = 'error'
            cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run, Fits_Files =Fits_Files,debug = Configuration['DEBUG'],exiting=registered_exception)
            if Configuration['OUTPUT_QUANTITY'] != 5:
                DHI = rf.get_DHI(Configuration,Model=Configuration['USED_FITTING'],debug=Configuration['DEBUG'])
                Totflux = rf.get_totflux(Configuration,f"/Finalmodel/Finalmodel_mom0.fits", debug=Configuration['DEBUG'])
                wf.basicinfo(Configuration, template=Tirific_Template,Tot_Flux = Totflux, DHI = [DHI,Configuration['BEAM'][0]*Configuration['RING_SIZE']],debug=Configuration['DEBUG'] )
                if Configuration['INSTALLATION_CHECK']:
                    cf.installation_check(Configuration,debug=Configuration['DEBUG'])
    except Exception as e:
        registered_exception = e
        Configuration['FINAL_COMMENT'] = e
        Configuration['OUTPUT_QUANTITY'] = 'error'
        cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'],exiting= registered_exception)

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
