# -*- coding: future_fstrings -*-

# This is the python version of FAT
import sys
import os
import copy
import numpy as np

from argparse import ArgumentParser
import traceback
import warnings

from datetime import datetime
from astropy.io import fits
from astropy.wcs import WCS
import pyFAT
import pyFAT.Support.read_functions as rf
import pyFAT.Support.support_functions as sf
# Functions that run external programs such as tirific and sofia
import pyFAT.Support.run_functions as runf
# function that keep things orderly and nicely
import pyFAT.Support.clean_functions as cf
# Functions that modify or produce fat fits file
import pyFAT.Support.fits_functions as ff
#functions that write files
import pyFAT.Support.write_functions as wf
#from pyFAT.Support.constants import initialize
from  pyFAT.Support.modify_template import write_new_to_template

class MissingProgramError(Exception):
    pass
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

def main(argv):


    try:
        #First we check for sofia and TiRiFiC
        if sf.is_available('sofia2'):
            print(f'''We have found SoFiA-2 as sofia2
''')
            sofia_available = "sofia2"
        elif sf.is_available('sofia'):
            sofia_available = "sofia"
            print(f'''We have found SoFiA-2 as sofia
''')
        else:
            print(f'''We can not find sofia. Please install Sofia-2
''')
            raise MissingProgramError('We can not find sofia.')

        if sf.is_available('tirific'):
            print(f'''We have found TiRiFiC as tirific
''')
        else:
            print(f'''We can not find tirific. Please install TiRiFiC
''')
            raise MissingProgramError('We can not find TiRiFiC.')
        #Get the directory we are running from, This is for the Installation Check
        start_dir = os.getcwd()
        #Then check the input options
        parser  = ArgumentParser()
        parser.add_argument('-c','--cf','--configuration_file', action ="store" ,dest = "configfile", default = 'No Default', help = 'Define the input configuration file.')
        parser.add_argument('-d','--debug', action ="store_true" ,dest = "debug", default = False, help = 'Print debug messages')
        parser.add_argument('-s','--sg','--single_galaxy', action ="store" ,dest = "single_cube", default=f'CataloguE', help = 'If set the locations in the config file are ignored and only this cube will be fitted in the directory it is located.')
        parser.add_argument('-i','--ic','--installation_check', action ="store_true" ,dest = "installation_check", default = False, help = 'Run the installation _check.')
        parser.add_argument('--LVT','--LVHIS_TEST', action ="store_true" ,dest = "lvhis_test", default = False, help = 'Run the LVHIS Test. Developer Only.')
        parser.add_argument('--PT','--PAPER_TEST', action ="store_true" ,dest = "paper_test", default = False, help = 'Run the PAPER Test. Developer Only.')
        parser.add_argument('--FD','--FULL_DATABASE', action ="store_true" ,dest = "full_test", default = False, help = 'Run the Full Database Test. Developer Only.')
        parser.add_argument('-p','--problems', action ="store_true" ,dest = "problems", default = False, help = 'Run the Problem test set. Developer Only.')
        parser.add_argument('-t','--timing', action ="store_true" ,dest = "timing", default = False, help = 'Create a file in the maindir that provides start and stop times for each galaxy.')
        parser.add_argument('-n','--ncpu', action ="store" ,dest = "ncpu", default = 6, help = 'Number of CPUs to use.')
        input_parameters = parser.parse_args()
        basic_info  = 'BasicInfo'
        if input_parameters.installation_check:
            input_parameters.configfile = 'ChecK.ConfiG'
        if input_parameters.lvhis_test:
            fat_main_test_dir = os.environ["FAT_TEST_DIR"]
            input_parameters.configfile=fat_main_test_dir+'/LVHIS-26_3/Input.config'
        if input_parameters.paper_test:
            fat_main_test_dir = os.environ["FAT_TEST_DIR"]
            input_parameters.configfile=fat_main_test_dir+'/SmallCat_Warps/Input_1.config'
        if input_parameters.full_test:
            fat_main_test_dir = os.environ["FAT_TEST_DIR"]
            input_parameters.configfile=fat_main_test_dir+'/Full_Database/FAT_INPUT.config'
        if input_parameters.problems:
            fat_main_test_dir = os.environ["FAT_TEST_DIR"]
            input_parameters.configfile=fat_main_test_dir+'/Problems/FAT_INPUT.config'

        try:
            Original_Configuration,log_write_config = rf.config_file(input_parameters,start_dir,debug=input_parameters.debug)
        except Exception as e:
            print(e)
            exit()
        # All Configuration parameters that are not set in the config file should be set here
        # Add the starting directory to the Configuration
        # Keys read in the config file
        #required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','FIX_SBR','FIX_VROT','HANNING',\
        #                               'STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT',\
        #                               'RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR',\
        #                                'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM',\
        #                                 'MAPS_OUTPUT','WARP_OUTPUT']
        Original_Configuration['START_DIR'] = start_dir
        # Also add the timing input and some other recurring parameters
        Original_Configuration['TIMING'] = input_parameters.timing
        Original_Configuration['DEBUG'] = input_parameters.debug
        Original_Configuration['NCPU'] = input_parameters.ncpu
        Original_Configuration['FINAL_COMMENT'] = "This fitting stopped with an unregistered exit."
        if Original_Configuration['DEBUG']:
            warnings.showwarning = warn_with_traceback
        # Keys that change depending on which type of fitting is run
        loop_counters=['RUN_COUNTER']
        timing_keys = ['PREP_END_TIME','START_TIME']
        fitting_status = ['OS_ACCEPTED']
        loop_counters.append('LOOPS')
        for key in timing_keys:
            Original_Configuration[key] = 'Not completed'
        for key in loop_counters:
            Original_Configuration[key] = 0
        for key in fitting_status:
            Original_Configuration[key] = False


        boolean_keys = ['OPTIMIZED', # Are we fitting an optimized cube
                        'TIRIFIC_RUNNING', # Is there a tirific initialized
                        'OUTER_RINGS_DOUBLED', #Do the outer rings have twice the size of the inner rings
                        'NEW_RING_SIZE',  #Have we update the size of the ring while not yet updating the template
                        'VEL_SMOOTH_EXTENDED', # Is the velocity smoothing extended ????
                        'EXCLUDE_CENTRAL', # Do we exclude the central part of the fitting due to blanks/an absorption source
                        'INSTALLATION_CHECK' # Are we running an installation check
                        ]
        #
        for key in boolean_keys:
            Original_Configuration[key] = False

        other_keys =  {'MINIMUM_WARP_SIZE': 3., # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc,set here.
                       'MINIMUM_RINGS': 3,  # we need at least this amount of rings (Including 0 and 1/5 beam), set here
                       'TOO_SMALL_GALAXY': 1., # if the number of beams across the major axis/2 is less than this we will not fit the galaxy, set here

                       'DISTANCE': 'Unset', # Distance to the galaxy, set from the catalogue at start of loop
                       'ID_NR': 'Unset', # ID of the galaxy in the catalogue , set from the catalogue at start of loop
                       'SUB_DIR': 'Unset', # Name of the directory in which galaxy resides, set from the catalogue at start of loop
                       'FITTING_DIR': 'Unset', # Full path of the directory in which the fitting takes place, set at start of loop
                       'SOFIA_BASENAME': 'Unset', #Basename of pre-processed sofia products, only set when provided in catalogue at start of loop
                       'BASENAME': 'Unset', #Basename for FAT products, typically {input_cube}_FAT, set at start of loop
                       'LOG_DIR': 'Unset', #Directory to put log files from run, set at start of loop

                       'CURRENT_STAGE': 'initial', #Current stage of the fitting process, set at switiching stages
                       'TIRIFIC_PID': 'Not Initialized', #Process ID of tirific that is running
                       'ACCEPTED': False, #Whether a fit is accepted or not

                       'MAX_SIZE_IN_BEAMS': 30, # The galaxy is not allowed to extend beyond this number of beams in radius, set in check_source
                       'MIN_SIZE_IN_BEAMS': 0., # Minimum allowed radius in number of beams of the galaxy, set in check_source
                       'SIZE_IN_BEAMS': 0, # The radius of the galaxy in number of beams, adapted after running Sofia
                       'NO_RINGS': 0., # The number of rings in the fit,
                       'LAST_RELIABLE_RINGS': [0.,0.], # Location of the rings where the SBR drops below the cutoff limits, adapted after every run. Should only be modified in check_size
                       'LIMIT_MODIFIER': [1.], #Modifier for the cutoff limits based on the inclination , adapted after every run.
                       'OLD_RINGS': [], # List to keep track of the ring sizes that have been fitted.

                       'NO_POINTSOURCES': 0. , # Number of point sources, set in run_tirific

                       'INNER_FIX': [4.,4.], #Number of rings that are fixed in the inner part for the INCL and PA, , adapted after every run in get_inner_fix in support_functions and for both sides
                       'WARP_SLOPE': [0.,0.], #Ring numbers from which outwards the warping should be fitted as a slope, set in get_warp_slope in modify_template
                       'OUTER_SLOPE_START': 1, # Ring number from where the RC is fitted as a slope
                       'RC_UNRELIABLE': 1, # Ring number from where the RC values are set flat. Should only be set in check_size

                       'NOISE': 0. , #Noise of the input cube, set in main
                       'BEAM_IN_PIXELS': [0.,0.,0.], #FWHM BMAJ, BMIN in pixels and total number of pixels in beam area, set in main
                       'BEAM': [0.,0.,0.], #  FWHM BMAJ, BMIN in arcsec and BPA, set in main
                       'NAXES': [0.,0.,0.], #  Size of the cube in pixels x,y,z arranged like sane people not python, set in main
                       'NAXES_LIMITS': [[0.,0.],[0.,0.],[0.,0.]], #  Size of the cube in degree and km/s,  x,y,z arranged like sane people not python, set in main updated in cut_cubes
                       'MAX_ERROR': {}, #The maximum allowed errors for the parameters, set in main derived from cube
                       'CHANNEL_WIDTH': 0., #Width of the channel in the cube in km/s, set in main derived from cube
                       'PIXEL_SIZE': 0., #'Size of the pixels in degree'
                       }

        for key in other_keys:
            Original_Configuration[key] = other_keys[key]

        # The parameters that need boundary limits are set here
        boundary_limit_keys = ['PA','INCL', 'SDIS', 'Z0','VSYS','XPOS','YPOS']
        for key in boundary_limit_keys:
            Original_Configuration[f"{key}_CURRENT_BOUNDARY"] = [[0.,0.],[0.,0.],[0.,0.]]
        #Then read the input Catalogue
        if input_parameters.single_cube != 'CataloguE':
            Full_Catalogue = sf.Proper_Dictionary({})
            Full_Catalogue['ENTRIES'] = ['ENTRIES','NUMBER','DISTANCE','DIRECTORYNAME','CUBENAME']
            Full_Catalogue['NUMBER'] = ['0']
            Full_Catalogue['DISTANCE'] = [-1]
            Full_Catalogue['DIRECTORYNAME'] = ['./']
            Full_Catalogue['CUBENAME'] = [f"{os.path.splitext(input_parameters.single_cube.split('/')[-1])[0]}"]
        else:
            Full_Catalogue = rf.catalogue(Original_Configuration['CATALOGUE'])
        stop_individual_errors = ['SmallSourceError','BadSourceError','SofiaFaintError','BadHeaderError','BadCubeError','BadMaskError','BadCatalogueError']
        # Get the longest directory name to format the output directory properly
        dirname = 'Directory Name'
        maximum_directory_length = len(dirname)
        for directory in Full_Catalogue['DIRECTORYNAME']:
            if directory == './':
                maximum_directory_length = len(Original_Configuration['MAINDIR'].split('/')[-2])
            if len(directory) > maximum_directory_length:
                maximum_directory_length = len(directory)

        # Create a file to write the results to if if required
        if Original_Configuration['OUTPUTCATALOGUE']:
            if not os.path.exists(Original_Configuration['OUTPUTCATALOGUE']) or Original_Configuration['NEW_OUTPUT']:
                with open(Original_Configuration['OUTPUTCATALOGUE'],'w') as output_catalogue:
                    comment = 'Comments on Fit Result'
                    AC1 = 'OS'
                    output_catalogue.write(f"{dirname:<{maximum_directory_length}s} {AC1:>6s} {comment}\n")

        if Original_Configuration['TIMING']:
            timing_result = open(Original_Configuration['MAINDIR']+'Timing_Result.txt','w')
            timing_result.write("This file contains the system start and end times for the fitting of each galaxy")
            timing_result.close()
        #if start_galaxy not negative then it is catalogue ID
        if -1 != Original_Configuration['STARTGALAXY']:
            Original_Configuration['STARTGALAXY'] = np.where(Original_Configuration['STARTGALAXY'] == Full_Catalogue['NUMBER'])[0][0]
        else:
            Original_Configuration['STARTGALAXY'] = 0
        # If the end galaxy is -1 fit the whole catalogue
        if Original_Configuration['ENDGALAXY'] == -1:
            Original_Configuration['ENDGALAXY'] = len(Full_Catalogue['NUMBER'])
            if Original_Configuration['ENDGALAXY'] == 0:
                Original_Configuration['ENDGALAXY'] = 1
        else:
            Original_Configuration['ENDGALAXY'] = np.where(Original_Configuration['ENDGALAXY'] == Full_Catalogue['NUMBER'])[0][0]
        # start the main fitting loop

        if float(Original_Configuration['STARTGALAXY']) > float(Original_Configuration['ENDGALAXY']):
            print(f''' Your starting galaxy (Line nr = {Original_Configuration['STARTGALAXY']}) is listed after your ending galaxy (Line nr = {Original_Configuration['ENDGALAXY']}), maybe you have double catalogue ids?''')
            exit()


        for current_galaxy_index in range(Original_Configuration['STARTGALAXY'],Original_Configuration['ENDGALAXY']):

            Configuration = copy.deepcopy(Original_Configuration)
            Configuration['START_TIME'] = datetime.now()
            # First check the starttime
            if input_parameters.installation_check:
                Configuration['INSTALLATION_CHECK'] = True
            Configuration['ID_NR'] = Full_Catalogue['NUMBER'][current_galaxy_index]
            Configuration['DISTANCE'] = Full_Catalogue['DISTANCE'][current_galaxy_index]
            Configuration['SUB_DIR'] = Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]
            if 'BASENAME' in Full_Catalogue['ENTRIES']:
                Configuration['SOFIA_BASENAME'] = Full_Catalogue['BASENAME'][current_galaxy_index]
            Configuration['BASE_NAME'] = Full_Catalogue['CUBENAME'][current_galaxy_index]+'_FAT'
            #Add our fitting directory to the Configuration
            #Maindir always ends in slash already
            if Full_Catalogue['DIRECTORYNAME'][current_galaxy_index] == './':
                Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}"
            else:
                Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}{Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]}/"
            if Configuration['FITTING_DIR'][-2:] == '//':
                Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'

            ini_mode_factor =25
            # We initially set the variations to fixed for all parameters
            #let's see what happens if we immediately

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}

            # If we have a fitting log we start writing
            log_statement = f'''This file is a log of the fitting process run at {Configuration ['START_TIME']}.
{"":8s}This is version {pyFAT.__version__} of the program.
'''


            # Adapt configuration to hold some specifics to this galaxy
            Configuration['LOG_DIR'] = f"{Configuration['FITTING_DIR']}Logs/"
            if Configuration['OUTPUTLOG']:
                if not os.path.isdir(Configuration['LOG_DIR']):
                    os.mkdir(Configuration['LOG_DIR'])
                Configuration['OUTPUTLOG'] = f"{Configuration['LOG_DIR']}{Configuration['OUTPUTLOG']}"
                #If it exists move the previous Log
                if os.path.exists(Configuration['OUTPUTLOG']):
                    os.rename(Configuration['OUTPUTLOG'],f"{Configuration['LOG_DIR']}/Previous_Log.txt")

            with open(Configuration['OUTPUTLOG'],'w') as log:
                if log_write_config != 'Empty':
                    log.write(log_write_config)
                    log_write_config = 'Empty'
                log.write(log_statement)



            #Make a dictionary for the fitsfiles we use
            Fits_Files = {'ORIGINAL_CUBE': f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"}
            Fits_Files['FITTING_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}_FAT.fits"
            Fits_Files['OPTIMIZED_CUBE'] = f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}_FAT_opt.fits"
            Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
            Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
            Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
            Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
            Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"

            # run cleanup
            cf.cleanup(Configuration,Fits_Files)

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
                with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'w') as file:
                    file.write("Creating a CPU RAM Log for analysis. \n")
            # Check if the input cube exists
            if not os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files['ORIGINAL_CUBE']}"):
                log_statement = f'''We cannot find the cube {Fits_Files['ORIGINAL_CUBE']} in the directory {Configuration['FITTING_DIR']}.
{'':8s}We skip this galaxy.
'''
                sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)
                Configuration['FINAL_COMMENT'] = "This galaxy has no fits cube to work with, it is skipped."
                cf.finish_galaxy(Configuration,maximum_directory_length)
                traceback.print_exc()
                continue



            # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
            if not os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"):
                try:
                    ff.create_fat_cube(Configuration, Fits_Files)
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['MAPS_OUTPUT'] = 5
                    else:
                        Configuration['MAPS_OUTPUT'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])
                    continue

            # We open the header of the fitting cube and get some parameters and make a header wcs structure
            cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}")
            Configuration['NOISE'] = cube_hdr['FATNOISE']
            Configuration['CHANNEL_WIDTH'] = cube_hdr['CDELT3']/1000.
            Configuration['PIXEL_SIZE'] = np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])])
            Configuration['NAXES'] = [cube_hdr['NAXIS1'],cube_hdr['NAXIS2'], cube_hdr['NAXIS3']]
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                coordinate_frame = WCS(cube_hdr)
                xlow,ylow,zlow = coordinate_frame.wcs_pix2world(1,1,1., 1.)
                xhigh,yhigh,zhigh = coordinate_frame.wcs_pix2world(*Configuration['NAXES'], 1.)
                Configuration['NAXES_LIMITS'] = [np.sort([xlow,xhigh]),np.sort([ylow,yhigh]),np.sort([zlow,zhigh])/1000.]
            # We write the pixels per beam info to Configuration such that it is easily accesible
            beamarea=(np.pi*abs(cube_hdr['BMAJ']*cube_hdr['BMIN']))/(4.*np.log(2.))
            Configuration['BEAM_IN_PIXELS'] = [cube_hdr['BMAJ']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                               cube_hdr['BMIN']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                               beamarea/(abs(cube_hdr['CDELT1'])*abs(cube_hdr['CDELT2']))]
            # Ad the major beam to configuration as we need it in many places
            Configuration['BEAM'] = [float(cube_hdr['BMAJ']*3600.), float(cube_hdr['BMIN']*3600.),float(cube_hdr['BPA'])]
            #Let's set some maximum errors based on the input cube
            Configuration['MAX_ERROR'] = {'VROT': Configuration['CHANNEL_WIDTH']*5., \
                                          'VSYS': Configuration['CHANNEL_WIDTH']*1.5, \
                                          'SBR': Configuration['NOISE']/beamarea*Configuration['CHANNEL_WIDTH']*3.,\
                                          'PA' : 15.,\
                                          'INCL': 15.,\
                                          'SDIS': Configuration['CHANNEL_WIDTH']*2.5,\
                                          'Z0' : Configuration['BEAM'][0],\
                                          'XPOS': cube_hdr['BMAJ'],\
                                          'YPOS': cube_hdr['BMAJ'],\
            }

            #If we have Sofia Preprocessed Output request make sure it all exists

            if Configuration['START_POINT'] >= 3:
                cf.copy_homemade_sofia(Configuration,Fits_Files,debug=debug)

            else:
                # Run sofia2
                try:
                    runf.sofia(Configuration, Fits_Files,sofia_to_use=sofia_available,debug=Configuration['DEBUG'])
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['MAPS_OUTPUT'] = 5
                    else:
                        Configuration['MAPS_OUTPUT'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])
                    continue

                    # We assume sofia is ran and created the proper files

            try:
                current_run = 'Not Initialized'
                # Process the found source in sofia to set up the proper fitting and make sure source can be fitted
                Initial_Parameters = runf.check_source(Configuration, Fits_Files,debug=Configuration['DEBUG'])
                sf.sofia_output_exists(Configuration,Fits_Files)

                sf.print_log(f'''The source is well defined and we will now setup the initial tirific file
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
                if Configuration['FINISHAFTER'] == 0:
                    Configuration['FINAL_COMMENT'] = 'You have chosen to end the fitting after preprocessing and sofia.'
                    cf.finish_galaxy(Configuration,maximum_directory_length,debug=Configuration['DEBUG'])
                    continue
                #Add your personal fitting types here
                if Configuration['FITTING_TYPE'].lower() == 'one_step_convergence' or Configuration['INSTALLATION_CHECK']:
                    current_run = runf.fitting_osc(Configuration,Fits_Files,Tirific_Template,Initial_Parameters)
                #cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
                #continue
                DHI = rf.get_DHI(Configuration,Model='One_Step_Convergence',debug=Configuration['DEBUG'])
                Configuration['FINAL_COMMENT'] = 'The fit has converged succesfully'
                if Configuration['MAPS_OUTPUT'] < 4:
                    Totflux = rf.get_totflux(Configuration,f"/Finalmodel/Finalmodel_mom0.fits", debug=Configuration['DEBUG'])
                else:
                    Totflux = [float('NaN'),float('NaN')]
                wf.basicinfo(Configuration, template=Tirific_Template,Tot_Flux = Totflux, DHI = DHI,debug=Configuration['DEBUG'] )
            except Exception as e:
                Configuration['FINAL_COMMENT'] = e
                if e.__class__.__name__ in stop_individual_errors:
                    Configuration['MAPS_OUTPUT'] = 5
                else:
                    Configuration['MAPS_OUTPUT'] = 'error'

            cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run, Fits_Files =Fits_Files,debug = Configuration['DEBUG'])
            if Configuration['MAPS_OUTPUT'] != 5:
                DHI = rf.get_DHI(Configuration,Model='One_Step_Convergence',debug=Configuration['DEBUG'])
                Totflux = rf.get_totflux(Configuration,f"/Finalmodel/Finalmodel_mom0.fits", debug=Configuration['DEBUG'])
                wf.basicinfo(Configuration, template=Tirific_Template,Tot_Flux = Totflux, DHI = DHI,debug=Configuration['DEBUG'] )
                if input_parameters.installation_check:
                    cf.installation_check(Configuration,debug=Configuration['DEBUG'])
    except Exception as e:
        Configuration['FINAL_COMMENT'] = e
        Configuration['MAPS_OUTPUT'] = 'error'
        cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])

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
     pyFAT -t --cf /home/your_computer/FAT_dir/FAT_INPUT.config'
'''
