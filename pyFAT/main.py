# -*- coding: future_fstrings -*-

# This is the python version of FAT
import sys
import os
import copy
import numpy as np
from optparse import OptionParser
import traceback
from datetime import datetime
from astropy.io import fits
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
from  pyFAT.Support.modify_template import write_new_to_template,flatten_the_curve

def main(argv):
    try:

        #Get the directory we are running from, This is for the Installation Check
        start_dir = os.getcwd()

        #Constants that are used in the code
        global H_0
        H_0 = 70. #km/s/Mpc #Hubble constant

        #Then check the input options
        parser  = OptionParser()
        parser.add_option('-c','--cf','--configuration_file', action ="store" ,dest = "configfile", default = 'FAT_INPUT.config', help = 'Define the input configuration file.',metavar='CONFIGURATION_FILE')
        parser.add_option('-d','--debug', action ="store_true" ,dest = "debug", default = False, help = 'Print debug messages',metavar = '')
        #parser.add_option('-s','--sd','--support_directory', action ="store" ,dest = "supportdir", default=f'{start_dir}/Support', help = 'location where the support files reside. Only required when FAT is not started in the directory where the Support dir resides.',metavar='SUPPORT_DIR')
        parser.add_option('-i','--ic','--installation_check', action ="store_true" ,dest = "installation_check", default = False, help = 'Run the installation _check.',metavar = '')
        parser.add_option('--LVT','--LVHIS_TEST', action ="store_true" ,dest = "lvhis_test", default = False, help = 'Run the LVHIS Test. Developer Only.')
        parser.add_option('--PT','--PAPER_TEST', action ="store_true" ,dest = "paper_test", default = False, help = 'Run the PAPER Test. Developer Only.')
        parser.add_option('--FD','--FULL_DATABASE', action ="store_true" ,dest = "full_test", default = False, help = 'Run the Full Database Test. Developer Only.')
        parser.add_option('-p','--problems', action ="store_true" ,dest = "problems", default = False, help = 'Run the Problem test set. Developer Only.')
        parser.add_option('-t','--timing', action ="store_true" ,dest = "timing", default = False, help = 'Create a file in the maindir that provides start and stop times for each galaxy.')
        parser.add_option('-n','--ncpu', action ="store" ,dest = "ncpu", default = 6, help = 'Number of CPUs to use.')
        input_parameters,args = parser.parse_args()

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

        #Add the support dir to the system path and read the functions from there
        #sys.path.insert(1, input_parameters.supportdir)
        # Functions that read files
        #import read_functions as rf
        # functions that are used often for menial tasks
        #import support_functions as sf
        # Functions that run external programs such as tirific and sofia
        #import run_functions as runf
        # function that keep things orderly and nicely
        #import clean_functions as cf
        # Functions that modify or produce fat fits file
        #import fits_functions as ff
        #functions that write files
        #import write_functions as wf
        #from modify_template import write_new_to_template,flatten_the_curve
        #Check the existence of the config file and read it

        try:
            Original_Configuration = rf.config_file(input_parameters,start_dir)
        except Exception as e:
            print(e)
            exit()
        # Add the starting directory to the Configuration
        Original_Configuration['START_DIR'] = start_dir
        # Also add the timing input and some other recurring parameters
        Original_Configuration['TIMING'] = input_parameters.timing
        Original_Configuration['DEBUG'] = input_parameters.debug
        Original_Configuration['NCPU'] = input_parameters.ncpu
        # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc
        Original_Configuration['MINIMUM_WARP_SIZE'] = 4.
        Original_Configuration['MINIMUM_RINGS'] = 3. # we need at least this amount of rings (Including 0 and 1/5 beam)
        # if the number of beams across the major axis/2 is less than this we will not fit the galaxy
        Original_Configuration['TOO_SMALL_GALAXY'] = 1.
        Original_Configuration['FINAL_COMMENT'] = "This fitting stopped with an unregistered exit."
        Original_Configuration['PREP_END_TIME'] = 'Not completed'
        Original_Configuration['CC_END_TIME'] = 'Not completed'
        Original_Configuration['EC_END_TIME'] = 'Not completed'
        Original_Configuration['CC_LOOPS'] = 0
        Original_Configuration['OS_LOOPS'] = 0
        Original_Configuration['EC_LOOPS'] = 0
        Original_Configuration['CC_ACCEPTED'] = False
        Original_Configuration['EC_ACCEPTED'] = False
        Original_Configuration['OS_ACCEPTED'] = False
        Original_Configuration['CURRENT_STAGE'] = 'initial'
        #Add some tracking paramaters
        Original_Configuration['MAX_RINGS'] = 0
        Original_Configuration['NEW_RING_SIZE'] = False
        Original_Configuration['OPTIMIZED'] = False
        Original_Configuration['OUTER_RINGS_DOUBLED'] = False
        Original_Configuration['TIRIFIC_RUNNING'] = False
        Original_Configuration['VEL_SMOOTH_EXTENDED'] = False
        Original_Configuration['TIRIFIC_PID'] = 'Not Initialized'
        Original_Configuration['RUN_COUNTER'] = 0
        Original_Configuration['LIMIT_MODIFIER'] = [1.]
        Original_Configuration['INNER_FIX'] = 3
        Original_Configuration['WARP_SLOPE'] = [0.,0.]
        Original_Configuration['OUTER_SLOPE'] = 1
        Original_Configuration['OLD_RINGS'] = []
        #Then read the input Catalogue
        Full_Catalogue = rf.catalogue(Original_Configuration['CATALOGUE'])
        stop_individual_errors = ['SmallSourceError','BadSourceError','SofiaFaintError','BadHeaderError','BadCubeError']
        # Get the longest directory name to format the output directory properly
        dirname = 'Directory Name'
        maximum_directory_length = len(dirname)
        for directory in Full_Catalogue['DIRECTORYNAME']:
            if directory == './':
                maximum_directory_length = len(Original_Configuration['MAINDIR'].split('/')[-2])
            if len(directory) > maximum_directory_length:
                maximum_directory_length = len(directory)

        # Create a file to write the results to if if required
        if not os.path.exists(Original_Configuration['OUTPUTCATALOGUE']) or Original_Configuration['NEW_OUTPUT']:
            with open(Original_Configuration['OUTPUTCATALOGUE'],'w') as output_catalogue:
                comment = 'Comments on Fit Result'
                if Original_Configuration['TWO_STEP']:
                    AC1 = 'CC'
                    AC2 = 'EC'
                    output_catalogue.write(f"{dirname:<{maximum_directory_length}s} {AC1:>6s} {AC2:>6s} {comment}\n")
                else:
                    AC1 = 'OS'
                    output_catalogue.write(f"{dirname:<{maximum_directory_length}s} {AC1:>6s} {comment}\n")

        if Original_Configuration['TIMING']:
            timing_result = open(Original_Configuration['MAINDIR']+'/Timing_Result.txt','w')
            timing_result.write("This file contains the system start and end times for the fitting of each galaxy")
            timing_result.close()
        #if start_galaxy not negative then it is catalogue ID
        if 0 <= Original_Configuration['STARTGALAXY']:
            Original_Configuration['STARTGALAXY'] = np.where(Original_Configuration['STARTGALAXY'] == Full_Catalogue['NUMBER'])[0][0]
        else:
            Original_Configuration['STARTGALAXY'] = 0
        # If the end galaxy is -1 fit the whole catalogue
        if Original_Configuration['ENDGALAXY'] == -1:
            Original_Configuration['ENDGALAXY'] = len(Full_Catalogue['NUMBER'])-1
            if Original_Configuration['ENDGALAXY'] == 0:
                Original_Configuration['ENDGALAXY'] = 1
        # start the main fitting loop
        for current_galaxy_index in range(Original_Configuration['STARTGALAXY'],Original_Configuration['ENDGALAXY']):
            Configuration = copy.deepcopy(Original_Configuration)
            Configuration['START_TIME'] = datetime.now()
            # First check the starttime
            Catalogue = {}
            for key in Full_Catalogue:
                if key != 'ENTRIES':
                    Catalogue[key] = Full_Catalogue[key][current_galaxy_index]
                else:
                    Catalogue[key] = Full_Catalogue[key]
            if Configuration['DEBUG']:
                    print(current_galaxy_index)
                    sf.print_log(f'''Catalogue:
    {'':8s}{Catalogue}
    ''',None,screen=True, debug = True )
            Configuration['DISTANCE'] = Catalogue['DISTANCE']
            doubled = 0
            ini_mode_factor =25
            # We initially set the variations to fixed for all parameters
            #let's see what happens if we immediately

            #Make a dictionary for the fitsfiles we use
            Fits_Files = {'ORIGINAL_CUBE': f"{Catalogue['CUBENAME']}.fits"}

            # If we have a fitting log we start writing
            log_statement = f'''This file is a log of the fitting process run at {Configuration ['START_TIME']}.
{"":8s}This is version {pyFAT.__version__} of the program.
'''


            # Adapt configuration to hold some specifics to this galaxy
            if Catalogue['DIRECTORYNAME'] == './':
                Configuration['LOG_DIR'] = f"{Configuration['MAINDIR']}Logs/"
            else:
                Configuration['LOG_DIR'] = f"{Configuration['MAINDIR']}{Catalogue['DIRECTORYNAME']}/Logs/"
            if Configuration['OUTPUTLOG']:
                if Catalogue['DIRECTORYNAME'] == './':
                    if not os.path.isdir(f"{Configuration['MAINDIR']}Logs/"):
                        os.mkdir(f"{Configuration['MAINDIR']}Logs/")
                    Configuration['OUTPUTLOG'] = f"{Configuration['MAINDIR']}Logs/{Configuration['OUTPUTLOG']}"
                else:
                    if not os.path.isdir(f"{Configuration['MAINDIR']}{Catalogue['DIRECTORYNAME']}/Logs/"):
                        os.mkdir(f"{Configuration['MAINDIR']}{Catalogue['DIRECTORYNAME']}/Logs/")
                    Configuration['OUTPUTLOG'] = f"{Configuration['MAINDIR']}{Catalogue['DIRECTORYNAME']}/Logs/{Configuration['OUTPUTLOG']}"
                #If it exists move the previous Log
                if os.path.exists(Configuration['OUTPUTLOG']):
                    os.rename(Configuration['OUTPUTLOG'],f"{Configuration['LOG_DIR']}/Previous_Log.txt")
            sf.print_log(log_statement,Configuration['OUTPUTLOG'])


            # Adapt configuration to hold some specifics to this galaxy
            #Never use the original cube only a fat modified one
            if 'BASENAME' in Catalogue['ENTRIES']:
                Configuration['BASE_NAME'] = Catalogue['BASENAME']
            else:
                Configuration['BASE_NAME'] = Catalogue['CUBENAME']+'_FAT'
            #Fits_Files['NOISEMAP'] = f"{Configuration['BASE_NAME']}_noisemap.fits"
            Fits_Files['FITTING_CUBE'] = f"{Catalogue['CUBENAME']}_FAT.fits"
            Fits_Files['OPTIMIZED_CUBE'] = f"{Catalogue['CUBENAME']}_FAT_opt.fits"
            Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
            Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
            Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
            Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
            Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"

            #Add our fitting directory to the Configuration
            if Catalogue['DIRECTORYNAME'] == './':
                Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/"
            else:
                Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/{Catalogue['DIRECTORYNAME']}/"
            if Configuration['FITTING_DIR'][-2:] == '//':
                Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'
            # run cleanup
            cf.cleanup(Configuration,Fits_Files)
            # cleanup removes any existing log file so now we add the start statement

            # then we want to read the template
            Tirific_Template = rf.tirific_template()

            log_statement = f'''We are in loop {current_galaxy_index}. This is catalogue number {Catalogue['NUMBER']} and the directory {Catalogue['DIRECTORYNAME']}.'''
            sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)
            #Run cleanup
            cf.cleanup(Configuration,Fits_Files)




            if Configuration['TIMING']:
                with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'w') as file:
                    file.write("Creating a CPU RAM Log for analysis. \n")
            # Check if the input cube exists
            if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['ORIGINAL_CUBE']}"):
                log_statement = f'''We cannot find the cube {Fits_Files['ORIGINAL_CUBE']} in the directory {Configuration['FITTING_DIR']}.
             We skip this galaxy.
    '''
                sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)
                Configuration['FINAL_COMMENT'] = "This galaxy has no fits cube to work with, it is skipped."
                cf.finish_galaxy(Configuration,maximum_directory_length)
                traceback.print_exc()
                continue



            # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
            if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}"):
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
            cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}")
            Configuration['NOISE'] = cube_hdr['FATNOISE']
            # We write the pixels per beam info to Configuration such that it is easily accesible
            beamarea=(np.pi*abs(cube_hdr['BMAJ']*cube_hdr['BMIN']))/(4.*np.log(2.))
            Configuration['PIX_PER_BEAM'] = beamarea/(abs(cube_hdr['CDELT1'])*abs(cube_hdr['CDELT2']))
            # Ad the major beam to configuration as we need it in many places
            Configuration['BMMAJ'] = float(cube_hdr['BMAJ']*3600.)
            #If we have Sofia Preprocessed Output request make sure it all exists

            if Configuration['START_POINT'] >= 3:
                try:
                    sf.sofia_output_exists(Configuration,Fits_Files)
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['MAPS_OUTPUT'] = 5
                    else:
                        Configuration['MAPS_OUTPUT'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])
                    continue
            else:
                # Run sofia2
                try:
                    runf.sofia(Configuration, Fits_Files,cube_hdr,debug=Configuration['DEBUG'])
                    # check that all is well
                    sf.sofia_output_exists(Configuration,Fits_Files,debug=Configuration['DEBUG'])
                except Exception as e:
                    Configuration['FINAL_COMMENT'] = e
                    if e.__class__.__name__ in stop_individual_errors:
                        Configuration['MAPS_OUTPUT'] = 5
                    else:
                        Configuration['MAPS_OUTPUT'] = 'error'
                    cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])
                    continue

                    # We assume sofia is ran and created the proper files
            allowed_loops = 15
            if input_parameters.installation_check:
                allowed_loops = 1
            try:
                current_run = 'Not Initialized'
                # Process the found source in sofia to set up the proper fitting and make sure source can be fitted
                Initial_Parameters = runf.check_source(Configuration, Fits_Files, Catalogue, cube_hdr,debug=Configuration['DEBUG'])

                sf.print_log(f'''The source is well defined and we will now setup the initial tirific file
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
                if Configuration['FINISHAFTER'] == 0:
                    Configuration['FINAL_COMMENT'] = 'You have chosen to end the fitting after preprocessing and sofia.'
                    cf.finish_galaxy(Configuration,maximum_directory_length,debug=Configuration['DEBUG'])
                    continue
                if not Configuration['TWO_STEP']:
                    if not os.path.isdir(Configuration['FITTING_DIR']+'/One_Step_Convergence'):
                        os.mkdir(Configuration['FITTING_DIR']+'/One_Step_Convergence')
                    wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                                            cube_hdr,Initial_Parameters= Initial_Parameters,fit_stage='One_Step_Convergence',debug=Configuration['DEBUG'])
                    sf.print_log(f'''The initial def file is written and we will now start fitting.
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
                    Configuration['PREP_END_TIME'] = datetime.now()
                    current_run = 'Not Initialized'
                        # If we have no directory to put the output we create it

                    while not Configuration['OS_ACCEPTED'] and Configuration['OS_LOOPS'] < allowed_loops:
                        Configuration['OS_LOOPS'] = Configuration['OS_LOOPS']+1
                        sf.print_log(f'''We are starting loop {Configuration['OS_LOOPS']} of trying to converge the center and extent.
''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                        # Run the step
                        current_run = runf.one_step_converge(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,debug = Configuration['DEBUG'],allowed_loops = allowed_loops)


                    if Configuration['OS_ACCEPTED']:
                        sf.print_log(f'''The model has converged in center and extent and we make a smoothed version.
''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                        current_run = runf.fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,stage = 'after_os', fit_stage = 'One_Step_Convergence',debug = Configuration['DEBUG'])
                        if Configuration['OPTIMIZED']:
                            runf.make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_stage = 'One_Step_Convergence',debug=Configuration['DEBUG'])
                    elif input_parameters.installation_check:
                        sf.print_log(f'''The Installation_check has run a fit suvccessfully.
''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                    else:
                        Configuration['FINAL_COMMENT'] = 'We could not converge on the extent or centre of the galaxy'
                        Configuration['MAPS_OUTPUT'] = 5
                        cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
                        continue
                    Configuration['OS_END_TIME'] = datetime.now()
                else:
                    if Configuration['START_POINT'] < 4:
                        #We first fixe the variations
                        Configuration['FIX_INCLINATION'][0] = True
                        Configuration['FIX_SDIS'][0] = True
                        Configuration['FIX_PA'][0] = True
                        Configuration['FIX_Z0'][0] = True
                        # setup the first def file to be used in the first loop
                        wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                                                cube_hdr,Initial_Parameters= Initial_Parameters,fit_stage='Centre_Convergence',debug=Configuration['DEBUG'])
                        sf.print_log(f'''The initial def file is written and we will now start fitting.
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
                        Configuration['PREP_END_TIME'] = datetime.now()
                        current_run = 'Not Initialized'
                        # If we have no directory to put the output we create it
                        if not os.path.isdir(Configuration['FITTING_DIR']+'/Centre_Convergence'):
                            os.mkdir(Configuration['FITTING_DIR']+'/Centre_Convergence')
                        #We skip the first fit atm
                        #Configuration['CC_ACCEPTED'] = True
                        #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}Cen_Conv.def", Tirific_Template)
                        #Upto here should be removed for real code
                        while not Configuration['CC_ACCEPTED'] and Configuration['CC_LOOPS'] < 10:
                            Configuration['CC_LOOPS'] = Configuration['CC_LOOPS']+1
                            sf.print_log(f'''We are starting loop {Configuration['CC_LOOPS']} of trying to converge the center.
        ''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                            current_run = runf.central_converge(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,Initial_Parameters, debug = Configuration['DEBUG'])

                        if Configuration['CC_ACCEPTED']:
                            sf.print_log(f''' The center has converged and we will adjust the smoothed profile and start to adjust the size of the galaxy.
        ''',Configuration['OUTPUTLOG'],screen =True,debug=Configuration['DEBUG'])
                        else:
                            sf.print_log(f''' We could not find a stable center for the the initial stages. We will now try while adapting the the size of the model.
        ''',Configuration['OUTPUTLOG'],screen =True,debug=Configuration['DEBUG'])

                        #Then we want to make a smoothed version that can be adapted
                        #current_run = runf.fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,stage = 'after_cc', fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])
                        #incl = rf.load_tirific(f"{Configuration['FITTING_DIR']}Centre_Convergence/Centre_Convergence.def",Variables = ['INCL'])
                        #sf.print_log(f'''BEFORE_CHECK_INCLINATION: CC_loops = {Configuration['CC_LOOPS']}
        #{'':8s} Incl = {incl}
        #{'':8s} Size in beams =  {Configuration['SIZE_IN_BEAMS']})
        #''',Configuration['OUTPUTLOG'],debug = Configuration['DEBUG'])
                        #if float(incl[0][0]) < 40.:
                            #If our fit stage is after cc we want to make sure we do an extra check on low inclinations or small Galaxies
                        #    runf.check_inclination(Configuration,Tirific_Template,Fits_Files,fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])

                        #if Configuration['OPTIMIZED']:
                        #    runf.make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])
                    else:
                        current_run = 'Not initialized'
                        write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}Centre_Convergence/Centre_Convergence.def", Tirific_Template, \
                                             Variables = ['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                                                             'RADI','INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2', 'CFLUX', 'CFLUX_2', 'NUR', 'CONDISP',\
                                                             'BMIN','BMAJ','RMS','BPA','NCORES','INIMODE', 'VARY','VARINDX','MODERATE','DELEND','DELSTART'\
                                                             ,'MINDELTA','PARMAX','PARMIN','DISTANCE','INSET'],debug=Configuration['DEBUG'])

                    Configuration['CC_END_TIME'] = datetime.now()


                    #If we only care about a centrally converged galaxy we stop here

                    # if our current run is not broken then we want to stop it
                    sf.finish_current_run(Configuration,current_run,debug= Configuration['DEBUG'])
                    # write the new values to the basic info file
                        #Write the info to the Basic info File

                    wf.basicinfo(Configuration,first_fit = True, template=Tirific_Template,Fits_Files=Fits_Files)



                    if Configuration['FINISHAFTER'] == 1:
                        Configuration['FINAL_COMMENT'] = 'You have chosen to end the fitting after preprocessing and sofia.'
                        cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run, Fits_Files =Fits_Files,debug=Configuration['DEBUG'])
                        continue
                    # Now set our variations to the original values but only if the galaxy is large enough
                    if Configuration['SIZE_IN_BEAMS'] > Configuration['MINIMUM_WARP_SIZE']:
                        Configuration['FIX_INCLINATION'][0] = Original_Configuration['FIX_INCLINATION'][0]
                        Configuration['FIX_SDIS'][0] = Original_Configuration['FIX_SDIS'][0]
                        Configuration['FIX_PA'][0] = Original_Configuration['FIX_PA'][0]
                        Configuration['FIX_Z0'][0] = Original_Configuration['FIX_Z0'][0]


                    #Then we want to setup for the next fit.
                    wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                                            cube_hdr,fit_stage='Extent_Convergence',debug=Configuration['DEBUG'])
                    if not os.path.isdir(Configuration['FITTING_DIR']+'/Extent_Convergence'):
                        os.mkdir(Configuration['FITTING_DIR']+'/Extent_Convergence')

                    while not Configuration['EC_ACCEPTED'] and Configuration['EC_LOOPS'] < allowed_loops:
                        Configuration['EC_LOOPS'] = Configuration['EC_LOOPS']+1
                        sf.print_log(f'''We are starting loop {Configuration['EC_LOOPS']} of trying to converge the extent.
        ''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                        if Configuration['DEBUG']:
                                sf.print_log(f'''Settings for the variations will be.
                    {'':8s} INCLINATION: Fixed = {Original_Configuration['FIX_INCLINATION']}
                    {'':8s} PA: Fixed = {Original_Configuration['FIX_PA']}
                    {'':8s} SDIS: Fixed = {Original_Configuration['FIX_SDIS']}
                    {'':8s} Z0: Fixed = {Original_Configuration['FIX_Z0']}
        ''',Configuration['OUTPUTLOG'],debug=Configuration['DEBUG'],screen =True)
                        if Configuration['SIZE_IN_BEAMS'] > Configuration['MINIMUM_WARP_SIZE']:
                            Configuration['FIX_INCLINATION'] = Original_Configuration['FIX_INCLINATION']
                            Configuration['FIX_SDIS'] = Original_Configuration['FIX_SDIS']
                            Configuration['FIX_PA'] = Original_Configuration['FIX_PA']
                            Configuration['FIX_Z0'] = Original_Configuration['FIX_Z0']
                        else:
                            Configuration['FIX_INCLINATION'] = True
                            Configuration['FIX_SDIS'] = True
                            Configuration['FIX_PA'] = True
                            Configuration['FIX_Z0'] = True
                            flatten_the_curve(Configuration,Tirific_Template)
                        current_run = runf.extent_converge(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,debug = Configuration['DEBUG'],allowed_loops = allowed_loops)


                    if Configuration['EC_ACCEPTED']:
                        sf.print_log(f'''The extent has converged and we make a smoothed version.
            ''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
                        current_run = runf.fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,stage = 'after_ec', fit_stage = 'Extent_Convergence',debug = Configuration['DEBUG'])
                        Configuration['FINAL_COMMENT'] = 'The galaxy has succesfully been fitted'
                        if Configuration['OPTIMIZED']:
                            runf.make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_stage = 'Extent_Convergence',debug=Configuration['DEBUG'])
                    else:
                        Configuration['FINAL_COMMENT'] = 'We could not converge on the extend of the galaxy'
                        Configuration['MAPS_OUTPUT'] = 5
                        cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
                        continue
                    Configuration['EC_END_TIME'] = datetime.now()
            except Exception as e:
                Configuration['FINAL_COMMENT'] = e
                if e.__class__.__name__ in stop_individual_errors:
                    Configuration['MAPS_OUTPUT'] = 5
                else:
                    Configuration['MAPS_OUTPUT'] = 'error'
                cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
                continue


            wf.basicinfo(Configuration,second_fit = True, template=Tirific_Template,Fits_Files=Fits_Files)
            cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run, Fits_Files =Fits_Files,debug = Configuration['DEBUG'])
            if input_parameters.installation_check:
                cf.installation_check(Configuration,debug=Configuration['DEBUG'])
    except Exception as e:
        Configuration['FINAL_COMMENT'] = e
        Configuration['MAPS_OUTPUT'] = 'error'
        cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run,debug=Configuration['DEBUG'])

main.__doc__ = '''
;+
; NAME:
;      FAT
; PURPOSE:
;      Fit Tilted Ring Models with Tirific in a fully automated manner
; CATEGORY:
;      Main for fitting galaxies. Tirific still requires interactive fitting this code attempts
;      to remedy that
;
; CALLING SEQUENCE:
;      FAT,support='supportdir',configuration_file='configfile'
;
; INPUTS:
;      -
; OPTIONAL INPUTS:
;      SUPPORT  = path to the directory where FAT's support
;      routines are located. The default location is ./Support/
;      CONFIGURATION_FILE = A configuration file for FAT. This file
;      should contain the locations of the galaxies to be fitted. See
;      readme for more detailed info.
;
; OPTIONAL INPUT KEYWORDS
;     /INSTALLATION_CHECK = Flag to run the Installation check.
; ---------------------------------------------------------------------------------
;     The following input keywords are only meant to be used by
;     developers. Except for the /debug flag they will not work for
;     the common user. If you want to know about these please contact Peter
;     Kamphuis.
; ---------------------------------------------------------------------------------
;     /DEBUG = Flag to print debugging information in several routines
;     /LVHIS_TEST = Flag to run the LVHIS Test.
;     /PAPER_TEST = Flag to run the paper artificial galaxies.
;     /RESOLUTION_TEST = Flag to run the additional resolution tests
; OPTIONAL KEYWORD OUTPUT:
;      -
;
; OUTPUTS:
;     See Readme or just run the code
;
; EXAMPLE:
;     python3 FAT.py --cf /home/your_computer/FAT_dir/FAT_INPUT.config'
'''
