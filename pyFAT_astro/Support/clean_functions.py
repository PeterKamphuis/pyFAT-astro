# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


import os,sys
import numpy as np
import traceback

from datetime import datetime
from astropy.io import fits
from make_moments.functions import moments
from pyFAT_astro.Support.modify_template import get_error,\
    set_fitting_parameters,get_ring_weights
from pyFAT_astro.Support.write_functions import make_overview_plot,plot_usage_stats,tirific
from pyFAT_astro.Support.log_functions import print_log,write_config
import pyFAT_astro.Support.support_functions as sf


class DummyLock():
    def __enter__(self):
        pass

    def __exit__(self, *args):
        pass

def check_legitimacy_osc(Configuration):
    print_log(f'''CHECK_LEGITIMACY_OSC: Start.
''',Configuration,case=['debug_start'])
    if Configuration['OUTPUT_QUANTITY'] == 'error':
        print_log(f'''CHECK_LEGITIMACY_OSC: An unspecified error is registered. The final message should reflect this.
''',Configuration)
        return
    elif Configuration['OUTPUT_QUANTITY'] == 5:
        print_log(f'''CHECK_LEGITIMACY_OSC: A FAT specific error is registered. The final message should reflect this.
''',Configuration)
        return
    else:
        fit_check=[True if 'fit_' in x.lower() else False for x in Configuration['FITTING_STAGES']]
        if any(fit_check):
            outfile = f"{Configuration['FITTING_DIR']}/{Configuration['USED_FITTING']}/{Configuration['USED_FITTING']}.def"
    inclination = sf.load_tirific(Configuration,outfile,Variables=['INCL'],array=True)

    if float(inclination[0]) < Configuration['UNRELIABLE_INCLINATION']:
        Configuration['ACCEPTED'] = False
        Configuration['FINAL_COMMENT'] = f"The final inclination is below {Configuration['UNRELIABLE_INCLINATION']}. FAT is not neccesarily reliable in this range."
        print_log(f'''CHECK_LEGITIMACY_OSC: The retrieved inclination {float(inclination[0])} is below {Configuration['UNRELIABLE_INCLINATION']} thus the fit is not accepted.
''',Configuration)
    if np.sum(Configuration['SIZE_IN_BEAMS']) < Configuration['UNRELIABLE_SIZE']:
        Configuration['ACCEPTED'] = False
        Configuration['FINAL_COMMENT'] = f"The final size is below {Configuration['UNRELIABLE_SIZE']}. FAT is not neccesarily reliable in this range."
        print_log(f'''CHECK_LEGITIMACY_OSC: The retrieved size {np.sum(Configuration['SIZE_IN_BEAMS'])} is below {Configuration['UNRELIABLE_SIZE']} thus the fit is not accepted.
''',Configuration)
    return

check_legitimacy_osc.__doc__ =f'''
 NAME:
    check_legitimacy_osc

 PURPOSE:
    Check the final output to see if it falls in a reliable rangesee if the output fall in the range where FAT is reliable.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    Updated final comment in the configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    !!!!!!! Until release these are set ridiculously low for testing purposes.
'''

def clean_before_sofia(Configuration):

    print_log(f'''CLEAN_BEFORE_SOFIA: Starting.
''',Configuration,case = ['debug_start'])
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits','_cat.txt']

    for file in files:
        try:
            os.remove(Configuration['BASE_NAME']+file)
        except FileNotFoundError:
            pass
        try:
            os.remove('Sofia_Output/'+Configuration['BASE_NAME']+file)
        except FileNotFoundError:
            pass
    try:
        os.remove('sofia_input.par')
    except FileNotFoundError:
        pass
    try:
        os.remove('Sofia_Output/sofia_input.par')
    except FileNotFoundError:
        pass

clean_before_sofia.__doc__ =f'''
 NAME:
    clean_before_sofia

 PURPOSE:
    Clean up the sofia out put files before running SoFiA

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    Previous sofia output from fat is removed

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def clean_after_sofia(Configuration):
    print_log(f'''CLEAN_AFTER_SOFIA: Starting clean
''',Configuration,case = ['debug_start'])
    files =['_mask.fits','_chan.fits','_cat.txt']

    for file in files:
        try:
            os.rename(Configuration['BASE_NAME']+file,'Sofia_Output/'+Configuration['BASE_NAME']+file )
        except FileNotFoundError:
            pass
    try:
        os.rename('sofia_input.par','Sofia_Output/sofia_input.par')
    except FileNotFoundError:
        pass

clean_after_sofia.__doc__ =f'''
 NAME:
    clean_after_sofia

 PURPOSE:
    Clean up the sofia output files by putting them in a dedicated directory.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    sofia output files are moved to a dedicated directory.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# cleanup dirty files before starting fitting
def cleanup(Configuration,Fits_Files):
    print_log(f'''CLEANUP: Starting clean
''',Configuration,case = ['debug_start'])
        #Move any existing output to the Log directory
    if os.path.exists(f"{Configuration['LOG_DIRECTORY']}ram_cpu.pdf"):
        if os.path.exists(f"{Configuration['LOG_DIRECTORY']}ram_cpu_prev.pdf"):
            os.remove(f"{Configuration['LOG_DIRECTORY']}ram_cpu_prev.pdf")
        os.rename( f"{Configuration['LOG_DIRECTORY']}ram_cpu.pdf",f"{Configuration['LOG_DIRECTORY']}ram_cpu_prev.pdf")
    #Move any existing Overview.png to the Log directory well
    if os.path.exists(f"{Configuration['FITTING_DIR']}Overview.png"):
        if os.path.exists(f"{Configuration['LOG_DIRECTORY']}Overview_prev.png"):
            os.remove(f"{Configuration['LOG_DIRECTORY']}Overview_prev.png")
            print_log(f'''CLEANUP: Removing an old Overview_prev.png from {Configuration['LOG_DIRECTORY']}
''',Configuration,case = ['debug_add'])
        os.rename( f"{Configuration['FITTING_DIR']}Overview.png",f"{Configuration['LOG_DIRECTORY']}Overview_prev.png")
        print_log(f'''CLEANUP: We moved an old Overview.png to {Configuration['LOG_DIRECTORY']}Overview_prev.png
''',Configuration,case = ['debug_add'])
    #clean the log directory of all files except those named Prev_ and not the Log as it is already moved if existing
    files_in_log = ['restart_One_Step_Convergence.txt','restart_Centre_Convergence.txt',f"restart_{Configuration['USED_FITTING']}.txt",\
                    'restart_Extent_Convergence.txt','Usage_Statistics.txt', 'clean_map_0.fits','clean_map_1.fits','clean_map.fits',\
                    'dep_map_0.fits','minimum_map_0.fits','rot_map_0.fits','dep_map.fits','minimum_map.fits','rot_map.fits',\
                    'dep_map_1.fits','minimum_map_1.fits','rot_map_1.fits','Convolved_Cube_FAT_opt.fits']

    #add all CFG files:
    for file in os.listdir(Configuration['LOG_DIRECTORY']):
        split_file = os.path.splitext(file)
        if split_file[-1] in ['.pkl', '.txt'] and 'CFG' in file:
            files_in_log.append(file)

    for file in files_in_log:
        try:
            os.remove(f"{Configuration['LOG_DIRECTORY']}{file}")
        except FileNotFoundError:
            pass
    #            !!!!!!!!!!!!!!! The Directories cleanup should be removed before release
    Directories = ['Extent_Convergence', 'Centre_Convergence','tmp_incl_check','One_Step_Convergence']
    for dir in Directories:
        try:
            for f in os.listdir(f"{Configuration['FITTING_DIR']}{dir}"):
                os.remove(os.path.join(f"{Configuration['FITTING_DIR']}{dir}", f))
        except FileNotFoundError:
            pass
        try:
            os.rmdir(f"{Configuration['FITTING_DIR']}{dir}")
        except FileNotFoundError:
            pass

    files_in_main = ['dep_map_0.0.fits','minimum_map_0.0.fits','clean_map_0.0.fits','rot_map_0.0.fits','tmp_incl_check_In.def']

    for file in files_in_main:
        try:
            os.remove(f"{Configuration['FITTING_DIR']}{file}")
        except FileNotFoundError:
            pass
    directories = []
    files = [Configuration['BASE_NAME']+'-Basic_Info.txt',Fits_Files['OPTIMIZED_CUBE']]
    if Configuration['USED_FITTING']:
        directories.append('Finalmodel')
        directories.append(Configuration['USED_FITTING'])
        files.append(f'{Configuration["USED_FITTING"]}_In.def')

    if 'create_fat_cube' in Configuration['FITTING_STAGES']:
        files.append(Fits_Files['FITTING_CUBE'])

    if 'run_sofia' in Configuration['FITTING_STAGES'] or 'external_sofia' in Configuration['FITTING_STAGES']:
        dir =f'{Configuration["FITTING_DIR"]}Sofia_Output/'
        file_ext=['_mask.fits','_mom0.fits','_mom1.fits','_mom2.fits','_chan.fits','_cat.txt','_sofia_xv.fits']
        print_log(f'''CLEANUP: We are cleaning the following files in the directory {dir}:
{"":8s}CLEANUP: sofia_input.par,{','.join([f'{Configuration["SOFIA_BASENAME"]}{x}' for x in file_ext])}
''',Configuration,case = ['verbose'])
        for extension in file_ext:
            try:
                os.remove(f'{dir}{Configuration["BASE_NAME"]}{extension}')
            except FileNotFoundError:
                pass
        try:
            os.remove(f'{dir}sofia_input.par')
        except FileNotFoundError:
            pass
    if 'tirshaker' in Configuration['FITTING_STAGES']:
        directories.append('Error_Shaker')
    # Existing_Sofia



    if directories:
        print_log(f'''CLEANUP: We are cleaning the following directories:
{"":8s}CLEANUP: {','.join(directories)}
{"":8s}CLEANUP: and the following files:
{"":8s}CLEANUP: {','.join(files)}
''',Configuration, case=['verbose'])


        ext=['.fits','_Prev.fits','.log','.ps','.def']
        moments = ['mom0','mom1','mom2', 'xv']
        #then specific files in the working directory
        #os.chdir(Configuration['FITTING_DIR'])
        #for configuration purposes we remove the old dirs


        for dir in directories:
            if os.path.isdir(f'{Configuration["FITTING_DIR"]}{dir}'):

                for fe in ext:
                    if dir == 'Finalmodel' and fe in ['.fits','.def']:
                        try:
                            os.unlink(f'{Configuration["FITTING_DIR"]}{dir}/{dir}{fe}')
                        except FileNotFoundError:
                            pass

                    elif dir == Configuration['USED_FITTING'] and fe in ['.def']:
                        target = sf.get_system_string(f"{Configuration['FITTING_DIR']}{dir}/{dir}*{fe}")
                        os.system(f'rm -f {target}')
                    else:
                        try:
                            os.remove(f'{Configuration["FITTING_DIR"]}{dir}/{dir}{fe}')
                        except FileNotFoundError:
                            pass

                for mom in moments:
                    try:
                        os.remove(f'{Configuration["FITTING_DIR"]}{dir}/{dir}_{mom}.fits')
                    except FileNotFoundError:
                        pass
                if dir == 'Finalmodel':
                    try:
                        os.remove(f'{Configuration["FITTING_DIR"]}{dir}/{Configuration["BASE_NAME"]}_final_xv.fits')
                    except FileNotFoundError:
                        pass


    for file in files:
        try:
            os.remove(f'{Configuration["FITTING_DIR"]}{file}')
        except FileNotFoundError:
            pass

cleanup.__doc__ =f'''
 NAME:
    cleanup

 PURPOSE:
    Clean up any existing files before fitting start, it will only remove specific
    files not directories or non-FAT files

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def cleanup_final(Configuration,Fits_Files):
    print_log(f'''CLEANUP_FINAL: Starting the final cleanup of the directory.
{'':8s} The request is {Configuration['OUTPUT_QUANTITY']}
''',Configuration,case = ['debug_start','verbose'] )

    if Configuration['USED_FITTING'] == 'Fit_Tirific_OSC':
        clean_files = [Fits_Files['OPTIMIZED_CUBE'],f"{Configuration['USED_FITTING']}_In.def",\
                        "clean_map_0.fits","dep_map_0.fits","minimum_map_0.fits","rot_map_0.fits",\
                        "clean_map_1.fits","dep_map_1.fits","minimum_map_1.fits","rot_map_1.fits",
                        "tmp_incl_check_In.def"\
                        ]
    else:
        clean_files = []
    for file in clean_files:
    # Not remove anything but cleanup all
        try:
            if Configuration['OUTPUT_QUANTITY'] >= 5 or Configuration['OUTPUT_QUANTITY'] == 0:
                os.rename(f"{Configuration['FITTING_DIR']}{file}",f"{Configuration['LOG_DIRECTORY']}{file}")
            else:
                os.remove(f"{Configuration['FITTING_DIR']}{file}")
        except FileNotFoundError:
            pass

    fit_directories = [Configuration['USED_FITTING']]
    delete_ext = ['.log']
    if Configuration['OUTPUT_QUANTITY'] == 4:
        delete_ext.append('.fits')

    for dir in fit_directories:
        try:
            files_in_dir = os.listdir(f"{Configuration['FITTING_DIR']}{dir}")
        except FileNotFoundError:
            files_in_dir = []
            pass
        for file in files_in_dir:
            name,extension = os.path.splitext(f"{Configuration['FITTING_DIR']}{dir}/{file}")
            if extension in delete_ext:
                try:
                    os.remove(f"{Configuration['FITTING_DIR']}{dir}/{file}")
                except FileNotFoundError:
                    pass
                except IsADirectoryError:
                    pass
            if (Configuration['OUTPUT_QUANTITY'] == 2 and extension != ".fits") \
                or (5 > Configuration['OUTPUT_QUANTITY'] >= 3):
                if file != f"{Configuration['USED_FITTING']}.def" and file != f"{Configuration['USED_FITTING']}.fits":
                    try:
                        if not (Configuration['DEBUG'] and extension == '.def'):
                            os.remove(f"{Configuration['FITTING_DIR']}{dir}/{file}")
                    except FileNotFoundError:
                        pass
    stage_dirs = []
    if os.path.isdir(f"{Configuration['FITTING_DIR']}tmp_incl_check"):
        stage_dirs.append('tmp_incl_check')
    if 'tirshaker' in Configuration['FITTING_STAGES']:
        stage_dirs.append('Error_Shaker')

    for dirs in stage_dirs:
        if 5 > Configuration['OUTPUT_QUANTITY'] >= 1:
            files_in_dir = os.listdir(f"{Configuration['FITTING_DIR']}{dirs}")
            for file in files_in_dir:
                try:
                    os.remove(f"{Configuration['FITTING_DIR']}{dirs}/{file}")
                except FileNotFoundError:
                    pass
            os.rmdir(f"{Configuration['FITTING_DIR']}{dirs}")
        else:
            # else move this directory to the LOG
            target = sf.get_system_string(f"{Configuration['LOG_DIRECTORY']}{dirs}")
            if  os.path.isdir(f"{Configuration['LOG_DIRECTORY']}{dirs}"):
                os.system(f"rm -Rf {target}")
            source = sf.get_system_string(f"{Configuration['FITTING_DIR']}{dirs}")
            os.system(f"mv {source} {target}")


cleanup_final.__doc__ =f'''
 NAME:
    cleanup_final

 PURPOSE:
    Clean up any files not requested by user after the fitting, it will only remove specific
    files not directories or non-FAT files

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''





def installation_check(Configuration):
    print_log(f'''INSTALLATION_CHECK: Starting to compare the output to what is expected.
''',Configuration, case = ['debug_start', 'main'])

    Model = sf.tirific_template(filename = 'Installation_Check')
    Variables_to_Compare = ['VROT','INCL','PA','SBR','SDIS','Z0','XPOS','YPOS','VSYS']
    Model_values = sf.load_tirific(Configuration,Model,\
        Variables = Variables_to_Compare,array=True)
    #Then the fitted file
    Fitted_values =sf.load_tirific(Configuration,\
        f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def",\
        Variables = Variables_to_Compare,array = True)
    succes = False
    diff = np.abs(Model_values-Fitted_values)

    print_log(f'''INSTALLATION_CHECK: the found differences
{'':8s}{diff}
''',Configuration, case = ['verbose'])
    too_much = np.array(np.where(diff > 1e-3),dtype=bool)

    print_log(f'''INSTALLATION_CHECK: at the locations
{'':8s}{too_much}{np.where(diff > 1e-3)}
{'':8s}{too_much.size}
''',Configuration, case = ['verbose'])

    if too_much.size == 0.:
        succes = True


    if succes:
        print_log(f'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! All parameters are fitted within the expected variance. !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        We think pyFAT is installed succesfully          !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration, case = ['main','screen'])
    else:
        print_log(f'''
!!!!---------------------------------------------!!!!!
!!!! FAT ran through the fitting process but the !!!!!
!!!! fitted values differ too much from their    !!!!!
!!!! expectations. Please update SoFiA and other !!!!!
!!!! dependencies. If you have done so and this !!!!!
!!!! message remains, the check is likely out of !!!!!
!!!! date. If you are unable to resolve the issue!!!!!
!!!! please file a bug report at:                !!!!!
!!!!                                             !!!!!
!!!! https://github.com/PeterKamphuis/pyFAT/issues !!!!!
!!!!                                             !!!!!
!!!! Please add the Log.txt file in the directory!!!!!
!!!!{Configuration['LOG_DIRECTORY']}!!!!!
!!!! and the Finalmodel.def.                     !!!!!
!!!!---------------------------------------------!!!!!
''',Configuration,case = ['main','screen'])

installation_check.__doc__ =f'''
 NAME:
    installation_check

 PURPOSE:
    Check the installation check run against the template.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    Message with succes or not to the screen and log

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def finish_galaxy(Configuration,current_run = 'Not initialized',\
        timing_lock= DummyLock(), catalogue_lock = DummyLock(),
        Fits_Files= {'ORIGINAL_CUBE': "Unset.fits"},exiting = None):
    Configuration['FULL_TIME'][1] = datetime.now()
    print_log(f'''FINISH_GALAXY: These fits files are used:
{'':8s} {Fits_Files}
''',Configuration,case = ['debug_start','verbose'])
    if Configuration['DEBUG']:
        write_config(f'{Configuration["LOG_DIRECTORY"]}CFG_At_Finish.txt',Configuration)

    #make sure we are not leaving stuff
    sf.finish_current_run(Configuration,current_run)
    # And place the correct velocity frame
    if Configuration['HDR_VELOCITY'] != 'VELO':
        print_log(f'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}Your velocity frame is not suitable for Tirific it is CTYPE3 = {Configuration['HDR_VELOCITY']} and should be VELO, FELO or FREQ 
{"":8s}If you want to continue to fit with tirific please use the cube marked _tirific , we have kept the correct conversion in the FAT cube
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration, case = ['main','screen'])
        



    # We need to check if the final output is legit
    if Configuration['USED_FITTING'] == 'Fit_Tirific_OSC':
        check_legitimacy_osc(Configuration)


    if Configuration['OUTPUT_QUANTITY'] == 'error':
        error_message = '''
            Your code has crashed for some reason. If this message completely baffles you then please submit the trace back as a bug report to: \n
            https://github.com/PeterKamphuis/pyFAT/issues \n
            If the error occured while fitting a galaxy, please attach your fitting log as well
'''
        log_statement = f'''------------When filing a bug report please copy all output  below this line------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT did not run the full fitting routines for catalog entry {Configuration['ID']}.
{"":8s}Which is the galaxy in directory {Configuration['FITTING_DIR']}.
{"":8s}Please check this log and output_catalogue carefully for what went wrong.
{"":8s}The detected exit reason is: "{Configuration['FINAL_COMMENT']}".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
        print_log(log_statement,Configuration, case = ['main','screen'])
        print_log(error_message,Configuration, case = ['main','screen'])

        if exiting != None:
            with open(Configuration['OUTPUTLOG'],'a') as log_file:
                traceback.print_exception(type(exiting),exiting,exiting.__traceback__,file=log_file)
            traceback.print_exception(type(exiting),exiting,exiting.__traceback__)
        
        if Configuration['MULTIPROCESSING']:
            Configuration['ACCEPTED'] = False
            Configuration['FINAL_COMMENT'] = f"The code crashed while fitting this galaxy please check it's log."
            Configuration['OUTPUT_QUANTITY'] = 5
        else:
            if exiting != None:
                sys.exit(1)
            else:
                Configuration['ACCEPTED'] = False
                Configuration['FINAL_COMMENT'] = f"The code crashed while fitting this galaxy please check it's log."
                Configuration['OUTPUT_QUANTITY'] = 5
    elif Configuration['OUTPUT_QUANTITY'] == 5:
        print_log(f'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT could not find an acceptable model for the galaxy in directory {Configuration['FITTING_DIR']}.
{"":8s}Please check the log in {Configuration['LOG_DIRECTORY']} and the output_catalogue carefully for more information.
{"":8s}The detected exit reason is: "{Configuration['FINAL_COMMENT']}".
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration, case = ['main','screen'])

    elif Configuration['OUTPUT_QUANTITY'] < 4:
        print_log( f'''Producing final output in {Configuration['FITTING_DIR']}.
''',Configuration)
        # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
        if any([True if 'fit_' in x else False for x in Configuration['FITTING_STAGES']]):
            if not 'Fit_Make_Your_Own' in  Configuration['USED_FITTING']:
                sf.create_directory('Finalmodel',Configuration['FITTING_DIR'])
                if 'tirshaker' not in Configuration['FITTING_STAGES']:
                    transfer_errors(Configuration,fit_type=Configuration['USED_FITTING'])
                    #Also transfer the last fitting settings
                    transfer_fitting_parameters(Configuration)    
                linkname = f"../{Configuration['USED_FITTING']}/{Configuration['USED_FITTING']}"
                os.symlink(f"{linkname}.fits",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
                os.symlink(f"{linkname}.def",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.def")

                # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
                if Fits_Files and os.path.exists(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits"):
                    # we do not mask the model as it is also important where the model extends and there is no data.  
                    messages = moments(filename =  f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits",\
                            overwrite = True, cube_velocity_unit='m/s',map_velocity_unit= 'km/s',\
                            debug = Configuration['DEBUG'], log=True, level = 0.25*Configuration['NOISE'],\
                            output_directory = f"{Configuration['FITTING_DIR']}/Finalmodel/",\
                            output_name = 'Finalmodel')
                    print_log(messages,Configuration,case=['verbose'])
                    #make_moments(Configuration,Fits_Files,fit_type = 'Generic_Final',vel_unit = 'm/s')
                    make_overview_plot(Configuration,Fits_Files)
    

    print_log(f'''Finished the fitting in {Configuration['FITTING_DIR']}.
''',Configuration, case = ['main','screen'])
    # Need to organize the fitting output orderly
    # Need to write date and Time to timing log
    if Configuration['TIMING']:
        plot_usage_stats(Configuration)
        with timing_lock:
            with open(Configuration['MAIN_DIRECTORY']+'/Timing_Result.txt','a') as timing_result:
               timing_result.write(f'''The galaxy in directory {Configuration['FITTING_DIR']} has the following timing results.
Initialization from {Configuration['INITIALIZATION_TIME'][0]} to {Configuration['INITIALIZATION_TIME'][1]}
Preparation time  from {Configuration['PREPARATION_TIME'][0]} to {Configuration['PREPARATION_TIME'][1]}
Time to run Sofia from {Configuration['SOFIA_TIME'][0]} to {Configuration['SOFIA_TIME'][1]}
Time to get Initial estimates {Configuration['INITIAL_GUESSES_TIME'][0]} to {Configuration['INITIAL_GUESSES_TIME'][1]}
Time for the fitting {Configuration['FIT_TIME'][0]} to {Configuration['FIT_TIME'][1]}
Time for TirShaking {Configuration['TIRSHAKER_TIME'][0]} to {Configuration['TIRSHAKER_TIME'][1]}
Full run time from {Configuration['FULL_TIME'][0]} to {Configuration['FULL_TIME'][1]} \n''')


            #    timing_result.write(f'''The galaxy in directory {Configuration['FITTING_DIR']} started at {Configuration['START_TIME']}.
#Finished preparations at {Configuration['PREP_END_TIME']}.
#Converged to a galaxy size at {Configuration['END_TIME']}.
#It finished the whole process at {datetime.now()}.
#''')

        print_log(f'''Finished timing statistics for the galaxy in {Configuration['FITTING_DIR']}.
''',Configuration)

    cleanup_final(Configuration,Fits_Files)
    # Need to write to results catalog
    catalogue_line = f"{Configuration['FITTING_DIR'].split('/')[-2]:{Configuration['MAXIMUM_DIRECTORY_LENGTH']}s} {str(Configuration['ACCEPTED']):>6s} {Configuration['FINAL_COMMENT']} \n"
    if Configuration['OUTPUT_CATALOGUE']:
        with catalogue_lock:
            with open(Configuration['OUTPUT_CATALOGUE'],'a') as output_catalogue:
                output_catalogue.write(catalogue_line)
    return catalogue_line

finish_galaxy.__doc__ =f'''
 NAME:
    finish_galaxy

 PURPOSE:
    Write the final logs, Produce an Overview plot if required and make sure the directory is clean and organized.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration


 OPTIONAL INPUTS:
    current_run = subproccess structure of the currently running tirific
    Fits_Files = Standard FAT dictionary with filenames

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This function should always be followed by a continue statement
'''

def mod_vrot_incl_err(profile,errors,incl,incl_errors):
    '''Error up and down have to be exactly the same'''
    vobs_up = np.sin(np.radians(incl+incl_errors))*\
                (profile+errors)
    errors = vobs_up/np.sin(np.radians(incl))-profile
    return errors

def transfer_errors(Configuration,fit_type='Undefined'):
    print_log(f'''TRANSFER_ERROR: Start.
''',Configuration, case = ['debug_start'])
    # Load the final file
    Tirific_Template = sf.tirific_template(filename = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def")
    variables = ['INCL','PA','VROT','SDIS','SBR','VSYS','XPOS','YPOS','Z0']
    weights = get_ring_weights(Configuration,Tirific_Template)
    incl_error = None
    incl = None
    for parameter in variables:
        print_log(f'''TRANSFER_ERROR: Creating errors for {parameter}.
''',Configuration, case = ['debug_add'])
        reg_profile = sf.load_tirific(Configuration,Tirific_Template,\
            [parameter,f"{parameter}_2"],array=True)
        profile = sf.load_tirific(Configuration,\
            f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Iteration_{Configuration['ITERATIONS']}.def",\
            Variables=[parameter,f"{parameter}_2"],array=True)
        
        # As we can cut rings away during the smoothing we have to ensure the profiles have the same length
        if len(profile[0,:]) != len(reg_profile[0,:]):
            tmp = [[],[]]    
            for i in [0,1]:
                tmp[i] = profile[i,:len(reg_profile[i,:]-1)]
            profile=np.array(tmp,dtype=float)
            #it is possible that the last iteration of the fitted smooth check
            #Crashes and is rerun with
        if parameter == 'VROT':
            apply_max= False
        else:
            apply_max = True
        errors = get_error(Configuration,profile,reg_profile,parameter,\
            min_error=Configuration['MIN_ERROR'][parameter],\
            apply_max_error = apply_max,weights = weights)
        if parameter == 'VROT':
            errors = mod_vrot_incl_err(reg_profile,errors,incl,incl_error)   
        if parameter == 'INCL':
            incl_error=errors
            incl = reg_profile

        format = sf.set_format(parameter)
        Tirific_Template.insert(parameter,f"# {parameter}_ERR",f"{' '.join([f'{x:{format}}' for x in errors[0]])}")
        Tirific_Template.insert(f"{parameter}_2",f"# {parameter}_2_ERR",f"{' '.join([f'{x:{format}}' for x in errors[1]])}")


    # write back to the File
    tirific(Configuration,Tirific_Template, name = f"{fit_type}/{fit_type}.def")

transfer_errors.__doc__ =f'''
 NAME:
    transfer_errors

 PURPOSE:
    Transfer the errors from the input file to the output as tirific removes the commented lines.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    fit_type = 'Undefined'

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified
'''

def transfer_fitting_parameters(Configuration):
    fit_type = Configuration['USED_FITTING']
    # first read the Final file
    print_log(f'''TRANSFER_FITTING_PARAMETERS: Start.
''',Configuration, case = ['debug_start'])
    # Load the final file
    Tirific_Template = sf.tirific_template(filename = \
        f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def")
    
    set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_os')
   
    # write back to the File
    tirific(Configuration,Tirific_Template, name = f"{fit_type}/{fit_type}.def")

