# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


import os,signal,sys
import numpy as np
import traceback
from datetime import datetime
from pyFAT.Support.support_functions import print_log,finish_current_run,set_format
from pyFAT.Support.fits_functions import make_moments
from pyFAT.Support.write_functions import make_overview_plot,plot_usage_stats,tirific
from pyFAT.Support.read_functions import tirific_template,load_tirific,load_template

def check_legitimacy(Configuration,debug=False):
    if debug:
        print_log(f'''CHECK_LEGITIMACY: Start.
''',Configuration['OUTPUTLOG'],debug =True)
    if Configuration['MAPS_OUTPUT'] == 'error':
        print_log(f'''CHECK_LEGITIMACY: An unspecified error is registered. The final message should reflect this.
''',Configuration['OUTPUTLOG'])
        return
    elif Configuration['MAPS_OUTPUT'] == 5:
        print_log(f'''CHECK_LEGITIMACY: A FAT specific error is registered. The final message should reflect this.
''',Configuration['OUTPUTLOG'])
        return
    else:
        if Configuration['FINISHAFTER'] > 0:
            outfile = f"{Configuration['FITTING_DIR']}/One_Step_Convergence/One_Step_Convergence.def"
    inclination = load_tirific(Configuration,outfile,Variables=['INCL'],debug=debug)[0]
    low_incl_limit = 10.
    low_beam_limit = 0.
    if float(inclination[0]) < low_incl_limit or 2.*Configuration['SIZE_IN_BEAMS'] < low_beam_limit:
        Configuration['ACCEPTED'] = False
        if float(inclination[0]) < low_incl_limit:
            Configuration['FINAL_COMMENT'] = f'The final inclination is below {low_incl_limit}. FAT is not neccesarily reliable in this range.'
        elif 2.*Configuration['SIZE_IN_BEAMS'] < low_beam_limit:
            Configuration['FINAL_COMMENT'] = f'The final size is below {low_beam_limit}. FAT is not neccesarily reliable in this range.'
    return

check_legitimacy.__doc__ =f'''
 NAME:
    check_legitimacy

 PURPOSE:
    Check the final output to see if it falls in a reliable rangesee if the output fall in the range where FAT is reliable.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Updated final comment in the configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    !!!!!!! Until release these are set ridiculously low for testing purposes.
'''

def clean_before_sofia(Configuration, debug = False):
    if debug:
        print_log(f'''CLEAN_BEFORE_SOFIA: Starting.
''',Configuration['OUTPUTLOG'],debug = True)
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
    debug = False

 OUTPUTS:
    Previous sofia output from fat is removed

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def clean_after_sofia(Configuration, debug = False):
    files =['_mask.fits','_chan.fits','_cat.txt']

    for file in files:
        try:
            os.rename(Configuration['BASE_NAME']+file,'Sofia_Output/'+Configuration['BASE_NAME']+file )
        except:
            pass
    try:
        os.rename('sofia_input.par','Sofia_Output/sofia_input.par')
    except:
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
    debug = False

 OUTPUTS:
    sofia output files are moved to a dedicated directory.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# cleanup dirty files before starting fitting
def cleanup(Configuration,Fits_Files, debug = False):
    #clean the log directory of all files except those named Prev_ and not the Log as it is already moved if existing
    files_in_log = ['restart_One_Step_Convergence.txt''restart_Centre_Convergence.txt',\
                    'restart_Extent_Convergence.txt','Usage_Statistics.txt', 'clean_map_0.fits','clean_map_1.fits','clean_map.fits',\
                    'dep_map_0.fits','minimum_map_0.fits','rot_map_0.fits','dep_map.fits','minimum_map.fits','rot_map.fits',\
                    'dep_map_1.fits','minimum_map_1.fits','rot_map_1.fits','Convolved_Cube_FAT_opt.fits']
    #            !!!!!!!!!!!!!!! The Directories cleanup should be removed before release
    Directories = ['Extent_Convergence', 'Centre_Convergence','tmp_incl_check']
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
    for file in files_in_log:
        try:
            os.remove(f"{Configuration['FITTING_DIR']}Logs/{file}")
        except FileNotFoundError:
            pass
    files_in_main = ['dep_map_0.0.fits','minimum_map_0.0.fits','clean_map_0.0.fits','rot_map_0.0.fits','tmp_incl_check_In.def']

    for file in files_in_main:
        try:
            os.remove(f"{Configuration['FITTING_DIR']}{file}")
        except FileNotFoundError:
            pass
    directories = ['Finalmodel','Sofia_Output','Def_Files']
    files = [Configuration['BASE_NAME']+'-Basic_Info.txt']
    directories.append('One_Step_Convergence')
    files.append('One_Step_Convergence_In.def')
    if Configuration['START_POINT'] == 1:
        files.append(Fits_Files['FITTING_CUBE'])
        files.append(Fits_Files['OPTIMIZED_CUBE'])
    elif Configuration['START_POINT'] == 2:
        pass
    elif Configuration['START_POINT'] == 4 and Configuration['FINISHAFTER'] > 1:
        files = ['Centre_Convergence_In.def',]
    else:
        directories= None
        files =  None
    if directories:
        print_log(f'''CLEANUP: We are cleaning the following directories:
{"":8s}CLEANUP: {','.join(directories)}
{"":8s}CLEANUP: and the following files:
{"":8s}CLEANUP: {','.join(files)}
''',Configuration['OUTPUTLOG'], screen =True,debug=debug)


        ext=['.fits','.log','.ps','.def']
        moments = ['mom0','mom1','mom2', 'xv']
        #then specific files in the working directory
        os.chdir(Configuration['FITTING_DIR'])
        #for configuration purposes we remove the old dirs
        if debug:
            try:
                os.system(f'rm -f Finalmodel/Convolved_Cube_FAT_final_xv.fits')
            except:
                pass

        for dir in directories:
            if os.path.isdir(dir):
                if dir == 'Sofia_Output':
                    os.system(f'rm -f {dir}/{Configuration["BASE_NAME"]}*_binmask* {dir}/{Configuration["BASE_NAME"]}*.ascii {dir}/{Configuration["BASE_NAME"]}*_sofia_xv.fits')
                else:
                    for fe in ext:
                        if dir == 'Finalmodel' and fe in ['.fits','.def']:
                            try:
                                os.unlink(f'{dir}/{dir}{fe}')
                            except FileNotFoundError:

                                pass
                        else:
                            try:
                                os.remove(f'{dir}/{dir}{fe}')
                            except FileNotFoundError:
                                pass
                            os.system(f'rm -f {dir}/{dir}_*{fe}')
                    for mom in moments:
                        try:
                            os.remove(f'{dir}/{dir}_{mom}.fits')
                        except FileNotFoundError:
                            pass


        for file in files:
            try:
                os.remove(file)
            except FileNotFoundError:
                pass
        #Change back to original dir
        os.chdir(Configuration['START_DIR'])

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
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def cleanup_final(Configuration,Fits_Files, debug =False):
    if debug:
         print_log(f'''Starting the final cleanup of the directory.
''',Configuration['OUTPUTLOG'],debug = True)
    clean_files = [Fits_Files['OPTIMIZED_CUBE'],'Centre_Convergence_In.def', 'Extent_Convergence_In.def',\
                    'clean_map_0.fits','dep_map_0.fits','minimum_map_0.fits','rot_map_0.fits',\
                    'clean_map_1.fits','dep_map_1.fits','minimum_map_1.fits','rot_map_1.fits',\
                    'One_Step_Convergence_In.def']
    dest_dir = ['Logs','Center_Convergence', 'Extent_Convergence', 'Logs', 'Logs', 'Logs', 'Logs', 'Logs', 'Logs', 'Logs', 'Logs','One_Step_Convergence']
    for file,dir in zip(clean_files,dest_dir):
    # Not remove anything but cleanup all
        try:
            if Configuration['MAPS_OUTPUT'] >= 5 or Configuration['MAPS_OUTPUT'] == 0:
                os.rename(f"{Configuration['FITTING_DIR']}/{file}",f"{Configuration['FITTING_DIR']}{dir}/{file}")
            else:
                os.remove(f"{Configuration['FITTING_DIR']}/{file}")
        except FileNotFoundError:
            pass

    #fit_directories = ['Centre_Convergence', 'Extent_Convergence']
    fit_directories = ['One_Step_Convergence']
    delete_ext = ['.log']
    if Configuration['MAPS_OUTPUT'] == 4:
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
            if (Configuration['MAPS_OUTPUT'] == 2 and extension != '.fits') or (5 >= Configuration['MAPS_OUTPUT'] >= 3) :
                if len(name.split('_')) > 3 and (file != 'One_Step_Convergence.def' and file != 'One_Step_Convergence.fits') :
                    try:
                        os.remove(f"{Configuration['FITTING_DIR']}{dir}/{file}")
                    except FileNotFoundError:
                        pass

    if 5 >= Configuration['MAPS_OUTPUT'] >= 1:
        if os.path.isdir(f"{Configuration['FITTING_DIR']}tmp_incl_check"):
            files_in_dir = os.listdir(f"{Configuration['FITTING_DIR']}tmp_incl_check")
            for file in files_in_dir:
                try:
                    os.remove(f"{Configuration['FITTING_DIR']}tmp_incl_check/{file}")
                except FileNotFoundError:
                    pass
            os.rmdir(f"{Configuration['FITTING_DIR']}tmp_incl_check")

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
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def copy_homemade_sofia(Configuration,Fits_Files,debug=False):
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits','_cat.txt']
    for file in files:
        try:
            os.system(f'''cp {Configuration['SOFIA_BASENAME']+file} {Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']+file}''')
        except:
            pass

copy_homemade_sofia.__doc__ =f'''
 NAME:
    copy_homemade_sofia

 PURPOSE:
    Copy user provided Sofia files to the FAT specified directory such that the original are a) kept in place, b) Never modified.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def installation_check(Configuration, debug =False):
    if debug:
         print_log(f'''INSTALLATION_CHECK: Starting to compare the output to what is expected.
''',Configuration['OUTPUTLOG'],debug = True)

    Model = tirific_template(filename = 'Installation_Check', debug= debug)
    Variables_to_Compare = ['VROT','INCL','PA','SBR','SDIS','Z0','XPOS','YPOS','VSYS']
    Model_values = load_template(Configuration,Model,Variables = Variables_to_Compare,unpack = False,debug=debug)
    #Then the fitted file
    Fitted_values =load_tirific(Configuration,f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def",Variables = Variables_to_Compare,unpack = False,debug=debug)
    succes = False

    try:
        diff = np.array(Model_values,dtype=float)-np.array(Fitted_values,dtype=float)
        if debug:
            print_log(f'''INSTALLATION_CHECK: the found differences
{'':8s}{diff}
''',Configuration['OUTPUTLOG'])
        too_much = np.array(np.where(diff > 1e-3),dtype=float)
        if debug:
            print_log(f'''INSTALLATION_CHECK: at the locations
{'':8s}{too_much}
{'':8s}{too_much.size}
''',Configuration['OUTPUTLOG'])

        if too_much.size == 0.:
            succes = True
    except:
        pass

    if succes:
        print_log(f'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!! All parameters are fitted within the expected variance. !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!        We think pyFAT is installed succesfully          !!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration['OUTPUTLOG'],screen =True)
    else:
        print_log(f'''
!!!!---------------------------------------------!!!!!
!!!! FAT ran through the fitting process but the !!!!!
!!!! fitted values differ too much from their    !!!!!
!!!! expectations. Please update SoFiA and other !!!!!
!!!! dependencies. If you are unable to resolve  !!!!!
!!!! resolve the issue please file a bug report  !!!!!
!!!! at:                                         !!!!!
!!!!                                             !!!!!
!!!! https://github.com/PeterKamphuis/pyFAT/issues !!!!!
!!!!                                             !!!!!
!!!! Please add the Log.txt file in the directory!!!!!
!!!! Installation_Check and the Finalmodel.def   !!!!!
!!!! in the Installation_Check/Finalmodel/       !!!!!
!!!! directory to your report.                   !!!!!
!!!!---------------------------------------------!!!!!
''',Configuration['OUTPUTLOG'],screen =True)

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
    debug = False

 OUTPUTS:
    Message with succes or not to the screen and log

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def finish_galaxy(Configuration,maximum_directory_length,current_run = 'Not initialized', Fits_Files= None, debug = False):
    Configuration['END_TIME'] = datetime.now()
    #make sure we are not leaving stuff
    finish_current_run(Configuration,current_run,debug=debug)
    # We need to check if the final output is legit
    check_legitimacy(Configuration,debug=debug)

    # Need to write to results catalog
    if Configuration['OUTPUTCATALOGUE']:
        with open(Configuration['OUTPUTCATALOGUE'],'a') as output_catalogue:
            output_catalogue.write(f"{Configuration['FITTING_DIR'].split('/')[-2]:{maximum_directory_length}s} {str(Configuration['ACCEPTED']):>6s} {Configuration['FINAL_COMMENT']} \n")

    if Configuration['MAPS_OUTPUT'] == 'error':
        error_message = '''
            Your code has crashed for some reason. If this message completely baffles you then please submit the trace back as a bug report to: \n
            https://github.com/PeterKamphuis/pyFAT/issues \n
            If the error occured while fitting a galaxy, please attach your fitting log as well
'''
        log_statement = f'''------------When filing a bug report please copy all output  below this line------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT did not run the full fitting routines for catalog nr {Configuration['ID_NR']}.
{"":8s}Which is the galaxy in directory {Configuration['FITTING_DIR']}.
{"":8s}Please check this log and output_catalogue carefully for what went wrong.
{"":8s}The detected exit reason is {Configuration['FINAL_COMMENT']}.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
        print_log(error_message,Configuration['OUTPUTLOG'], screen = True)
        with open(Configuration['OUTPUTLOG'],'a') as log_file:
                    traceback.print_exc(file=log_file)
        traceback.print_exc()
        sys.exit(1)
    elif Configuration['MAPS_OUTPUT'] == 5:
        log_statement = f'''
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT did not run the full fitting routines for the galaxy in directory {Configuration['FITTING_DIR']}.
{"":8s}Please check this log and output_catalogue carefully for what went wrong.
{"":8s}The detected exit reason is {Configuration['FINAL_COMMENT']}.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
        with open(Configuration['OUTPUTLOG'],'a') as log_file:
                    traceback.print_exc(file=log_file)
        traceback.print_exc()
        print_log(log_statement,Configuration['OUTPUTLOG'])
    elif Configuration['MAPS_OUTPUT'] < 4:
        log_statement = f'''Producing final output in {Configuration['FITTING_DIR']}.
'''
        #
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
        # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
        if Configuration['FINISHAFTER'] > 0:
            if not os.path.isdir(Configuration['FITTING_DIR']+'/Finalmodel'):
                os.mkdir(Configuration['FITTING_DIR']+'/Finalmodel')


            transfer_errors(Configuration,fit_type='One_Step_Convergence')
            linkname = f"../One_Step_Convergence/One_Step_Convergence"
            os.symlink(f"{linkname}.fits",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
            os.symlink(f"{linkname}.def",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.def")

            # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
            if Fits_Files:
                if Configuration['FINISHAFTER'] == 1 and Configuration['START_POINT'] == 4:
                    print_log("Moments should exist",Configuration['OUTPUTLOG'],screen =True)
                else:
                    make_moments(Configuration,Fits_Files,fit_type = 'Generic_Final',vel_unit = 'm/s',debug=debug)
                make_overview_plot(Configuration,Fits_Files,debug=debug)

    log_statement = f'''Finished final output in {Configuration['FITTING_DIR']}.
'''
    print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
    # Need to organize the fitting output orderly
    # Need to write date and Time to timing log
    if Configuration['TIMING']:

        plot_usage_stats(Configuration,debug = debug)
        timing_result = open(Configuration['MAINDIR']+'/Timing_Result.txt','a')
        timing_result.write(f'''The galaxy in directory {Configuration['FITTING_DIR']} started at {Configuration['START_TIME']}.
Finished preparations at {Configuration['PREP_END_TIME']} \n''')
        timing_result.write(f'''Converged to a galaxy size at {Configuration['END_TIME']}. \n''')
        timing_result.write(f'''It finished the whole process at {datetime.now()}
''')
        timing_result.close()
        log_statement = f'''Finished timing statistics for the galaxy in {Configuration['FITTING_DIR']}.
'''
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)

    cleanup_final(Configuration,Fits_Files, debug=debug)

finish_galaxy.__doc__ =f'''
 NAME:
    finish_galaxy

 PURPOSE:
    Write the final logs, Produce an Overview plot if required and make sure the directory is clean and organized.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration
    maximum_directory_length = the maximum string length of the input directory names

 OPTIONAL INPUTS:
    debug = False
    current_run = subproccess structure of the currently running tirific
    Fits_Files = Standard FAT dictionary with filenames

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This function should always be followed by a continue statement
'''

def transfer_errors(Configuration,fit_type='Undefined',debug = False):
    # Load the final file
    Tirific_Template = tirific_template(filename = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",debug= debug)
    # Get the errors from the input
    errors_to_transfer= ['VROT_ERR','VROT_2_ERR','INCL_ERR','INCL_2_ERR','PA_ERR','PA_2_ERR','SDIS_ERR','SDIS_2_ERR','Z0_ERR','Z0_2_ERR']
    FAT_Model = load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}_In.def",Variables=errors_to_transfer,unpack=False,debug=debug)
    # add to the templatethere
    Tirific_Template['GR_CONT']=' '
    Tirific_Template.insert('GR_CONT','RESTARTID','0')
    Tirific_Template.insert('MINDELTA','DISTANCE',Configuration['DISTANCE'])
    for key in errors_to_transfer:
        format = set_format(key[:-4])
        Tirific_Template.insert(key[:-4],f"# {key}",f"{' '.join([f'{x:{format}}' for x in FAT_Model[:,errors_to_transfer.index(key)]])}")
    # write back to the File
    tirific(Configuration,Tirific_Template, name = f"{fit_type}/{fit_type}.def", debug = debug)

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
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified
'''
