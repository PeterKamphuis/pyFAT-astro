#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


import os,signal
import traceback
from datetime import datetime
from support_functions import print_log,finish_current_run
from fits_functions import make_moments
from write_functions import make_overview_plot,plot_usage_stats,tirific
from read_functions import tirific_template,load_tirific

def clean_before_sofia(Configuration, debug = False):
    if debug:
        print_log(f'''CLEAN_BEFORE_SOFIA: Starting.
''',Configuration['OUTPUTLOG'],debug = debug)
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

clean_before_sofia.__doc__ = '''
Clean up the sofia out put files before running SoFiA
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



clean_after_sofia.__doc__ = '''
Clean up the sofia output files by putting them in a dedicated directory.
'''
# cleanup dirty files before starting fitting
def cleanup(Configuration,Fits_Files, debug = False):

    if Configuration['START_POINT'] == 1:
        directories=['Finalmodel','Sofia_Output','Centre_Convergence','Def_Files','Extent_Convergence']
        files = [Fits_Files['FITTING_CUBE'],Fits_Files['OPTIMIZED_CUBE'], \
                'Centre_Convergence_In.def',Configuration['BASE_NAME']+'-BasicInfo.txt','Centre_Convergence/Centre_Convergence_In_Before_Smooth.def',\
                'Extent_Convergence_In.def','Extent_Convergence/Extent_Convergence_In_Before_Smooth.def']
    elif Configuration['START_POINT'] == 2:
        directories=['Finalmodel','Sofia_Output','Centre_Convergence','Def_Files','Extent_Convergence']
        files = ['Centre_Convergence_In.def',Configuration['BASE_NAME']+'-BasicInfo.txt','Centre_Convergence/Centre_Convergence_In_Before_Smooth.def']
    elif Configuration['START_POINT'] == 3:
        directories=['Finalmodel','Centre_Convergence','Def_Files','Extent_Convergence']
        files = ['Centre_Convergence_In.def','Centre_Convergence/Centre_Convergence_In_Before_Smooth.def']
    elif Configuration['START_POINT'] == 4 and Configuration['FINISHAFTER'] > 1:
        directories=['Finalmodel','Def_Files','Extent_Convergence']
        files = ['Centre_Convergence_In.def',]
    else:
        directories= None
        file =  None
    if directories:
        print_log(f'''CLEANUP: We are cleaning the following directories:
{"":8s}CLEANUP: {','.join(directories)}
{"":8s}CLEANUP: and the following files:
{"":8s}CLEANUP: {','.join(files)}
''',Configuration['OUTPUTLOG'], screen =True)


        ext=['.fits','.log','.ps','.def']
        moments = ['mom0','mom1','mom2', 'xv']
        #then specific files in the working directory
        os.chdir(Configuration['FITTING_DIR'])

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
                            try:
                                os.remove(f'{dir}/{dir}_Prev{fe}')
                            except FileNotFoundError:
                                pass
                            try:
                                os.remove(f'{dir}/{dir}_Prev2{fe}')
                            except FileNotFoundError:
                                pass
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


cleanup.__doc__ = '''
    ;+
; NAME:
;       cleanup
;
; PURPOSE:
;       Clean up any existing files before fitting start, it will only remove specific
;       files not directories or non-FAT files
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       cleanup(Configuration)
;
; INPUTS:
;     Configuration = Structure that has the current fitting directory and expected steps
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;       os
;
; EXAMPLE:
;
'''



def finish_galaxy(Configuration,maximum_directory_length,current_run = 'Not initialized', Fits_Files= None, debug = False):
    #make sure we are not leaving stuff
    finish_current_run(Configuration,current_run,debug=debug)

    # Need to write to results catalog
    output_catalogue = open(Configuration['OUTPUTCATALOGUE'],'a')
    output_catalogue.write(f"{Configuration['FITTING_DIR'].split('/')[-2]:{maximum_directory_length}s} {Configuration['CC_ACCEPTED']} {Configuration['EC_ACCEPTED']}   {Configuration['FINAL_COMMENT']} \n")
    output_catalogue.close()

    if Configuration['MAPS_OUTPUT'] == 5:
        log_statement = f'''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT did not run the full fitting routines for the galaxy in directory {Configuration['FITTING_DIR']}.
{"":8s}Please check this log and output_catalogue carefully for what went wrong.
{"":8s}The detected exit reason is {Configuration['FINAL_COMMENT']}.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
        sys.exit(1)
    elif Configuration['MAPS_OUTPUT'] == 4:
        print("We just want def files but this aint working yet")
    else:
        log_statement = f'''Producing final output in {Configuration['FITTING_DIR']}.
'''
        #
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
        # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
        if Configuration['FINISHAFTER'] > 0:
            if not os.path.isdir(Configuration['FITTING_DIR']+'/Finalmodel'):
                os.mkdir(Configuration['FITTING_DIR']+'/Finalmodel')
            if Configuration['FINISHAFTER'] == 1 and Configuration['START_POINT'] < 4:
                # As tirific does not transfer the errors we have to do This
                transfer_errors(Configuration,fit_stage='Centre_Convergence')
                os.symlink(f"{Configuration['FITTING_DIR']}/Centre_Convergence/Centre_Convergence.fits",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
                os.symlink(f"{Configuration['FITTING_DIR']}/Centre_Convergence/Centre_Convergence.def",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.def")
            elif Configuration['FINISHAFTER'] == 2:
                # As tirific does not transfer the errors we have to do This
                transfer_errors(Configuration,fit_stage='Extent_Convergence')
                os.symlink(f"{Configuration['FITTING_DIR']}/Extent_Convergence/Extent_Convergence.fits",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
                os.symlink(f"{Configuration['FITTING_DIR']}/Extent_Convergence/Extent_Convergence.def",f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.def")
            # We need to produce a FinalModel Directory with moment maps and an XV-Diagram of the model.
            if Fits_Files:
                if Configuration['FINISHAFTER'] == 1 and Configuration['START_POINT'] == 4:
                    print_log("Moments should exist",Configuration['OUTPUTLOG'],screen =True)
                else:
                    make_moments(filename = f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits",\
                             basename = 'Finalmodel', directory = f"{Configuration['FITTING_DIR']}/Finalmodel/",\
                             mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}",vel_unit = 'm/s',debug=debug)
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
Finished preparations at {Configuration['PREP_END_TIME']}
Converged to a central position at {Configuration['CC_END_TIME']}.
Converged to a galaxy size at {Configuration['EC_END_TIME']}.
It finished the whole process at {datetime.now()}
''')
        timing_result.close()
        log_statement = f'''Finished timing statistics for the galaxy in {Configuration['FITTING_DIR']}.
'''
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)
finish_galaxy.__doc__ = '''
    ;+
; NAME:
;       finish_galaxy
;
; PURPOSE:
;       Write the final logs, Produce an Overview plot if required and make sure the directory is clean and organized.
;
; CATEGORY:
;       Clean
;
; CALLING SEQUENCE:
;       finish_galaxy(Configuration)
;
; INPUTS:
;     Configuration = Dictionary that has the current fitting directory and expected steps
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;       datetime
: NOTES:
;  This function should always be followed by a continue statement
;
; EXAMPLE:
;
'''
def transfer_errors(Configuration,fit_stage='Not Initialized',debug = False):
    # Load the final file
    Tirific_Template = tirific_template(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",debug= debug)
    # Get the errors from the input
    errors_to_transfer= ['VROT_ERR','VROT_2_ERR','INCL_ERR','INCL_2_ERR','PA_ERR','PA_2_ERR','SDIS_ERR','SDIS_2_ERR']
    FAT_Model = load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}_In.def",Variables=errors_to_transfer,unpack=False,debug=debug)
    # add to the templatethere
    Tirific_Template.insert('GR_CONT','RESTARTID','0')
    for key in errors_to_transfer:
        if key[:-4] in ['SBR','SBR_2']:
            format = '.2e'
        else:
            format= '.2f'
        Tirific_Template.insert(key[:-4],f"# {key}",f"{' '.join([f'{x:{format}}' for x in FAT_Model[:,errors_to_transfer.index(key)]])}")
    # write back to the File
    tirific(Configuration,Tirific_Template, name = f"{fit_stage}/{fit_stage}.def", debug = debug)
