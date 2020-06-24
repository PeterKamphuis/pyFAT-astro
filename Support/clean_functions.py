#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


import os
from datetime import datetime
from support_functions import print_log

def clean_before_sofia(Configuration):
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

def clean_after_sofia(Configuration):
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits','_cat.txt']

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
def cleanup(Configuration):

    if Configuration['START_POINT'] == 1:
        directories=['Optimized','Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output','Def_Files']
        files = [Configuration['BASE_NAME']+'.fits',Configuration['BASE_NAME']+'_opt.fits', 'tirific.def',Configuration['BASE_NAME']+'-BasicInfo.txt']
    elif Configuration['START_POINT'] == 2:
        directories=['Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Sofia_Output','Def_Files']
        files = ['tirific.def',Configuration['BASE_NAME']+'-BasicInfo.txt']
    elif Configuration['START_POINT'] == 3:
        directories=['Intermediate','Finalmodel','No_Warp','Moments','PV-Diagrams','Def_Files']
        files = ['tirific.def',Configuration['BASE_NAME']+'-BasicInfo.txt']

    print_log(f'''CLEANUP: We are cleaning the following directories:
{"":8s}CLEANUP: {','.join(directories)}
{"":8s}CLEANUP: and the following files:
{"":8s}CLEANUP: {','.join(files)}
''',Configuration['OUTPUTLOG'], screen =True)


    ext=['.fits','.log','.ps','.def']
    moments = ['mom0','mom1','mom2']
    #then specific files in the working directory
    os.chdir(Configuration['FITTING_DIR'])

    for dir in directories:
        if os.path.isdir(dir):
            if dir == 'Optimized' or dir == 'Intermediate':
                os.system(f'rm -f {dir}/Finalmodel* {dir}/Center_Convergence*')
            elif dir == 'Finalmodel':
                for fe in ext:
                    try:
                        os.remove(f'{dir}/{dir}.{fe}')
                    except FileNotFoundError:
                        pass
            elif dir == 'Moments':
                for moments in moments:
                    try:
                        os.remove(f'{dir}/Finalmodel_{moment}.fits')
                    except FileNotFoundError:
                        pass
                    try:
                        os.remove(f'{dir}/{Configuration["BASE_NAME"]}_{moment}.fits')
                    except FileNotFoundError:
                        pass
            elif dir == 'PV-Diagrams':
                try:
                    os.remove(f'{dir}/Finalmodel_xv.fits')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(f'{dir}/{Configuration["BASE_NAME"]}_sofia_xv.fits')
                except FileNotFoundError:
                    pass
                try:
                    os.remove(f'{dir}/{Configuration["BASE_NAME"]}_cc_xv.fits')
                except FileNotFoundError:
                    pass
            elif dir == 'Sofia_Output':
                os.system(f'rm -f {dir}/{Configuration["BASE_NAME"]}*_binmask* {dir}/{Configuration["BASE_NAME"]}*.ascii')
            elif dir == 'Def_Files':
                os.system(f'rm -f {dir}/*.def')

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

def finish_galaxy(Configuration,maximum_directory_length):
    # Need to write to results catalog
    output_catalogue = open(Configuration['OUTPUTCATALOGUE'],'a')
    output_catalogue.write(f"{Configuration['FITTING_DIR'].split('/')[-2]:{maximum_directory_length}s} {Configuration['AC1']:6d} {Configuration['AC2']:6d}   {Configuration['FINAL_COMMENT']} \n")
    output_catalogue.close()

    if Configuration['MAPS_OUTPUT'] == 5:
        log_statement = f'''!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
{"":8s}FAT did not run the full fitting routines on this galaxy.
{"":8s}Please check this log and output_catalogue carefully for what went wrong.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
'''
        print_log(log_statement,Configuration['OUTPUTLOG'], screen = True)

    # Need to organize the fitting output orderly
    # Need to write date and Time to timing log
    if Configuration['TIMING']:
        timing_result = open(Configuration['MAINDIR']+'/Timing_Result.txt','w')
        timing_result.write(f'''The galaxy in directory {Configuration['FITTING_DIR']} started at {Configuration['START_TIME']}.
Converged to a central position at {Configuration['CONVERGENCE_1']}.
It finished the whole process at {datetime.now()}
''')
        timing_result.close()
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
