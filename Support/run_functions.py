#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to run external programs



class SofiaRunError(Exception):
    pass
class TirificRunError(Exception):
    pass

import read_functions as rf
import write_functions as wf
from support_functions import print_log
from clean_functions import clean_before_sofia,clean_after_sofia
import os
import subprocess

def sofia(Configuration, Fits_Files, hdr, supportdir):

    sofia_template = rf.sofia_template(supportdir+'/sofia_template.par')
    if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
        os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')
    os.chdir(Configuration['FITTING_DIR'])
    threshold = 5.
    counter = 3
    sofia_template['input.data'] = Fits_Files['FITTING_CUBE']
    beam_in_pixels = int(round(hdr['BMAJ']/((abs(hdr['CDELT1'])+abs(hdr['CDELT2']))/2.)))
    spatial_kernels = [0,beam_in_pixels,beam_in_pixels*2]
    if (hdr['NAXIS1']+hdr['NAXIS2'])/(2.*beam_in_pixels)  > 30:
        spatial_kernels.append(beam_in_pixels*3)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube is more than 30 beams across.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'])

    velocity_kernels = [0,3,6,12]
    if hdr['NAXIS3'] > 52:
        velocity_kernels.append(16)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'])
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par')
        print_log("RUN_SOFIA: Running SoFiA.",Configuration['OUTPUTLOG'],screen = True)
        sfrun = subprocess.Popen(["sofia2",'sofia_input.par'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sofia_run, sofia_warnings_are_annoying = sfrun.communicate()
        print_log(sofia_run.decode("utf-8"), Configuration['OUTPUTLOG'])

        if sfrun.returncode == 8:
            if threshold > 3.:
                log_statement = f'''RUN_SOFIA: We did not find a source at a threshold of {threshold}
{"":8s} RUN_SOFIA: Lowering the threshold and trying again."
'''
                print_log(log_statement,Configuration['OUTPUTLOG'])
                threshold -= 1
            else:
                clean_after_sofia(Configuration)
                log_statement = f'''RUN_SOFIA: We did not find a source above a threshold of {threshold}.
{"":8s}RUN_SOFIA: We cannot lower the threshold lower as the risk of fitting noise becomes too high.
{"":8s}Continuing to the next galaxy.
'''
                print_log(log_statement,Configuration['OUTPUTLOG'])
                raise SofiaRunError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'])
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration)

    os.chdir(Configuration['START_DIR'])


sofia.__doc__ ='''
;+
; NAME:
;       SOFIA
;
; PURPOSE:
;       Run SoFiA to create a mask Moment maps and a catalogue
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;        sofia(Configuration, Fits_Files)
;       RUN_SOFIA,allnew,new_dir,currentfitcube,catcatalogname,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log
;
;
; INPUTS:
;      Configuration, Fits_files hdr and supportdir
; OPTIONAL INPUTS:
;       LOG = name of the tracing log
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       FILE_TEST(),READ_FITS(),READ_TEMPLATE,SPAWN
;
'''
