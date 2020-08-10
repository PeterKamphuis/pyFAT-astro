#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to run external programs

class BadSourceError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class TirificRunError(Exception):
    pass

from support_functions import print_log, convert_type,set_limits,sbr_limits,rename_fit_products,set_ring_size,get_usage_statistics,get_inner_fix
from clean_functions import clean_before_sofia,clean_after_sofia
from fits_functions import cut_cubes,extract_pv,make_moments
from modify_template import write_new_to_template,smooth_profile,set_cflux,check_sbr,regularise_profile,set_fitting_parameters,check_size
from astropy.wcs import WCS
from astropy.io import fits

import read_functions as rf
import write_functions as wf
import os
import time
import subprocess
import numpy as np
import traceback
import warnings
import re
from datetime import datetime
def central_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr,Initial_Parameters, debug = False):
    #This function will run while it does not return 1, this means when it is called it is not in acceptance

    #First we actually run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = 'run_cc', fit_stage = 'Centre_Convergence')

    accepted = check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = 'Centre_Convergence')
    set_cflux(Configuration,Tirific_Template,debug = debug)
    check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
    if accepted:
        Configuration['CC_ACCEPTED'] = True
    else:
        print_log('''CENTRAL_CONVERGE: Tirific ran the maximum amount of loops which means the fit is not accepted and we retry.
''',Configuration['OUTPUTLOG'],debug = debug)
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',hdr,debug=debug)
        smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr,debug=debug)
        #after smoothing the sbr we should check it
        check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
        xpos,ypos,vsys,pa,incl = rf.load_template(Tirific_Template,Variables = ['XPOS','YPOS','VSYS','PA','INCL'])
        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = 'run_cc',\
                               systemic = [vsys[0],Initial_Parameters['VSYS'][1]], \
                               inclination = [incl[0],Initial_Parameters['INCL'][1]],
                               pa = [pa[0],Initial_Parameters['PA'][1]],
                               rotation = [np.max(smoothed_vrot),Initial_Parameters['VROT'][1]],
                               ra = [xpos[0],Initial_Parameters['RA'][1]], dec = [ypos[0],Initial_Parameters['DEC'][1]], debug = debug)

        wf.tirific(Configuration,Tirific_Template,name = 'Centre_Convergence_In.def')

    return current_run


def fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,hdr, stage = 'initial',fit_stage='Undefined_Stage', debug = False):
    #if we have only a few rings we only smooth. else we fit a polynomial to the RC and smooth the SBR
    smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr,debug = debug)
    check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
    if Configuration['NO_RINGS'] < 6:
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',hdr,debug = debug)
    else:
        smoothed_vrot = regularise_profile(Configuration,Tirific_Template,'VROT',hdr,min_error = hdr['CDELT3']/1000.,debug = debug)
    if stage == 'after_ec':
        min_error = []
        pars_to_smooth = []
        if not Configuration['FIX_INCLINATION']:
            pars_to_smooth.append('INCL')
            min_error.append(2.*Configuration['LIMIT_MODIFIER'])
        if not Configuration['FIX_Z0']:
            pars_to_smooth.append('Z0')
            min_error.append(convertskyangle(0.1,Configuration['DISTANCE'],physical= True))
        if not Configuration['FIX_PA']:
            pars_to_smooth.append('PA')
            min_error.append(1.)
        if not Configuration['FIX_SDIS']:
            pars_to_smooth.append('SDIS')
            min_error.append(hdr['CDELT3']/2000.*Configuration['LIMIT_MODIFIER'])
        for key,min_err in zip(pars_to_smooth,min_error):
            if Configuration['NO_RINGS'] < 6:
                smoothed = smooth_profile(Configuration,Tirific_Template,key,hdr,debug = debug)
            else:
                smoothed = regularise_profile(Configuration,Tirific_Template,key,hdr,min_error = min_err,debug = debug)

    xpos,ypos,vsys,pa,incl = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = ['XPOS','YPOS','VSYS','PA','INCL'])
    set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = stage,\
                           systemic = [vsys[0], Configuration['CHANNEL_WIDTH']], \
                           inclination = [incl[0],10.],
                           pa = [pa[0],5.],
                           rotation = [np.max(smoothed_vrot),np.max(smoothed_vrot)*0.1],
                           ra = [xpos[0],abs(hdr['CDELT1'])], dec = [ypos[0],abs(hdr['CDELT2'])], debug = debug)
    os.system(f"cp {Configuration['FITTING_DIR']}{fit_stage}_In.def {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_In_Before_Smooth.def")
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_stage}_In.def')
    accepted,current_run = run_tirific(Configuration,current_run, stage = stage, fit_stage=fit_stage)


    return current_run

def check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = 'Undefined_Stage', debug = False):

    new_xpos,new_ypos,new_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = ['XPOS','YPOS','VSYS'])
    old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}_In.def",Variables = ['XPOS','YPOS','VSYS'])
    if Configuration['OUTER_RINGS_DOUBLED']:
        shift_beam_frac =0.3
    else:
        shift_beam_frac =0.15
    ra_lim = set_limits(shift_beam_frac*hdr['BMAJ'],np.max([abs(hdr['CDELT1']),1./3600.]),hdr['BMAJ'] )
    dec_lim =set_limits(shift_beam_frac*hdr['BMAJ'],np.max([abs(hdr['CDELT2']),1./3600.]),hdr['BMAJ'] )
    sys_lim = set_limits(0.5*hdr['CDELT3']/1000.,2.5, 2.*hdr['CDELT3']/1000. )
    if abs(new_xpos[0] - old_xpos[0]) > ra_lim or \
       abs(new_ypos[0] - old_ypos[0]) > dec_lim or \
       abs(new_vsys[0] - old_vsys[0]) > sys_lim:
        if Configuration['SIZE_IN_BEAMS'] < 25:
            apply_limit = 2*hdr['BMAJ']
        else:
            apply_limit = hdr['BMAJ']*Configuration['SIZE_IN_BEAMS']*0.08
        if abs(new_xpos[0] - old_xpos[0]) > apply_limit or\
           abs(new_ypos[0] - old_ypos[0]) > apply_limit:
           print_log(f'''CHECK_CONVERGENCE: The center shifted more than {apply_limit/hdr['BMAJ']} FWHM.
{"":8s}CHECK_CONVERGENCE: Not applying this shift
''', Configuration['OUTPUTLOG'],debug = debug)
           write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template, Variables = ['VROT',
                         'Z0', 'SBR', 'INCL','PA','SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','SDIS_2'])
           return 0
        else:
            try:
                old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev2.def",Variables = ['XPOS','YPOS','VSYS'])
                if abs(new_xpos[0] - old_xpos[0]) < ra_lim and \
                   abs(new_ypos[0] - old_ypos[0]) < dec_lim and \
                   abs(new_vsys[0] - old_vsys[0]) < sys_lim:
                   print_log(f'''CHECK_CONVERGENCE: The center shifted back to the old position. Moving on to the next stage.
''', Configuration['OUTPUTLOG'],debug = debug)
                   write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template)
                   return accepted
            except:
                pass
            print_log(f'''CHECK_CONVERGENCE: The center shifted too much trying again with new center.
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
            write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template)
            return 0
    else:
        print_log(f'''CHECK_CONVERGENCE: The center is accepted. The shift is:
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
        write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template)
        return accepted

def check_source(Configuration, Fits_Files, Catalogue, header, debug = False):


    name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa, \
        w50,err_f_sum, err_x,err_y,err_z= rf.sofia_catalogue(Configuration,header=header)

    x_min,x_max,y_min,y_max,z_min,z_max = convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int')
    x,y,z,ra,dec,v_app,f_sum,kin_pa,f_sum_err , err_x,err_y,err_z= convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa,err_f_sum, err_x,err_y,err_z])
    #How does sofia 2 deal with the fully flagged channels?
    # Need to check for that here if NaNs are included
    if f_sum < 0.:
            print_log(f'''CHECK_SOURCE: This galaxy has negative total flux. That will not work. Aborting.
    ''',Configuration['OUTPUTLOG'],debug = debug)
            raise BadSourceError('We found an initial negative total flux.')
    galaxy_box = [[z_min,z_max],[y_min,y_max],[x_min,x_max]]

    # If the provided distance  = -1 we assume a Hubble follow
    if Catalogue['DISTANCE'] == -1:
        Catalogue['DISTANCE'] == v_app/(1000.*H_0)
    if Catalogue['DISTANCE'] < 0.5:
        Catalogue['DISTANCE'] == 0.5


    #Check whether the cube is very large, if so cut it down

    new_box = cut_cubes(Configuration, Fits_Files, galaxy_box, header)
    #update our pixel values to match the new sizes
    for i in range(len(new_box)):
        shift = new_box[i,0]
        if i == 0:
            z -= shift; z_min -= shift; z_max -= shift
        elif i == 1:
            y -= shift; y_min -= shift; y_max -= shift
        elif i == 2:
            x -= shift; x_min -= shift; x_max -= shift

    Cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'],uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    data = Cube[0].data
    header = Cube[0].header

    Configuration['CHANNEL_WIDTH'] = header['CDELT3']/1000.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cube_wcs = WCS(header)
    # convert the boundaries to real coordinates
    ralow,declow,vellow = cube_wcs.wcs_pix2world(x_min,y_min,z_min,0.)
    rahigh,dechigh,velhigh = cube_wcs.wcs_pix2world(x_max,y_max,z_max,0.)
    DECboun = np.sort([float(declow),float(dechigh)])
    RAboun = np.sort([float(ralow),float(rahigh)])
    VELboun = np.sort([float(vellow),float(velhigh)])
    # We write the results of the cut cube to the log
    print_log(f'''CHECK_SOURCE: The source finder found the following center in pixels.
{"":8s}CHECK_SOURCE: RA center = {x} with boundaries {x_min}, {x_max}
{"":8s}CHECK_SOURCE: DEC center = {y} with boundaries {y_min}, {y_max}
{"":8s}CHECK_SOURCE: V_sys center = {z} with boundaries {z_min}, {z_max}
{"":8s}CHECK_SOURCE: This should correspond to the WCS coordinates:
{"":8s}CHECK_SOURCE: RA center = {ra} with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: V_sys center = {v_app/1000.:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'],debug = debug)

    #There is a factor of two missing here but this is necessary otherwise the maxima are far to small
    Configuration['MAX_SIZE_IN_BEAMS'] = int(round(np.sqrt(((x_max-x_min)/2.)**2+((y_max-y_min)/2.)**2) \
                /((header['BMAJ'])/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])]))\
                +3.))

    # Determine whether the centre is blanked or not
    Central_Flux = data[int(round(z)),int(round(y)),int(round(x))]
    print_log(f'''CHECK_SOURCE: In the center we find {Central_Flux} Jy/beam at the location:
{"":8s}CHECK_SOURCE: x,y,z = {int(round(x))}, {int(round(y))}, {int(round(z))}.
''',Configuration['OUTPUTLOG'],debug = debug)

    if not np.isfinite(Central_Flux):
        Configuration['EXCLUDE_CENTRAL'] = True
        print_log(f'''CHECK_SOURCE: The flux in the central part is blanked. We exclude the central rings.
''',Configuration['OUTPUTLOG'],debug = debug)
    else:
        Configuration['EXCLUDE_CENTRAL'] = False

    #Check that the source is bright enough
    Max_SNR = np.nanmax(data)/Configuration['NOISE']
    if Max_SNR < 2.5:
        log_statement = f'''CHECK_SOURCE: The maximum Signal to Noise in this cube is {Max_SNR} that is not enough for a fit.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)
        raise BadSourceError(log_statement)
    #Let's get the initial estimates of the PA and inclination from the axis ratios

    try:
        pa, inclination, SBR_initial, maj_extent = rf.guess_orientation(Configuration,Fits_Files, center = [x,y])
    except:
        print_log(f'''CHECK_SOURCE: We could not establish proper initial estimates from the moment maps.
{"":8s} {traceback.print_exc()}
''', Configuration['OUTPUTLOG'], screen = True,debug = debug)
        raise BadSourceError("No initial estimates. Likely the source is too faint.")
    if abs(pa[0]-kin_pa) < 25:
        pa[0] = (pa[0]/pa[1]+kin_pa/10.)/(1./pa[1]+1./10.)
    else:
        if np.isfinite(kin_pa):
            pa[0] = kin_pa

    # Size of the galaxy in beams
    Configuration['SIZE_IN_BEAMS'] = set_limits(int(round(maj_extent/(header['BMAJ']))),3,Configuration['MAX_SIZE_IN_BEAMS'])
    Configuration['NO_RINGS'] = (Configuration['SIZE_IN_BEAMS']/Configuration['RING_SIZE']+2)
    print_log(f'''CHECK_SOURCE: From the original Configuration and SoFiA we find:
{"":8s}CHECK_SOURCE: The maximum diameter we will allow  is  {2.*Configuration['MAX_SIZE_IN_BEAMS']} beams.
{"":8s}CHECK_SOURCE: We start with a diameter of {2*Configuration['SIZE_IN_BEAMS']} beams in the model.
{"":8s}CHECK_SOURCE: SoFiA found a PA of {kin_pa:.2f} and we use a PA = {pa[0]:.2f} +/- {pa[1]:.2f}
{"":8s}CHECK_SOURCE: We start with an inclination of {inclination[0]:.2f} +/- {inclination[1]:.2f}{"":8s}
{"":8s}CHECK_SOURCE: SoFiA found a W50 of {w50/1000.:.2f} km/s
{"":8s}CHECK_SOURCE: We will now check whether the size of the galaxy is sufficient.
''',Configuration['OUTPUTLOG'],debug = debug)

    set_ring_size(Configuration, debug = debug)



    if Configuration['HANNING']:
        vres = header['CDELT3']*2.
    else:
        vres = header['CDELT3']

    if inclination[0] < 40:
        max_vrot=w50/2./np.sin(np.radians(abs(inclination[0]+5.)))
    else:
        max_vrot=w50/2./np.sin(np.radians(abs(inclination[0])))
    max_vrot_dev=set_limits((VELboun[1]-VELboun[0])/4./np.sin(np.radians(inclination[0])),4.*vres,20*vres)
    #Write the info to the Basic info File
    wf.basicinfo(Configuration,initialize = True,
              RA=[ra,abs(err_x*header['CDELT1'])],
              DEC=[dec,abs(err_y*header['CDELT2'])],
              VSYS =[v_app,abs(err_z*header['CDELT3'])],
              PA=pa, Inclination = inclination, Max_Vrot = [max_vrot,max_vrot_dev], Tot_Flux = [f_sum,f_sum_err], V_mask = [VELboun[1]-VELboun[0],vres],
              Distance = Configuration['DISTANCE'] , DHI = maj_extent*3600.)


    # extract a PV-Diagram
    if Configuration['START_POINT'] < 3 or not os.path.exists(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits"):
        PV =extract_pv(Cube, pa[0], center=[ra,dec,v_app], convert = 1000.,
                       finalsize = [int(round(maj_extent/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])*1.25+header['NAXIS1']*0.2)),
                                    int(round(z_max-z_min)+10.)])
        if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
                os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')

        fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits",PV[0].data,PV[0].header)
    Cube.close()
    Initial_Parameters = {}
    Initial_Parameters['RA'] = [ra,abs(err_x*header['CDELT1'])]
    Initial_Parameters['DEC'] = [dec,abs(err_y*header['CDELT2'])]
    Initial_Parameters['VSYS'] =[v_app,abs(err_z*header['CDELT3'])]
    Initial_Parameters['SBR'] = SBR_initial
    Initial_Parameters['VROT'] = [max_vrot/1000.,max_vrot_dev/1000.]
    Initial_Parameters['PA'] = pa
    Initial_Parameters['INCL'] = inclination
    Initial_Parameters['FLUX'] = [f_sum,f_sum_err]
    return Initial_Parameters
check_source.__doc__='''

; NAME:
;       check_source(Configuration, Fits_Files, Tirific_Template, Catalogue)
;
; PURPOSE:
;       Check that the source found by sofia can be fitted and the cube is suitable for fitting
;
; CATEGORY:
;       run function
;
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      set_limits, print_log
;
; EXAMPLE:
;
;
'''


def extent_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr, debug = False,allowed_loops = 10.):
    fit_stage = 'Extent_Convergence'
    accepted,current_run = run_tirific(Configuration,current_run,stage = 'run_ec', fit_stage = fit_stage, debug= debug)
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def" , Tirific_Template, debug = debug)
    accepted_size = check_size(Configuration,Tirific_Template,hdr, debug=debug)
    check_sbr(Configuration,Tirific_Template,hdr,debug = debug,stage = 'run_ec')
    if accepted and accepted_size:
            Configuration['EC_ACCEPTED'] = True
    else:
        if not accepted:
            print_log('''EXTENT_CONVERGE: Tirific ran the maximum amount of loops which means the fit is not accepted and we smooth and retry.
''',Configuration['OUTPUTLOG'],debug = debug)
        else:
            print_log(f'''EXTENT_CONVERGE: FAT adjusted the rings. Refitting with new settings after smoothing them.
''',Configuration['OUTPUTLOG'],debug = debug)
        Configuration['EC_ACCEPTED'] = False
        if Configuration['EC_LOOPS'] > allowed_loops:
            print_log(f'''EXTENT_CONVERGE: We have ran the convergence more than {allowed_loops} times aborting the fit.
''',Configuration['OUTPUTLOG'],debug = debug)
            return current_run

        Configuration['INNER_FIX'] = get_inner_fix(Configuration, Tirific_Template,debug=debug)
        set_cflux(Configuration,Tirific_Template,debug = debug)
        keys_to_smooth =['INCL','PA','SDIS','Z0','VROT']
        for key in keys_to_smooth:
            smoothed = smooth_profile(Configuration,Tirific_Template,key,hdr,debug=debug)

        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = 'run_ec',\
                             debug = debug)

        #After smoothing the SBR we should check it again?
        check_sbr(Configuration,Tirific_Template,hdr,debug = debug,stage = 'run_ec')
        wf.tirific(Configuration,Tirific_Template,name = 'Extent_Convergence_In.def')


    return current_run


def make_full_resolution(Configuration,Tirific_Template,Fits_Files,fit_stage = 'Undefined_Stage', debug = False):
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template)
    Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Tirific_Template['LOOPS'] = "0"
    Tirific_Template['INIMODE'] = "0"
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_stage}_In.def')
    accepted,current_run = run_tirific(Configuration,'Not Initialized', stage = 'full_res', fit_stage=fit_stage)
    try:
        current_run.kill()
        Configuration['TIRIFIC_RUNNING'] = False
    except AttributeError:
        pass

def run_tirific(Configuration, current_run, stage = 'initial',fit_stage = 'Undefined_Stage', debug = False):
    # First move the previous fits
    rename_fit_products(Configuration,fit_stage = fit_stage)
    # Then if already running change restart file
    if Configuration['TIRIFIC_RUNNING']:
        print_log('''RUN_TIRIFIC: We are using an initialized tirific
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
    else:
        try:
            current_run.kill()
            Configuration['TIRIFIC_RUNNING'] = False
        except AttributeError:
            pass
        print_log('''RUN_TIRIFIC: We are starting a new TiRiFiC
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'w') as file:
            file.write("Initialized a new run")
        current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_stage}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
                                    stderr = subprocess.PIPE,cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        Configuration['TIRIFIC_RUNNING'] = True
        Configuration['TIRIFIC_PID'] = current_run.pid
    currentloop =1
    max_loop = 0
    counter = 0
    if Configuration['TIMING']:
        time.sleep(0.1)
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write(f"# Started Tirific at stage = {fit_stage}\n")
            CPU,mem = get_usage_statistics(current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
    print(f"RUN_TIRIFIC: Starting loop 1")
    triggered = False
    for tir_out_line in current_run.stdout:
        tmp = re.split(r"[/: ]+",tir_out_line.strip())
        counter += 1
        if (counter % 50) == 0:
            if Configuration['TIMING']:
                with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
                    if tmp[0] == 'L' and not triggered:
                        if tmp[1] == '1':
                            file.write("# Started the actual fitting \n")
                            triggered = True
                    CPU,mem = get_usage_statistics(current_run.pid)
                    file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
        if tmp[0] == 'L':
            if int(tmp[1]) != currentloop:
                print(f"RUN_TIRIFIC: Starting loop {tmp[1]} out of a maximum {max_loop}")
            currentloop  = int(tmp[1])
            if max_loop == 0:
                max_loop = int(tmp[2])
            Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
        if tmp[0].strip() == 'Finished':
            break
    if Configuration['TIMING']:
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write("# Finished this run \n")
            CPU,mem = get_usage_statistics(current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb\n")
    print(f"RUN_TIRIFIC: Finished the current tirific run.")
    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    while not os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def"):
        time.sleep(0.1)
    if currentloop != max_loop:
        return 1,current_run
    else:
        return 0,current_run

run_tirific.__doc__= '''

; NAME:
;       tirific(Configuration)
;
; PURPOSE:
;       Check whether we have an initialized tirific if not initialize and reun else run.
;
; CATEGORY:
;       run function
;
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      set_limits, print_log
;
; EXAMPLE:
;
;
'''

def sofia(Configuration, Fits_Files, hdr, supportdir, debug = False):

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
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)

    velocity_kernels = [0,3,6,12]
    if hdr['NAXIS3'] > 52:
        velocity_kernels.append(16)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par')
        print_log("RUN_SOFIA: Running SoFiA. \n",Configuration['OUTPUTLOG'],screen = True,debug = debug)
        sfrun = subprocess.Popen(["sofia2",'sofia_input.par'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sofia_run, sofia_warnings_are_annoying = sfrun.communicate()
        print_log(sofia_run.decode("utf-8"), Configuration['OUTPUTLOG'],debug = debug)

        if sfrun.returncode == 8:
            if threshold > 3.:
                log_statement = f'''RUN_SOFIA: We did not find a source at a threshold of {threshold}
{"":8s} RUN_SOFIA: Lowering the threshold and trying again."
'''
                print_log(log_statement,Configuration['OUTPUTLOG'],debug = debug)
                threshold -= 1
            else:
                clean_after_sofia(Configuration)
                log_statement = f'''RUN_SOFIA: We did not find a source above a threshold of {threshold}.
{"":8s}RUN_SOFIA: We cannot lower the threshold lower as the risk of fitting noise becomes too high.
{"":8s}Continuing to the next galaxy.
'''
                print_log(log_statement,Configuration['OUTPUTLOG'],debug = debug)
                raise SofiaRunError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'],debug = debug)
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration)

    os.chdir(Configuration['START_DIR'])
# to ensure compatible units and calculations with th models we make the maps ourselves
    mask = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}")
    del mask[0].header['C*3']
    mask[0].data[mask[0].data> 0.5] = 1.
    chan_map = np.nansum(mask[0].data,axis=0)
    mask[0].header['BITPIX'] = -32
    fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_chan.fits",chan_map,mask[0].header)
    mask.close()
    make_moments(filename = f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}",\
                 basename = f"{Configuration['BASE_NAME']}", directory = f"{Configuration['FITTING_DIR']}/Sofia_Output/",\
                 mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}",vel_unit = 'm/s')

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
