# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to run external programs

class BadSourceError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class SofiaFaintError(Exception):
    pass

from pyFAT.Support.support_functions import print_log, convert_type,set_limits,rename_fit_products,\
                              set_ring_size,calc_rings,get_usage_statistics,get_inner_fix,convertskyangle,\
                              finish_current_run, remove_inhomogeneities,get_from_template,set_format, \
                              set_rings, convertRADEC, is_available
from pyFAT.Support.clean_functions import clean_before_sofia,clean_after_sofia
from pyFAT.Support.fits_functions import cut_cubes,extract_pv,make_moments
from pyFAT.Support.modify_template import write_new_to_template,smooth_profile,set_cflux,fix_sbr, \
                                          regularise_profile,set_fitting_parameters,check_size, \
                                          no_declining_vrot, set_new_size,set_errors
from pyFAT.Support.constants import H_0
from astropy.wcs import WCS
from astropy.io import fits

import pyFAT.Support.read_functions as rf
import pyFAT.Support.write_functions as wf
import os
import sys
import time
import copy
import subprocess
import numpy as np
import traceback
import warnings
import re
from datetime import datetime

def fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,hdr, stage = 'initial',fit_type='Undefined', debug = False):
    if debug:
        print_log(f'''FIT_SMOOTHED_CHECK: Starting stage {stage} and fit_type {fit_type}.
''',Configuration['OUTPUTLOG'],debug = debug)
    #if we have only a few rings we only smooth. else we fit a polynomial to the RC and smooth the SBR
    #smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr, min_error= np.max([float(Tirific_Template['CFLUX']),float(Tirific_Template['CFLUX_2'])]),debug = debug)
    fix_sbr(Configuration,Tirific_Template,smooth = True, debug = debug)
    if stage == 'after_cc':
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',min_error = float(hdr['CDELT3']/1000.),debug = debug)
    else:
        smoothed_vrot = regularise_profile(Configuration,Tirific_Template,'VROT',min_error = float(hdr['CDELT3']/1000.),debug = debug)
    no_declining_vrot(Configuration, Tirific_Template, debug = debug)
    if stage in ['after_ec', 'after_os']:
        min_error = []
        pars_to_smooth = []
        not_to_smooth = []
        fixed_errors = []
        if not Configuration['FIX_INCLINATION'][0]:
            pars_to_smooth.append('INCL')
            min_error.append(set_limits(3.*np.mean(Configuration['LIMIT_MODIFIER']),2,5))
        else:
            not_to_smooth.append('INCL')
            fixed_errors.append(set_limits(2.*np.mean(Configuration['LIMIT_MODIFIER']),2,5))
        if not Configuration['FIX_Z0'][0]:
            pars_to_smooth.append('Z0')
            min_error.append(convertskyangle(Configuration,0.1,Configuration['DISTANCE'],physical= True))
        else:
            not_to_smooth.append('Z0')
            fixed_errors.append(convertskyangle(Configuration,0.1,Configuration['DISTANCE'],physical= True))
        if not Configuration['FIX_PA'][0]:
            pars_to_smooth.append('PA')
            min_error.append(2.)
        else:
            not_to_smooth.append('PA')
            fixed_errors.append(1.)

        if not Configuration['FIX_SDIS'][0]:
            pars_to_smooth.append('SDIS')
            min_error.append(hdr['CDELT3']/2000.*Configuration['LIMIT_MODIFIER'])
        else:
            not_to_smooth.append('SDIS')
            fixed_errors.append(hdr['CDELT3']/2000.*Configuration['LIMIT_MODIFIER'])

        for key,min_err in zip(pars_to_smooth,min_error):
            smoothed = regularise_profile(Configuration,Tirific_Template,key,min_error = min_err,debug = debug)

        for key,min_err in zip(not_to_smooth,fixed_errors):
            set_errors(Configuration,Tirific_Template,key,min_error = min_err,debug = debug)
    #If our fit stage is after cc we want to make sure we do an extra check on low inclinations or small Galaxies
    #if stage =='after_cc' and (Configuration['SIZE_IN_BEAMS'] < 5 or incl[0] < 50.):
    #    check_inclination(Configuration,Tirific_Template,Fits_Files,debug=debug)
    if stage == 'after_cc':

        Tirific_Template['LOOPS'] = 1.
        if Configuration['OPTIMIZED']:
            Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
            Tirific_Template['INIMODE'] = "0"
            finish_current_run(Configuration,current_run,debug= debug)
    else:
        Tirific_Template['LOOPS'] = 10.
        xpos,ypos,vsys,pa,incl = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",Variables = ['XPOS','YPOS','VSYS','PA','INCL'])
        parameters = {'VSYS': [vsys[0], Configuration['CHANNEL_WIDTH']], \
                      'XPOS': [xpos[0],abs(hdr['CDELT1'])],
                      'YPOS': [ypos[0],abs(hdr['CDELT2'])],
                      'INCL':[incl[0],10.],
                      'PA':  [pa[0],5.],
                      'VROT': [np.max(smoothed_vrot),np.max(smoothed_vrot)*0.1],  }

        set_fitting_parameters(Configuration, Tirific_Template,stage = stage,\
                          initial_estimates = parameters, debug = debug)
    os.system(f"cp {Configuration['FITTING_DIR']}{fit_type}_In.def {Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Last_Unsmoothed_Input.def")
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def',debug=debug)
    accepted,current_run = run_tirific(Configuration,current_run, stage = stage, fit_type=fit_type,debug=debug)


    return current_run

def check_inclination(Configuration,Tirific_Template,Fits_Files, fit_type = 'Unstaged', debug =False):
    if debug:
        print_log(f'''CHECK_INCLINATION: estimating whether our inclination estimate is decent.
''',Configuration['OUTPUTLOG'],debug = True)
    incl_to_check = np.linspace(5.,50.,20)
    # Let's make a directory to put things in
    tmp_stage = 'tmp_incl_check'
    try:
        os.mkdir(f"{Configuration['FITTING_DIR']}/{tmp_stage}")
    except:
        pass
    other_run = Configuration['TIRIFIC_RUNNING']
    Configuration['TIRIFIC_RUNNING'] = False
    #and a copy of the tirific template
    Check_Template = copy.deepcopy(Tirific_Template)
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Check_Template,debug = debug)
    Check_Template['LOOPS'] = '0'
    Check_Template['INIMODE'] = '0'
    Check_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Check_Template['RESTARTNAME'] = f"Logs/restart_{tmp_stage}.txt"
    out_keys = ['LOGNAME','OUTSET','TIRDEF']
    out_extensions = ['log','fits','def']
    incl_run= 'Not Initialized'
    for i,key in enumerate(out_keys):
        Check_Template[key] = f"{tmp_stage}/{tmp_stage}.{out_extensions[i]}"
    current = rf.load_template(Check_Template,Variables= ['VROT','INCL','PA','XPOS','YPOS'])


    vobs = [x*np.sin(np.radians(y)) for x,y in zip(current[0][:],current[1][:])]
    if debug:
        print_log(f'''CHECK_INCLINATION: These are the values we get as input
{'':8s}Inclination = {current[1][:]}
{'':8s}Vrot = {current[0][:]}
{'':8s}Vobs = {vobs}
''',Configuration['OUTPUTLOG'],debug = False)
    mom_chi = []
    model_mom0 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT0']}")
    model_mom0 = remove_inhomogeneities(Configuration,model_mom0,inclination=float(current[1][0]), pa = float(current[2][0]), center = [current[3][0],current[4][0]],debug=debug)
    chan_map = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['CHANNEL_MAP']}")
    noisemap = np.sqrt(chan_map[0].data)*Configuration['NOISE']/np.nanmax(model_mom0[0].data)
    model_mom0[0].data = model_mom0[0].data/np.nanmax(model_mom0[0].data)
    for incl in incl_to_check:
        #print(f'We are doing this inclination {incl}')
        vrot = [x/np.sin(np.radians(incl)) for x in vobs]
        for key in ['INCL','INCL_2']:
            Check_Template[key]= f"{incl:.2f}"
        for key in ['VROT','VROT_2']:
            Check_Template[key]= f"{' '.join([f'{x:.2f}' for x in vrot])}"
        wf.tirific(Configuration,Check_Template,name = f'{tmp_stage}_In.def',debug=debug)
        accepted,incl_run = run_tirific(Configuration,incl_run, stage = 'incl_check', fit_type=tmp_stage,debug=debug)
        make_moments(filename = f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}.fits", basename = 'tmp_incl', directory = f"{Configuration['FITTING_DIR']}{tmp_stage}/",\
                     moments = [0],mask_cube = f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['MASK']}", \
                     overwrite = True, log= Configuration['OUTPUTLOG'], vel_unit = 'm/s',debug=debug)
        #make_moments(filename = f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}.fits", basename = 'tmp_incl', directory = f"{Configuration['FITTING_DIR']}{tmp_stage}/",\
        #             moments = [0],level = 3.*Configuration['NOISE'], \
        #             overwrite = True, log= Configuration['OUTPUTLOG'], vel_unit = 'm/s',debug=debug)
        incl_mom0 = fits.open(f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_incl_mom0.fits")
        if debug:
            try:
                os.remove(f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
            except:
                pass
            os.rename(f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_incl_mom0.fits",f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
        chi = np.nansum((model_mom0[0].data[noisemap > 0.]/np.nanmax(model_mom0[0].data)-incl_mom0[0].data[noisemap > 0.]/np.nanmax(incl_mom0[0].data))**2/noisemap[noisemap > 0.]**2)
        mom_chi.append(chi)
    finish_current_run(Configuration, incl_run)
    Configuration['TIRIFIC_RUNNING'] = other_run
    low= np.where(mom_chi == np.nanmin(mom_chi))[0]
    new_incl = incl_to_check[low]
    if debug:
        print_log(f'''CHECK_INCLINATION: This is the new inclination {new_incl} it was {current[1][0]}.
{'':8s}mom_chi = {mom_chi}
{'':8s}low = {low}
''',Configuration['OUTPUTLOG'],debug = False)


    #exit()
    if new_incl < current[1][0]:
        if debug:
            print_log(f'''CHECK_INCLINATION: The inclination has decreased writing to file.
''',Configuration['OUTPUTLOG'],debug = False)
        for key in ['INCL','INCL_2']:
                Check_Template[key]= f"{new_incl:.2f}"
        vrot = [x/np.sin(np.radians(new_incl)) for x in vobs]
        for key in ['VROT','VROT_2']:
                Check_Template[key]= f"{' '.join([f'{x:.2f}' for x in vrot])}"
        Tirific_Template = copy.deepcopy(Check_Template)
        Check_Template = []
        wf.tirific(Configuration,Tirific_Template,name = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",debug=debug)


def check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_type = 'Undefined_Stage', debug = False):
    if debug:
        print_log(f'''CHECK_CENTRAL_CONVERGE: Starting stage {fit_type}.
''',Configuration['OUTPUTLOG'],debug = debug)
    new_xpos,new_ypos,new_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_type}_In.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    if Configuration['OUTER_RINGS_DOUBLED']:
        shift_beam_frac =Configuration['SIZE_IN_BEAMS']*0.05
    else:
        shift_beam_frac =0.15
    ra_lim = set_limits(shift_beam_frac*Configuration['BEAM'][0]/3600.,np.max([abs(hdr['CDELT1']),1./3600.]),Configuration['BEAM'][0]/3600.,debug = debug )
    dec_lim =set_limits(shift_beam_frac*Configuration['BEAM'][0]/3600.,np.max([abs(hdr['CDELT2']),1./3600.]),Configuration['BEAM'][0]/3600.,debug = debug )
    sys_lim = set_limits(0.5*Configuration['CHANNEL_WIDTH'],2.5, 2.*Configuration['CHANNEL_WIDTH'],debug = debug )
    if abs(new_xpos[0] - old_xpos[0]) > ra_lim or \
       abs(new_ypos[0] - old_ypos[0]) > dec_lim or \
       abs(new_vsys[0] - old_vsys[0]) > sys_lim:
        if Configuration['SIZE_IN_BEAMS'] < 20:
            apply_limit = 2*Configuration['BEAM'][0]/3600.
        else:
            apply_limit = Configuration['BEAM'][0]/3600.*Configuration['SIZE_IN_BEAMS']*0.2
        if abs(new_xpos[0] - old_xpos[0]) > apply_limit or\
           abs(new_ypos[0] - old_ypos[0]) > apply_limit:
           print_log(f'''CHECK_CONVERGENCE: The center shifted more than {apply_limit/Configuration['BEAM'][0]/3600.} FWHM.
{"":8s}CHECK_CONVERGENCE: Not applying this shift
''', Configuration['OUTPUTLOG'],debug = debug)
           write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template, Variables = ['VROT',
                         'Z0', 'SBR', 'INCL','PA','SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','SDIS_2'],debug=debug)
           return 0
        else:
            try:
                old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev2.def",Variables = ['XPOS','YPOS','VSYS'],debug=debug)
                if abs(new_xpos[0] - old_xpos[0]) < ra_lim and \
                   abs(new_ypos[0] - old_ypos[0]) < dec_lim and \
                   abs(new_vsys[0] - old_vsys[0]) < sys_lim:
                   print_log(f'''CHECK_CONVERGENCE: The center shifted back to the old position. Moving on to the next stage.
''', Configuration['OUTPUTLOG'],debug = debug)
                   write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
                   return accepted
            except:
                pass
            print_log(f'''CHECK_CONVERGENCE: The center shifted too much trying again with new center.
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
            write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
            return 0
    else:
        print_log(f'''CHECK_CONVERGENCE: The center is accepted. The shift is:
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
        write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
        return accepted

def check_source(Configuration, Fits_Files, debug = False):
    if debug:
        print_log(f'''CHECK_SOURCE: Starting.
''',Configuration['OUTPUTLOG'],debug = debug)

    name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa, \
        w50,err_f_sum, err_x,err_y,err_z= rf.sofia_catalogue(Configuration,Fits_Files,debug=debug)

    x_min,x_max,y_min,y_max,z_min,z_max = convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int')
    x,y,z,ra,dec,v_app,f_sum,kin_pa,f_sum_err , err_x,err_y,err_z= convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa,err_f_sum, err_x,err_y,err_z])
    #How does sofia 2 deal with the fully flagged channels?
    # Need to check for that here if NaNs are included
    if f_sum < 0.:
        print_log(f'''CHECK_SOURCE: This galaxy has negative total flux. That will not work. Aborting.
''',Configuration['OUTPUTLOG'],debug = debug)
        raise BadSourceError('We found an initial negative total flux.')
    galaxy_box = [[z_min,z_max],[y_min,y_max],[x_min,x_max]]
    if debug:
        print_log(f'''CHECK_SOURCE:  From the catalogue we got {Configuration['DISTANCE']}
''',Configuration['OUTPUTLOG'],debug = False,screen=True)
    # If the provided distance  = -1 we assume a Hubble follow
    if float(Configuration['DISTANCE']) == -1.:
        Configuration['DISTANCE'] = v_app/(1000.*H_0)
    if float(Configuration['DISTANCE']) < 0.5:
        Configuration['DISTANCE'] = 0.5
    if debug:
        print_log(f'''CHECK_SOURCE: We use a distance of {Configuration['DISTANCE']}.
''',Configuration['OUTPUTLOG'],debug = False,screen=True)
    #Check whether the cube is very large, if so cut it down

    new_box = cut_cubes(Configuration, Fits_Files, galaxy_box,debug=debug)
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


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cube_wcs = WCS(header)
    # convert the boundaries to real coordinates
    ralow,declow,vellow = cube_wcs.wcs_pix2world(x_min,y_min,z_min,1)
    rahigh,dechigh,velhigh = cube_wcs.wcs_pix2world(x_max,y_max,z_max,1)
    DECboun = np.sort([float(declow),float(dechigh)])
    RAboun = np.sort([float(ralow),float(rahigh)])
    VELboun = np.sort([float(vellow),float(velhigh)])
    rahr,dechr = convertRADEC(Configuration,ra,dec,debug=debug)
    # We write the results of the cut cube to the log
    print_log(f'''CHECK_SOURCE: The source finder found the following center in pixels.
{"":8s}CHECK_SOURCE: RA center = {x} with boundaries {x_min}, {x_max}
{"":8s}CHECK_SOURCE: DEC center = {y} with boundaries {y_min}, {y_max}
{"":8s}CHECK_SOURCE: V_sys center = {z} with boundaries {z_min}, {z_max}
{"":8s}CHECK_SOURCE: This should correspond to the WCS coordinates:
{"":8s}CHECK_SOURCE: RA center = {ra} deg or {rahr}  with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} deg  or {dechr}  with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: V_sys center = {v_app/1000.:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'],debug = debug)

    #There is a factor of two missing here but this is necessary otherwise the maxima are far to small
    Configuration['MAX_SIZE_IN_BEAMS'] = int(round(np.sqrt(((x_max-x_min)/2.)**2+((y_max-y_min)/2.)**2) \
                /(Configuration['BEAM_IN_PIXELS'][0])+5.))

    pa, inclination, SBR_initial, maj_extent,x_new,y_new,VROT_initial = rf.guess_orientation(Configuration,Fits_Files, center = [x,y],debug=debug)

    if x_new != x or y_new != y:
        x=x_new
        y=y_new
        ra,dec,v_app = cube_wcs.wcs_pix2world(x,y,z,1)
        print_log(f'''CHECK_SOURCE: The center is updated to.
{"":8s}CHECK_SOURCE: RA center = {ra} with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: V_sys center = {v_app/1000.:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'],debug = debug)
    if np.sum(pa) == 0. or any(np.isnan(pa)) or \
        np.sum(inclination) == 0. or any(np.isnan(inclination)) or \
        np.sum(maj_extent) == 0. or np.isnan(maj_extent) or \
        np.sum(SBR_initial) == 0. or all(np.isnan(SBR_initial)) or \
        np.sum(VROT_initial) == 0. or all(np.isnan(VROT_initial)):
        print_log(f'''CHECK_SOURCE: We could not establish proper initial estimates from the moment maps.
''', Configuration['OUTPUTLOG'], screen = True,debug = debug)
        raise BadSourceError("No initial estimates. Likely the source is too faint.")
     # Determine whether the centre is blanked or not
    Central_Flux = np.mean(data[int(round(z-1)):int(round(z+1)),\
                                int(round(y-Configuration['BEAM_IN_PIXELS'][0]/2.)):int(round(y+Configuration['BEAM_IN_PIXELS'][0]/2.)),\
                                int(round(x-Configuration['BEAM_IN_PIXELS'][0]/2.)):int(round(x+Configuration['BEAM_IN_PIXELS'][0]/2.))])
    print_log(f'''CHECK_SOURCE: In the center we find an average flux of  {Central_Flux} Jy/beam around the location:
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


    if abs(pa[0]-kin_pa) < 25:
        pa[0] = (pa[0]/pa[1]+kin_pa/10.)/(1./pa[1]+1./10.)
    else:
        if np.isfinite(kin_pa):
            pa[0] = kin_pa

    # Size of the galaxy in beams
    Configuration['SIZE_IN_BEAMS'] = set_limits(maj_extent/(Configuration['BEAM'][0]/3600.),1.0,Configuration['MAX_SIZE_IN_BEAMS'])
    if Configuration['SIZE_IN_BEAMS'] <= Configuration['TOO_SMALL_GALAXY']:
        print_log(f'''CHECK_SOURCE: This galaxy has an estimated size of  {2*Configuration['SIZE_IN_BEAMS']} beams in diameter.
{'':8s}This is not large enough too fit. We will exit this fit.
''',Configuration['OUTPUTLOG'],debug = debug)
        raise BadSourceError('The extracted source is too small')
    set_ring_size(Configuration, debug = debug)
    new_radii = set_rings(Configuration, debug = debug)
    old_radii = np.linspace(0.,Configuration['BEAM'][0]*(len(SBR_initial)-1),len(SBR_initial))
    new_radii = set_rings(Configuration, debug = debug)
    SBR_initial = np.interp(new_radii,old_radii, SBR_initial)
    old_radii = np.linspace(0.,Configuration['BEAM'][0]*(len(VROT_initial)-1),len(VROT_initial))
    VROT_initial = np.interp(new_radii,old_radii, VROT_initial)
    Configuration['OUTER_SLOPE_START'] = Configuration['NO_RINGS']
    if debug:
        print_log(f'''CHECK_SOURCE: Interpolating the SBR and VROT estimates to these radi.
{'':8s} {new_radii}
{'':8s} We got SBR = {SBR_initial}, VROT = {VROT_initial}
''',Configuration['OUTPUTLOG'],debug = False)


    # The extent is fairly well determined and the maximum should be no more than +3 beams and a minimum no less than 4
    # Swithcing here from doubled outer rings causes problems though

    Configuration['MAX_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']+2.,1,Configuration['MAX_SIZE_IN_BEAMS'])
    Configuration['MIN_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']-4.,Configuration['TOO_SMALL_GALAXY'],Configuration['SIZE_IN_BEAMS'])
    Configuration['NO_RINGS'] = calc_rings(Configuration,debug=debug)
    Configuration['LAST_RELIABLE_RINGS'] = [Configuration['NO_RINGS'],Configuration['NO_RINGS']]
    print_log(f'''CHECK_SOURCE: From the original Configuration and SoFiA we find:
{"":8s}CHECK_SOURCE: The maximum diameter we will allow  is  {2.*Configuration['MAX_SIZE_IN_BEAMS']} beams.
{"":8s}CHECK_SOURCE: We start with a diameter of {2*Configuration['SIZE_IN_BEAMS']} beams in the model.
{"":8s}CHECK_SOURCE: SoFiA found a PA of {kin_pa:.2f} and we use a PA = {pa[0]:.2f} +/- {pa[1]:.2f}
{"":8s}CHECK_SOURCE: We start with an inclination of {inclination[0]:.2f} +/- {inclination[1]:.2f}{"":8s}
{"":8s}CHECK_SOURCE: SoFiA found a W50 of {w50:.2f} km/s
{"":8s}CHECK_SOURCE: We will use {2.*Configuration['NO_RINGS']} for the model with a ring size of {Configuration['RING_SIZE']}.
''',Configuration['OUTPUTLOG'],debug = debug)




    if Configuration['HANNING']:
        vres = (Configuration['CHANNEL_WIDTH']*2)/(2.*np.sqrt(2.*np.log(2)))
    else:
        vres = (Configuration['CHANNEL_WIDTH']*1.2)/(2.*np.sqrt(2.*np.log(2)))
    #if inclination[0] < 40:
    #    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0]+5.)))
    #else:
    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0])))
    max_vrot_dev=set_limits((VELboun[1]-VELboun[0])/4000./np.sin(np.radians(inclination[0])),4.*vres,20*vres,debug=debug)


    #Write the info to the Basic info File
    wf.basicinfo(Configuration,initialize = True,
              RA=[ra,abs(err_x*header['CDELT1'])],
              DEC=[dec,abs(err_y*header['CDELT2'])],
              VSYS =[v_app,abs(err_z*header['CDELT3'])],
              PA=pa, Inclination = inclination, Max_Vrot = [max_vrot,max_vrot_dev], Tot_Flux = [f_sum,f_sum_err], V_mask = [VELboun[1]-VELboun[0],vres],
              Distance = Configuration['DISTANCE'] , DHI = maj_extent*3600.,debug=debug)


    # extract a PV-Diagram
    if Configuration['START_POINT'] < 3 or not os.path.exists(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits"):
        PV =extract_pv(Configuration,Cube, pa[0], center=[ra,dec,v_app], convert = 1000.,
                       finalsize = [int(round(maj_extent/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])*1.25+header['NAXIS1']*0.2)),
                                    int(round(z_max-z_min)+10.)],debug=debug)
        if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
                os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')

        fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits",PV[0].data,PV[0].header)
    Cube.close()
    Initial_Parameters = {}
    Initial_Parameters['XPOS'] = [ra,set_limits(abs(err_x*header['CDELT1']),0.1/3600.*Configuration['BEAM'][0],3./3600.*Configuration['BEAM'][0] )]
    Initial_Parameters['YPOS'] = [dec,set_limits(abs(err_y*header['CDELT2']),0.1/3600.*Configuration['BEAM'][0],3./3600.*Configuration['BEAM'][0] )]
    Initial_Parameters['VSYS'] =[v_app,set_limits(abs(err_z*header['CDELT3']),abs(header['CDELT3']),5.*abs(header['CDELT3']))]
    Initial_Parameters['SBR_profile'] = SBR_initial
    #Initial_Parameters['VROT'] = [max_vrot/1000.,max_vrot_dev/1000.]
    Initial_Parameters['VROT_profile'] = VROT_initial
    Initial_Parameters['VROT'] = [np.max(VROT_initial) ,(np.max(VROT_initial)-np.min(VROT_initial))]
    Initial_Parameters['PA'] = pa
    Initial_Parameters['INCL'] = inclination
    Initial_Parameters['FLUX'] = [f_sum,f_sum_err]
    return Initial_Parameters

check_source.__doc__='''
 NAME:
    check_source

 PURPOSE:
    Check that the source found by sofia can be fitted and the cube is suitable for fitting

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Locations of the Fits files used by FAT
    header = header of the input cube

 OPTIONAL INPUTS:
    debug = False

 KEYWORD PARAMETERS:

 OUTPUTS:
    Initial_Parameters = Dictionary with all initial values obtained from the Sofia source finding and processing.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 EXAMPLE:
'''

def make_full_resolution(Configuration,Tirific_Template,Fits_Files,fit_type = 'Undefined', current_run = 'Not zed', debug = False):
    if debug:
        print_log(f'''MAKE_FULL_RESOLUTION: creating full resolution for stage {fit_type}.
''',Configuration['OUTPUTLOG'],debug = debug)
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug = debug)
    Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Tirific_Template['LOOPS'] = "0"
    Tirific_Template['INIMODE'] = "0"
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def')
    finish_current_run(Configuration,current_run,debug= debug)
    accepted,current_run = run_tirific(Configuration,'Not zed', stage = 'full_res', fit_type=fit_type,debug = debug)
    finish_current_run(Configuration,current_run,debug= debug)

def one_step_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr, debug = False,allowed_loops = 10.):
    if debug:
        print_log(f'''ONE_STEP_CONVERGENCE: Starting with loop {Configuration['OS_LOOPS']} out of maximum {allowed_loops}.
''',Configuration['OUTPUTLOG'],debug = True)
    fit_type = 'One_Step_Convergence'
    stage = 'run_os'
    #First we run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = stage, fit_type = fit_type, debug= debug)
    #Then we load the produced output into our template
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def" , Tirific_Template, debug = debug)
    #Check that the centre does not deviate too much
    accepted_central = check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_type = fit_type,debug=debug)
    # Check whether we have the correct sizes,
    accepted_size = check_size(Configuration,Tirific_Template, fit_type = fit_type, stage = stage, current_run = current_run, debug=debug,Fits_Files=Fits_Files)
    if accepted and accepted_size and accepted_central:
        Configuration['OS_ACCEPTED'] = True
    else:
        Configuration['OS_ACCEPTED'] = False
        if Configuration['OS_LOOPS'] > allowed_loops:
                print_log(f'''ONE_STEP_CONVERGENCE: We have ran the convergence more than {allowed_loops} times aborting the fit.
    ''',Configuration['OUTPUTLOG'])
                return current_run
        if not accepted:
            print_log(f'''ONE_STEP_CONVERGENCE: Tirific ran the maximum amount of loops hence we do not accept and we smooth and retry.
''',Configuration['OUTPUTLOG'])
        if not accepted_central:
            print_log(f'''ONE_STEP_CONVERGENCE: The center varied too much hence we do not accept and we smooth and retry.
''',Configuration['OUTPUTLOG'])
        if not accepted_size:
            print_log(f'''ONE_STEP_CONVERGENCE: FAT adjusted the rings. Refitting with new settings after smoothing them.
''',Configuration['OUTPUTLOG'])
                    # First we fix the SBR we are left with also to set the reliable ring to configuration.
        fix_sbr(Configuration,Tirific_Template,debug = debug)    # Then we determine the inner rings that should remain fixed

        Configuration['INNER_FIX'] = get_inner_fix(Configuration, Tirific_Template,debug=debug)
        set_cflux(Configuration,Tirific_Template,debug = debug)
        keys_to_smooth =['INCL','PA','SDIS','Z0','VROT']
        min_errors = [3.*np.mean(Configuration['LIMIT_MODIFIER']),2.,hdr['CDELT3']/(2000.*Configuration['LIMIT_MODIFIER']), \
                        convertskyangle(Configuration,0.1,Configuration['DISTANCE'],physical= True)/Configuration['LIMIT_MODIFIER'],\
                        hdr['CDELT3']/(1000.*Configuration['LIMIT_MODIFIER'])]
        for j,key in enumerate(keys_to_smooth):
            #Smoothing the profile also fixes it
            smoothed = smooth_profile(Configuration,Tirific_Template,key,debug=debug,min_error=min_errors[j])

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_os',\
                             debug = debug)

        wf.tirific(Configuration,Tirific_Template,name = 'One_Step_Convergence_In.def',debug = debug)
    return current_run

def run_tirific(Configuration, current_run, stage = 'initial',fit_type = 'Undefined', debug = False):
    if debug:
        print_log(f'''RUN_TIRIFIC: Starting a new run in stage {stage} and fit_type {fit_type}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    # First move the previous fits
    rename_fit_products(Configuration,fit_type = fit_type, stage=stage, debug = debug)
    # Then if already running change restart file
    if Configuration['TIRIFIC_RUNNING']:
        print_log(f'''RUN_TIRIFIC: We are using an initialized tirific in {Configuration['FITTING_DIR']}
''',Configuration['OUTPUTLOG'], screen = True)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_type}.txt",'a') as file:
            file.write("Restarting from previous run")
    else:
        print_log(f'''RUN_TIRIFIC: We are starting a new TiRiFiC in {Configuration['FITTING_DIR']}
''',Configuration['OUTPUTLOG'], screen = True)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_type}.txt",'w') as file:
            file.write("Initialized a new run")
        current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_type}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
                                    stderr = subprocess.PIPE,cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        Configuration['TIRIFIC_RUNNING'] = True
        Configuration['TIRIFIC_PID'] = current_run.pid
    currentloop =1
    max_loop = 0
    counter = 0
    if Configuration['TIMING']:
        time.sleep(0.1)
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write(f"# Started Tirific at stage = {fit_type}\n")
            CPU,mem = get_usage_statistics(Configuration,current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
    else:
        time.sleep(0.1)
    print(f"\r RUN_TIRIFIC: 0 % Completed", end =" ",flush = True)
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
                    CPU,mem = get_usage_statistics(Configuration,current_run.pid)
                    file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
        if tmp[0] == 'L':
            if int(tmp[1]) != currentloop:
                print(f"\r RUN_TIRIFIC: {float(tmp[1])/float(max_loop)*100.:.1f} % Completed", end =" ",flush = True)
            currentloop  = int(tmp[1])
            if max_loop == 0:
                max_loop = int(tmp[2])
            try:
                Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
            except:
                #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                pass
        if tmp[0].strip() == 'Finished':
            break
        if tmp[0].strip() == 'Abort':
            break
    print(f'\n')
    if Configuration['TIMING']:
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write("# Finished this run \n")
            CPU,mem = get_usage_statistics(Configuration,current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb\n")
    print(f"RUN_TIRIFIC: Finished the current tirific run.")
    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    time.sleep(1.0)
    wait_counter = 0
    while not os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.fits") and wait_counter < 1000.:
        print(f"\r Waiting ", end = "", flush = True)
        time.sleep(0.5)
        wait_counter += 1

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

def sofia(Configuration, Fits_Files,sofia_to_use = 'sofia', debug = False):
    if debug:
        print_log(f'''RUN_SOFIA: starting sofia run from the template.
''',Configuration['OUTPUTLOG'],debug = debug)
    sofia_template = rf.sofia_template(debug=debug)
    if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
        os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')
    os.chdir(Configuration['FITTING_DIR'])
    threshold = 5.
    counter = 3
    sofia_template['input.data'] = Fits_Files['FITTING_CUBE']
    spatial_kernels = [0,int(round(Configuration['BEAM_IN_PIXELS'][0])),int(round(Configuration['BEAM_IN_PIXELS'][0]))*2]
    if np.sum(Configuration['NAXES'][:2])/(2.*int(round(Configuration['BEAM_IN_PIXELS'][0])))  > 30:
        spatial_kernels.append(int(round(Configuration['BEAM_IN_PIXELS'][0]))*3)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube is more than 30 beams across.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = False)

    velocity_kernels = [0,3,6,12]
    if Configuration['NAXES'][2] > 52:
        velocity_kernels.append(16)
        Configuration['VEL_SMOOTH_EXTENDED'] = True
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = False)
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration,debug=debug)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par')
        print_log("RUN_SOFIA: Running SoFiA. \n",Configuration['OUTPUTLOG'],screen = True,debug = False)
        # Check which sofia to start
        sfrun = subprocess.Popen([sofia_to_use,'sofia_input.par'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sofia_run, sofia_warnings_are_annoying = sfrun.communicate()
        print_log(sofia_run.decode("utf-8"), Configuration['OUTPUTLOG'],debug = False)

        if sfrun.returncode == 8:
            if threshold > 3.:
                log_statement = f'''RUN_SOFIA: We did not find a source at a threshold of {threshold}
{"":8s} RUN_SOFIA: Lowering the threshold and trying again."
'''
                print_log(log_statement,Configuration['OUTPUTLOG'],debug = False)
                threshold -= 1
            else:
                clean_after_sofia(Configuration)
                log_statement = f'''RUN_SOFIA: We did not find a source above a threshold of {threshold}.
{"":8s}RUN_SOFIA: We cannot lower the threshold lower as the risk of fitting noise becomes too high.
{"":8s}Continuing to the next galaxy.
'''
                print_log(log_statement,Configuration['OUTPUTLOG'],debug = False)
                raise SofiaFaintError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'],debug = False)
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration,debug=debug)

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
