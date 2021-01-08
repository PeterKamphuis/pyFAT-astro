# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to run external programs

class BadSourceError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class SofiaFaintError(Exception):
    pass

from pyFAT.Support.support_functions import print_log, convert_type,set_limits,sbr_limits,rename_fit_products,\
                              set_ring_size,calc_rings,get_usage_statistics,get_inner_fix,convertskyangle,\
                              finish_current_run, rotateImage, remove_inhomogeneities,the_bends,get_from_template,set_format, \
                              wait_for_tirific,set_rings
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
def central_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr,Initial_Parameters, debug = False):
    #This function will run while it does not return 1, this means when it is called it is not in acceptance
    if debug:
        print_log('''CENTRAL_CONVERGE: Starting.
''',Configuration['OUTPUTLOG'],debug = debug)
    #First we actually run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = 'run_cc', fit_stage = 'Centre_Convergence',debug=debug)
    #accepted,current_run = run_parameterized_tirific(Configuration,Tirific_Template,current_run,loops = 10, stage = 'run_cc', fit_stage = 'Centre_Convergence',debug=debug)

    accepted = check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = 'Centre_Convergence',debug=debug)
    set_cflux(Configuration,Tirific_Template,debug = debug)
    check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
    if accepted:
        Configuration['CC_ACCEPTED'] = True
    else:
        print_log('''CENTRAL_CONVERGE: Tirific ran the maximum amount of loops which means the fit is not accepted and we retry.
''',Configuration['OUTPUTLOG'],debug = debug)
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',hdr,min_error = float(Configuration['CHANNEL_WIDTH']), debug=debug)
        smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr,min_error= np.max([float(Tirific_Template['CFLUX']),float(Tirific_Template['CFLUX_2'])]) ,debug=debug)
        #after smoothing the sbr we should check it
        check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
        no_declining_vrot(Configuration,Tirific_Template,debug = debug)
        for key in ['XPOS','YPOS','VSYS','PA','INCL']:
            Initial_Parameters[key][0] = rf.load_template(Tirific_Template,Variables = [key])[0][0]
        Initial_Parameters['VROT'] =  [np.nanmax(smoothed_vrot),set_limits(np.nanstd(smoothed_vrot[1:]),np.nanmax(smoothed_vrot)-np.nanmin(smoothed_vrot[1:]),np.nanmax(smoothed_vrot))]

        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = 'run_cc',\
                               initial_estimates = Initial_Parameters, debug = debug)
        Tirific_Template['SIZE'] = f"{float(Tirific_Template['SIZE']) * 2.}"
        wf.tirific(Configuration,Tirific_Template,name = 'Centre_Convergence_In.def',debug=debug)

    return current_run


def fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,hdr, stage = 'initial',fit_stage='Undefined_Stage', debug = False):
    if debug:
        print_log(f'''FIT_SMOOTHED_CHECK: Starting stage {stage} and fit_stage {fit_stage}.
''',Configuration['OUTPUTLOG'],debug = debug)
    #if we have only a few rings we only smooth. else we fit a polynomial to the RC and smooth the SBR
    smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr, min_error= np.max([float(Tirific_Template['CFLUX']),float(Tirific_Template['CFLUX_2'])]),debug = debug)
    fix_sbr(Configuration,Tirific_Template,hdr,debug = debug)
    if stage == 'after_cc':
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',hdr,min_error = float(hdr['CDELT3']/1000.),debug = debug)
    else:
        smoothed_vrot = regularise_profile(Configuration,Tirific_Template,'VROT',hdr,min_error = float(hdr['CDELT3']/1000.),debug = debug)
    no_declining_vrot(Configuration, Tirific_Template, debug = debug)
    if stage in ['after_ec', 'after_os']:
        min_error = []
        pars_to_smooth = []
        not_to_smooth = []
        fixed_errors = []
        if not Configuration['FIX_INCLINATION'][0]:
            pars_to_smooth.append('INCL')
            min_error.append(2.*Configuration['LIMIT_MODIFIER'])
        else:
            not_to_smooth.append('INCL')
            fixed_errors.append(2.*Configuration['LIMIT_MODIFIER'])
        if not Configuration['FIX_Z0'][0]:
            pars_to_smooth.append('Z0')
            min_error.append(convertskyangle(0.1,Configuration['DISTANCE'],physical= True))
        else:
            not_to_smooth.append('Z0')
            fixed_errors.append(convertskyangle(0.1,Configuration['DISTANCE'],physical= True))
        if not Configuration['FIX_PA'][0]:
            pars_to_smooth.append('PA')
            min_error.append(1.)
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
            smoothed = regularise_profile(Configuration,Tirific_Template,key,hdr,min_error = min_err,debug = debug)

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
        xpos,ypos,vsys,pa,incl = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = ['XPOS','YPOS','VSYS','PA','INCL'])
        parameters = {'VSYS': [vsys[0], Configuration['CHANNEL_WIDTH']], \
                      'XPOS': [xpos[0],abs(hdr['CDELT1'])],
                      'YPOS': [ypos[0],abs(hdr['CDELT2'])],
                      'INCL':[incl[0],10.],
                      'PA':  [pa[0],5.],
                      'VROT': [np.max(smoothed_vrot),np.max(smoothed_vrot)*0.1],  }

        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = stage,\
                          initial_estimates = parameters, debug = debug)
    os.system(f"cp {Configuration['FITTING_DIR']}{fit_stage}_In.def {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Last_Unsmoothed_Input.def")
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_stage}_In.def',debug=debug)
    accepted,current_run = run_tirific(Configuration,current_run, stage = stage, fit_stage=fit_stage,debug=debug)


    return current_run

def check_inclination(Configuration,Tirific_Template,Fits_Files, fit_stage = 'Unstaged', debug =False):
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
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Check_Template,debug = debug)
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
    noisemap = np.sqrt(chan_map[0].data)*model_mom0[0].header['FATNOISE']/np.nanmax(model_mom0[0].data)
    model_mom0[0].data = model_mom0[0].data/np.nanmax(model_mom0[0].data)
    for incl in incl_to_check:
        #print(f'We are doing this inclination {incl}')
        vrot = [x/np.sin(np.radians(incl)) for x in vobs]
        for key in ['INCL','INCL_2']:
            Check_Template[key]= f"{incl:.2f}"
        for key in ['VROT','VROT_2']:
            Check_Template[key]= f"{' '.join([f'{x:.2f}' for x in vrot])}"
        wf.tirific(Configuration,Check_Template,name = f'{tmp_stage}_In.def',debug=debug)
        accepted,incl_run = run_tirific(Configuration,incl_run, stage = 'incl_check', fit_stage=tmp_stage,debug=debug)
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
    print(mom_chi)
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
        wf.tirific(Configuration,Tirific_Template,name = f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",debug=debug)


def check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = 'Undefined_Stage', debug = False):
    if debug:
        print_log(f'''CHECK_CENTRAL_CONVERGE: Starting stage {fit_stage}.
''',Configuration['OUTPUTLOG'],debug = debug)
    new_xpos,new_ypos,new_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}_In.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    if Configuration['OUTER_RINGS_DOUBLED']:
        shift_beam_frac =Configuration['SIZE_IN_BEAMS']*0.05
    else:
        shift_beam_frac =0.15
    ra_lim = set_limits(shift_beam_frac*hdr['BMAJ'],np.max([abs(hdr['CDELT1']),1./3600.]),hdr['BMAJ'],debug = debug )
    dec_lim =set_limits(shift_beam_frac*hdr['BMAJ'],np.max([abs(hdr['CDELT2']),1./3600.]),hdr['BMAJ'],debug = debug )
    sys_lim = set_limits(0.5*hdr['CDELT3']/1000.,2.5, 2.*hdr['CDELT3']/1000.,debug = debug )
    if abs(new_xpos[0] - old_xpos[0]) > ra_lim or \
       abs(new_ypos[0] - old_ypos[0]) > dec_lim or \
       abs(new_vsys[0] - old_vsys[0]) > sys_lim:
        if Configuration['SIZE_IN_BEAMS'] < 20:
            apply_limit = 2*hdr['BMAJ']
        else:
            apply_limit = hdr['BMAJ']*Configuration['SIZE_IN_BEAMS']*0.2
        if abs(new_xpos[0] - old_xpos[0]) > apply_limit or\
           abs(new_ypos[0] - old_ypos[0]) > apply_limit:
           print_log(f'''CHECK_CONVERGENCE: The center shifted more than {apply_limit/hdr['BMAJ']} FWHM.
{"":8s}CHECK_CONVERGENCE: Not applying this shift
''', Configuration['OUTPUTLOG'],debug = debug)
           write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template, Variables = ['VROT',
                         'Z0', 'SBR', 'INCL','PA','SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','SDIS_2'],debug=debug)
           return 0
        else:
            try:
                old_xpos,old_ypos,old_vsys = rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev2.def",Variables = ['XPOS','YPOS','VSYS'],debug=debug)
                if abs(new_xpos[0] - old_xpos[0]) < ra_lim and \
                   abs(new_ypos[0] - old_ypos[0]) < dec_lim and \
                   abs(new_vsys[0] - old_vsys[0]) < sys_lim:
                   print_log(f'''CHECK_CONVERGENCE: The center shifted back to the old position. Moving on to the next stage.
''', Configuration['OUTPUTLOG'],debug = debug)
                   write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template,debug=debug)
                   return accepted
            except:
                pass
            print_log(f'''CHECK_CONVERGENCE: The center shifted too much trying again with new center.
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
            write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template,debug=debug)
            return 0
    else:
        print_log(f'''CHECK_CONVERGENCE: The center is accepted. The shift is:
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'],debug = debug)
        write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template,debug=debug)
        return accepted

def check_source(Configuration, Fits_Files, Catalogue, header, debug = False):
    if debug:
        print_log(f'''CHECK_SOURCE: Starting.
''',Configuration['OUTPUTLOG'],debug = debug)

    name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa, \
        w50,err_f_sum, err_x,err_y,err_z= rf.sofia_catalogue(Configuration,header=header,debug=debug)

    x_min,x_max,y_min,y_max,z_min,z_max = convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int',debug=debug)
    x,y,z,ra,dec,v_app,f_sum,kin_pa,f_sum_err , err_x,err_y,err_z= convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa,err_f_sum, err_x,err_y,err_z],debug=debug)
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

    new_box = cut_cubes(Configuration, Fits_Files, galaxy_box, header,debug=debug)
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
    bmmaj_in_pixels = header['BMAJ']/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])
    #There is a factor of two missing here but this is necessary otherwise the maxima are far to small
    Configuration['MAX_SIZE_IN_BEAMS'] = int(round(np.sqrt(((x_max-x_min)/2.)**2+((y_max-y_min)/2.)**2) \
                /(bmmaj_in_pixels)+5.))

    pa, inclination, SBR_initial, maj_extent,x_new,y_new,VROT_initial = rf.guess_orientation(Configuration,Fits_Files, center = [x,y],debug=debug)

    if x_new != x or y_new != y:
        x=x_new
        y=y_new
        ra,dec,v_app = cube_wcs.wcs_pix2world(x,y,z,0.)
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
    Central_Flux = np.mean(data[int(round(z-1)):int(round(z+1)),int(round(y-bmmaj_in_pixels/2.)):int(round(y+bmmaj_in_pixels/2.)),int(round(x-bmmaj_in_pixels/2.)):int(round(x+bmmaj_in_pixels/2.))])
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


    if abs(pa[0]-kin_pa) < 25:
        pa[0] = (pa[0]/pa[1]+kin_pa/10.)/(1./pa[1]+1./10.)
    else:
        if np.isfinite(kin_pa):
            pa[0] = kin_pa

    # Size of the galaxy in beams
    Configuration['SIZE_IN_BEAMS'] = set_limits(maj_extent/(header['BMAJ']),1.0,Configuration['MAX_SIZE_IN_BEAMS'])
    if Configuration['SIZE_IN_BEAMS'] <= Configuration['TOO_SMALL_GALAXY']:
        print_log(f'''CHECK_SOURCE: This galaxy has an estimated size of  {2*Configuration['SIZE_IN_BEAMS']} beams in diameter.
{'':8s}This is not large enough too fit. We will exit this fit.
''',Configuration['OUTPUTLOG'],debug = debug)
        raise BadSourceError('We found an initial negative total flux.')
    set_ring_size(Configuration, debug = debug)
    new_radii = set_rings(Configuration, debug = debug)
    old_radii = np.linspace(0.,Configuration['BMMAJ']*(len(SBR_initial)-1),len(SBR_initial))
    new_radii = set_rings(Configuration, debug = debug)
    SBR_initial = np.interp(new_radii,old_radii, SBR_initial)
    old_radii = np.linspace(0.,Configuration['BMMAJ']*(len(VROT_initial)-1),len(VROT_initial))
    VROT_initial = np.interp(new_radii,old_radii, VROT_initial)
    Configuration['OUTER_SLOPE_START'] = Configuration['NO_RINGS']
    if debug:
        print_log(f'''CHECK_SOURCE: Interpolating the SBR and VROT estimates to these radi.
{'':8s} {new_radii}
{'':8s} We got SBR = {SBR_initial}, VROT = {VROT_initial}
''',Configuration['OUTPUTLOG'],debug = debug)


    # The extent is fairly well determined and the maximum should be no more than +3 beams and a minimum no less than 4
    Configuration['MAX_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']+2.,1.0,Configuration['MAX_SIZE_IN_BEAMS'])
    Configuration['MIN_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']-4.,Configuration['TOO_SMALL_GALAXY'],Configuration['SIZE_IN_BEAMS'])
    Configuration['NO_RINGS'] = calc_rings(Configuration,debug=debug)
    Configuration['LAST_RELIABLE_RINGS'] = [Configuration['NO_RINGS'],Configuration['NO_RINGS']]
    print_log(f'''CHECK_SOURCE: From the original Configuration and SoFiA we find:
{"":8s}CHECK_SOURCE: The maximum diameter we will allow  is  {2.*Configuration['MAX_SIZE_IN_BEAMS']} beams.
{"":8s}CHECK_SOURCE: We start with a diameter of {2*Configuration['SIZE_IN_BEAMS']} beams in the model.
{"":8s}CHECK_SOURCE: SoFiA found a PA of {kin_pa:.2f} and we use a PA = {pa[0]:.2f} +/- {pa[1]:.2f}
{"":8s}CHECK_SOURCE: We start with an inclination of {inclination[0]:.2f} +/- {inclination[1]:.2f}{"":8s}
{"":8s}CHECK_SOURCE: SoFiA found a W50 of {w50/1000.:.2f} km/s
{"":8s}CHECK_SOURCE: We will use {2.*Configuration['NO_RINGS']} for the model with a ring size of {Configuration['RING_SIZE']}.
''',Configuration['OUTPUTLOG'],debug = debug)



    if Configuration['HANNING']:
        vres = header['CDELT3']*2.
    else:
        vres = header['CDELT3']

    #if inclination[0] < 40:
    #    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0]+5.)))
    #else:
    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0])))
    max_vrot_dev=set_limits((VELboun[1]-VELboun[0])/4./np.sin(np.radians(inclination[0])),4.*vres,20*vres,debug=debug)


    #Write the info to the Basic info File
    wf.basicinfo(Configuration,initialize = True,
              RA=[ra,abs(err_x*header['CDELT1'])],
              DEC=[dec,abs(err_y*header['CDELT2'])],
              VSYS =[v_app,abs(err_z*header['CDELT3'])],
              PA=pa, Inclination = inclination, Max_Vrot = [max_vrot,max_vrot_dev], Tot_Flux = [f_sum,f_sum_err], V_mask = [VELboun[1]-VELboun[0],vres],
              Distance = Configuration['DISTANCE'] , DHI = maj_extent*3600.,debug=debug)


    # extract a PV-Diagram
    if Configuration['START_POINT'] < 3 or not os.path.exists(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits"):
        PV =extract_pv(Cube, pa[0], center=[ra,dec,v_app], convert = 1000.,
                       finalsize = [int(round(maj_extent/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])*1.25+header['NAXIS1']*0.2)),
                                    int(round(z_max-z_min)+10.)],debug=debug)
        if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
                os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')

        fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits",PV[0].data,PV[0].header)
    Cube.close()
    Initial_Parameters = {}
    Initial_Parameters['XPOS'] = [ra,set_limits(abs(err_x*header['CDELT1']),0.1*header['BMAJ'],3.*header['BMAJ'] )]
    Initial_Parameters['YPOS'] = [dec,set_limits(abs(err_y*header['CDELT2']),0.1*header['BMAJ'],3.*header['BMAJ'] )]
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
    if debug:
        print_log(f'''EXTENT_CONVERGENCE: Starting with loop {Configuration['EC_LOOPS']} out of maximum {allowed_loops}.
''',Configuration['OUTPUTLOG'],debug = debug)
    fit_stage = 'Extent_Convergence'
    stage = 'run_ec'
    accepted,current_run = run_tirific(Configuration,current_run,stage = stage, fit_stage = fit_stage, debug= debug)
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def" , Tirific_Template, debug = debug)
    accepted_size = check_size(Configuration,Tirific_Template,hdr,  fit_stage = fit_stage, stage = stage,current_run = current_run, debug=debug,Fits_Files=Fits_Files)
    check_sbr(Configuration,Tirific_Template,hdr,stage = 'run_ec',debug=debug)
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
        min_errors = [2.*Configuration['LIMIT_MODIFIER'],1.,hdr['CDELT3']/(2000.*Configuration['LIMIT_MODIFIER']), \
                        convertskyangle(0.1,Configuration['DISTANCE'],physical= True)/Configuration['LIMIT_MODIFIER'],\
                        hdr['CDELT3']/(1000.*Configuration['LIMIT_MODIFIER'])]
        for j,key in enumerate(keys_to_smooth):
            smoothed = smooth_profile(Configuration,Tirific_Template,key,hdr,debug=debug,min_error=min_errors[j])

        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = 'run_ec',\
                             debug = debug)

        #After smoothing the SBR we should check it again?
        check_sbr(Configuration,Tirific_Template,hdr,stage = 'run_ec',debug = debug)
        wf.tirific(Configuration,Tirific_Template,name = 'Extent_Convergence_In.def',debug = debug)


    return current_run


def make_full_resolution(Configuration,Tirific_Template,Fits_Files,fit_stage = 'Undefined_Stage', current_run = 'Not zed', debug = False):
    if debug:
        print_log(f'''MAKE_FULL_RESOLUTION: creating full resolution for stage {fit_stage}.
''',Configuration['OUTPUTLOG'],debug = debug)
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def", Tirific_Template,debug = debug)
    Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Tirific_Template['LOOPS'] = "0"
    Tirific_Template['INIMODE'] = "0"
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_stage}_In.def')
    finish_current_run(Configuration,current_run,debug= debug)
    accepted,current_run = run_tirific(Configuration,'Not zed', stage = 'full_res', fit_stage=fit_stage,debug = debug)
    finish_current_run(Configuration,current_run,debug= debug)



def one_step_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr, debug = False,allowed_loops = 10.):
    if debug:
        print_log(f'''ONE_STEP_CONVERGENCE: Starting with loop {Configuration['OS_LOOPS']} out of maximum {allowed_loops}.
''',Configuration['OUTPUTLOG'],debug = True)
    fit_stage = 'One_Step_Convergence'
    stage = 'run_os'
    #First we run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = stage, fit_stage = fit_stage, debug= debug)
    #Then we load the produced output into our template
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def" , Tirific_Template, debug = debug)
    #Check that the centre does not deviate too much
    accepted_central = check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = fit_stage,debug=debug)
    # Check whether we have the correct sizes,
    accepted_size = check_size(Configuration,Tirific_Template,hdr, fit_stage = fit_stage, stage = stage, current_run = current_run, debug=debug,Fits_Files=Fits_Files)
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
        fix_sbr(Configuration,Tirific_Template,hdr,debug = debug)    # Then we determine the inner rings that should remain fixed

        Configuration['INNER_FIX'] = get_inner_fix(Configuration, Tirific_Template,debug=debug)
        set_cflux(Configuration,Tirific_Template,debug = debug)
        keys_to_smooth =['INCL','PA','SDIS','Z0','VROT','SBR']
        min_errors = [2.*Configuration['LIMIT_MODIFIER'],1.,hdr['CDELT3']/(2000.*Configuration['LIMIT_MODIFIER']), \
                        convertskyangle(0.1,Configuration['DISTANCE'],physical= True)/Configuration['LIMIT_MODIFIER'],\
                        hdr['CDELT3']/(1000.*Configuration['LIMIT_MODIFIER']),np.max([float(Tirific_Template['CFLUX']),\
                        float(Tirific_Template['CFLUX_2'])])]
        for j,key in enumerate(keys_to_smooth):
            #Smoothing the profile also fixes it
            smoothed = smooth_profile(Configuration,Tirific_Template,key,hdr,debug=debug,min_error=min_errors[j])

        set_fitting_parameters(Configuration, Tirific_Template, hdr = hdr,stage = 'run_os',\
                             debug = debug)

        wf.tirific(Configuration,Tirific_Template,name = 'One_Step_Convergence_In.def',debug = debug)


    return current_run

def run_tirific(Configuration, current_run, stage = 'initial',fit_stage = 'Undefined_Stage', debug = False):
    if debug:
        print_log(f'''RUN_TIRIFIC: Starting a new run in stage {stage} and fit_stage {fit_stage}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    # First move the previous fits
    rename_fit_products(Configuration,fit_stage = fit_stage, stage=stage, debug = debug)
    # Then if already running change restart file
    if Configuration['TIRIFIC_RUNNING']:
        print_log(f'''RUN_TIRIFIC: We are using an initialized tirific in {Configuration['FITTING_DIR']}
''',Configuration['OUTPUTLOG'], screen = True)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
    else:
        print_log(f'''RUN_TIRIFIC: We are starting a new TiRiFiC in {Configuration['FITTING_DIR']}
''',Configuration['OUTPUTLOG'], screen = True)
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
                    CPU,mem = get_usage_statistics(current_run.pid)
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
            CPU,mem = get_usage_statistics(current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb\n")
    print(f"RUN_TIRIFIC: Finished the current tirific run.")
    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    time.sleep(1.0)
    wait_counter = 0
    while not os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.fits") and wait_counter < 1000.:
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

def run_parameterized_tirific(Configuration, Tirific_Template, current_run, loops = 5, stage = 'initial',fit_stage = 'Undefined_Stage', debug = False):
    if debug:
        print_log(f'''RUN_PARAMETIRIZED_TIRIFIC: Starting a new paremetirized run in stage {stage} and fit_stage {fit_stage}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    # First scale the inclination and the vrot

    Tirific_Template['LOOPS']  = 1
    wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = debug)
    if Configuration['TIRIFIC_RUNNING']:
        print_log('''RUN_TIRIFIC: We are using an initialized tirific
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
    else:
        print_log('''RUN_TIRIFIC: We are starting a new TiRiFiC
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'w') as file:
            file.write("Initialized a new run")
        #current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_stage}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
        #                            stderr = subprocess.PIPE,cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_stage}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
                                    cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        Configuration['TIRIFIC_RUNNING'] = True
        Configuration['TIRIFIC_PID'] = current_run.pid
    #rename_fit_products(Configuration,fit_stage = fit_stage, stage=stage, debug = debug)
    # Then if already running change restart file
    sys.stdout.flush()
    currentloop =1
    max_loop = 0
    counter = 0
    Chi_Current = 0.
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
            try:
                Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
            except:
                #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                pass
        if tmp[0].strip() == 'Finished':
            Chi_Current = (float(tmp[tmp.index('C')+1]))
            break
        if tmp[0].strip() == 'Abort':
            break
    time.sleep(1.)
    # Now do the actual fitting loops
    for i in range(loops):
        print(f"RUN_Parameterized TIRIFIC: Starting loop {i+1}")
        resid= Tirific_Template['RESTARTID']
        resname =  Tirific_Template['RESTARTNAME']
        Tirific_Template = rf.tirific_template(f"{Configuration['FITTING_DIR']}{Tirific_Template['TIRDEF']}")
        if i == 0.:
            prev_bends = [Tirific_Template['INCL'],Tirific_Template['INCL_2'],Tirific_Template['PA'],Tirific_Template['PA_2']]
        Tirific_Template['RESTARTID'] = resid
        Tirific_Template['RESTARTNAME'] = resname
    # whatever has been fitted we will now attempt to modify by fitting a warped __version__
        current_run,bends = fit_warp(Configuration,Tirific_Template,current_run, fit_stage=fit_stage,debug=debug)
        Tirific_Template['LOOPS'] = 1
        wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = False)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
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
                try:
                    Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
                except:
                    #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                    pass
            if tmp[0].strip() == 'Finished':
                Chi = float(tmp[tmp.index('C')+1])
                break
            if tmp[0].strip() == 'Abort':
                break
        comp_loop = i+1
        print(f"Chi = {Chi}, Chi_Current = {Chi_Current}, loop = {comp_loop}" )
        if Chi < Chi_Current:
            Chi_Current = Chi
            prev_bends = bends
        else:
            break
    restore_bends(Configuration,Tirific_Template,prev_bends,debug=debug)
    Tirific_Template['LOOPS'] = 0
    wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = False)
    with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
        file.write("Restarting from previous run")
    if Configuration['TIMING']:
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write("# Finished this run \n")
            CPU,mem = get_usage_statistics(current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb\n")
    print(f"RUN_TIRIFIC: Finished the current tirific run.")
    time.sleep(1.0)
    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    wait_counter = 0
    while not os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def") and wait_counter < 1000.:
        time.sleep(0.1)
        wait_counter += 1

    if comp_loop != loops:
        return 1,current_run
    else:
        return 0,current_run

run_parameterized_tirific.__doc__= '''

; NAME:
;       tirific(Configuration)
;
; PURPOSE:
;       Fit a parametirized functions in tirific
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
def restore_bends(Configuration,Tirific_Template,bends,debug=False):
    if debug:
        print_log(f'''RESTORE_BENDS: setting the best fit bends to be output.
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    bend_keys = ['INCL','INCL_2','PA','PA_2']
    for bend,key in zip(bends,bend_keys):
        Tirific_Template[key] = bend

def fit_warp(Configuration,Tirific_Template,current_run, fit_stage = 'None',debug=False):
    if debug:
        print_log(f'''FIT_WARP: Starting a new paremeterized run in fit_stage {fit_stage}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)

    Tirific_Template['LOOPS'] = 0
    incl_bend = planet_telex(Configuration, Tirific_Template,current_run,'INCL',fit_stage=fit_stage,debug=debug)
    incl_bend_2 = planet_telex(Configuration, Tirific_Template,current_run,'INCL_2',fit_stage=fit_stage,debug=debug)
    pa_bend = planet_telex(Configuration, Tirific_Template,current_run,'PA',max_bend = 50, fit_stage=fit_stage,debug=debug)
    pa_bend_2 = planet_telex(Configuration, Tirific_Template,current_run,'PA_2',max_bend = 50, fit_stage=fit_stage,debug=debug)

    print(f"PA = {Tirific_Template['PA']}" )
    print(f"PA_2 = {Tirific_Template['PA_2']}" )
    print(f"INCL = {Tirific_Template['INCL']}" )
    print(f"INCL_2 = {Tirific_Template['INCL_2']}" )

    return current_run,[Tirific_Template['INCL'],Tirific_Template['INCL_2'],Tirific_Template['PA'],Tirific_Template['PA_2']]



def planet_telex(Configuration, Tirific_Template,current_run, key, black_stars = 0.,max_bend =30 ,step_bend = 5.,fit_stage= 'None', debug=False):
    Chi = []
    thebends = []
    if key[-2:] == '_2':
        SBR_key = 'SBR_2'
        other_key = key[:-2]
        fitkey = other_key
    else:
        SBR_key = 'SBR'
        other_key = f"{key}_2"
        fitkey = key
    initial_profile = np.array(get_from_template(Tirific_Template,['RADI',key,other_key]))
    if black_stars == 0.:
        black_stars = initial_profile[0][-1]/np.linspace(3,1.1,10)
        black_index = np.where(initial_profile[0] <= black_stars[0])[0][-1]

    # The start outermost bent of inclination
    SBRformat = set_format(SBR_key)
    format = set_format(key)
    SBR = np.array(get_from_template(Tirific_Template,[SBR_key])[0],dtype=float)
    new_SBR = []
    for i,x in enumerate(SBR):
        if i < black_index:
            new_SBR.append(x)
        else:
            new_SBR.append(SBR[black_index])

    Tirific_Template[SBR_key] = f"{' '.join([f'{x:{SBRformat}}' for x in new_SBR ])}"
    Tirific_Template['LOOPS'] = 1.
    parameters = {key : [np.mean(initial_profile[1:2]), 10.]}
    set_fitting_parameters(Configuration, Tirific_Template,stage = 'parameterized',\
                            parameters_to_adjust = [fitkey],initial_estimates = parameters, debug = debug)
    Chi = []
    all_profs = []
    total_fits = len(black_stars)*2*max_bend/step_bend
    counter = 0.
    for black_star in black_stars:
        Chis,profiles,counter = fake_plastic_trees(Configuration,Tirific_Template,current_run,\
                                            key,other_key,initial_profile,max_bend,step_bend,\
                                            black_star= black_star,fit_stage= fit_stage, \
                                            counter = counter, total_fits = total_fits, debug=debug)
        Chi.append(Chis)
        all_profs.append(profiles)
    Chi = np.array(Chi,dtype=float)
    all_profs = np.array(all_profs,dtype=float)
    print(f"Chi = {Chi}")
    #print(np.where(Chi == np.min(Chi))[0][0],np.min(Chi))
    high_and_dry = np.where(Chi == np.min(Chi))
    print(f"This is the bend value we find {high_and_dry[0]}, {high_and_dry[1]}for {key}")
    new_profile = all_profs[high_and_dry][0][0]
    print(f" This is the profile {new_profile}")
    #Restore the brightness profiles
    Tirific_Template[SBR_key] = f"{' '.join([f'{x:{SBRformat}}' for x in SBR])}"
    diff = float(initial_profile[2][0]-new_profile[0])
    Tirific_Template[key] =  f"{' '.join([f'{x:{format}}' for x in new_profile])}"
    Tirific_Template[other_key] =  f"{' '.join([f'{x-diff:{format}}' for x in initial_profile[2]])}"

    return new_profile


def fake_plastic_trees(Configuration,Tirific_Template,current_run,key,other_key,\
                        initial_profile,max_bend,step_bend,black_star = 0., \
                        counter = 0, total_fits = 1.,fit_stage= 'unknown',debug=False):
    prev_bend = -1*max_bend
    Chi = []
    profiles = []
    format = set_format(key)
    while prev_bend < max_bend:
        #first side 1
        bend = float(prev_bend+step_bend)
        new_profile = the_bends(initial_profile[1],initial_profile[0],black_star,bend)
        diff = float(initial_profile[2][0]-new_profile[0])
        Tirific_Template[key] =  f"{' '.join([f'{x:{format}}' for x in new_profile])}"
        Tirific_Template[other_key] =  f"{' '.join([f'{x-diff:{format}}' for x in initial_profile[2]])}"
        wf.tirific(Configuration,Tirific_Template,name =  f"{fit_stage}_In.def",debug = False)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write(f"Restarting for bend = {bend} and RI = {Tirific_Template['RESTARTID']} \n")
        Chi_here = wait_for_tirific(Configuration,Tirific_Template['RESTARTID'],current_run, counter= counter, total_fits= total_fits)
        if ~np.isnan(Chi_here):
            Chi.append(Chi_here)
            profiles.append(rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = [key]))
            prev_bend = copy.deepcopy(bend)
        else:
            prev_bend = copy.deepcopy(bend-step_bend)
        counter += 1
    return Chi,profiles,counter

def sofia(Configuration, Fits_Files, hdr, debug = False):
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
        Configuration['VEL_SMOOTH_EXTENDED'] = True
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration,debug=debug)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par',debug=debug)
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
                raise SofiaFaintError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'],debug = debug)
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration,debug=debug)

    os.chdir(Configuration['START_DIR'])
# to ensure compatible units and calculations with th models we make the maps ourselves
    mask = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}")
    del mask[0].header['C*3']
    mask[0].data[mask[0].data> 0.5] = 1.
    chan_map = np.array(np.nansum(mask[0].data,axis=0),dtype =float)
    mask[0].header['BITPIX'] = -32
    fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_chan.fits",chan_map,header=mask[0].header,overwrite = True)
    mask.close()
    make_moments(filename = f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}",\
                 basename = f"{Configuration['BASE_NAME']}", directory = f"{Configuration['FITTING_DIR']}/Sofia_Output/",\
                 mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}",vel_unit = 'm/s',debug=debug)

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
