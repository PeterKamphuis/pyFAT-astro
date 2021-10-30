# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to run external programs



from pyFAT_astro.Support.support_functions import print_log, convert_type,set_limits,rename_fit_products,\
                              set_ring_size,calc_rings,get_inner_fix,convertskyangle,\
                              finish_current_run, remove_inhomogeneities,get_from_template,set_format, \
                              set_rings, convertRADEC,sbr_limits, create_directory,get_system_string,\
                              get_fit_groups,run_tirific,update_statistic
from pyFAT_astro.Support.clean_functions import clean_before_sofia,clean_after_sofia
from pyFAT_astro.Support.fits_functions import cut_cubes,extract_pv,make_moments
from pyFAT_astro.Support.read_functions import load_template,tirific_template

from pyFAT_astro.Support.modify_template import write_new_to_template,smooth_profile,set_cflux,fix_sbr, \
                                          regularise_profile,set_fitting_parameters,check_size, \
                                          no_declining_vrot,set_errors,get_warp_slope,check_angles,write_center
from pyFAT_astro.Support.constants import H_0
from pyFAT_astro.Support.fat_errors import SofiaFaintError,BadConfigurationError,\
                                              InclinationRunError,SofiaRunError,BadSourceError
from pyFAT_astro.Support.tirshaker import tirshaker

from astropy.wcs import WCS
from astropy.io import fits

import pyFAT_astro.Support.read_functions as rf
import pyFAT_astro.Support.write_functions as wf
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


def check_central_convergence(Configuration,Tirific_Template, fit_type = 'Undefined', debug = False):
    #The new values are already loaded into the Tirific_Template so if we do accept we have to reset the values
    if debug:
        print_log(f'''CHECK_CENTRAL_CONVERGE: Starting stage {fit_type}.
''',Configuration['OUTPUTLOG'],debug = True)
    update_statistic(Configuration, message= "Starting the central convergence run", debug=debug)

    new_xpos,new_ypos,new_vsys = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    old_xpos,old_ypos,old_vsys = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}_In.def",Variables = ['XPOS','YPOS','VSYS'],debug = debug)
    if Configuration['OUTER_RINGS_DOUBLED']:
        shift_beam_frac =Configuration['SIZE_IN_BEAMS']*0.05
    else:
        shift_beam_frac =0.15
    ra_lim = set_limits(shift_beam_frac*Configuration['BEAM'][0]/3600.,np.max([Configuration['PIXEL_SIZE'],1./3600.]),Configuration['BEAM'][0]/3600.,debug = debug )
    dec_lim =set_limits(shift_beam_frac*Configuration['BEAM'][0]/3600.,np.max([Configuration['PIXEL_SIZE'],1./3600.]),Configuration['BEAM'][0]/3600.,debug = debug )
    sys_lim = set_limits(0.5*Configuration['CHANNEL_WIDTH'],2.5, 2.*Configuration['CHANNEL_WIDTH'],debug = debug )
    if abs(new_xpos[0] - old_xpos[0]) > ra_lim or \
       abs(new_ypos[0] - old_ypos[0]) > dec_lim or \
       abs(new_vsys[0] - old_vsys[0]) > sys_lim:
        if Configuration['SIZE_IN_BEAMS'] < 20:
            apply_limit = 2.*Configuration['BEAM'][0]/3600.
        else:
            apply_limit = Configuration['BEAM'][0]/3600.*Configuration['SIZE_IN_BEAMS']*0.2
        if abs(new_xpos[0] - old_xpos[0]) > apply_limit or\
           abs(new_ypos[0] - old_ypos[0]) > apply_limit:
           print_log(f'''CHECK_CONVERGENCE: The center shifted more than {apply_limit/(Configuration['BEAM'][0]/3600.)} FWHM.
{"":8s}CHECK_CONVERGENCE: Not applying this shift
''', Configuration['OUTPUTLOG'])
           write_center(Configuration,Tirific_Template, [old_xpos[0],old_ypos[0],old_vsys[0]])

           #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template, Variables = ['VROT',
            #             'Z0', 'SBR', 'INCL','PA','SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','SDIS_2'],debug=debug)
           return False
        else:
            try:
                old_xpos,old_ypos,old_vsys = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev2.def",Variables = ['XPOS','YPOS','VSYS'],debug=debug)
                if abs(new_xpos[0] - old_xpos[0]) < ra_lim and \
                   abs(new_ypos[0] - old_ypos[0]) < dec_lim and \
                   abs(new_vsys[0] - old_vsys[0]) < sys_lim:
                   print_log(f'''CHECK_CONVERGENCE: The center shifted back to the old position. Moving on to the next stage.
''', Configuration['OUTPUTLOG'])
                   #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
                   return True
            except:
                pass
            print_log(f'''CHECK_CONVERGENCE: The center shifted too much trying again with new center.
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'])
            #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
            return False
    else:
        print_log(f'''CHECK_CONVERGENCE: The center is accepted. The shift is:
{"":8s}The RA has shifted from {old_xpos[0]} to  {new_xpos[0]} which is a difference of {abs(new_xpos[0] - old_xpos[0])} needed = {ra_lim}.
{"":8s}The DEC has shifted from {old_ypos[0]} to  {new_ypos[0]} which is a difference of {abs(new_ypos[0] - old_ypos[0])} needed = {dec_lim}.
{"":8s}The VSYS has shifted from {old_vsys[0]} to  {new_vsys[0]} which is a difference of {abs(new_vsys[0] - old_vsys[0])} needed = {sys_lim}.
''', Configuration['OUTPUTLOG'])
        #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug=debug)
        return True
check_central_convergence.__doc__ =f'''
 NAME:
    check_central_convergence

 PURPOSE:
    Check whether the center coodinates have converged to within the limits.

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting being done.

 OUTPUTS:
    Returns a boolean that is true when the change is within the limits false when not

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_inclination(Configuration,Tirific_Template,Fits_Files, fit_type = 'Undefined', debug =False):
    if debug:
        print_log(f'''CHECK_INCLINATION: estimating whether our inclination estimate is decent.
''',Configuration['OUTPUTLOG'],debug = True)
    update_statistic(Configuration, message= "Starting the the inclination run", debug=debug)

    to_extract=['VROT','INCL','INCL_2','PA','XPOS','YPOS']
    current = rf.load_template(Configuration,Tirific_Template,
                Variables= to_extract)


    inclination = float(current[to_extract.index('INCL')][0])
    if debug:
        print_log(f'''CHECK_INCLINATION: This is the initial inclination {inclination}
''',Configuration['OUTPUTLOG'],debug = True)

    if Configuration['SIZE_IN_BEAMS'] < 5:
        incl_to_check = np.linspace(-15.,+15.,20)
    else:
        max=inclination-10
        min=inclination-50
        incl_to_check = np.linspace(min,max,20)
    # Let's make a directory to put things in
    tmp_stage = 'tmp_incl_check'
    create_directory(tmp_stage,f"{Configuration['FITTING_DIR']}")

    other_run = [Configuration['TIRIFIC_RUNNING'],Configuration['TIRIFIC_PID']]
    try:
        Configuration['TIRIFIC_RUNNING'] = False
        #and a copy of the tirific template
        if debug:
                print_log(f'''CHECK_INCLINATION: python is so stupid
    {'':8s}PA = {Tirific_Template['PA']}
    ''',Configuration['OUTPUTLOG'])
        Check_Template = copy.deepcopy(Tirific_Template)
        #write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Check_Template,debug = debug)
        Check_Template['LOOPS'] = '0'
        Check_Template['INIMODE'] = '0'
        Check_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
        current_cwd = os.getcwd()
        short_log = Configuration['LOG_DIRECTORY'].replace(current_cwd,'.')
        Check_Template['RESTARTNAME'] = f"{short_log}restart_{tmp_stage}.txt"
        #Check_Template['RESTARTNAME'] = get_system_string(f"{Configuration['LOG_DIRECTORY']}restart_{tmp_stage}.txt")
        #Check_Template['RESTARTNAME'] = f"./Logs/06-09-2021/restart_tmp_incl_check.txt"

        #These are small galaxies make sure the VARINDX is not meesing things up
        Check_Template['VARINDX'] = ''

        out_keys = ['LOGNAME','OUTSET','TIRDEF']
        out_extensions = ['log','fits','def']
        incl_run= 'Not Initialized'
        for i,key in enumerate(out_keys):
            Check_Template[key] = f"{tmp_stage}/{tmp_stage}.{out_extensions[i]}"


        vobs = [x*np.sin(np.radians(np.mean([float(y),float(z)]))) for x,y,z in \
                zip(current[to_extract.index('VROT')][:],current[to_extract.index('INCL')][:],current[to_extract.index('INCL_2')][:])]
        if debug:
            print_log(f'''CHECK_INCLINATION: These are the values we get as input
    {'':8s}Inclination = {current[to_extract.index('INCL')][:]}, {current[to_extract.index('INCL_2')][:]}
    {'':8s}Vrot = {current[to_extract.index('VROT')][:]}
    {'':8s}Vobs = {vobs}
    {'':8s}PA = {Check_Template['PA']}
    ''',Configuration['OUTPUTLOG'])
        mom_chi = []
        model_mom0 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT0']}")
        #model_mom0 = remove_inhomogeneities(Configuration,model_mom0,inclination=float(current[1][0]), pa = float(current[2][0]), center = [current[3][0],current[4][0]],debug=debug)
        chan_map = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['CHANNEL_MAP']}")
        noisemap = np.sqrt(chan_map[0].data)*Configuration['NOISE']/np.nanmax(model_mom0[0].data)
        max_in_moment = np.nanmax(model_mom0[0].data)
        model_mom0[0].data = model_mom0[0].data
        #/max_in_moment
        for incl in incl_to_check:
            #print(f'We are doing this inclination {incl}')
            vrot = [x/np.sin(np.radians(inclination+incl)) for x in vobs]
            format = set_format('INCL')
            for key in ['INCL','INCL_2']:
                Check_Template[key]= f"{' '.join([f'{x+incl:{format}}' for x in current[to_extract.index(key)][:]])}"
            format = set_format('VROT')
            for key in ['VROT','VROT_2']:
                Check_Template[key]= f"{' '.join([f'{x:{format}}' for x in vrot])}"
            wf.tirific(Configuration,Check_Template,name = f'{tmp_stage}_In.def',debug=True)
            accepted,incl_run = run_tirific(Configuration,incl_run, stage = 'incl_check', fit_type=tmp_stage,debug=False)
            make_moments(Configuration,Fits_Files,fit_type=tmp_stage,\
                         moments = [0], \
                         overwrite = True, vel_unit = 'm/s',debug=False)
            #make_moments(filename = f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}.fits", basename = 'tmp_incl', directory = f"{Configuration['FITTING_DIR']}{tmp_stage}/",\
            #             moments = [0],level = 3.*Configuration['NOISE'], \
            #             overwrite = True, log= Configuration['OUTPUTLOG'], vel_unit = 'm/s',debug=debug)
            incl_mom0 = fits.open(f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}_mom0.fits")
            if debug:
                try:
                    os.remove(f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
                except:
                    pass
                os.rename(f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}_mom0.fits",f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
            chi = np.nansum((model_mom0[0].data[noisemap > 0.]-incl_mom0[0].data[noisemap > 0.])**2/noisemap[noisemap > 0.]**2)
            mom_chi.append(abs(chi))
            incl_mom0.close()
        Check_Template = []
        chan_map.close()
        model_mom0.close()
        finish_current_run(Configuration, incl_run)
    except:
        finish_current_run(Configuration, incl_run)
        Configuration['TIRIFIC_RUNNING'] = other_run[0]
        Configuration['TIRIFIC_PID'] =  other_run[1]
        raise InclinationRunError('Something went wrong while estimating the inclination. This should not happen.')
    Configuration['TIRIFIC_RUNNING'] = other_run[0]
    Configuration['TIRIFIC_PID'] =  other_run[1]
    low= np.where(mom_chi == np.nanmin(mom_chi))[0]
    if low.size > 1:
        low = low[0]
    new_incl = float(incl_to_check[low])
    if debug:
        print_log(f'''CHECK_INCLINATION: This is the new inclination {new_incl} it was {current[1][0]}.
{'':8s}mom_chi = {mom_chi}
{'':8s}low = {low}
''',Configuration['OUTPUTLOG'])


    #exit()
    incl_err = 5.
    #np.mean(np.array([current[to_extract.index('INCL_ERR')],current[to_extract.index('INCL_2_ERR')]],dtype=float))
    if incl_err < 5.:
        incl_err = 5.
    if not current[1][0]-incl_err < new_incl < current[1][0]+incl_err:
        if debug:
            print_log(f'''CHECK_INCLINATION: The inclination has changed, writing to file.
''',Configuration['OUTPUTLOG'])

        vrot = [x/np.sin(np.radians(np.mean([float(y),float(z)])+new_incl)) for x,y,z in \
            zip(vobs,current[to_extract.index('INCL')][:],current[to_extract.index('INCL_2')][:])]
        format = set_format(key)
        for key in ['INCL','INCL_2']:
            Tirific_Template[key]= f"{' '.join([f'{x+new_incl:{format}}' for x in current[to_extract.index(key)][:]])}"
        format = set_format(key)
        for key in ['VROT','VROT_2']:
            Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in vrot])}"
            if debug:
                print_log(f'''CHECK_INCLINATION: This has gone to the template
{'':8s} vrot = {Tirific_Template['VROT']}
{'':8s} incl = {Tirific_Template['INCL']}
{'':8s} incl_2 = {Tirific_Template['INCL_2']}
''',Configuration['OUTPUTLOG'])

        #Tirific_Template = copy.deepcopy(Check_Template)
        wf.tirific(Configuration,Tirific_Template,name = f"{fit_type}/{fit_type}.def",debug=debug)
    #return Tirific_Template
check_inclination.__doc__ =f'''
 NAME:
    check_inclination

 PURPOSE:
    For low inclinations do a better check on the inclinations. Comparing the Chi^2 of the models over different inclinations

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting

 OUTPUTS:
    A new file with modified inclination if it has decreased.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This should efinitely be further explored.
'''

def check_source(Configuration, Fits_Files, debug = False):
    if debug:
        print_log(f'''CHECK_SOURCE: Starting.
''',Configuration['OUTPUTLOG'],debug = True)
    update_statistic(Configuration, message= "Starting the check source run", debug=debug)


    name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa, \
        w50,err_f_sum, err_x,err_y,err_z= rf.sofia_catalogue(Configuration,Fits_Files,debug=debug)

    x_min,x_max,y_min,y_max,z_min,z_max = convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int')
    x,y,z,ra,dec,v_app,f_sum,kin_pa,f_sum_err , err_x,err_y,err_z= convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa,err_f_sum, err_x,err_y,err_z])
    #How does sofia 2 deal with the fully flagged channels?
    v_app = v_app/1000.


    # Need to check for that here if NaNs are included
    if f_sum < 0.:
        print_log(f'''CHECK_SOURCE: This galaxy has negative total flux. That will not work. Aborting.
''',Configuration['OUTPUTLOG'],screen = True)
        raise BadSourceError('We found an initial negative total flux.')
    galaxy_box = [[z_min,z_max],[y_min,y_max],[x_min,x_max]]
    if debug:
        print_log(f'''CHECK_SOURCE:  From the catalogue we got {Configuration['DISTANCE']}
''',Configuration['OUTPUTLOG'])
    # If the provided distance  = -1 we assume a Hubble follow
    if float(Configuration['DISTANCE']) == -1.:
        Configuration['DISTANCE'] = v_app/H_0
    if float(Configuration['DISTANCE']) < 0.5:
        Configuration['DISTANCE'] = 0.5
    if debug:
        print_log(f'''CHECK_SOURCE: We use a distance of {Configuration['DISTANCE']}.
''',Configuration['OUTPUTLOG'])
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
    vsys_error= np.mean(np.abs(np.array(VELboun,dtype=float)-v_app))*0.05
    rahr,dechr = convertRADEC(Configuration,ra,dec,debug=debug)
    # We write the results of the cut cube to the log
    print_log(f'''CHECK_SOURCE: The source finder found the following center in pixels.
{"":8s}CHECK_SOURCE: RA center = {x} with boundaries {x_min}, {x_max}
{"":8s}CHECK_SOURCE: DEC center = {y} with boundaries {y_min}, {y_max}
{"":8s}CHECK_SOURCE: V_sys center = {z} with boundaries {z_min}, {z_max}
{"":8s}CHECK_SOURCE: This should correspond to the WCS coordinates:
{"":8s}CHECK_SOURCE: RA center = {ra} deg or {rahr}  with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} deg  or {dechr}  with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: Sofia V_sys center = {v_app:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'])

    #There is a factor of two missing here but this is necessary otherwise the maxima are far to small
    Configuration['MAX_SIZE_IN_BEAMS'] = int(round(np.sqrt(((x_max-x_min)/2.)**2+((y_max-y_min)/2.)**2) \
                /(Configuration['BEAM_IN_PIXELS'][0])+5.))


    pa, inclination, SBR_initial, maj_extent,x_new,y_new,new_vsys,VROT_initial = rf.guess_orientation(Configuration,Fits_Files, center = [x,y],debug=debug)

    if x_new != x or y_new != y or new_vsys != v_app:

        x,y,z_new=cube_wcs.wcs_world2pix(ra,dec,new_vsys*1000.,1)
        x=float(x_new)
        y=float(y_new)
        z=float(z_new)
        ra,dec,v_app = cube_wcs.wcs_pix2world(x,y,z,1)
        v_app = v_app/1000.
        print_log(f'''CHECK_SOURCE: The center is updated to.
{"":8s}CHECK_SOURCE: RA center = {ra} with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: V_sys center = {v_app:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'])
    if np.sum(pa) == 0. or any(np.isnan(pa)) or \
        np.sum(inclination) == 0. or any(np.isnan(inclination)) or \
        np.sum(maj_extent) == 0. or np.isnan(maj_extent) or \
        np.sum(SBR_initial) == 0. or all(np.isnan(SBR_initial)) or \
        np.sum(VROT_initial) == 0. or all(np.isnan(VROT_initial)):
        print_log(f'''CHECK_SOURCE: We could not establish proper initial estimates from the moment maps. These are what we got
{"":8s}CHECK_SOURCE: pa = {pa}
{"":8s}CHECK_SOURCE: inclination = {inclination}
{"":8s}CHECK_SOURCE: maj_extent = {maj_extent}
{"":8s}CHECK_SOURCE: SBR_initial = {SBR_initial}
{"":8s}CHECK_SOURCE: VROT_initial = {VROT_initial}
''', Configuration['OUTPUTLOG'])

        raise BadSourceError("No initial estimates. Likely the source is too faint.")
     # Determine whether the centre is blanked or not


    #if debug:
    #    print_log(f'''CHECK_SOURCE: In the center we find the vsys {Central_Velocity} km/s around the location:
#{"":8s}CHECK_SOURCE: x,y,z = {int(round(x))}, {int(round(y))}, {int(round(z))}.
#{'':8s}CHECK_SOURCE: we will use a systemic velocity of {v_app}
#''',Configuration['OUTPUTLOG'])
    Central_Flux = np.mean(data[int(round(z-1)):int(round(z+1)),\
                                int(round(y-Configuration['BEAM_IN_PIXELS'][0]/2.)):int(round(y+Configuration['BEAM_IN_PIXELS'][0]/2.)),\
                                int(round(x-Configuration['BEAM_IN_PIXELS'][0]/2.)):int(round(x+Configuration['BEAM_IN_PIXELS'][0]/2.))])
    if debug:
        print_log(f'''CHECK_SOURCE: In the center we find an average flux of  {Central_Flux} Jy/beam around the location:
{"":8s}CHECK_SOURCE: x,y,z = {int(round(x))}, {int(round(y))}, {int(round(z))}.
''',Configuration['OUTPUTLOG'])

    if not np.isfinite(Central_Flux):
        Configuration['EXCLUDE_CENTRAL'] = True
        print_log(f'''CHECK_SOURCE: The flux in the central part is blanked. We exclude the central rings.
''',Configuration['OUTPUTLOG'])
    else:
        Configuration['EXCLUDE_CENTRAL'] = False

    #Check that the source is bright enough
    Max_SNR = np.nanmax(data)/Configuration['NOISE']
    if Max_SNR < 2.5:
        log_statement = f'''CHECK_SOURCE: The maximum Signal to Noise in this cube is {Max_SNR} that is not enough for a fit.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'],screen = True)
        raise BadSourceError(log_statement)

    # Size of the galaxy in beams
    Configuration['SIZE_IN_BEAMS'] = set_limits(maj_extent/(Configuration['BEAM'][0]/3600.),1.0,Configuration['MAX_SIZE_IN_BEAMS'])
    if Configuration['SIZE_IN_BEAMS'] <= Configuration['TOO_SMALL_GALAXY']:
        print_log(f'''CHECK_SOURCE: This galaxy has an estimated size of  {2*Configuration['SIZE_IN_BEAMS']} beams in diameter.
{'':8s}This is not large enough too fit. We will exit this fit.
''',Configuration['OUTPUTLOG'])
        raise BadSourceError('The extracted source is too small')
    set_ring_size(Configuration, debug = debug)
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
''',Configuration['OUTPUTLOG'])


    # The extent is fairly well determined and the maximum should be no more than +3 beams and a minimum no less than 4
    # Swithcing here from doubled outer rings causes problems though
    SNR_range=set_limits(Max_SNR,1.9, 2.6)

    Configuration['MAX_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']*3.125/SNR_range+3.,1.,Configuration['MAX_SIZE_IN_BEAMS'])
    Configuration['MIN_SIZE_IN_BEAMS'] = set_limits(Configuration['SIZE_IN_BEAMS']-3.,Configuration['TOO_SMALL_GALAXY'],Configuration['SIZE_IN_BEAMS'])
    Configuration['NO_RINGS'] = calc_rings(Configuration,debug=debug)
    Configuration['LAST_RELIABLE_RINGS'] = [Configuration['NO_RINGS'],Configuration['NO_RINGS']]
    print_log(f'''CHECK_SOURCE: From the original Configuration and SoFiA we find:
{"":8s}CHECK_SOURCE: The maximum diameter we will allow  is  {2.*Configuration['MAX_SIZE_IN_BEAMS']} beams. This is based on a SNR range of {SNR_range}
{"":8s}CHECK_SOURCE: The minimum diameter we will allow  is  {2.*Configuration['MIN_SIZE_IN_BEAMS']} beams.
{"":8s}CHECK_SOURCE: We start with a diameter of {2*Configuration['SIZE_IN_BEAMS']} beams in the model.
{"":8s}CHECK_SOURCE: SoFiA found a PA of {kin_pa:.2f} and we use a PA = {pa[0]:.2f} +/- {pa[1]:.2f}
{"":8s}CHECK_SOURCE: We start with an inclination of {inclination[0]:.2f} +/- {inclination[1]:.2f}{"":8s}
{"":8s}CHECK_SOURCE: SoFiA found a W50 of {w50:.2f} km/s
{"":8s}CHECK_SOURCE: We will use {2.*Configuration['NO_RINGS']} for the model with a ring size of {Configuration['RING_SIZE']}.
''',Configuration['OUTPUTLOG'])




    if Configuration['CHANNEL_DEPENDENCY'].lower() == 'hanning':
        vres = (Configuration['CHANNEL_WIDTH']*2)/(2.*np.sqrt(2.*np.log(2)))
    elif Configuration['CHANNEL_DEPENDENCY'].lower() == 'sinusoidal' :
        vres = (Configuration['CHANNEL_WIDTH']*1.2)/(2.*np.sqrt(2.*np.log(2)))
    elif Configuration['CHANNEL_DEPENDENCY'].lower() == 'independent':
        vres = Configuration['CHANNEL_WIDTH']
    else:
        raise BadConfigurationError('Something went wrong in the Configuration setup')

    #if inclination[0] < 40:
    #    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0]+5.)))
    #else:
    max_vrot=w50/2./np.sin(np.radians(abs(inclination[0])))
    max_vrot_dev=set_limits((VELboun[1]-VELboun[0])/4000./np.sin(np.radians(inclination[0])),4.*vres,20*vres,debug=debug)


    #Write the info to the Basic info File
    wf.basicinfo(Configuration,initialize = True,
              RA=[ra,abs(err_x*header['CDELT1'])],
              DEC=[dec,abs(err_y*header['CDELT2'])],
              VSYS =[v_app,abs(err_z*header['CDELT3']/1000.)],
              PA=pa, Inclination = inclination, Max_Vrot = [max_vrot,max_vrot_dev], Tot_Flux = [f_sum,f_sum_err], V_mask = [VELboun[1]-VELboun[0],vres],
              Distance = Configuration['DISTANCE'] , DHI = [2*maj_extent*3600.,Configuration['BEAM'][0]*Configuration['RING_SIZE']],debug=debug)


    # extract a PV-Diagram
    if not os.path.exists(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['SOFIA_BASENAME']}_sofia_xv.fits"):
        PV =extract_pv(Configuration,Cube, pa[0], center=[ra,dec,v_app*1000.], convert = 1000.,
                       finalsize = [int(round(maj_extent/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])*1.25+header['NAXIS1']*0.2)),
                                    int(round(z_max-z_min)+10.)],debug=debug)
        fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_sofia_xv.fits",PV[0].data,PV[0].header)
    Cube.close()
    Initial_Parameters = {}
    Initial_Parameters['XPOS'] = [ra,set_limits(abs(err_x*header['CDELT1']),0.1/3600.*Configuration['BEAM'][0],3./3600.*Configuration['BEAM'][0] )]
    Initial_Parameters['YPOS'] = [dec,set_limits(abs(err_y*header['CDELT2']),0.1/3600.*Configuration['BEAM'][0],3./3600.*Configuration['BEAM'][0] )]
    Initial_Parameters['VSYS'] =[v_app*1000.,vsys_error]
    Initial_Parameters['SBR_profile'] = SBR_initial
    #Initial_Parameters['VROT'] = [max_vrot/1000.,max_vrot_dev/1000.]
    Initial_Parameters['VROT_profile'] = VROT_initial
    Initial_Parameters['VROT'] = [np.max(VROT_initial) ,(np.max(VROT_initial)-np.min(VROT_initial))]
    Initial_Parameters['PA'] = pa
    Initial_Parameters['INCL'] = inclination
    Initial_Parameters['FLUX'] = [f_sum,f_sum_err]
    new_errors = np.array([pa[1],inclination[1],vsys_error])*0.1
    para = ['PA','INCL','VSYS']
    for err,parameter in zip(new_errors,para):
        Configuration['MIN_ERROR'][parameter] = [set_limits(err,Configuration['MIN_ERROR'][parameter][0],Configuration['MAX_ERROR'][parameter][0]/2.)]

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

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Initial_Parameters = Dictionary with all initial values obtained from the Sofia source finding and processing.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fitting_osc(Configuration,Fits_Files,Tirific_Template,Initial_Parameters):
    current_run = 'Not Initialized'
    create_directory(Configuration['USED_FITTING'],Configuration['FITTING_DIR'])
    wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                           Initial_Parameters= Initial_Parameters, \
                           fit_type=Configuration['USED_FITTING'],\
                           debug=Configuration['DEBUG'])
    print_log(f'''FITTING_OSC: The initial def file is written and we will now start fitting.
''' ,Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
    Configuration['PREP_END_TIME'] = datetime.now()
        # If we have no directory to put the output we create it

    while not Configuration['ACCEPTED'] and Configuration['ITERATIONS'] <  Configuration['MAX_ITERATIONS']:
        Configuration['ITERATIONS'] = Configuration['ITERATIONS']+1
        # Run the step
        current_run = one_step_converge(Configuration, Fits_Files,Tirific_Template,current_run,debug = Configuration['DEBUG'])

        if (Configuration['ITERATIONS'] == 1 or Configuration['ACCEPTED']) and Configuration['SIZE_IN_BEAMS'] < 4:
            if Configuration['DEBUG']:
                        print_log(f'''FITTING_OSC: Checking the inclination due to small galaxy size.
{'':8s}PA = {Tirific_Template['PA']}
{'':8s}INCL = {Tirific_Template['INCL']}
''',Configuration['OUTPUTLOG'])
            check_inclination(Configuration,Tirific_Template,Fits_Files, fit_type =Configuration['USED_FITTING'], debug = Configuration['DEBUG'])

    if Configuration['ACCEPTED']:
        print_log(f'''FITTING_OSC: The model has converged in center and extent and we make a smoothed version.
''',Configuration['OUTPUTLOG'],screen =True)
        check_angles(Configuration, Tirific_Template,debug=Configuration['DEBUG'] )
        current_run = fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,stage = 'after_os', fit_type = Configuration['USED_FITTING'],debug = Configuration['DEBUG'])
        if Configuration['OPTIMIZED']:
            make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_type = Configuration['USED_FITTING'],debug=Configuration['DEBUG'])
    elif Configuration['INSTALLATION_CHECK']:
        print_log(f'''FITTING_OSC: The Installation_check has run a fit successfully.
''',Configuration['OUTPUTLOG'],screen =True)
    else:
        Configuration['FINAL_COMMENT'] = 'We could not converge on the extent or centre of the galaxy'
        Configuration['OUTPUT_QUANTITY'] = 5
    return current_run
fitting_osc.__doc__ =f'''
 NAME:
    fitting_osc

 PURPOSE:
    Do the full one step convergence routines.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template
    Initial_Parameters = The initial guesses obtained from sofia

 OPTIONAL INPUTS:

 OUTPUTS:
    the tirific subprocess structure
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This is the main fitting routine if you want to make your own start by looking here. And then add it in main.py
       This is the fitting type done for the installation check.
'''

def fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run, stage = 'initial',fit_type='Undefined', debug = False):
    update_statistic(Configuration, message= "Starting the smoothed check run", debug=debug)
    if debug:
        print_log(f'''FIT_SMOOTHED_CHECK: Starting stage {stage} and fit_type {fit_type}.
''',Configuration['OUTPUTLOG'],debug = True)
    #if we have only a few rings we only smooth. else we fit a polynomial to the RC and smooth the SBR
    #smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',hdr, min_error= np.max([float(Tirific_Template['CFLUX']),float(Tirific_Template['CFLUX_2'])]),debug = debug)
    fix_sbr(Configuration,Tirific_Template,smooth = True, debug = debug)
    if stage == 'after_cc':
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',min_error = Configuration['CHANNEL_WIDTH'],debug = debug)
    else:
        min_error = []
        pars_to_smooth = []
        not_to_smooth = []
        fixed_errors = []
        for parameter in ['VROT','INCL','Z0','SDIS','PA','XPOS','YPOS','VSYS']:
            if parameter in Configuration['FIXED_PARAMETERS'][0]:
                not_to_smooth.append(parameter)
            else:
                pars_to_smooth.append(parameter)
            min_error.append(Configuration['MIN_ERROR'][parameter])

        for key,min_err in zip(pars_to_smooth,min_error):
                smoothed = regularise_profile(Configuration,Tirific_Template,key,min_error = Configuration['MIN_ERROR'][parameter],debug = debug)
                if key == 'VROT':
                    smoothed_vrot=copy.deepcopy(smoothed)

        for key,min_err in zip(not_to_smooth,fixed_errors):
            set_errors(Configuration,Tirific_Template,key,min_error = Configuration['MIN_ERROR'][parameter],debug = debug)

    #if stage == 'after_cc':
    #    smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',min_error = Configuration['CHANNEL_WIDTH'],debug = debug)
    #else:
    #    smoothed_vrot = regularise_profile(Configuration,Tirific_Template,'VROT',min_error = Configuration['CHANNEL_WIDTH'],debug = debug)

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
        Tirific_Template['LOOPS'] = f"{Configuration['LOOPS']}"
        xpos,ypos,vsys,pa,incl = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",Variables = ['XPOS','YPOS','VSYS','PA','INCL'])
        parameters = {'VSYS': [vsys[0], Configuration['CHANNEL_WIDTH']], \
                      'XPOS': [xpos[0], Configuration['PIXEL_SIZE']],
                      'YPOS': [ypos[0], Configuration['PIXEL_SIZE']],
                      'INCL':[incl[0],10./np.sin(np.radians(incl[0]))],
                      'PA':  [pa[0],5.],
                      'VROT': [np.max(smoothed_vrot),np.max(smoothed_vrot)*0.05],  }

        if incl[0] < 40 and Configuration['NO_RINGS'] > 5.:
                parameters_to_adjust = ['VSYS','SBR','XPOS','YPOS','PA','SDIS','VROT']
        else:
            parameters_to_adjust = ['NO_ADJUSTMENT'] #This triggers the default settings in set_fitting_parameters

        set_fitting_parameters(Configuration, Tirific_Template,stage = stage,\
                          initial_estimates = parameters,parameters_to_adjust = parameters_to_adjust, debug = debug)


    target = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Last_Unsmoothed_Input.def"
    source = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}_In.def")
    target = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Last_Unsmoothed_Input.def")
    os.system(f"cp {source} {target}")
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def',debug=debug)
    accepted,current_run = run_tirific(Configuration,current_run, stage = stage, fit_type=fit_type,debug=debug)

    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug = debug)
    if Configuration['NO_RINGS'] > 5.:
        smoothed_vrot = regularise_profile(Configuration,Tirific_Template,'VROT',min_error = Configuration['CHANNEL_WIDTH'],debug = debug)
        set_fitting_parameters(Configuration, Tirific_Template,stage = 'final_os',\
                          initial_estimates = parameters,parameters_to_adjust  = ['VROT'], debug = debug)
        source = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}_In.def")
        target = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_First_Smoothed_Input.def")
        os.system(f"cp {source} {target}")
        wf.tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def',debug=debug)
        accepted,current_run = run_tirific(Configuration,current_run, stage = 'final_os', fit_type=fit_type,debug=debug)

    return current_run
fit_smoothed_check.__doc__ =f'''
 NAME:
    fit_smoothed_check

 PURPOSE:
    fine tune the fit after smoothing all variables.

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template
    current_run = subprocces structure of active tirific

 OPTIONAL INPUTS:
    debug = False

    stage = 'initial'
    stage of the fitting

    fit_type = 'Undefined'
    type of the fitting

 OUTPUTS:
    tirific subprocess call

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def make_full_resolution(Configuration,Tirific_Template,Fits_Files,fit_type = 'Undefined', current_run = 'Not zed', debug = False):
    update_statistic(Configuration, message= "Starting to make a full resolution model run", debug=debug)
    if debug:
        print_log(f'''MAKE_FULL_RESOLUTION: creating full resolution for stage {fit_type}.
''',Configuration['OUTPUTLOG'],debug = True)
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Tirific_Template,debug = debug)
    Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Tirific_Template['LOOPS'] = "0"
    Tirific_Template['INIMODE'] = "0"
    wf.tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def')
    finish_current_run(Configuration,current_run,debug= debug)
    accepted,current_run = run_tirific(Configuration,'Not zed', stage = 'full_res', fit_type=fit_type,debug = debug)
    finish_current_run(Configuration,current_run,debug= debug)
make_full_resolution.__doc__ =f'''
 NAME:
    make_full_resolution

 PURPOSE:
    If we used an optimized cube we want to make the final fit at the original resolution.

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting

    current_run = 'Not zed'
    subprocess structure of tirific

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def one_step_converge(Configuration, Fits_Files,Tirific_Template,current_run, debug = False):

    print_log(f'''ONE_STEP_CONVERGENCE: Starting with loop {Configuration['ITERATIONS']} out of maximum {Configuration['MAX_ITERATIONS']}.
''',Configuration['OUTPUTLOG'],screen=True,debug = debug)
    fit_type = 'Fit_Tirific_OSC'
    stage = 'run_os'
    #First we run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = stage, fit_type = fit_type, debug= debug)
    #Then we load the produced output into our template
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def" , Tirific_Template, debug = debug)
    #Check that the centre does not deviate too much
    print_log(f'''CHECK THAT WE CAHANGE !!!!!!
{Tirific_Template['VROT']}
{Tirific_Template['PA']}
{Tirific_Template['XPOS']}
''',Configuration['OUTPUTLOG'])
    accepted_central = check_central_convergence(Configuration,Tirific_Template, fit_type = fit_type,debug=debug)
    print_log(f'''CHECK THAT WE CAHANGE !!!!!!
{Tirific_Template['VROT']}
{Tirific_Template['PA']}
{Tirific_Template['XPOS']}
''',Configuration['OUTPUTLOG'])

    # Check whether we have the correct sizes,
    accepted_size = check_size(Configuration,Tirific_Template, fit_type = fit_type, stage = stage, current_run = current_run, debug=debug,Fits_Files=Fits_Files)
    if accepted and accepted_size and accepted_central:
        Configuration['ACCEPTED'] = True
    else:
        Configuration['ACCEPTED'] = False
        if Configuration['ITERATIONS'] > Configuration['MAX_ITERATIONS']:
                print_log(f'''ONE_STEP_CONVERGENCE: We have ran the convergence more than {Configuration['MAX_ITERATIONS']} times aborting the fit.
    ''',Configuration['OUTPUTLOG'])
                if debug:
                    Configuration['ACCEPTED'] = True
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

        check_angles(Configuration,Tirific_Template,debug = debug)

                    # First we fix the SBR we are left with also to set the reliable ring to configuration.
        fix_sbr(Configuration,Tirific_Template,debug = debug)    # Then we determine the inner rings that should remain fixed

        get_inner_fix(Configuration, Tirific_Template,debug=debug)

        if all([True if x  in Configuration['FIXED_PARAMETERS'][0] else False for x in ['PA','INCL','Z0','SDIS']] ):
        #if all([Configuration['FIX_INCLINATION'][0],Configuration['FIX_PA'][0],Configuration['FIX_Z0'][0],Configuration['FIX_SDIS'][0]]):
            Configuration['WARP_SLOPE'] = [0.,0.]
        else:
            get_warp_slope(Configuration,Tirific_Template, debug=debug)

        set_cflux(Configuration,Tirific_Template,debug = debug)
        keys_to_smooth =['INCL','PA','SDIS','Z0','VROT']
        min_errors = [3.*np.mean(Configuration['LIMIT_MODIFIER']),2.,Configuration['CHANNEL_WIDTH']/(2.*Configuration['LIMIT_MODIFIER']), \
                        convertskyangle(Configuration,0.1,Configuration['DISTANCE'],physical= True)/Configuration['LIMIT_MODIFIER'],\
                        Configuration['CHANNEL_WIDTH']/(Configuration['LIMIT_MODIFIER'])]
        for j,key in enumerate(keys_to_smooth):
            #Smoothing the profile also fixes it
            smoothed = smooth_profile(Configuration,Tirific_Template,key,debug=debug,min_error=min_errors[j])

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_os',\
                             debug = debug)
        print_log(f'''CHECK THAT WE CAHANGE !!!!!!
        {Tirific_Template['VROT']}
        {Tirific_Template['PA']}
        {Tirific_Template['XPOS']}
        ''',Configuration['OUTPUTLOG'])
        wf.tirific(Configuration,Tirific_Template,name = f"{Configuration['USED_FITTING']}_In.def",debug = debug)
    return current_run
one_step_converge.__doc__ =f'''
 NAME:
    one_step_converge

 PURPOSE:
    Converge on a best fitmodel without smoothing.

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template
    current_run = subprocess structure of tirific

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    current_run = subprocess structure with the current tirific call.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This routine sets the Configuration['ACCEPTED']
'''

def sofia(Configuration, Fits_Files, debug = False):
    if debug:
        print_log(f'''RUN_SOFIA: starting sofia run from the template.
''',Configuration['OUTPUTLOG'],debug = True)
    update_statistic(Configuration, message= "Starting the SoFiA run", debug=debug)

    sofia_template = rf.sofia_template(debug=debug)
    create_directory('Sofia_Output',Configuration['FITTING_DIR'])
    os.chdir(Configuration['FITTING_DIR'])
    threshold = 5.
    counter = 3
    sofia_template['input.data'] = Fits_Files['FITTING_CUBE']
    spatial_kernels = [0,int(round(Configuration['BEAM_IN_PIXELS'][0])),int(round(Configuration['BEAM_IN_PIXELS'][0]))*2]
    if np.sum(Configuration['NAXES'][:2])/(2.*int(round(Configuration['BEAM_IN_PIXELS'][0])))  > 30:
        spatial_kernels.append(int(round(Configuration['BEAM_IN_PIXELS'][0]))*3)
        print_log(f'''RUN_SOFIA: Adding an extra kernel scale as the cube is more than 30 beams across.
''', Configuration['OUTPUTLOG'])

    velocity_kernels = [0,3,6,12]
    if Configuration['NAXES'][2] > 52:
        velocity_kernels.append(16)
        Configuration['VEL_SMOOTH_EXTENDED'] = True
        print_log(f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
''', Configuration['OUTPUTLOG'])
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration,debug=debug)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par')
        print_log("RUN_SOFIA: Running SoFiA. \n",Configuration['OUTPUTLOG'],screen = True)
        # Check which sofia to start
        sfrun = subprocess.Popen([Configuration['SOFIA2'],'sofia_input.par'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
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
                raise SofiaFaintError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'])
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration,debug=debug)
    Configuration['SOFIA_RAN'] = True
    os.chdir(Configuration['START_DIRECTORY'])
sofia.__doc__ =f'''
 NAME:
    sofia

 PURPOSE:
    run sofia on the cube to get a source mask and initial estimates.

 CATEGORY:
    run_functions

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

def tirshaker_call(Configuration,debug = False):
    # First we make a directory to keep all contained
    update_statistic(Configuration, message= "Starting the Tirshaker call run", debug=debug)
    if not os.path.isdir(f"{Configuration['FITTING_DIR']}/Error_Shaker/"):
        os.mkdir(f"{Configuration['FITTING_DIR']}/Error_Shaker/")


    # Then we open the final file
    filename = f"{Configuration['FITTING_DIR']}Error_Shaker/Error_Shaker_In.def"

    final_FAT_file= f"{Configuration['FITTING_DIR']}{Configuration['USED_FITTING']}/{Configuration['USED_FITTING']}.def"

    Tirific_Template = tirific_template(filename = final_FAT_file \
                    , debug = debug)
    #Change the name and run only 2 LOOPS
    Tirific_Template['RESTARTNAME']= f"{Configuration['LOG_DIRECTORY']}restart_Error_Shaker.txt"
    Tirific_Template['INSET'] = f"../{Tirific_Template['INSET']}"
    Tirific_Template['TIRDEF']= f"Error_Shaker_Out.def"
    Tirific_Template['LOOPS'] = '1'

    outfilename = 'Error_Shaker.def'
    outfileprefix = 'Error_Shaker'

    #This we need to remove from the actual run
    #Configuration['NO_RINGS'] = int(Tirific_Template['NUR'])
    #Configuration['RING_SIZE'] = 1.
    #Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']-2
    set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_os',\
                         debug = debug)

    #Write it back to the file
    wf.tirific(Configuration,Tirific_Template, name =f"Error_Shaker/Error_Shaker_In.def" , debug = debug)

    #Determine the error block from the last fit settings.
    parameter_groups,rings,block,variation,variation_type = get_fit_groups(Configuration,Tirific_Template,debug = debug)
    if debug:
        print_log(f'''TIRSHAKER_CALL: We are shaking with the following parameters:
{'':8s}groups = {parameter_groups}
{'':8s}rings = {rings}
{'':8s}block = {block}
{'':8s}variation = {variation}
{'':8s}variation_type = {variation_type}
''',Configuration['OUTPUTLOG'])

    iterations = Configuration['SHAKER_ITERATIONS']
    random_seed = 2
    os.chdir(f"{Configuration['FITTING_DIR']}/Error_Shaker/")
    tirshaker(Configuration,Tirific_Template, outfilename = outfilename,\
                outfileprefix = outfileprefix, parameter_groups = parameter_groups, \
                rings = rings, block = block, variation = variation,\
                 variation_type = variation_type, iterations = iterations,
                 random_seed = random_seed, mode = 'mad',debug=debug)
    os.chdir(f"{Configuration['START_DIRECTORY']}")
    '''
    try:
        os.remove(f"{Configuration['FITTING_DIR']}Error_Shaker/blatirshin")
    except:
        pass
    try:
        os.remove(f"{Configuration['FITTING_DIR']}Error_Shaker/blatirlog")
    except:
        pass
    #remove the duplicate line in the shaker output
    # make sure all duplicate are removed

    outshaker = open(f"{Configuration['FITTING_DIR']}Error_Shaker/{outfilename}", 'r').readlines()
    lines_set = set(outshaker)

    with open(f"{Configuration['FITTING_DIR']}Error_Shaker/{outfilename}", 'w') as out:
        for line in lines_set:
            out.write(line)

    #Let's load the   the tirshaker output
    err_var =[]
    for group in parameter_groups:
        for par in group:
            if f'ERR_{par[:-1]}' not in err_var:
                err_var.append(f'ERR_{par[:-1]}')
    errors = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}Error_Shaker/{outfilename}",Variables = err_var,
                     unpack = False , debug = debug )

    for para in err_var:
        key=para[4:]
        FAT_err = f'# {key}_ERR'
        format = set_format(key)
        if FAT_err in Tirific_Template:
            Tirific_Template[FAT_err]= f"{' '.join([f'{x:{format}}' for x in errors[:int(Configuration['NO_RINGS']),err_var.index(para)]])}"
        else:
            Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in errors[:int(Configuration['NO_RINGS']),err_var.index(para)]])}")
'''
    wf.tirific(Configuration,Tirific_Template,name=f"{Configuration['USED_FITTING']}/{Configuration['USED_FITTING']}.def", debug = debug)

tirshaker_call.__doc__ =f'''
 NAME:
    tirshaker

 PURPOSE:
    function to setup the right input for tirshaker and call it and afterwards process the results

 CATEGORY:
    run_functions

 INPUTS:
    Configuration

 OPTIONAL INPUTS:
    debug = False


 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:


'''
