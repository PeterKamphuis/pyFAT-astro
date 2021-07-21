# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to Modify the Tirific_Template



import copy
from pyFAT_astro.Support.support_functions import set_rings,convertskyangle,sbr_limits,set_limits,print_log,set_limit_modifier,\
                              set_ring_size,calc_rings,finish_current_run,set_format,get_from_template,gaussian_function,fit_gaussian,\
                              get_ring_weights

import numpy as np
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline,Akima1DInterpolator

#Define some errors
class InitializeError(Exception):
    pass
class CfluxError(Exception):
    pass
class FunctionCallError(Exception):
    pass




def arc_tan_function(axis,center,length,amplitude,mean):
    # to prevent instant turnover
    c = axis[-1]*0.1
    #and the turnover has to be beyon 20*of the inner part
    c2 = set_limits(axis[-1]*0.2,axis[2],axis[int(len(axis)/1.5)])
    return -1*np.arctan((axis-(c2+abs(center)))/(c+abs(length)))*amplitude/np.pi + mean
arc_tan_function.__doc__ =f'''
 NAME:
    arc_tan_function

 PURPOSE:
    arc tangent function

 CATEGORY:
    modify_template

 INPUTS:
    axis = xaxis
    center = 0 point of the acrtan
    length = turnover length
    amplitude = height of the arctan
    mean = amplitude zero point

 OPTIONAL INPUTS:

 OUTPUTS:
    arctan function

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_angles(Configuration,Tirific_Template, debug = False):
    incl = get_from_template(Configuration,Tirific_Template, ['INCL','INCL_2'],debug=debug)

    if incl[0][0] > 90. and incl[1][0] > 90.:
        Tirific_Template['INCL'] = f"{' '.join(f'{180.-x:.2e}' for x in incl[0])}"
        Tirific_Template['INCL_2'] = f"{' '.join(f'{180.-x:.2e}' for x in incl[1])}"
    pa = get_from_template(Configuration,Tirific_Template, ['PA','PA_2'],debug=debug)

    if pa[0][0] > 360. and pa[1][0] > 360.:
        Tirific_Template['PA'] = f"{' '.join(f'{x-360.:.2e}' for x in pa[0])}"
        Tirific_Template['PA_2'] = f"{' '.join(f'{x-360.:.2e}' for x in pa[1])}"

    if pa[0][0] < 0. and pa[1][0] < 0.:
        Tirific_Template['PA'] = f"{' '.join(f'{360.-x:.2e}' for x in pa[0])}"
        Tirific_Template['PA_2'] = f"{' '.join(f'{360.-x:.2e}' for x in pa[1])}"
check_angles.__doc__=f'''
 NAME:
    check_angles

 PURPOSE:
       Check whether PA and INCLination are in the range 0-360. and < 90 in the center. If not modify all angles such that they are.

 CATEGORY:
       modify_template

 INPUTS:
    Configuration  = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    No output Tirific_Template is modified

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      get_from_template, join

 NOTE:
'''


def check_flat(Configuration,profile,error,key, last_reliable_ring = -1,inner_fix = 4, debug = False):
    if last_reliable_ring == -1:
        last_reliable_ring = len(profile)-1
    if debug:
        print_log(f'''CHECK_FLAT: checking flatness
{'':8s}profile = {profile}
{'':8s}error = {error}
{'':8s}last_reliable_ring = {last_reliable_ring}
''',Configuration['OUTPUTLOG'],debug = True)

    inner = np.mean(profile[:inner_fix])
    mean_error = np.mean(error[:last_reliable_ring])
    if inner_fix+1 < last_reliable_ring:
        outer = np.mean(profile[inner_fix+1:last_reliable_ring])
        outer_std = np.std(profile[inner_fix+1:last_reliable_ring])
    else:
        outer = inner
        outer_std = np.mean(error)

    if key in ['SDIS']:
        outer_std = 2.*mean_error
    if abs(outer-inner) < mean_error or outer_std  < mean_error:
        if debug:
            print_log(f'''CHECK_FLAT: If  {abs(outer-inner)} less than {mean_error} or
{'':8s} The outer variation {outer_std} less than the median error {mean_error} we break and set flat.
''',Configuration['OUTPUTLOG'])
        return True
    if key in ['INCL','INCL_2']:
        if outer < 40. or inner < 40.:
            if debug:
                print_log(f'''CHECK_FLAT: If  the outer profile is {outer} hence we set this to flat.
''',Configuration['OUTPUTLOG'])
            return True

    for e,x,y in zip(error[1:last_reliable_ring],profile[1:last_reliable_ring],profile[0:last_reliable_ring]):
        if debug:
            print_log(f'''CHECK_FLAT: x = {x}, y = {y}, e = {e}
''',Configuration['OUTPUTLOG'])
        if not x-e/2. < y < x+e/2.:
            if debug:
                print_log(f'''CHECK_FLAT: This taco is bend
''',Configuration['OUTPUTLOG'])
            return False
    if debug:
            print_log(f'''CHECK_FLAT: All values were within the error of eachother
''',Configuration['OUTPUTLOG'])
    return True
check_flat.__doc__ = '''
 NAME:
    check_flat

 PURPOSE:
       Check whether within its errors a routine is varying compared to the prvious rings

 CATEGORY:
       modify_template

 INPUTS:
    Configuration  = Standard FAT configuration
    profile = the profile to examine
    error = The accompanying error

 OPTIONAL INPUTS:
    debug = False

    last_reliable_ring = -1
    the last ring in the profile that can be trusted

 OUTPUTS:
       True if no variation is found false if variation is found

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      zip(),np.std,np.mean,abs,print_log

 NOTE:
'''

def check_size(Configuration,Tirific_Template, fit_type = 'Undefined', stage = 'initial', Fits_Files= 'No Files' ,debug = False,current_run='Not Initialized'):
    if debug:
        print_log(f'''CHECK_SIZE: Starting a new Check_size with the following parameters:
{'':8s}CHECK_SIZE: Rings = {Configuration['NO_RINGS']}
{'':8s}CHECK_SIZE: Size in Beams = {Configuration['SIZE_IN_BEAMS']}
''',Configuration['OUTPUTLOG'],debug=True)


    radii, sbr_ring_limits = sbr_limits(Configuration,systemic = float(Tirific_Template['VSYS'].split()[0]),debug=debug)
    #get the sbr profiles

    sbr = np.array(get_from_template(Configuration,Tirific_Template, ['SBR','SBR_2'],debug=debug),dtype = float)
    if debug:
        print_log(f'''CHECK_SIZE: This is the sizes
{'':8s}CHECK_SIZE: SBR = {len(sbr[0])},{len(sbr[1])}
{'':8s}CHECK_SIZE: limits = {len(sbr_ring_limits)}
{'':8s}CHECK_SIZE: No. Rings  = {Configuration['NO_RINGS']}
''',Configuration['OUTPUTLOG'])
    if len(sbr[0]) != len(sbr_ring_limits):
        #Check what the appropriate size should be
        if debug:
            print_log(f'''CHECK_SIZE: Equalizing the sizes
''',Configuration['OUTPUTLOG'])
        if len(sbr[0]) != Configuration['NO_RINGS']:
            if debug:
                print_log(f'''CHECK_SIZE: Interpolating SBR
''',Configuration['OUTPUTLOG'])
            old_radii = np.array(get_from_template(Configuration,Tirific_Template, ['RADI'],debug=debug),dtype = float)
            for i in [0,1]:
                sbr[i] = np.interp(np.array(radii,dtype=float),np.array(old_radii[0],dtype=float),np.array(sbr[i],dtype=float))

    if debug:
        print_log(f'''CHECK_SIZE: This after correcting.
{'':8s}CHECK_SIZE: SBR = {len(sbr[0])}
{'':8s}CHECK_SIZE: limits = {len(sbr_ring_limits)}
{'':8s}CHECK_SIZE: No. Rings  = {Configuration['NO_RINGS']}
''',Configuration['OUTPUTLOG'])
    #sbr_ring_limits = 1.25*np.array([sbr_ring_limits,sbr_ring_limits])
    sbr_ring_limits = np.array([sbr_ring_limits,sbr_ring_limits],dtype=float)
    new_rings = get_number_of_rings(Configuration,sbr,sbr_ring_limits, debug=debug)
    # if all of the last points are  below the limit we start checking how far to cut
    #the lower the inclination the sooner the RC becomes unreliable
    limit_factor = set_limits(2.5*np.mean(Configuration['LIMIT_MODIFIER']) ,2.,4.)
    if debug:
        print_log(f'''CHECK_SIZE: Using a limit factor for reliable RC of  {limit_factor}
''',Configuration['OUTPUTLOG'])
    Configuration['RC_UNRELIABLE'] = get_number_of_rings(Configuration,sbr,limit_factor*sbr_ring_limits, debug=debug)-1
    if Configuration['RC_UNRELIABLE'] == Configuration['NO_RINGS']:
        Configuration['RC_UNRELIABLE'] -= 1
    for i in [0,1]:
        corr_val = np.where(sbr[i,2:] > sbr_ring_limits[i,2:])[0]+2
        if corr_val.size > 0:
            Configuration['LAST_RELIABLE_RINGS'][i] = corr_val[-1]+1
        else:
            Configuration['LAST_RELIABLE_RINGS'][i] = Configuration['NO_RINGS']
    if debug:
        print_log(f'''CHECK_SIZE: We set these as the last reliable rings {Configuration['LAST_RELIABLE_RINGS']}
''',Configuration['OUTPUTLOG'])
    #if we haven't subtracted we check if we should add
    if int(new_rings) == int(Configuration['NO_RINGS']):
        if (np.any(sbr[:,-2] > sbr_ring_limits[:,-2]*7.) and np.any(sbr[:,-1] > sbr_ring_limits[:,-1]*3.)) or \
            np.any(sbr[:,-1] > sbr_ring_limits[:,-1]*5.):
            if debug:
                print_log(f'''CHECK_SIZE: The last rings were found to be:
{'':8s}{sbr[:,-2:]}
{'':8s}and the limits:
{'':8s}{sbr_ring_limits[:,-2:]}
{'':8s}Thus we add a ring.
''', Configuration['OUTPUTLOG'])
            new_rings += 1
        else:
            if debug:
                print_log(f'''CHECK_SIZE: The last rings were found to be:
{'':8s}{sbr[:,-2:]}
{'':8s}and the limits:
{'':8s}{sbr_ring_limits[:,-2:]}
{'':8s}Thus we keep the ring size.
''', Configuration['OUTPUTLOG'])
    else:
        if debug:
            print_log(f'''CHECK_SIZE: The last rings were found to be:
{'':8s}{sbr[:,-2:]}
{'':8s}and the limits:
{'':8s}{sbr_ring_limits[:,-2:]}
{'':8s}Thus we have subtracted a set of rings.
''', Configuration['OUTPUTLOG'])

    if new_rings <= Configuration['NO_RINGS']:
        size_in_beams = (radii[new_rings-1]-1./5.*Configuration['BEAM'][0])/Configuration['BEAM'][0]
    else:
        size_in_beams = (radii[-1]-1./5.*Configuration['BEAM'][0])/Configuration['BEAM'][0]+Configuration['RING_SIZE']
    if debug:
        print_log(f'''CHECK_SIZE: Before checking against the minimum and maximum the size in beams = {size_in_beams}
''', Configuration['OUTPUTLOG'],debug=True)
    size_in_beams = set_limits(size_in_beams, Configuration['MIN_SIZE_IN_BEAMS'], Configuration['MAX_SIZE_IN_BEAMS'])
    # limit between 3 and the maximum allowed from the sofia estimate
    size_in_beams,ring_size,number_of_rings = set_ring_size(Configuration, size_in_beams=size_in_beams,check_set_rings = True, debug=debug)

    print_log(f'''CHECK_SIZE: We find the following size in beams {size_in_beams:.1f} with size {ring_size:.1f}.
{'':8s}CHECK_SIZE: The previous iteration had a size of {Configuration['SIZE_IN_BEAMS']:.1f} with rings  {Configuration['RING_SIZE']:.1f} times the beam. .
{'':8s}CHECK_SIZE: This results in {int(number_of_rings):.1f} rings in the model compared to {int(Configuration['NO_RINGS'])} previously.
''', Configuration['OUTPUTLOG'],screen=True)
    if f"{ring_size:.1f}" != f"{Configuration['RING_SIZE']:.1f}":
        Configuration['NEW_RING_SIZE'] = True
    else:
        Configuration['NEW_RING_SIZE'] = False
    if debug:
        print_log(f'''CHECK_SIZE: The previous rings were {Configuration['OLD_RINGS']}
''',Configuration['OUTPUTLOG'])

    if f"{size_in_beams:.1f}" == Configuration['OLD_RINGS'][-1]:
        return True
    elif f"{size_in_beams:.1f}" in Configuration['OLD_RINGS']:
        print_log(f'''CHECK_SIZE: We have processed this size before.
''', Configuration['OUTPUTLOG'])
        ind = Configuration['OLD_RINGS'].index(f"{size_in_beams:.1f}")
        if Configuration['OLD_RINGS'][ind+1] > Configuration['OLD_RINGS'][ind]:
            print_log(f'''CHECK_SIZE: After which we increased the size.
''', Configuration['OUTPUTLOG'])
            if Configuration['OLD_RINGS'][ind+1] == Configuration['OLD_RINGS'][-1]:
                print_log(f'''CHECK_SIZE: Which is the current fit so that's ok.
''', Configuration['OUTPUTLOG'])
                return True
            else:
                print_log(f'''CHECK_SIZE: Which is not the current fit so we refit this size.
''', Configuration['OUTPUTLOG'])
            Configuration['OLD_RINGS'].append(f"{size_in_beams:.1f}")
        else:
            print_log(f'''CHECK_SIZE: After which we decreased so we allow this addition. But no more.
''', Configuration['OUTPUTLOG'])
            Configuration['OLD_RINGS'].append(f"{size_in_beams:.1f}")
    else:
        if debug:
            print_log(f'''CHECK_SIZE: Adding the new ring size to OLD_RINGS.
''', Configuration['OUTPUTLOG'])
        Configuration['OLD_RINGS'].append(f"{size_in_beams:.1f}")

    Configuration['RING_SIZE'] = ring_size
    Configuration['SIZE_IN_BEAMS'] = size_in_beams
    Configuration['NO_RINGS'] = calc_rings(Configuration,debug=debug)
    print_log(f'''CHECK_SIZE: We need to modify the number of rings in the model.
''', Configuration['OUTPUTLOG'],screen=True)

    if Configuration['RC_UNRELIABLE'] > Configuration['NO_RINGS']:
        Configuration['RC_UNRELIABLE'] = Configuration['NO_RINGS']
    for i in [0,1]:
        if Configuration['LAST_RELIABLE_RINGS'][i] > Configuration['NO_RINGS']:
            Configuration['LAST_RELIABLE_RINGS'][i] = Configuration['NO_RINGS']
    if debug:
        print_log(f'''CHECK_SIZE: We trust the RC upto ring {Configuration['RC_UNRELIABLE']}
{'':8s} and the rings in general upto {Configuration['LAST_RELIABLE_RINGS']}
''',Configuration['OUTPUTLOG'])
    if Fits_Files == 'No Files':
        raise InitializeError('CHECK_SIZE: Trying to adapt the model size but the fits files were not provided.')
    else:
        # Do not move this from here else other routines such as sbr_limits are messed up
        set_new_size(Configuration,Tirific_Template,Fits_Files,fit_type= fit_type, stage = stage ,debug = debug,current_run = current_run)
    return False
check_size.__doc__  =f'''
 NAME:
    check_size

 PURPOSE:
    Check and update the size of the model. Update RC_UNRELIABLE and LAST_RELIABLE_RINGS

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting

    stage = 'initial'
    Stage of the fitting process

    Fits_Files= 'No Files'
    Standard FAT dictionary with filenames

    current_run='Not Initialized'
    subproccess to be stopped for the tirific run when rings change

 OUTPUTS:
    Boolean which is True when rings are not updated, False if they are

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Configuration['RC_UNRELIABLE'] = modified and updated here.
'''


def fit_arc(Configuration,radii,sm_profile,error, debug = False ):

    c2 = set_limits(radii[-1]*0.2,radii[2],radii[int(len(radii)/1.5)])
    est_center = radii[-1]/2.-c2
    est_length = radii[-1]*0.1
    est_amp = abs(np.max(sm_profile)-np.min(sm_profile))
    est_mean = np.mean(sm_profile)

    if not error.any():
        error = np.full(len(y),1.)
        absolute_sigma = False
    else:
        absolute_sigma = True

    arc_par,arc_cov  =  curve_fit(arc_tan_function, radii, sm_profile,p0=[est_center,est_length,est_amp,est_mean]\
                                ,sigma=error,absolute_sigma=absolute_sigma)

    new_profile = arc_tan_function(radii,*arc_par)

    new_profile[:3] = np.mean(new_profile[:3])

    return new_profile#,new_error
fit_arc.__doc__ =f'''
 NAME:
    fit_arc
 PURPOSE:
    Fit an arc tangent function
 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    radii = the horizontal axis
    sm_profile = the  profile to examine
    error = The accompanying error

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    the fitted profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fit_polynomial(Configuration,radii,profile,sm_profile,error, key, Tirific_Template,inner_fix = 4,min_error =0.,boundary_limits = [0,0.], debug = False ):
    if debug:
        print_log(f'''FIT_POLYNOMIAL: starting to fit the polynomial with the following input:
{'':8s} key = {key}
{'':8s} radii = {radii}
{'':8s} profile = {profile}
{'':8s} sm_profile = {sm_profile}
{'':8s} error = {error}
''',Configuration['OUTPUTLOG'], debug=debug)
    only_inner = False
    if key in ['PA','INCL','Z0']:
        fixed = inner_fix
        only_inner =True
        error[:fixed] = error[:fixed]/10.
    elif key in ['VROT']:
        #fixed =len(radii)-Configuration['OUTER_SLOPE_START']
        fixed =set_limits(len(radii)-np.min(Configuration['LAST_RELIABLE_RINGS']),1,len(radii))
        error[np.min(Configuration['LAST_RELIABLE_RINGS']):] = Configuration['CHANNEL_WIDTH']
        error[0] = Configuration['CHANNEL_WIDTH']
        error[1] = error[1]*3.
    else:
        fixed = 1
    if len(radii) > 15.:
        start_order = int(len(radii)/5)
    else:
        start_order = 0

    st_fit = int(0)
    if key in ['VROT']:
        if np.mean(profile[1:3]) > 120.:
            st_fit = int(1)

        #This needs another -1 because the 0 and 1/5. ring are more or less 1 ring
        max_order = set_limits(len(radii)-fixed-2,3,8)
        #The rotation curve varies a lot so the lower limit should be as high as possible
        #But at least 3 less than max order and maximally 4
        if len(radii)-fixed-2 <= 6:
            lower_limit=set_limits(3,3,max_order-2)
        elif len(radii)-fixed-2 <= 10:
            lower_limit=set_limits(4,3,max_order-2)
        else:
            lower_limit=set_limits(5,3,max_order-2)
        start_order = set_limits(start_order,lower_limit,max_order)

    else:
        if key in ['PA','INCL','Z0']:
            max_order = set_limits(len(radii)-fixed,3,7)
        else:
            max_order = set_limits(len(radii)-1,3,7)


    if start_order >= max_order:
        max_order = max_order+1
    if debug:
        print_log(f'''FIT_POLYNOMIAL: For {key} we start at {start_order} because we have {len(radii)} rings of which {fixed} are fixed
{'':8s} this gves us a maximum order of {max_order}
''',Configuration['OUTPUTLOG'])

    reduced_chi = []
    order = range(start_order,max_order+1)
    if debug:
        print_log(f'''FIT_POLYNOMIAL: We will fit the following radii.
{'':8s}{radii[st_fit:]}
{'':8s} and the following profile:
{'':8s}{profile[st_fit:]}
{'':8s} weights = {1./error[st_fit:]}
''',Configuration['OUTPUTLOG'])

    #make sure there are no 0. in the errors
    zero_locations = np.where(error[st_fit:] == 0.)[0]
    if zero_locations.size > 0.:
        error[zero_locations+st_fit] = 1./np.nanmax(1./error)

    for ord in order:
        fit_prof = np.poly1d(np.polyfit(radii[st_fit:],profile[st_fit:],ord,w=1./error[st_fit:]))
        if st_fit > 0.:
            fit_profile = np.concatenate(([sm_profile[0]],[e for e in fit_prof(radii[st_fit:])]))
        else:
            fit_profile = fit_prof(radii)
        #fit_profile = fit_prof(radii)
        if key != 'SBR':
            fit_profile = fix_profile(Configuration, key, fit_profile, Tirific_Template,inner_fix=inner_fix, singular = True,only_inner =only_inner)
        red_chi = np.sum((profile[st_fit:]-fit_profile[st_fit:])**2/error[st_fit:])/(len(radii[st_fit:])-ord)
        #We penailze profiles that go outside the boundaries

        if np.sum(np.absolute(np.array(boundary_limits,dtype=float))) != 0.:
                diff = np.sum(np.array([abs(x-set_limits(x,\
                                                  boundary_limits[0],\
                                                  boundary_limits[1])) \
                                        for x in  fit_profile[st_fit:]],dtype = float))
                if diff > 1.:
                    red_chi = red_chi*(diff)

        reduced_chi.append(red_chi)
        #if key in ['VROT'] and Configuration['NO_RINGS'] < 2.5*max_order:
        #    reduced_chi[-1] = reduced_chi[-1]*(ord/Configuration['NO_RINGS'])**2.5
    if debug:
        print_log(f'''FIT_POLYNOMIAL: We have fitted these:
{'':8s} order = {[x for x in order]}
{'':8s} reducuced chi = {reduced_chi}
''',Configuration['OUTPUTLOG'])
    reduced_chi = np.array(reduced_chi,dtype = float)
    final_order = order[np.where(np.min(reduced_chi ) == reduced_chi )[0][0]]

    print_log(f'''FIT_POLYNOMIAL: We have regularised {key} with a polynomial of order {final_order}.
''',Configuration['OUTPUTLOG'])
    fit_profile = np.poly1d(np.polyfit(radii[st_fit:],profile[st_fit:],final_order,w=1./error[st_fit:]))
    if st_fit > 0.:
        new_profile = np.concatenate(([sm_profile[0]],[e for e in fit_profile(radii[st_fit:])]))
    else:
        new_profile = fit_profile(radii)
    #if key in ['VROT'] and profile[1] < profile[2]:
    #    new_profile[1] = profile[1]
    if key != 'SBR':
        new_profile = fix_profile(Configuration, key, new_profile, Tirific_Template,debug =debug,inner_fix=inner_fix,singular = True,only_inner =only_inner)

    return new_profile#,new_error
fit_polynomial.__doc__ =f'''
 NAME:
    fit_polynomial

 PURPOSE:
    Fit a polynomial between 3 and 8 order and determine the one with the optimal reduced Chi^2 and return the regularised profile

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    radii = rings in arcsec of profiles to regularise
    profile = profile to be regularised
    sm_profile = smoothed profile
    error = errors on the profile
    key = parameter that is being regularised
    Tirific_Template = standard tirific template

 OPTIONAL INPUTS:
    debug = False

    min_error =0.
    the error should always be large than this value

 OUTPUTS:
    polynomial fitted profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fix_outer_rotation(Configuration,profile, debug = False):
    if debug:
        print_log(f'''FIX_OUTER_ROTATION: adjust last rings of VROT profile:
{'':8s}{profile}
''',Configuration['OUTPUTLOG'],debug = True)
    profile = np.array(profile,dtype=float)
    inner_slope = Configuration['RC_UNRELIABLE']
    if debug:
        print_log(f'''FIX_OUTER_ROTATION: this is the inner slope {inner_slope}
''',Configuration['OUTPUTLOG'])
    NUR = Configuration['NO_RINGS']
    #inner_slope = int(round(set_limits(NUR*(4.-Configuration['LIMIT_MODIFIER'][0])/4.,round(NUR/2.),NUR-2)))
    if inner_slope != NUR-1 and np.mean(profile[1:3]) > 180.:
        profile[inner_slope:] = profile[inner_slope-1]

    for i in range(int(Configuration['NO_RINGS']*3./4),Configuration['NO_RINGS']-1):
        if profile[i+1] > profile[i]*1.3:
            profile[i+1] = profile[i]*1.3

    if debug:
        print_log(f'''FIX_OUTER_ROTATION: this is corrected profile:
{profile}
''',Configuration['OUTPUTLOG'])

    return profile
fix_outer_rotation.__doc__ =f'''
 NAME:
    fix_outer_rotation

 PURPOSE:
    Fix the outer parts of the rotation curve against upward outliers

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    profile = the RC

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    profile = the corrected profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Declining rotation curves are dealt with in no_declining_vrot
'''

def fix_profile(Configuration, key, profile, Tirific_Template, debug= False, inner_fix = [4,4], singular = False,only_inner = False ):

    if isinstance(inner_fix,int):
        inner_fix = [inner_fix]
    if debug:
        print_log(f'''FIX_PROFILE: Starting to fix {key} with the input values:
{'':8s}{profile}
''', Configuration['OUTPUTLOG'],debug=True)

    if key == 'SBR':
        print_log(f'''FIX_PROFILE: To fix sbr profiles use FIX_SBR.
''', Configuration['OUTPUTLOG'],screen=True)
        raise FunctionCallError("FIX_PROFILE: To fix sbr profiles use check SBR.")
    if singular:
        indexes = [0]
        profile = np.array([profile,profile])
        inner_mean = np.nanmean([profile[0,:inner_fix[0]]])
    else:
        indexes = [0,1]
        if np.sum(inner_fix) != 0.:
            inner_mean = np.nanmean(np.concatenate((profile[0,:inner_fix[0]],profile[1,:inner_fix[1]])))


    profile = np.array(profile,dtype=float)


    if key in ['VROT']:
        indexes = [0]
        inner_mean = 0.
    else:
        for i in indexes:
            profile[i,:inner_fix[i]] = inner_mean
        if debug:
            print_log(f'''FIX_PROFILE: the  {inner_fix} inner rings are fixed for the profile:
{'':8s} profile = {profile[i,:]}
''', Configuration['OUTPUTLOG'])
    if debug:
        print_log(f'''FIX_PROFILE: the  inner mean is {inner_mean}.
''', Configuration['OUTPUTLOG'])

    for i in indexes:
        if Configuration['LAST_RELIABLE_RINGS'][i] < len(profile[i,:]) and not only_inner:
            if debug:
                print_log(f'''FIX_PROFILE: From ring {Configuration['LAST_RELIABLE_RINGS'][i]} on we do not trust these rings.
''', Configuration['OUTPUTLOG'])
            #profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:] = profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:]*0.25+profile[i,Configuration['LAST_RELIABLE_RINGS'][i]-1]*0.75
            profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:] = profile[i,Configuration['LAST_RELIABLE_RINGS'][i]-1]
        if debug:
            print_log(f'''FIX_PROFILE:After fixing the last reliable rings
{'':8s} profile = {profile[i,:]}
''', Configuration['OUTPUTLOG'])
        if key == 'VROT':
            profile[i] =fix_outer_rotation(Configuration,profile[i],debug= debug)
        if key in ['PA','INCL','Z0']:
            xrange = set_limits((int(round(len(profile[0])-5.)/4.)),1,4)
        # need to make sure this connects smoothly
            for x in range(0,xrange):
                if inner_fix[i]+x < len(profile[i,:]):
                    profile[i,inner_fix[i]+x] = 1/(x+4./xrange)*inner_mean+ (1-1/(x+4./xrange))*profile[i,inner_fix[i]+x]

        if debug:
            print_log(f'''FIX_PROFILE:After smoothe transition
{'':8s} profile = {profile[i,:]}
''', Configuration['OUTPUTLOG'])
        #profile[:,:Configuration['INNER_FIX']] = np.nanmean(profile[:,:Configuration['INNER_FIX']])
        if key in ['SDIS']:
            inner_max = np.nanmax(profile[i,:int(len(profile[i,:])/2.)])
            if inner_mean < inner_max:
                ind = np.where(profile[i,:] == inner_max)[0]
                if ind.size > 1:
                    ind = int(ind[0])
                else:
                    ind = int(ind)
                profile[i,:ind] = inner_max
            profile[i] =np.hstack([[profile[i,0]],[y if y <= x else x*0.95 for x,y in zip(profile[i,:],profile[i,1:])]])
    if key == 'VROT':
        tmp  = no_declining_vrot(Configuration,Tirific_Template,profile = profile[0],debug=debug)
        profile[0] = tmp
        profile[1] = tmp
    if singular:
        profile = profile[0]


    if debug:
        print_log(f'''FIX_PROFILE: The final profile for {key} is:
{'':8s}{profile}
''', Configuration['OUTPUTLOG'])
    return profile
fix_profile.__doc__ =f'''
 NAME:
    fix_profile

 PURPOSE:
    Modify a fitted profile to be stable against outliers

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    key = Parameter to be fixed
    profile = profile to be fixed # Maybe we can read this from the template?
    Tirific_Template = Standard Tirific template

 OPTIONAL INPUTS:
    debug = False

    singular = False
    Profile is a single side

    only_inner = False
    Only fix the inner part not the outer parts

 OUTPUTS:
    profile
    the modified profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fix_sbr(Configuration,Tirific_Template, smooth = False, debug=False):
    if debug:
        print_log(f'''FIX_SBR: Starting a SBR check.
''',Configuration['OUTPUTLOG'],debug=True)

    # get the cutoff limits
    vsys = float(Tirific_Template['VSYS'].split()[0])
    radii,cutoff_limits = sbr_limits(Configuration, systemic=vsys)
    cutoff_limits = np.array([cutoff_limits,cutoff_limits],dtype=float)
    # Then get the profile from the template
    sbr = np.array(get_from_template(Configuration,Tirific_Template,['SBR','SBR_2']),dtype=float)
    if debug:
        print_log(f'''FIX_SBR: Before modify.
{'':8s}sbr from template = {sbr}
''',Configuration['OUTPUTLOG'])
    # First make a correction on the inner 2 values
    sbr = inner_sbr_fix(Configuration,sbr,cutoff_limits,debug=debug)
    # Let's use a smoothed profile for the fittings
    sm_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',
                            min_error= cutoff_limits,no_apply = True,
                            fix_sbr_call = True,profile_in = sbr ,debug=debug)

    if debug:
        print_log(f'''FIX_SBR: Before modify.
{'':8s}sbr  = {sbr}
{'':8s}sm_sbr  = {sm_sbr}
''',Configuration['OUTPUTLOG'])

    # We interpolate negative values as well as values below the limits inner part with a cubi
    errors = get_error(Configuration,sbr,sm_sbr,'SBR',min_error=cutoff_limits,debug=debug)
    if debug:
        print_log(f'''FIX_SBR: retrieved errors.
{'':8s}errors  = {errors}
''',Configuration['OUTPUTLOG'])
    error_weights = cutoff_limits/sm_sbr*1./np.nanmin(cutoff_limits[:,2:]/sm_sbr[:,2:])
    error_weights[:,0] = 3.
    errors =errors*error_weights
    if debug:
        print_log(f'''FIX_SBR: weighed errors errors.
{'':8s}sm_sbr = {sm_sbr}
{'':8s}cutoff_limits  = {cutoff_limits}
{'':8s}errors  = {errors}
{'':8s}weights  = {error_weights}
''',Configuration['OUTPUTLOG'])


    store_gaussian = []
    for i in [0,1]:
        corr_val = np.where(sbr[i,2:] > cutoff_limits[i,2:])[0]+2
        # If we have enough safe values in the profile we attempt to fit it
        if corr_val.size > 3.:
            fit_sbr = sm_sbr[i,corr_val]
            if debug:
                print_log(f'''FIX_SBR: The values used for fitting are {fit_sbr}.
''',Configuration['OUTPUTLOG'])
            try:
                if 'SBR' in Configuration['FIXED_PARAMETERS'][0]:
                    vals = fit_gaussian(Configuration,radii[corr_val],fit_sbr,errors=errors[i,corr_val],debug=debug)
                    gaussian = gaussian_function(radii,*vals)
                else:
                    gaussian = fit_polynomial(Configuration,radii,sbr[i,:],sm_sbr[i,:],errors[i,:],'SBR', Tirific_Template,\
                                             min_error=cutoff_limits[i,:],debug= debug)
                # if the peak of this gaussian is in the inner two points replace it with the smoothed profile
                if np.any(np.where(np.max(gaussian) == gaussian)[0] < 2):
                    if debug:
                        print_log(f'''FIX_SBR: We are trying to replace the inner gaussian.
{'':8s} gaussian = {gaussian[[0,1]]}
{'':8s} sm_sbr = {sm_sbr[i,[0,1]]}
''',Configuration['OUTPUTLOG'])
                    gaussian[[0,1]] = sm_sbr[i,[0,1]]
            except RuntimeError:
                # If we fail we try a CubicSpline interpolation
                try:
                    tmp = CubicSpline(radii[corr_val],fit_sbr,extrapolate = True,bc_type ='natural')
                    gaussian = tmp(radii)
                except:
                    # and if that fails we just let the values be what they are in the smoothed profile
                    gaussian= sm_sbr[i,:]
            store_gaussian.append(gaussian)
        else:
            #if there are not enough points we simply use the smoothed profile
            store_gaussian.append(sm_sbr[i,:])

    store_gaussian = np.array(store_gaussian,dtype=float)
    if smooth:
        #If we smooth we take the fit in total
        sbr = store_gaussian
    else:
        sbr[np.where(sbr<cutoff_limits)] = store_gaussian[np.where(sbr<cutoff_limits)]
        sbr[:,[0,1,-1]] = store_gaussian[:,[0,1,-1]]

    # Need to make sure there are no nans
    sbr[np.isnan(sbr)] = 2.*cutoff_limits[np.isnan(sbr)]
    # and where we are lower than the cutoff we replace with the 1.2 *cutoff unless we smoothed
    if not smooth:
        sbr[np.where(sbr<cutoff_limits)] = 1.2*cutoff_limits[np.where(sbr<cutoff_limits)]
    else:
        sbr[np.where(sbr<cutoff_limits/2.)] = 1e-16
        #no rising outer end profiles
        for i in [0,1]:
            if sbr[i,-2] < sbr[i,-1]:
                last = sbr[i,-1]
                second = sbr[i,-2]
                sbr[i,-2] = last
                sbr[i,-1] = second

    # and ensure that both sides the inner two rings are the same
    sbr[:,[0,1]] = np.mean(sbr[:,[0,1]])

    if debug:
        print_log(f'''FIX_SBR: After modify.
{'':8s}{sbr}
''',Configuration['OUTPUTLOG'])
    Tirific_Template['SBR'] = f"{' '.join([f'{x:.2e}' for x in sbr[0]])}"
    Tirific_Template['SBR_2'] = f"{' '.join([f'{x:.2e}' for x in sbr[1]])}"
    print_log(f'''FIX_SBR: We checked the surface brightness profiles.
''',Configuration['OUTPUTLOG'])
fix_sbr.__doc__ =f'''
 NAME:
    fix_sbr

 PURPOSE:
    Correct the surface brightness profile against outliers.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = standard tirific template

 OPTIONAL INPUTS:
    debug = False

    smooth= False
    Normally only bad points (Not bright enough) and problematic inner points are corrected.
    When smooth is set the profile is either fitted  with a Gaussian function When FIX_SBR = True
    or fitted with a polynomial of n degree. When these fail the profile is
    interpolated with a cubic spline or smoothed with a savgol kernel.

 OUTPUTS:
    The template is corrected

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fix_vrot_for_incl_change(Configuration,Tirific_Template, incl_original,incl_modified,debug=False):
    vrot = np.array(get_from_template(Configuration,Tirific_Template,["VROT","VROT_2"]),dtype = float)
    new_vrot = copy.deepcopy(vrot)
    for i in [0,1]:
        vobs = np.array([x*np.sin(np.radians(y)) for x,y in zip(vrot[i],incl_original[i])],dtype=float)
        new_vrot[i] = np.array([x/np.sin(np.radians(y)) for x,y in zip(vobs,incl_modified[i])],dtype=float)
    vrot = np.array([np.mean([x,y]) for x,y in zip(new_vrot[0],new_vrot[1])],dtype=float)
    format = set_format("VROT")
    Tirific_Template["VROT"]= f"{' '.join([f'{x:{format}}' for x in vrot[:int(Configuration['NO_RINGS'])]])}"
    Tirific_Template[f"VROT_2"]= f"{' '.join([f'{x:{format}}' for x in vrot[:int(Configuration['NO_RINGS'])]])}"

fix_vrot_for_incl_change.__doc__ =f'''
 NAME:
    fix_vrot_for_incl_change

 PURPOSE:
    Correct the rotation curve when modifying the inclination

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = standard tirific template
    incl_original = profile before modification
    incl_modified = profile after modification

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    The template is corrected

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def flatten_the_curve(Configuration,Tirific_Template,debug = False):
    to_flatten = ['INCL','Z0','PA','SDIS']
    for key in to_flatten:
        profile = get_from_template(Configuration,Tirific_Template, [key,f'{key}_2'] ,debug=debug)
        new_profile = [np.mean(profile) for x in profile[0]]
        Tirific_Template[key] =f" {' '.join([f'{x:.2f}' for x in new_profile])}"
        Tirific_Template[f'{key}_2'] =f" {' '.join([f'{x:.2f}' for x in new_profile])}"
flatten_the_curve.__doc__ =f'''
 NAME:
    flatten_the_curve

 PURPOSE:
    Return radially varying profiles back to a mean value when the model becomes to small.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:
    the flattened profiles for PA, INLC, SDIS and Z0 are written to the template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_error(Configuration,profile,sm_profile,key,min_error = [0.],singular = False,weights= [1.],apply_max_error = False,debug=False):

    try:
        size= len(min_error)
        min_error = np.array(min_error,dtype=float)
    except TypeError:
        min_error = np.array([min_error],dtype=float)

    if debug:
        print_log(f'''GET_ERROR: starting;
{'':8s}original profile = {profile}
{'':8s}new profile = {sm_profile}
{'':8s}weights = {weights}
{'':8s}singular = {singular}
{'':8s}min_error = {min_error}
{'':8s}max_error = {apply_max_error}
''',Configuration['OUTPUTLOG'],debug = True)

    if singular:
        if len(weights) == 1:
            weights = np.full((len(profile)),weights[0])
        profile = [profile]
        sm_profile = [sm_profile]
        error =[[]]
        sides = [0]

    else:
        error = [[],[]]
        if len(weights) == 1:
            weights = np.full((len(profile[0]),len(profile[1])),1.)
        sides =[0,1]
    if debug:
        print_log(f'''GET_ERROR: using these weights =
{'':8s}{weights}
''',Configuration['OUTPUTLOG'])
    for i in sides:
        error[i] = abs(profile[i]-sm_profile[i])/2.
        error[i]= error[i]/weights[i]
        if len(min_error.shape) == 2:
            error[i] = [np.max([y,x]) for x,y in zip(error[i],min_error[i])]
        elif len(min_error) == len(error[i]):
            error[i] = [np.max([y,x]) for x,y in zip(error[i],min_error)]
        else:
            error[i] = [np.max([x,min_error[0]]) for x in error[i]]
        if apply_max_error:
            if len(Configuration['MAX_ERROR'][key]) == len(error[i]):
                error[i] = [np.nanmin([x,y]) for x,y in zip(error[i],Configuration['MAX_ERROR'][key])]
            else:
                error[i] = [np.nanmin([x,float(Configuration['MAX_ERROR'][key][0])]) for x in error[i]]
        if key in ['PA','INCL','Z0']:
            error[i][:Configuration['INNER_FIX'][i]+1] = [np.min(min_error) for x in error[i][:Configuration['INNER_FIX'][i]+1]]
    if singular:
        error = np.array(error[0],dtype=float)
    else:
        error = np.array(error,dtype=float)
    error[error == 0.] = np.min(error[error > 0.])

    if debug:
        print_log(f'''GET_ERROR: error =
{'':8s}{error}
''',Configuration['OUTPUTLOG'])
    return error
get_error.__doc__ =f'''
 NAME:
    get_error

 PURPOSE:
    get the errors associated with a profile

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    profile = the profile to consider
    sm_profile = the smoothed version of the profile
    key = the parameter to consider

 OPTIONAL INPUTS:
    debug = False

    min_error = [0.]
    Errors should always be bigger than this value

    singular = False
    If true a single parameter is expected not both sides

    weights= [1.]
    Weights to apply to the error withe most imporatant heighest, i.e. the error is divided by the weights

    apply_max_error = False
    aplly the maximum error values

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: the maximum errors are defined in main.py
'''

def get_number_of_rings(Configuration,sbr,sbr_ring_limits, debug=False):
    new_rings = Configuration['NO_RINGS']
    difference_with_limit = np.array(sbr-sbr_ring_limits,dtype=float)
    if np.all(difference_with_limit[:,-1] < 0.):
        if debug:
            print_log(f'''GET_NUMBER_OF_RINGS: both last rings are below the limit
''',Configuration['OUTPUTLOG'],debug=True)
        for i in range(len(difference_with_limit[0,:])-1,int(new_rings/2.),-1):
            if debug:
                print_log(f'''GET_NUMBER_OF_RINGS: Checking ring {i}
''',Configuration['OUTPUTLOG'])
            if np.all(difference_with_limit[:,i] < 0.):
                #check that 1 any of the lesser rings are bright enough
                if np.any(sbr[:,i-1] > 1.5 *sbr_ring_limits[:,i-1]):
                    new_rings = i+1
                    if debug:
                        print_log(f'''GET_NUMBER_OF_RINGS: we find that the previous rings are bright enough, rings = {new_rings}
''',Configuration['OUTPUTLOG'])
                    break
                else:
                    if debug:
                        print_log(f'''GET_NUMBER_OF_RINGS: the previous rings are not bright enough so we reduce 1, old_rings = {new_rings}, new_rings = {i}
''',Configuration['OUTPUTLOG'])
                    new_rings = i
            else:
                #if not both values are below than this is the extend we want
                new_rings = i+1
                if debug:
                    print_log(f'''GET_NUMBER_OF_RINGS: Not both rings warrant cutting, rings = {new_rings}
''',Configuration['OUTPUTLOG'])
                break
    else:
        if debug:
            print_log(f'''GET_NUMBER_OF_RINGS: Not both last rings are below the limit
''',Configuration['OUTPUTLOG'])
        # if they are not we first check wether both second to last rings are
        if ((difference_with_limit[0,-2] < 0.) and (sbr[0,-1] < 2*sbr_ring_limits[0,-1]) and (difference_with_limit[1,-1] < 0.)) or \
            ((difference_with_limit[1,-2] < 0.) and (sbr[1,-1] < 2*sbr_ring_limits[1,-1]) and (difference_with_limit[0,-1] < 0.)) or\
            ((difference_with_limit[0,-2] < 0.) and (difference_with_limit[1,-2] < 0.) and (sbr[0,-1] < 3*sbr_ring_limits[0,-1]) and (sbr[1,-1] < 3*sbr_ring_limits[1,-1])):
            new_rings -= 1
            if debug:
                print_log(f'''GET_NUMBER_OF_RINGS: A second ring is too faint, rings = {new_rings}
''',Configuration['OUTPUTLOG'])
        elif np.all(difference_with_limit[:,-2] < 0.) and np.all(sbr[:,-1] < 5*sbr_ring_limits[:,-1]):
            new_rings -= 1
            if debug:
                print_log(f'''GET_NUMBER_OF_RINGS: Both second rings are too faint, rings = {new_rings}
''',Configuration['OUTPUTLOG'])
        else:
            if debug:
                print_log(f'''GET_NUMBER_OF_RINGS: The second rings are too bright and do not allow for cutting.
''',Configuration['OUTPUTLOG'])
    return new_rings
get_number_of_rings.__doc__ =f'''
 NAME:
    get_number_of_rings

 PURPOSE:
    Determine whether the amount of rings is good for the limits or not should change or not

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    sbr = sbr profiles
    sbr_ring_limits = the limits to evaluate

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    new_rings = the required number of rings

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_warp_slope(Configuration,Tirific_Template, debug = False):
    if debug:
        print_log(f'''GET_WARP_SLOPE: We have {Tirific_Template['NUR']} rings in the template. and this should be {Configuration['NO_RINGS']}
''', Configuration['OUTPUTLOG'],debug = True)
    radii, sbr_ring_limits = sbr_limits(Configuration,systemic = float(Tirific_Template['VSYS'].split()[0]),debug=debug)
    #get the sbr profiles
    sbr = np.array(get_from_template(Configuration,Tirific_Template, ['SBR','SBR_2'],debug=debug),dtype = float)
    if debug:
        print_log(f'''GET_WARP_SLOPE: We have {len(sbr_ring_limits)} rings in our limits.
{'':8s}GET_WARP_SLOPE: And we have {len(sbr[0])} rings in our profiles.
''', Configuration['OUTPUTLOG'])
    warp_slope = [Tirific_Template['NUR'],Tirific_Template['NUR']]
    sbr_ring_limits = 2.*np.array([sbr_ring_limits,sbr_ring_limits],dtype=float)
    difference_with_limit = np.array(sbr-sbr_ring_limits,dtype=float)
    for i in [0,1]:
        slope = difference_with_limit[i]
        final = slope[slope < 0.]
        if len(slope) > len(final) > 0.:
            not_found = True
            counter = len(slope)-1
            while not_found:
                if not slope[counter] < 0.:
                    not_found = False
                    final = counter+1
                else:
                    counter -= 1
        elif len(final) == len(slope):
            final = 1.
            for parameter in ['INCL',"PA",'SDIS','Z0']:
                if parameter not in Configuration['FIXED_PARAMETERS'][0]:
                    Configuration['FIXED_PARAMETERS'][0].append(parameter)
        else:
            final = len(slope)
        if final > Configuration['LAST_RELIABLE_RINGS'][i]:
            final = Configuration['LAST_RELIABLE_RINGS'][i]
        warp_slope[i] = final
    if debug:
        print_log(f'''GET_WARP_SLOPE: We find a slope of {warp_slope}.
''', Configuration['OUTPUTLOG'])
    Configuration['WARP_SLOPE'] = warp_slope
    incl = np.array(get_from_template(Configuration,Tirific_Template, ['INCL','INCL_2'],debug=debug),dtype = float)
    if np.mean(incl[:,:int(Configuration['NO_RINGS']/2.)]) < 35. :
        if 'INCL' not in Configuration['FIXED_PARAMETERS'][0]:
            Configuration['FIXED_PARAMETERS'][0].append('INCL')

get_warp_slope.__doc__ =f'''
 NAME:
    get_warp_slope

 PURPOSE:
    Get the rings where the warp should be sloped

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def inner_sbr_fix(Configuration,sbr,cutoff_limits,debug=False):
    if debug:
        print_log(f'''INNER_SBR_FIX: Checking the SBR inner points for runaway values
{'':8s} sbr in  = {sbr}
''',Configuration['OUTPUTLOG'], debug = True)

    if np.all(sbr[:,0] > 2*sbr[:,2]) or np.all(sbr[:,1] > 2*sbr[:,2]):
        if debug:
            print_log(f'''INNER_SBR_FIX: We need to correct
{'':8s} sbr 0  = {sbr[:,0]} sbr 1  = {sbr[:,1]} sbr 2  = {sbr[:,2]}
''',Configuration['OUTPUTLOG'])
        if np.mean(sbr[:,2]) > cutoff_limits[0,2]:
            sbr[:,[0,1]] = np.mean(sbr[:,2])
            if debug:
                print_log(f'''INNER_SBR_FIX: We need to correct with mean
{'':8s} mean = {np.mean(sbr[:,2])}
''',Configuration['OUTPUTLOG'])
        else:
            if debug:
                print_log(f'''INNER_SBR_FIX: We need to correct with cut_off_limits
{'':8s} limit = {1.5*cutoff_limits[0,2]}
''',Configuration['OUTPUTLOG'])
            sbr[:,[0,1,2]] = 1.5*cutoff_limits[0,2]
    if np.any(sbr[:,0] > sbr[:,1]):
        if debug:
                print_log(f'''INNER_SBR_FIX: We correct 0 point
{'':8s} mean 1 = {np.mean(sbr[:,1])}
''',Configuration['OUTPUTLOG'])
        sbr[:,0] = np.mean(sbr[:,1])

    for i in [0,1]:
        if np.any(sbr[:,i] < cutoff_limits[:,2]):
            if debug:
                print_log(f'''INNER_SBR_FIX: correcting ring {i}
''',Configuration['OUTPUTLOG'])
            sbr[:,i] = 1.5*cutoff_limits[0,2]

    if debug:
                print_log(f'''INNER_SBR_FIX: the fixed sbr {sbr}
''',Configuration['OUTPUTLOG'])

    return sbr
inner_sbr_fix.__doc__ =f'''
 NAME:
    inner_sbr_fix

 PURPOSE:
    Make sure the inner two points of the SBR are not run away

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    sbr = SBR profile for bothe sides
    cutoff_limits = The reliability limits of the fit

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    sbr = profile with the modified inner points if required

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def modify_flat(Configuration,profile,original_profile,errors,key,inner_fix = [4,4], debug=False):

    if debug:
         print_log(f'''MODIFY_FLAT: These {key} profiles are checked to be flat.
{'':8s} profile = {profile}
{'':8s} original_profile = {original_profile}
{'':8s} errors = {errors}
''',Configuration['OUTPUTLOG'],debug = True)

    flatness = []


    for side in [0,1]:
         flatness.append(check_flat(Configuration,profile[side],errors[side],key,inner_fix=inner_fix[side]\
                                    ,last_reliable_ring= Configuration['LAST_RELIABLE_RINGS'][side],debug=debug))
    if debug:
         print_log(f'''MODIFY_FLAT: Side 0 is flat = {flatness[0]}
{'':8s} Side 1 is flat = {flatness[1]}
''',Configuration['OUTPUTLOG'])

    if all(flatness):
        if key in ['PA','INCL']:
            profile[:] = profile[0,0]
        #elif key in ['INCL'] and 30. < np.nanmedian(original_profile[:,-4:]) < 50.:
        #    profile[:] = np.nanmedian(original_profile[:,:])
        #elif key in ['INCL'] and np.nanmedian(original_profile[:,-4:]) <= 30.:
        #    profile[:] = np.nanmedian(original_profile[:,:])
        else:
            profile[:] = np.nanmedian(original_profile[:,:round(len(original_profile)/2.)])
        errors = get_error(Configuration,original_profile,profile,key,apply_max_error = True,min_error =np.nanmin(errors) , debug=debug)
    else:
        if any(flatness):
            if key not in ['SDIS']:
                for side in [0,1]:
                    if flatness[side]:
                        if key in ['PA','INCL']:
                            profile[side,:] = profile[side,0]
                        #elif key in ['INCL'] and 30.< np.nanmedian(original_profile[side,round(len(original_profile)/2.):]) < 50:
                        #    profile[side,:] = np.nanmedian(original_profile[side,:])
                        #elif key in ['INCL'] and 30.< np.nanmedian(original_profile[side,round(len(original_profile)/2.):]) <= 30.:
                        #    profile[side,:] = np.nanmedian(original_profile[side,:])
                        else:
                            profile[side,:]  = np.median(original_profile[side,:round(len(original_profile)/2.)])

                        flat_val = profile[side,0]
                        errors[side] = get_error(Configuration,original_profile[side],profile[side],key,apply_max_error = True,min_error =np.nanmin(errors[side]),singular = True ,debug=debug)
                profile[:,0:3] = flat_val
                profile[:,4] = (flat_val+profile[:,4])/2.
            else:
                profile[:] = np.nanmedian(original_profile[:,:round(len(original_profile)/2.)])
                errors = get_error(Configuration,original_profile,profile,key,apply_max_error = True,min_error =np.nanmin(errors) , debug=debug)

    if debug:
        print_log(f'''MODIFY_FLAT: Returning:
{'':8s} profile = {profile}
{'':8s} errors = {errors}
''',Configuration['OUTPUTLOG'])
    return profile,errors
modify_flat.__doc__ =f'''
 NAME:
    modify_flat

 PURPOSE:
    Check if a profile should be flat within its errors and if so make it flat

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    profile = the profile to check
    original_profile = the original unmodified profile
    errors = the errors on the profile
    key = the parameter to be checked

 OPTIONAL INPUTS:
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:
    the final profile and errors

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def no_declining_vrot(Configuration, Tirific_Template, profile = None, debug = False):
    if debug:
        print_log(f'''NO_DECLINING_VROT: make RC flat from highest point on.
{'':8s}NO_DECLINING_VROT: But only for low value RCs
''',Configuration['OUTPUTLOG'],debug = True)
    no_input = False
    if profile is None:
        no_input = True
        profile = np.array(get_from_template(Configuration,Tirific_Template,['VROT'], debug = debug)[0],dtype = float)

    RCval = np.mean(profile[2:])
    RCmax = np.where(profile == np.max(profile))[0]
    if len(RCmax) > 1:
        RCmax = RCmax[0]
    if debug:
            print_log(f'''NO_DECLINING_VROT: We find the maximum at ring {RCmax}
{'':8s}NO_DECLINING_VROT: And a mean value of {RCval}.
''',Configuration['OUTPUTLOG'])
    Configuration['OUTER_SLOPE_START'] = Configuration['NO_RINGS']-1
    if RCmax < len(profile)/2. or RCval > 180.:
        if debug:
            print_log(f'''NO_DECLINING_VROT: We shan't adapt the RC
''',Configuration['OUTPUTLOG'])
    else:
        for i in range(int(len(profile)/2.),len(profile)-1):
            if profile[i+1] < profile[i]:
                profile[i:] =profile[i]
                if debug:
                    print_log(f'''NO_DECLINING_VROT: Flattening from ring {i} on.)
    ''',Configuration['OUTPUTLOG'])
                Configuration['OUTER_SLOPE_START'] = i+2
                break

    #and we check that the last parts are not declining too much in anycase
    # For galaxies with more than 10 rings let's make sure the last quarter is not declinining steeply
    if Configuration['NO_RINGS'] > 10:
        for i in range(int(Configuration['NO_RINGS']*3./4),Configuration['NO_RINGS']-1):
            if profile[i+1] < profile[i]*0.85:
                profile[i+1] = profile[i]*0.85
    else:
        if profile[-1] < profile[-2]*0.8:
            profile[-1] = profile[-2]*0.8


    if Configuration['OUTER_SLOPE_START'] > Configuration['NO_RINGS']:
        Configuration['OUTER_SLOPE_START'] = Configuration['NO_RINGS']
    format = set_format('VROT')
    if no_input:
        Tirific_Template['VROT'] = f"{' '.join([f'{x:{format}}' for x in profile])}"
        Tirific_Template['VROT_2'] = f"{' '.join([f'{x:{format}}' for x in profile])}"
    else:
        return profile
no_declining_vrot.__doc__ =f'''
 NAME:
    no_declining_vrot

 PURPOSE:
    Ensure that the RC is not declining in an unphysical manner, i.e. if the maximum lies in the outer rings it should be the last ring.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard tirific template

 OPTIONAL INPUTS:
    debug = False

    profile = None
    Normally the RC is read from the Template however if profile is set it is taken from there.
    This allows for the profile to be modified before it is checked.

 OUTPUTS:
    If profile was set the modified profile is returned else the new profile is written to the template which return modified.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def regularise_profile(Configuration,Tirific_Template, key,min_error= [0.],debug = False, no_apply =False):
    # We start by getting an estimate for the errors
    min_error=np.array(min_error,dtype=float)
    profile = np.array(get_from_template(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype=float)
    weights = get_ring_weights(Configuration,Tirific_Template,debug=debug)


    #First if we have an RC we flatten the curve
    if debug:
        print_log(f'''REGULARISE_PROFILE: profile of {key} before regularistion
{'':8s}{profile[0]}
{'':8s}{profile[1]}
{'':8s}The minmal error is
{'':8s}{min_error}
''',Configuration['OUTPUTLOG'],debug = True)
    # get a smoothed profiles
    sm_profile = smooth_profile(Configuration,Tirific_Template, key ,min_error=min_error,debug=debug,no_apply=True)

    error = get_error(Configuration,profile,sm_profile,key,min_error=min_error,weights = weights,debug=debug)

    #Check that we have two profiles
    diff = np.sum(profile[0]-profile[1])
    if key in ['SDIS','VROT']:
        diff = False
    if diff:
        if debug:
            print_log(f'''REGULARISE_PROFILE: Treating both sides independently.
''',Configuration['OUTPUTLOG'])
        sides = [0,1]
    else:
        if debug:
            print_log(f'''REGULARISE_PROFILE: Found symmetric profiles.
''',Configuration['OUTPUTLOG'])
        sides = [0]
        error[0] = np.array([np.mean([x,y]) for x,y in zip(error[0],error[1])],dtype=float)

    radii =set_rings(Configuration,debug=debug)
    for i in sides:

        if key == 'SDIS':
            try:
                fit_profile = fit_arc(Configuration,radii,sm_profile[i],error[i],debug= debug)
            except:
                fit_profile = np.full(len(sm_profile[i]), np.mean(sm_profile[i]))


        else:
            if f"{key}_CURRENT_BOUNDARY" in Configuration:
                boundary = Configuration[f"{key}_CURRENT_BOUNDARY"][i+1]
            else:
                boundary = [0.,0.]
            fit_profile = fit_polynomial(Configuration,radii,profile[i],sm_profile[i],error[i],key, Tirific_Template,\
                                         inner_fix = Configuration['INNER_FIX'][i],min_error=min_error,\
                                         boundary_limits= boundary,debug= debug)
        profile[i] = fit_profile

    if not diff:
        profile[1] = profile[0]

    original = np.array(get_from_template(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype=float)
    error = get_error(Configuration,original,profile,key,weights=weights,apply_max_error = True,min_error=min_error,debug=debug)
    if debug:
            print_log(f'''REGULARISE_PROFILE: This the fitted profile without corrections:
{'':8s}{profile}
''',Configuration['OUTPUTLOG'])
#then we want to fit the profiles with a polynomial
    if key not in ['SBR','VROT','SDIS']:
        #We should not fix the profile again as the fitted profile is fixed should be good

        profile,error = modify_flat(Configuration, profile, original, error,key,inner_fix= Configuration['INNER_FIX'],debug=debug)


    format = set_format(key)

    if not no_apply:
        Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in profile[0,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:{format}}' for x in profile[1,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in error[0,:int(Configuration['NO_RINGS'])]])}")
        Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in error[1,:int(Configuration['NO_RINGS'])]])}")
        #if key in ['INCL'] and np.mean( profile[:,int(Configuration['NO_RINGS']/2.):int(Configuration['NO_RINGS'])]) < 40.:
        #    fix_vrot_for_incl_change(Configuration,Tirific_Template,original,profile,debug=debug)

        if debug:
            print_log(f'''REGULARISE_PROFILE: And this has gone to the template.
{'':8s}{key} = {Tirific_Template[key]}
{'':8s}{key}_2 ={Tirific_Template[f"{key}_2"]}
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration['OUTPUTLOG'])
    return profile
regularise_profile.__doc__ =f'''
 NAME:
    regularise_profile

 PURPOSE:
    Regularise a parameter profile with a polynomial or a arctan when it is the SDIS

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    key = parameter to fix

 OPTIONAL INPUTS:
    debug = False

    no_apply = false
    if true do not apply the regularised profile to the template

    min_error = [0.]
    error should alway be larger than this

 OUTPUTS:
    the regularised profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Errors are not returned but they are written to the template
'''


def set_boundary_limits(Configuration,Tirific_Template,key, tolerance = 0.01, values = [10,1],upper_bracket = [10.,100.], lower_bracket=[0., 50.], fixed = False,increase=10., debug = False):
    if debug:
        print_log(f'''SET_BOUNDARY_LIMITS: checking limits for {key},
''',Configuration['OUTPUTLOG'],debug = True)
    profile = np.array(get_from_template(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype = float)

    current_boundaries = Configuration[f"{key}_CURRENT_BOUNDARY"]
    if debug:
        print_log(f'''SET_BOUNDARY_LIMITS: We have found the following limits,
{'':8s} current Boundaries = {current_boundaries}
{'':8s} Applying to the profiles {profile}
''',Configuration['OUTPUTLOG'])

    if np.sum(current_boundaries) == 0.:
        current_boundaries =[[set_limits(values[0]-values[1]*5.,*lower_bracket),\
                    set_limits(values[0]+values[1]*5.,*upper_bracket)] for x in range(3)]
        if debug:
            print_log(f'''SET_BOUNDARY_LIMITS: We have set the boundaries to the following as all were 0.
{'':8s} current Boundaries = {current_boundaries}
''',Configuration['OUTPUTLOG'])

    if fixed:
        range_to_check = [0]
    else:
        range_to_check = [0,1,2]
    for i in range_to_check:
        buffer = float(current_boundaries[i][1]-current_boundaries[i][0]) * tolerance
        if debug:
            print_log(f'''SET_BOUNDARY_LIMITS: Using a buffer of {buffer}.
''',Configuration['OUTPUTLOG'])
        if i == 0:
            profile_part = profile[0,:int(np.mean(Configuration['INNER_FIX']))+1]
        else:
            profile_part = profile[i-1,Configuration['INNER_FIX'][i-1]:]
        if debug:
            print_log(f'''SET_BOUNDARY_LIMITS: Checking {profile_part}.
''',Configuration['OUTPUTLOG'])
        #check the upper bounderies
        on_boundary = np.where(profile_part > float(current_boundaries[i][1])-buffer)[0]
        if debug:
            print_log(f'''SET_BOUNDARY_LIMITS: Found the following on the upper {on_boundary}.
''',Configuration['OUTPUTLOG'])
        if len(on_boundary) > 0:
            if on_boundary[0] != len(profile[0])-1:
                current_boundaries[i][1] = set_limits(current_boundaries[i][1] + buffer*increase,*upper_bracket)
        #check the lower boundaries.
        on_boundary = np.where(profile_part < float(current_boundaries[i][0])+buffer)[0]
        if debug:
            print_log(f'''SET_BOUNDARY_LIMITS: Found the following on the lower {on_boundary}.
''',Configuration['OUTPUTLOG'])
        if len(on_boundary) > 0:
            if on_boundary[0] != len(profile[0])-1:
                current_boundaries[i][0] = set_limits(current_boundaries[i][0] - buffer*increase,*lower_bracket)
    Configuration[f"{key}_CURRENT_BOUNDARY"] = current_boundaries
    if debug:
        print_log(f'''SET_BOUNDARY_LIMITS: We have adjusted the boundaries to  {Configuration[f"{key}_CURRENT_BOUNDARY"]}.
''',Configuration['OUTPUTLOG'])
    return current_boundaries
set_boundary_limits.__doc__ =f'''
 NAME:
    set_boundary_limits

 PURPOSE:
    Update the boundary limits of parameters when too many rings are on the boundary.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template
    key = parameter to check

 OPTIONAL INPUTS:
    debug = False

    tolerance = 0.01,
    The minimum difference with the boundary

    values = [10,1],upper_bracket = [10.,100.], lower_bracket=[0., 50.]
    Initial definition of the boundary limits

    fixed = False
    If true it is assumed there is only one boundary value to check else three

    increase=10.
    how much to widen the boundary by

 OUTPUTS:
    The new Boundaries, they are also updated in Configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_cflux(Configuration,Tirific_Template,debug = False):

    if any(np.isnan(Configuration['NO_POINTSOURCES'])):
        print_log(f'''SET_CFLUX: We detected an infinite number of model point sources.
{"":8s}SET_CFLUX: This must be an error. Exiting the fitting.
''',Configuration['OUTPUTLOG'],screen = True)
        raise CfluxError('The model had infinite point sources')
    if Configuration['SIZE_IN_BEAMS'] < 15:
        factor = 1.
    else:
        factor=(Configuration['SIZE_IN_BEAMS']/15.)**1.5
    triggered = 0
    if not 0.5e6 < Configuration['NO_POINTSOURCES'][0] < 2.2e6:
        new_cflux = set_limits(float(Tirific_Template['CFLUX'])*Configuration['NO_POINTSOURCES'][0]/(factor*1e6),1e-7,5e-3)
        print_log(f'''SET_CFLUX: CFLUX is adapted from {Tirific_Template['CFLUX']} to {new_cflux:.2e}
''',Configuration['OUTPUTLOG'])
        Tirific_Template['CFLUX'] = f"{new_cflux:.2e}"
        triggered = 1
    if not 0.5e6 < Configuration['NO_POINTSOURCES'][1] < 2.2e6:
        new_cflux = set_limits(float(Tirific_Template['CFLUX_2'])*Configuration['NO_POINTSOURCES'][1]/(factor*1e6),1e-7,5e-3)
        print_log(f'''SET_CFLUX: CFLUX_2 is adapted from {Tirific_Template['CFLUX_2']} to {new_cflux:.2e}
''',Configuration['OUTPUTLOG'])
        Tirific_Template['CFLUX_2'] = f"{new_cflux:.2e}"
        triggered = 1
    if not triggered:
        print_log(f'''SET_CFLUX: CFLUXES are within the required limits.
''',Configuration['OUTPUTLOG'])
set_cflux.__doc__ =f'''
 NAME:
    set_cflux

 PURPOSE:
    Check CFLUX values and make sure they are in the right order for the amount of point sources

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Updated CFLUX values in template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_errors(Configuration,Tirific_Template,key,min_error = 0.,debug = False):
    error = np.full((2,int(Configuration['NO_RINGS'])),min_error)
    format=set_format(key)
    Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in error[0,:int(Configuration['NO_RINGS'])]])}")
    Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in error[1,:int(Configuration['NO_RINGS'])]])}")
    if debug:
        print_log(f'''SET_ERRORS: This has gone to the template.
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration['OUTPUTLOG'],debug = True)
set_errors.__doc__ =f'''
 NAME:
    set_errors

 PURPOSE:
    Write the errors for flat profiles to the template

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard tirific template
    key = parameter to write the errors for

 OPTIONAL INPUTS:
    debug = False

    min_error =0 .
    the error to be used for all rings

 OUTPUTS:
    Updated Tirific Template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_fitting_parameters(Configuration, Tirific_Template, parameters_to_adjust  = ['NO_ADJUSTMENT'], modifiers = ['EMPTY'], stage = 'initial', initial_estimates = ['EMPTY'],debug = False):
    if debug:
        print_log(f'''SET_FITTING_PARAMETERS: We are starting with these modifiers.
{'':8s} {modifiers}
''',Configuration['OUTPUTLOG'],debug=True)
    try:
        if modifiers[0] == 'EMPTY':
            modifiers = {}
    except KeyError:
        pass
    try:
        if initial_estimates[0] == 'EMPTY':
            initial_estimates = {}
    except KeyError:
        pass
    fitting_settings = {}
    fitting_keys = ['VARY','VARINDX','MODERATE','DELEND','DELSTART','MINDELTA','PARMAX','PARMIN']

    if 'INCL' not in initial_estimates:
        profile = np.array([np.mean([x,y]) for x,y in \
                    zip(get_from_template(Configuration,Tirific_Template, ['INCL']),\
                    get_from_template(Configuration,Tirific_Template, [f"INCL_2"]) )],dtype=float)
        diff = abs(np.max(profile)-np.min(profile))/10.
        initial_estimates['INCL'] = [profile[0],set_limits(diff,1,5)/np.sin(np.radians(profile[0]))]

    if parameters_to_adjust[0] == 'NO_ADJUSTMENT':
        if stage in ['initial','run_cc','after_cc']:
            if initial_estimates['INCL'][0] < 30.:
                parameters_to_adjust = ['VSYS','SBR','XPOS','YPOS','PA','SDIS','INCL','VROT']
            elif initial_estimates['INCL'][0] < 50.:
                parameters_to_adjust = ['VSYS','XPOS','YPOS','SBR','PA','SDIS','VROT','INCL']
            elif initial_estimates['INCL'][0] > 75.:
                parameters_to_adjust = ['VSYS','XPOS','YPOS','VROT','SBR','PA','INCL','SDIS','Z0']
            else:
                parameters_to_adjust = ['VSYS','XPOS','YPOS','SBR','VROT','PA','INCL','SDIS']
        elif stage in ['initialize_ec','run_ec','after_ec']:
            parameters_to_adjust = ['INCL','PA','VROT','SDIS','SBR','Z0','XPOS','YPOS','VSYS']
        elif stage in  ['initialize_os','run_os','after_os']:
            parameters_to_adjust = ['VSYS','XPOS','YPOS','INCL','PA','SBR','VROT','SDIS','Z0']
        else:
            if debug:
                print_log(f'''SET_FITTING_PARAMETERS: No default adjustment for unknown stage.
''',Configuration['OUTPUTLOG'])
            raise InitializeError('No default adjustment for unknown stage. ')

    for key in parameters_to_adjust:
        if key not in initial_estimates:
            profile = np.array([np.mean([x,y]) for x,y in \
                                zip(get_from_template(Configuration,Tirific_Template, [key])[0],\
                                get_from_template(Configuration,Tirific_Template, [f"{key}_2"])[0] )],dtype=float)
            diff = abs(np.max(profile)-np.min(profile))/10.
            if key == 'PA':
                initial_estimates['PA'] = [profile[0],set_limits(diff,0.5,10)]
            elif key == 'VROT':
                initial_estimates['VROT'] = [np.nanmax(profile),set_limits(np.nanstd(profile[1:]),np.nanmax(profile)-np.nanmin(profile[1:]),np.nanmax(profile))]
            elif key == 'XPOS':
                initial_estimates['XPOS'] = [profile[0],Configuration['BEAM_IN_PIXELS'][1]*Configuration['PIXEL_SIZE']]
            elif key == 'YPOS':
                initial_estimates['YPOS'] = [profile[0],Configuration['BEAM_IN_PIXELS'][1]*Configuration['PIXEL_SIZE']]
            elif key == 'VSYS':
                initial_estimates['VSYS'] = [profile[0],Configuration['CHANNEL_WIDTH']/2.]
            elif key == 'SDIS':
                initial_estimates['SDIS'] = [np.mean(profile),Configuration['CHANNEL_WIDTH']]
            elif key == 'Z0':
                initial_estimates['Z0'] = convertskyangle(Configuration,[0.2,0.05],Configuration['DISTANCE'], physical = True)
        if key not in modifiers:
            if debug:
                print_log(f'''SET_FITTING_PARAMETERS: Adding {key} to the modifiers
''',Configuration['OUTPUTLOG'])
            if key == 'Z0': modifiers['Z0'] = [0.5,0.5,2.]
            elif key in ['XPOS','YPOS']: modifiers[key] = [1.,1.,1.]
            elif key == 'VSYS': modifiers[key] = [2.,0.5,0.1]
            elif stage in  ['initial','run_cc','after_cc','after_ec','after_os','final_os']:
                if key == 'INCL': modifiers['INCL'] = [1.,1.,1.]
                elif key == 'PA': modifiers['PA'] = [1.,1.,1.]
                elif key == 'SDIS': modifiers['SDIS'] =  [1.,1.,2.]
            else:
                if key == 'INCL': modifiers['INCL'] = [2.0,0.5,0.5]
                elif key == 'PA': modifiers['PA'] =[1.0,1.0,2.0]
                elif key == 'SDIS': modifiers['SDIS'] =  [1.,1.,0.5]
    if debug:
        for key in modifiers:
            print_log(f'''SET_FITTING_PARAMETERS: This {key} is in modifiers
{'':8s} With these values {modifiers[key]}
''',Configuration['OUTPUTLOG'])
    if stage not in ['final_os']:
        if initial_estimates['INCL'][0] < 30.:
            if 'Z0' in modifiers:
                modifiers['Z0'][0:1] = np.array(modifiers['Z0'][0:1],dtype=float)*(0.2/0.5)
                modifiers['Z0'][2] = float(modifiers['Z0'][2])*1.5
            if 'INCL' in modifiers:
                modifiers['INCL'][0:1] =np.array(modifiers['INCL'][0:1],dtype=float)*(0.1/1.0)
                modifiers['INCL'][2] = float(modifiers['INCL'][2])*2.
            if 'SDIS' in modifiers:
                modifiers['SDIS'][0:1] = np.array(modifiers['SDIS'][0:1],dtype=float)*1.5
                modifiers['SDIS'][2] = float(modifiers['SDIS'][2])*0.5
            if debug:
                print_log(f'''SET_FITTING_PARAMETERS: These are the  modifiers after correcting < 30.
{'':8s} {modifiers}
''',Configuration['OUTPUTLOG'])
        elif initial_estimates['INCL'][0] < 50.:
            if 'Z0' in modifiers:
                modifiers['Z0'][0:1] = np.array(modifiers['Z0'][0:1],dtype=float)*(0.4/0.5)
                modifiers['Z0'][2] = float(modifiers['Z0'][2])*1.2
            if 'INCL' in modifiers:
                modifiers['INCL'][0:1] =np.array( modifiers['INCL'][0:1],dtype=float)*(0.5/1.0)
                modifiers['INCL'][2] = float(modifiers['INCL'][2])*1.5
            if debug:
                print_log(f'''SET_FITTING_PARAMETERS: These are the  modifiers after correcting < 50.
{'':8s} {modifiers}
''',Configuration['OUTPUTLOG'])
        elif initial_estimates['INCL'][0] > 75.:
            if 'Z0' in modifiers:
                modifiers['Z0'][0:1] = np.array(modifiers['Z0'][0:1],dtype=float)*(1.25)
                modifiers['Z0'][2] = float(modifiers['Z0'][2])*0.5
            if 'INCL' in modifiers:
                modifiers['INCL'][0:1] = np.array(modifiers['INCL'][0:1],dtype=float)*(1.5/1.0)
                modifiers['INCL'][2] = float(modifiers['INCL'][2])*0.25
            if 'SDIS' in modifiers:
                modifiers['SDIS'][0:1] = np.array(modifiers['SDIS'][0:1],dtype=float)*0.5
                modifiers['SDIS'][2] = float(modifiers['SDIS'][2])*1.5
            if debug:
                print_log(f'''SET_FITTING_PARAMETERS: These are the modifiers after correcting > 75.
{'':8s} {modifiers}
''',Configuration['OUTPUTLOG'])
        else:
            pass
    else:
        if initial_estimates['INCL'][0] < 40.:
            if 'INCL' in modifiers:
                modifiers['INCL'] = np.array(modifiers['INCL'],dtype=float)*(0.2)

    for key in parameters_to_adjust:
        if key == 'VROT':
            fitting_settings['VROT'] = set_vrot_fitting(Configuration,stage = stage, rotation = initial_estimates['VROT'], debug = debug )
        elif key == 'SBR':
            fitting_settings['SBR'] = set_sbr_fitting(Configuration,stage = stage, systemic = initial_estimates['VSYS'][0], debug = debug)
        else:
            if key in Configuration['FIXED_PARAMETERS'][0]:
                fixed=True
            else:
                fixed =False
            if key == 'INCL':
                brackets = [[60.,90.],[5.,50.]]
            elif key == 'PA':
                brackets = [[190.,370.],[-10,170]]
            elif key == 'Z0':
                brackets = [convertskyangle(Configuration,[0.2,2.5],Configuration['DISTANCE'], physical = True), \
                            convertskyangle(Configuration,[0.05,0.2],Configuration['DISTANCE'], physical = True)]
            elif key in ['XPOS','YPOS','VSYS']:
                if key == 'XPOS': i = 0
                elif key == 'YPOS': i = 1
                elif key == 'VSYS': i = 2
                brackets = [Configuration['NAXES_LIMITS'][i],Configuration['NAXES_LIMITS'][i]]
            elif key == 'SDIS':
                if stage in ['initial','run_cc','after_cc','after_ec','after_os']:
                    limits = [[Configuration['CHANNEL_WIDTH'], 16.],\
                              [Configuration['CHANNEL_WIDTH']/4., 16.], \
                              [Configuration['CHANNEL_WIDTH']/4., 16.]]
                else:
                    limits = [[Configuration['CHANNEL_WIDTH'], set_limits(initial_estimates['SDIS'][0]*2.,Configuration['CHANNEL_WIDTH'],25.)], \
                             [Configuration['CHANNEL_WIDTH']/4., set_limits(initial_estimates['SDIS'][0]*2.,Configuration['CHANNEL_WIDTH']/2.,25.)],\
                             [Configuration['CHANNEL_WIDTH']/4., set_limits(initial_estimates['SDIS'][0]*2.,Configuration['CHANNEL_WIDTH']/2.,25.)]]
            if key in ['PA','INCL','Z0']:
                slope = Configuration['WARP_SLOPE']
                if stage in ['initialize_os']:
                    inner = int(set_limits(Configuration['NO_RINGS']*1./3., 3,Configuration['NO_RINGS']-2 ))
                else:
                    inner =  Configuration['INNER_FIX']
            elif key in ['SDIS']:
                flat_slope = True

                inner = int(set_limits(Configuration['NO_RINGS']*0.33,4,Configuration['NO_RINGS']*0.5))
                slope = [int(Configuration['NO_RINGS']*0.66),int(Configuration['NO_RINGS']*0.66)]
            else:
                inner = 4
                slope = [0.,0.]

            if key in ['SDIS','INCL']:
                fact = set_limits(float(initial_estimates['INCL'][0])/20.-2.5,1.,2.)
                try:
                    inner[0] = int(set_limits(inner[0]*fact,4,Configuration['NO_RINGS']/2.))
                    inner[1] = int(set_limits(inner[1]*fact,4,Configuration['NO_RINGS']/2.))
                except:
                    inner = int(set_limits(inner*fact,4,Configuration['NO_RINGS']/2.))

            if key != 'SDIS':
                limits = set_boundary_limits(Configuration,Tirific_Template,key, tolerance = 0.1, values = initial_estimates[key],\
                                upper_bracket = brackets[0],lower_bracket = brackets[1],fixed = fixed,debug=debug)
                flat_slope = False
                symmetric = False
            else:
                symmetric = True


            fitting_settings[key] =  set_generic_fitting(Configuration,key,stage = stage, values = initial_estimates[key],\
                                                        debug = debug, limits=limits,slope= slope, flat_slope = flat_slope,\
                                                         fixed =fixed, flat_inner = inner, step_modifier = modifiers[key], \
                                                         symmetric = symmetric)

    # Reset the fitting parameters
    for fit_key in fitting_keys:
        Tirific_Template[fit_key]= ''
    #write the new parameters
    for key in parameters_to_adjust:
        if key in fitting_settings:
            for fit_key in fitting_keys:
                if  fit_key in fitting_settings[key]:

                    if fit_key in ['DELSTART','DELEND','MINDELTA']:
                    #if fit_key in ['DELEND','MINDELTA']:
                        # These should never be 0.
                        format = set_format(key)
                        for i,x in enumerate(fitting_settings[key][fit_key]):
                            while float(f'{fitting_settings[key][fit_key][i]:{format}}') == 0.:
                                if float(fitting_settings[key][fit_key][i]) == 0.:
                                    fitting_settings[key][fit_key][i] += 0.01
                                else:
                                    fitting_settings[key][fit_key][i] *= 2.

                    if fit_key == 'VARY':
                        if len(Tirific_Template[fit_key]) == 0:
                            Tirific_Template[fit_key] = ', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])
                        else:
                            Tirific_Template[fit_key] = f"{Tirific_Template[fit_key]}, {', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])}"
                    else:
                        if fit_key == 'VARINDX':
                            format = '<10s'
                        else:
                            format = set_format(key)
                        Tirific_Template[fit_key] = f"{Tirific_Template[fit_key]} {' '.join([f'{x:{format}}' for x in fitting_settings[key][fit_key]])}"
set_fitting_parameters.__doc__ = '''
 NAME:
    set_fitting_parameters

 PURPOSE:
    Set the parameters that control the fitting in the Tirific template

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = The tirific template to be modified

 OPTIONAL INPUTS:
    debug = False

    parameters_to_adjust  = ['NO_ADJUSTMENT']
    a list of parameters that have to be set

    modifiers = ['EMPTY']
    modifies for startdelt,enddelts,mindelt

    stage = 'initial'
    fitting stage

    initial_estimates = ['EMPTY']
    initial_estimates for the various values,
    if not provided they are guessed from the Template and cube specifics

 OUTPUTS:
    Nothing is returned but the Tirific Template is updated

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_generic_fitting(Configuration, key , stage = 'initial', values = [60,5.], \
                        limits = [[0.,0.],[0.,0.],[0.,0.]],debug = False, slope = [0, 0], flat_slope = False , symmetric = False,\
                        upper_bracket = [10.,100.], lower_bracket=[0., 50.], fixed = True, moderate = 3, step_modifier = [1.,1.,1.],\
                        flat_inner = 3):
    if debug:
        print_log(f'''SET_GENERIC_FITTING: We are processing {key}.
''', Configuration['OUTPUTLOG'] ,debug = True)
    if isinstance(flat_inner,int):
        flat_inner = [flat_inner,flat_inner]
    NUR = Configuration['NO_RINGS']
    if all(x == 0. for x in np.array(slope,dtype=float)):
        slope = [NUR,NUR]
    if all(x == 0. for x in np.array(limits,dtype=float).reshape(6)):
        if debug:
            print_log(f'''SET_GENERIC_FITTING: Implementing limits
''', Configuration['OUTPUTLOG'])

        limits = [[set_limits(values[0]-values[1]*5.,*lower_bracket),\
                    set_limits(values[0]+values[1]*5.,*upper_bracket)] for x in limits]



    input= {}
    if debug:
            print_log(f'''SET_GENERIC_FITTING: flat is {fixed}
''', Configuration['OUTPUTLOG'])
    if (stage in ['after_os','final_os','after_cc','after_ec','parameterized']) or fixed:
        if debug:
            print_log(f'''SET_GENERIC_FITTING: Fitting all as 1.
''', Configuration['OUTPUTLOG'])
        input['VARY'] =  np.array([f"{key} 1:{NUR} {key}_2 1:{NUR}"],dtype=str)
        input['PARMAX'] = np.array([limits[0][1]],dtype=float)
        input['PARMIN'] = np.array([limits[0][0]],dtype=float)
        input['MODERATE'] = np.array([moderate],dtype=int) #How many steps from del start to del end
        input['DELSTART'] = np.array([values[1]*step_modifier[0]],dtype=float) # Starting step
        input['DELEND'] = np.array([0.1*values[1]*step_modifier[1]],dtype=float) #Ending step
        input['MINDELTA'] = np.array([0.05*values[1]*step_modifier[2]],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
    else:
        if not symmetric:
            if debug:
                print_log(f'''SET_GENERIC_FITTING: implementing a varying non-symmetric profile.
{'':8s} step_modifier = {step_modifier}
{'':8s} values = {values}
{'':8s} limits = {limits}
''', Configuration['OUTPUTLOG'])
            input['VARY'] = []
            end = []
            add = ''
            for i in[0,1]:
                if i == 1: add='_2'
                if flat_inner[i]+1 >= NUR:
                    input['VARY'].append(f"!{key}{add} {NUR}")
                    end.append(NUR-1)
                else:
                    input['VARY'].append(f"!{key}{add} {NUR}:{flat_inner[i]+1}")
                    end.append(flat_inner[i])

            input['VARY'].append(f"{key} 1:{end[0]} {key}_2 1:{end[1]}")
            input['VARY'] = np.array(input['VARY'],dtype=str)
            input['PARMAX'] = np.concatenate((np.array([limits[1][1]],dtype=float),\
                                              np.array([limits[2][1]],dtype=float),\
                                              np.array([limits[0][1]],dtype=float)))
            input['PARMIN'] = np.concatenate((np.array([limits[1][0]],dtype=float),\
                                              np.array([limits[2][0]],dtype=float),\
                                              np.array([limits[0][0]],dtype=float)))
            input['MODERATE'] =np.array([moderate,moderate,moderate],dtype=float) #How many steps from del start to del end
            input['DELSTART'] =np.array([2.,2.,0.5],dtype=float)*step_modifier[0]*values[1]# Starting step
            input['DELEND'] = np.array([0.1,0.1,0.05],dtype=float)*step_modifier[1]*values[1] #Ending step
            input['MINDELTA'] = np.array([0.1,0.1,0.075],dtype=float)*step_modifier[2]*values[1] #saturation criterum when /SIZE SIZE should be 10 troughout the code
        else:
            if debug:
                print_log(f'''SET_GENERIC_FITTING: implementing a varying symmetric profile
{'':8s} step_modifier = {step_modifier}
{'':8s} values = {values}
{'':8s} limits = {limits}
''', Configuration['OUTPUTLOG'])
            flat_inner = int(np.min(flat_inner))
            if flat_inner+1 >= NUR:
                input['VARY'] =  np.concatenate((np.array([f"!{key} {NUR} {key}_2 {NUR}"],dtype=str),\
                                                 np.array([f"{key} 1:{NUR-1} {key}_2 1:{NUR-1}"],dtype=str)))
            else:
                input['VARY'] =  np.concatenate((np.array([f"!{key} {NUR}:{flat_inner+1} {key}_2 {NUR}:{flat_inner+1}"],dtype=str),\
                                                 np.array([f"{key} 1:{flat_inner} {key}_2 1:{flat_inner}"],dtype=str)))
            input['PARMAX'] = np.concatenate((np.array([limits[1][1]],dtype=float),\
                                              np.array([limits[0][1]],dtype=float)))
            input['PARMIN'] = np.concatenate((np.array([limits[1][0]],dtype=float),\
                                              np.array([limits[0][0]],dtype=float)))
            input['MODERATE'] =np.array([moderate,moderate],dtype=float) #How many steps from del start to del end
            input['DELSTART'] =np.array([2.,0.5],dtype=float)*step_modifier[0]*values[1] # Starting step
            input['DELEND'] = np.array([0.1,0.05],dtype=float)*step_modifier[1]*values[1] #Ending step
            input['MINDELTA'] = np.array([0.1,0.075],dtype=float)*step_modifier[2]*values[1] #saturation criterum when /SIZE SIZE should be 10 troughout the code
        # then we need to set the warp slope

        forvarindex = ''
        if Configuration['NO_RINGS'] > 5:
            implement_keys = [key, f"{key}_2"]
            for i,cur_key in enumerate(implement_keys):
                if slope[i] < NUR:
                    if flat_slope:
                        if slope[i]+1 == NUR:
                            forvarindex = forvarindex+f"{cur_key} {NUR} "
                        else:
                            forvarindex = forvarindex+f"{cur_key} {NUR}:{slope[i]+1} "
                    else:
                        if slope[i]+1 >= NUR-1:
                            forvarindex = forvarindex+f"{cur_key} {NUR-1} "
                        else:
                            forvarindex = forvarindex+f"{cur_key} {NUR-1}:{slope[i]+1} "

        input['VARINDX'] = np.array([forvarindex],dtype=str)


    return input
set_generic_fitting.__doc__ =f'''
 NAME:
    set_generic_fitting

 PURPOSE:
    Generic routine for setting fitting parameters, SBR, VROT are separate

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    key = parameter to set the fitting for

 OPTIONAL INPUTS:
    debug = False

    stage = 'initial'
    STAGE OF THE FITTING WE ARE IN

    values = [60,5.]
    mean values of the parameter

    limits = [[0.,0.],[0.,0.],[0.,0.]]
    fitting limits

    slope = [0, 0]
    outer rings to be fitted as slope

    flat_slope = False
    if true slope should include the last ring

    symmetric = False
    Both sides are fitted as one

    upper_bracket = [10.,100.], lower_bracket=[0., 50.]
    limiting values

    fixed = True
    no radial variations

    moderate = 3
    moderate parameter

    step_modifier = [1.,1.,1.]
    array that modifies the fitting steps corresponding to [STARTDELTA, ENDDELTA, MINDELTA]
    for flat disks the standard values are [error,0.1*error,0.05 error]
    for varying disk [1.,0.1,1.] for the varying part and [error/2., 0.05*error, 0.1*error]  for the flat part

    flat_inner = 3
    amount of inner rings to fix

 OUTPUTS:
    A directory with the specified fitting values

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Function
def set_model_parameters(Configuration, Tirific_Template,Model_Values, stage = 'initial',debug = False):
    parameters_to_set = ['RADI','VROT_profile','Z0','SBR_profile','INCL','PA','XPOS','YPOS','VSYS','SDIS']


    check_parameters = []
    if 'VSYS' in Model_Values:
        vsys =Model_Values['VSYS'][0]/1000.
    else:
        vsys=100.
    scramble = np.zeros(len(parameters_to_set))
    # If the inclination is low we want to throw the inclination and VROT out of Vobs hence we scramble with random values
    #if stage == 'initialize_def_file':
    #        if 'INCL' in parameters_to_set:
    #            if Model_Values['INCL'][0] < 40.:
    #                scramble[parameters_to_set.index('INCL')] = np.random.uniform(-5,5,1)[0]
                    #if 'VROT_profile' in parameters_to_set:
                    #    Model_Values['VROT_profile'] = [set_limits(x+np.random.uniform(-2.*Configuration['CHANNEL_WIDTH'],2.*Configuration['CHANNEL_WIDTH'],1)[0], 0,600.) for x in Model_Values['VROT_profile'] ]
                    #    scramble[parameters_to_set.index('VROT_profile')] = np.mean([x*np.sin(np.radians(Model_Values['INCL'][0]))/np.sin(np.radians(Model_Values['INCL'][0]))])

    for key in parameters_to_set:
        if key in Model_Values:
            if key in ['VROT_profile','SBR_profile']:
                key_to_set = key.split('_')[0]
            else:
                key_to_set = key
            # if 2 long we have a value and error
            format = set_format(key_to_set)
            if len(Model_Values[key]) == 2:
                Tirific_Template[key_to_set]= f"{Model_Values[key][0]+scramble[parameters_to_set.index(key)]:{format}}"
            else:
                Tirific_Template[key_to_set]= f"{' '.join([f'{x+scramble[parameters_to_set.index(key)]:{format}}' for x in Model_Values[key][:int(Configuration['NO_RINGS'])]])}"
            if key_to_set != 'RADI':
                key_write = f"{key_to_set}_2"
                if f"{key}_2" in Model_Values:
                    key = f"{key}_2"
                if len(Model_Values[key]) == 2:
                    Tirific_Template[key_write]= f"{Model_Values[key][0]+scramble[parameters_to_set.index(key)]:{format}}"
                else:
                    Tirific_Template[key_write]= f"{' '.join([f'{x+scramble[parameters_to_set.index(key)]:{format}}' for x in Model_Values[key][:int(Configuration['NO_RINGS'])]])}"

            check_parameters.append(key)
        else:
            if key == 'RADI':
                rad = set_rings(Configuration,debug=debug)
                if len(Configuration['OLD_RINGS']) == 0:
                    size_in_beams = (rad[-1]-1./5.*Configuration['BEAM'][0])/Configuration['BEAM'][0]
                    Configuration['OLD_RINGS'].append(f"{size_in_beams:.1f}")
                Tirific_Template['RADI']= f"{' '.join([f'{x:.2f}' for x in rad])}"
                Tirific_Template['NUR']=str(len(rad))
                check_parameters.append('RADI')
            elif key == 'Z0':
                check_parameters.append('Z0')
                if Model_Values['INCL'][0] > 80:
                    Tirific_Template['Z0'] = f"{np.max([convertskyangle(Configuration,0.2,distance=Configuration['DISTANCE'],physical= True),Configuration['BEAM'][0]/4.]):.3f}"

                else:
                    Tirific_Template['Z0'] = f"{convertskyangle(Configuration,0.2,distance=Configuration['DISTANCE'],physical= True):.3f}"


                Tirific_Template['Z0_2'] = Tirific_Template['Z0']

            elif key == 'SDIS':
                check_parameters.append('SDIS')
                Tirific_Template['SDIS'] = '8.'
                Tirific_Template['SDIS_2'] = '8.'
            elif key == 'XPOS':
                if 'RA' in Model_Values:
                    if len(Model_Values['RA']) == 2:
                        Tirific_Template[key]= f"{Model_Values['RA'][0]:.8e}"
                        Tirific_Template[f"{key}_2"]= f"{Model_Values['RA'][0]:.8e}"
                    else:
                        Tirific_Template[key]= f"{' '.join([f'{x:.8e}' for x in Model_Values['RA'][:int(Configuration['NO_RINGS'])]])}"
                        Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:.8e}' for x in Model_Values['RA'][:int(Configuration['NO_RINGS'])]])}"
                    check_parameters.append('XPOS')
            elif key == 'YPOS':
                if 'DEC' in Model_Values:
                    if len(Model_Values['DEC']) == 2:
                        Tirific_Template[key]= f"{Model_Values['DEC'][0]:.8e}"
                        Tirific_Template[f"{key}_2"]= f"{Model_Values['DEC'][0]:.8e}"
                    else:
                        Tirific_Template[key]= f"{' '.join([f'{x:.8e}' for x in Model_Values['DEC'][:int(Configuration['NO_RINGS'])]])}"
                        Tirific_Template[f"{key}_2"]=f"{' '.join([f'{x:.8e}' for x in Model_Values['DEC'][:int(Configuration['NO_RINGS'])]])}"
                    check_parameters.append('YPOS')



    #if we are in the initial stage check that all parameters are set
    if stage == 'initial':
        for key in parameters_to_set:
            if not key in check_parameters:
                raise InitializeError(f"The parameter {key} is not set in the initialization")

    #make sure that the first value in VROT = 0
    vrot = Tirific_Template['VROT'].split()
    if float(vrot[0]) != 0.:
        Tirific_Template['VROT']=f" 0. {' '.join([f'{x}' for x in vrot[:int(Configuration['NO_RINGS']-1)]])}"
        Tirific_Template['VROT_2']=f" 0. {' '.join([f'{x}' for x in vrot[:int(Configuration['NO_RINGS']-1)]])}"
    no_declining_vrot(Configuration, Tirific_Template, debug = debug)
set_model_parameters.__doc__ =f'''
 NAME:
    set_model_parameters

 PURPOSE:
    Set the model values parameters in the tirific file that are singular values that apply to all of the fitting such as names and other issues

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template
    Model_Values = the values for setting

 OPTIONAL INPUTS:
    debug = False

    stage = 'initial'
    stage of the fitting process

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
#function to check that all parameters in template have the proper length.
def set_new_size(Configuration,Tirific_Template, Fits_Files, fit_type = 'Undefined', stage = 'initial',
                    current_run='Not Initialized', debug = False, Variables =
                    ['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                     'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2', 'AZ1P', 'AZ1W' ,'AZ1P_2','AZ1W_2', 'RADI']):
    if debug:
        print_log(f'''SET_NEW_SIZE: Starting the adaptation of the size in the template
''',Configuration['OUTPUTLOG'],debug =True)
    for key in Variables:
        if not key in Tirific_Template:
            Variables.remove(key)
    parameters = get_from_template(Configuration,Tirific_Template,Variables,debug = debug)
    # do a check on whether we need to modify all
    radii = set_rings(Configuration,debug=debug)
    if not Configuration['NEW_RING_SIZE']:
        print_log(f'''SET_NEW_SIZE: The rings size is stable and we are not interpolating
''',Configuration['OUTPUTLOG'])
        interpolate = False
    else:
        print_log(f'''SET_NEW_SIZE: The rings size is updated and we are interpolating the old values.
''',Configuration['OUTPUTLOG'])
        old_radii = parameters[Variables.index('RADI')]
        Configuration['NEW_RING_SIZE'] = False
        interpolate = True
    #if we add a ring we need to make sure that the radius gets extended
    if len(radii) > len(parameters[Variables.index('RADI')]):
        parameters[Variables.index('RADI')] = radii
        #and the last SBR values are far too high, set them to zero such that fix_sbr will correct them
        for i in [-3,-2,-1]:
            parameters[Variables.index('SBR')][i] = 0.
            parameters[Variables.index('SBR_2')][i] = 0.


    for i,key in enumerate(Variables):

        if debug:
            print_log(f'''SET_NEW_SIZE: We are processing {key}
{'':8s}SET_NEW_SIZE: We have a parameter of length {len(parameters[i])}.
{'':8s}SET_NEW_SIZE: Our current number of rings in the model is {Configuration['NO_RINGS']}.
{'':8s}SET_NEW_SIZE: {parameters[i]}
''',Configuration['OUTPUTLOG'])
        if key == 'RADI':
            Tirific_Template[key] = f" {' '.join([f'{x:.2f}' for x in radii])}"
        else:
            if interpolate:
                if debug:
                    print_log(f'''SET_NEW_SIZE: We are interpolating par = {parameters[i]} old radii={old_radii} new radii={radii}
    ''',Configuration['OUTPUTLOG'])
                if len(parameters[i]) > len(old_radii):
                    if debug:
                        print_log(f'''SET_NEW_SIZE: The parameters have more values than the radii. Cutting the end.
    ''',Configuration['OUTPUTLOG'])
                        parameters[i] = parameters[i][:len(old_radii)-1]
                elif len(parameters[i]) < len(old_radii):
                    if debug:
                        print_log(f'''SET_NEW_SIZE: The parameters have less values than the radii. Adding the last value until match.
    ''',Configuration['OUTPUTLOG'])
                        while len(parameters[i]) < len(old_radii):
                            parameters[i].append(parameters[i][-1])

                parameters[i] = list(np.interp(np.array(radii,dtype=float),np.array(old_radii,dtype=float),np.array(parameters[i],dtype=float)))

            format = set_format(key)

            if len(parameters[i]) > Configuration['NO_RINGS']-1:
                # if we are cutting a ring it is likely the outer ring have done weird stuff so we flatten the curve

                Tirific_Template[key] =f" {' '.join([f'{x:{format}}' for x in parameters[i][:int(Configuration['NO_RINGS'])]])}"
            elif len(parameters[i]) <= Configuration['NO_RINGS']-1:
                # We simply add the last value
                counter = len(parameters[i])
                while counter !=  Configuration['NO_RINGS']:
                    Tirific_Template[key] = Tirific_Template[key] + f" {parameters[i][-1]:{format}}"
                    counter += 1


        if debug:
            print_log(f'''SET_NEW_SIZE: We wrote the following line {Tirific_Template[key]}
''',Configuration['OUTPUTLOG'])

    #Replace the old ring numbers in VARY and VARINDX
    old_rings = calc_rings(Configuration,size_in_beams = int(round(float(Configuration['OLD_RINGS'][-2]))), ring_size  = 0.,debug=debug)
    current_rings = calc_rings(Configuration,debug=debug)
    #This could lead to replacing a value smaller than the other side
    Tirific_Template['VARY'] = Tirific_Template['VARY'].replace(f"{old_rings}",f"{current_rings}")
    Tirific_Template['VARINDX'] = Tirific_Template['VARINDX'].replace(f"{old_rings-1}",f"{current_rings-1}")
    Tirific_Template['NUR'] = f"{current_rings}"
    # if we cut we want to flatten things
    #if current_rings < old_rings:
    #    flatten_the_curve(Configuration,Tirific_Template,debug=debug)
    # Check whether the galaxy has become to small for variations
    if Configuration['SIZE_IN_BEAMS'] > Configuration['MINIMUM_WARP_SIZE']:
        Configuration['FIXED_PARAMETERS'][0] = Configuration['FIXED_PARAMETERS'][1]
    else:
        for parameter in ['INCL','Z0','SDIS','PA']:
            if parameter not in Configuration['FIXED_PARAMETERS'][0]:
                Configuration['FIXED_PARAMETERS'][0].append(parameter)
        flatten_the_curve(Configuration,Tirific_Template,debug=debug)
    if debug:
        print_log(f'''The following parameters will be fixed'
{'':8s} {Configuration['FIXED_PARAMETERS'][0]}
''',Configuration['OUTPUTLOG'])

    #update the limit_modifier
    Inclination = np.array([(float(x)+float(y))/2. for x,y in zip(Tirific_Template['INCL'].split(),Tirific_Template['INCL_2'].split())],dtype=float)
    set_limit_modifier(Configuration,Inclination,debug=debug)
    # Maybe increase the amount of loops in smaller galaxies
    set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,fit_type=fit_type, debug=debug,stage=stage)
    # if we change the radii we need to restart tirific
    finish_current_run(Configuration,current_run,debug=debug)
set_new_size.__doc__ =f'''
 NAME:
    set_new_size

 PURPOSE:
    Set the parameters in the template when the rings need to be adjusted.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting

    stage = 'initial'
    stage of the fitting

    current_run='Not Initialized'
    subprocess structure running tirific

    Variables =['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                     'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2', 'AZ1P', 'AZ1W' ,'AZ1P_2','AZ1W_2', 'RADI']
    Variables to be updated

 OUTPUTS:
    No output but the tirific process is stopped and the template updated

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_overall_parameters(Configuration, Fits_Files,Tirific_Template,stage = 'initial',fit_type='Undefined', flux = None,debug = False):

            if Configuration['OPTIMIZED']:
                Tirific_Template['INSET'] = f"{Fits_Files['OPTIMIZED_CUBE']}"
            else:
                Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
            #if stage in ['run_ec','initialize_ec']:
            #    if Configuration['NO_RINGS'] < 5:
            #        Tirific_Template['INIMODE'] = '1'
            #    elif Configuration['NO_RINGS'] < 15:
            #        Tirific_Template['INIMODE'] = '2'
            #    else:
            #        Tirific_Template['INIMODE'] = '3'
            #else:
            if Configuration['NO_RINGS'] < 20:
                Tirific_Template['INIMODE'] = '1'
            elif Configuration['NO_RINGS'] < 30:
                Tirific_Template['INIMODE'] = '2'
            else:
                Tirific_Template['INIMODE'] = '3'

            if Configuration['NO_RINGS'] > 3.:

                Tirific_Template['INTY'] = 0
                #Tirific_Template['INDINTY'] = 0
            else:
                #This should not be used with 8 <  rings.
                Tirific_Template['INTY'] = 1
                #preferably we'd use the akima spline but there is an issue with that where the final ring does not get modified
                #Tirific_Template['INDINTY'] = 2
            Tirific_Template['NUR'] = f"{Configuration['NO_RINGS']}"

            Tirific_Template['RESTARTNAME'] = f"{Configuration['LOG_DIRECTORY']}restart_{fit_type}.txt"
            #this could be fancier
            '''
            if Configuration['NO_RINGS'] < 3:
                Tirific_Template['NCORES'] = '2'
            elif Configuration['NO_RINGS'] < 6:
                Tirific_Template['NCORES'] = '3'
            elif Configuration['NO_RINGS'] < 12:
                Tirific_Template['NCORES'] = '4'
            else:
                Tirific_Template['NCORES'] = '6'
            '''

            Tirific_Template['NCORES'] = Configuration['NCPU']

            Tirific_Template['LOOPS'] = f"{Configuration['LOOPS']}"
            Tirific_Template['DISTANCE'] = f"{Configuration['DISTANCE']}"
            out_keys = ['LOGNAME','OUTSET','TIRDEF']
            out_extensions = ['log','fits','def']
            for i,key in enumerate(out_keys):
                Tirific_Template[key] = f"{fit_type}/{fit_type}.{out_extensions[i]}"
            #some things we only set if a header is provided

            Tirific_Template['BMAJ'] = f"{Configuration['BEAM'][0]:.2f}"
            Tirific_Template['BMIN'] = f"{Configuration['BEAM'][1]:.2f}"
            Tirific_Template['BPA'] = f"{Configuration['BEAM'][2]:.2f}"
            Tirific_Template['RMS'] = f"{Configuration['NOISE']:.2e}"

            if Configuration['HANNING_SMOOTHED']:
                instrumental_vres = (Configuration['CHANNEL_WIDTH']*2)/(2.*np.sqrt(2.*np.log(2)))
            else:
                instrumental_vres = (Configuration['CHANNEL_WIDTH']*1.2)/(2.*np.sqrt(2.*np.log(2)))
            Tirific_Template['CONDISP'] = f"{instrumental_vres:.2f}"
            if flux:
                Tirific_Template['CFLUX'] = f"{set_limits(flux/1.5e7,1e-6,1e-3):.2e}"
                Tirific_Template['CFLUX_2'] = f"{set_limits(flux/1.5e7,1e-6,1e-3):.2e}"
set_overall_parameters.__doc__ =f'''
 NAME:
    set_overall_parameters

 PURPOSE:
    Set the parameters in the tirific file that are singular values that apply to all of the fitting such as names and other issues

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    stage = 'initial'
    stage of fitting

    fit_type='Undefined_Stage'
    type of fitting

    flux = None,
    Total flux in the model

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 '''

def set_sbr_fitting(Configuration,systemic = 100., stage = 'no_stage',debug = False):
    if debug:
        print_log(f'''SET_SBR_FITTING: We are setting the SBR limits.
{'':8s} No_Rings = {Configuration['NO_RINGS']}
''',Configuration['OUTPUTLOG'],debug = True)
    sbr_input = {}
    inner_ring = 2
    if stage in ['initial','run_cc','initialize_ec','run_ec','initialize_os','run_os']:
        radii,sbr_ring_limits = sbr_limits(Configuration,systemic = systemic, debug = debug)
        if stage in ['run_ec','run_os']:
            sbr_ring_limits[-4:]=[x/5 for x in sbr_ring_limits]
        if debug:
            print_log(f'''SET_SBR_FITTING: Using these SBR limits.
{'':8s} limits = {sbr_ring_limits}
{'':8s} no of limits = {len(sbr_ring_limits)}
''',Configuration['OUTPUTLOG'])
        if stage in ['initial','run_cc']:
            max_size = 4
        elif stage in ['initialize_ec','run_ec','initialize_os','run_os']:
            max_size = 2.

        if Configuration['SIZE_IN_BEAMS'] < max_size:
            sbr_input['VARY'] =  np.array([f"SBR {x+1} SBR_2 {x+1}" for x in range(len(radii)-1,inner_ring-1,-1)],dtype=str)
            sbr_input['PARMAX'] = np.array([1 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float)
            #if stage in ['initial','run_cc']:
            #    sbr_input['PARMIN'] = np.array([sbr_ring_limits[x]/2. if x <= (3./4.)*len(radii) else 0 for x in range(len(radii)-1,inner_ring-1,-1)])
            #elif stage in ['initialize_ec','run_ec']:
            sbr_input['PARMIN'] = np.array([sbr_ring_limits[x]/1.5 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float)
            sbr_input['MODERATE'] = np.array([5 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float) #How many steps from del start to del end
            sbr_input['DELSTART'] = np.array([1e-4 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float) # Starting step
            sbr_input['DELEND'] = np.array([2.5e-6 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float) #Ending step
            sbr_input['MINDELTA'] = np.array([5e-6 for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
        else:
            sbr_input['VARY'] =  np.array([[f"SBR {x+1}",f"SBR_2 {x+1}"] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=str).reshape((len(radii)-inner_ring)*2)
            sbr_input['PARMAX'] = np.array([[1,1] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
        #if stage in ['initial','run_cc']:
            #    sbr_input['PARMIN'] = np.array([[sbr_ring_limits[x]/2.,sbr_ring_limits[x]/2.] if x <= (3./4.)*len(radii) else [0.,0.] for x in range(len(radii)-1,inner_ring-1,-1)]).reshape((len(radii)-inner_ring)*2)
            #elif stage in ['initialize_ec','run_ec']:
            sbr_input['PARMIN'] = np.array([[sbr_ring_limits[x]/4.,sbr_ring_limits[x]/4.] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            sbr_input['MODERATE'] = np.array([[5,5] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2) #How many steps from del start to del end
            sbr_input['DELSTART'] = np.array([[1e-4,1e-4] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2) # Starting step
            sbr_input['DELEND'] = np.array([[sbr_ring_limits[x]/20.,sbr_ring_limits[x]/20.] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            sbr_input['MINDELTA'] = np.array([[sbr_ring_limits[x]/20.,sbr_ring_limits[x]/20.] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)

        sbr_input['VARY'] = np.concatenate((sbr_input['VARY'],[f"SBR {' '.join([str(int(x)) for x in range(1,inner_ring+1)])} SBR_2 {' '.join([str(int(x)) for x in range(1,inner_ring+1)])}"]),axis=0)
        sbr_input['PARMAX'] = np.concatenate((sbr_input['PARMAX'],[2e-3]))
        sbr_input['PARMIN'] = np.concatenate((sbr_input['PARMIN'],[np.min(sbr_ring_limits)]))
        sbr_input['MODERATE'] = np.concatenate((sbr_input['MODERATE'],[5]))
        sbr_input['DELSTART'] = np.concatenate((sbr_input['DELSTART'],[1e-5]))
        sbr_input['DELEND'] = np.concatenate((sbr_input['DELEND'],[1e-6]))
        sbr_input['MINDELTA'] = np.concatenate((sbr_input['MINDELTA'],[2e-6]))
    elif stage in ['after_cc','after_ec','after_os']:
        #Used in Fit_Smoothed_Check
        sbr_input['VARY'] = [f"SBR 3:{Configuration['NO_RINGS']}, SBR_2 3:{Configuration['NO_RINGS']}"]
        sbr_input['PARMAX'] = np.concatenate(([2e-3],[2e-3]))
        sbr_input['PARMIN'] = np.concatenate(([0],[0]))
        sbr_input['MODERATE'] = np.concatenate(([5],[5]))
        sbr_input['DELSTART'] = np.concatenate(([1e-5],[1e-5]))
        sbr_input['DELEND'] = np.concatenate(([1e-6],[1e-6]))
        sbr_input['MINDELTA'] = np.concatenate(([2e-6],[2e-6]))

    return sbr_input
set_sbr_fitting.__doc__ =f'''
 NAME:
    set_sbr_fitting

 PURPOSE:
    set the sbr fitting parameters

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration


 OPTIONAL INPUTS:
    debug = False

    systemic = 100.
    systemic velocity of the source

    stage = 'no_stage'
    stage of the fitting

 OUTPUTS:
    a fitting dictionary

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_vrot_fitting(Configuration, stage = 'initial', rotation = [100,5.], debug = False):
    NUR = Configuration['NO_RINGS']
    vrot_input = {}
    vrot_limits = [set_limits(rotation[0]-rotation[1]-10,Configuration['CHANNEL_WIDTH'],360.), \
                   set_limits(rotation[0]+rotation[1]+10,80.,600.)]
    if debug:
        print_log(f'''SET_VROT_FITTING: We are setting the VROT limits.
{'':8s} No_Rings = {Configuration['NO_RINGS']}
{'':8s} Limits = {vrot_limits}
''',Configuration['OUTPUTLOG'],debug = True)
    if stage in ['final_os']:
        modifier= [Configuration['LIMIT_MODIFIER'][0]/2.,Configuration['LIMIT_MODIFIER'][0],Configuration['LIMIT_MODIFIER'][0]/4.]
    elif stage in ['after_os']:
        modifier= [Configuration['LIMIT_MODIFIER'][0]/3.,Configuration['LIMIT_MODIFIER'][0]/2.,Configuration['LIMIT_MODIFIER'][0]/5.]
    else:
        modifier= [Configuration['LIMIT_MODIFIER'][0]/2,Configuration['LIMIT_MODIFIER'][0],Configuration['LIMIT_MODIFIER'][0]/3.]



    if stage in ['after_cc', 'after_ec', 'final_os']:
        vrot_input['VARY'] =  np.array([f"VROT {NUR}:2 VROT_2 {NUR}:2"],dtype=str)
    else:
        vrot_input['VARY'] =  np.array([f"!VROT {NUR}:2 VROT_2 {NUR}:2"],dtype=str)
    vrot_input['PARMAX'] = np.array([vrot_limits[1]],dtype=float)
    vrot_input['PARMIN'] = np.array([vrot_limits[0]],dtype=float)
    vrot_input['MODERATE'] = np.array([5],dtype=float) #How many steps from del start to del end
    vrot_input['DELSTART'] = np.array([2.*Configuration['CHANNEL_WIDTH']*modifier[0]],dtype=float) # Starting step
    #These were lower in the original fat
    vrot_input['DELEND'] = np.array([0.1*Configuration['CHANNEL_WIDTH']*modifier[1]],dtype=float) #Ending step
    vrot_input['MINDELTA'] = np.array([0.1*Configuration['CHANNEL_WIDTH']*modifier[2]],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
    #if there is not values in the center we connect the inner ring to the next ring
    forvarindex = ''
    if Configuration['NO_RINGS'] > 5:
        if Configuration['EXCLUDE_CENTRAL'] or rotation[0] > 150.:
            forvarindex = 'VROT 2 VROT_2 2 '
        if Configuration['OUTER_SLOPE_START'] == NUR:
            if Configuration['NO_RINGS'] > 5:
                inner_slope =  int(round(set_limits(Configuration['RC_UNRELIABLE'],round(NUR/2.),NUR-1)))
            else:
                inner_slope = NUR
            if Configuration['NO_RINGS'] > 15 and inner_slope > int(Configuration['NO_RINGS']*4./5.) :
                inner_slope = int(Configuration['NO_RINGS']*4./5.)
            Configuration['OUTER_SLOPE_START'] = inner_slope
        else:
            inner_slope = Configuration['OUTER_SLOPE_START']


        if inner_slope >= NUR-1:
            if rotation[0] > 180.:
                forvarindex = forvarindex+f"VROT {NUR}:{NUR-1} VROT_2 {NUR}:{NUR-1} "
            else:
                forvarindex = forvarindex+f"VROT {NUR-1} VROT_2 {NUR-1} "
        else:
            if rotation[0] > 180.:
                forvarindex = forvarindex+f"VROT {NUR}:{inner_slope} VROT_2 {NUR}:{inner_slope} "
            else:
                forvarindex = forvarindex+f"VROT {NUR-1}:{inner_slope} VROT_2 {NUR-1}:{inner_slope} "
    vrot_input['VARINDX'] = np.array([forvarindex],dtype=str)


    return vrot_input
set_vrot_fitting.__doc__ =f'''
 NAME:
    set_vrot_fitting

 PURPOSE:
    Set the fitting parameters for VROT

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

    stage = 'initial'
    stage of the fitting

    rotation = [100,5.]
    estimate of mean rotation curve and error

 OUTPUTS:
    a fitting dictionary

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def smooth_profile(Configuration,Tirific_Template,key,min_error = 0.,debug=False ,profile_in = None, no_apply =False,fix_sbr_call = False):

    if key == 'SBR' and not fix_sbr_call:
        error_message = f'''SMOOTH_PROFILE: Do not use smooth_profile for the SBR, SBR is regularised in fix_sbr'''
        print_log(error_message,Configuration['OUTPUTLOG'],screen=True,debug = debug)
        raise FunctionCallError(error_message)

    if debug:
        print_log(f'''SMOOTH_PROFILE: Starting to smooth the {key} profile.
''',Configuration['OUTPUTLOG'],debug = True)

    if profile_in is None:
        profile = np.array(get_from_template(Configuration,Tirific_Template,[key,f"{key}_2"]),dtype = float)
    else:
        profile= copy.deepcopy(profile_in)

    original_profile = copy.deepcopy(profile)
    min_error = np.array(min_error,dtype=float)
    if min_error.size == 1:
        min_error = np.full(profile.size,min_error)
    if key in ['INCL','Z0', 'PA']:
        inner_fixed = Configuration['INNER_FIX']
    elif key in ['SDIS']:
        inner_fixed = [4,4]
    else:
        inner_fixed = [0,0]

    #he sbr profile is already fixed before geting to the smoothing
    if not fix_sbr_call:
        profile =fix_profile(Configuration, key, profile, Tirific_Template,inner_fix = inner_fixed,debug=debug)

    if key == 'VROT':
        #if profile[0,1] > profile[0,2] or np.mean(profile[1:3]) > 120.:
        if np.mean(profile[1:3]) > 120.:
            shortened =True
            profile = np.delete(profile, 0, axis = 1)
        else:
            shortened = False
    if debug:
        print_log(f'''SMOOTH_PROFILE: retrieved profile.
{'':8s}{profile}
''',Configuration['OUTPUTLOG'])
    # savgol filters do not work for small array
    for i in [0,1]:
        if len(profile[i]) < 8:
            #In this case we are using a cubic spline so no smoothing required
            pass
        elif len(profile[i]) < 15:
            profile[i] = savgol_filter(profile[i], 3, 1)
        elif len(profile[i]) < 20:
            profile[i] = savgol_filter(profile[i], 5, 2)
        elif len(profile[i]) < 25:
            profile[i] = savgol_filter(profile[i], 7, 3)
        else:
            profile[i] = savgol_filter(profile[i], 9, 4)
        if fix_sbr_call:
            below_zero = np.where(profile[i,:int(len(profile)/2.)] < 0.)[0]
            if below_zero.size > 0.:
                profile[i,below_zero] = min_error[i,below_zero]
            below_zero = np.where(profile[i,int(len(profile)/2.):] < 0.)[0]+int(len(profile)/2.)
            if below_zero.size > 0.:
                profile[i,below_zero] = profile[i,below_zero-1]/2.

    # Fix the settings
    format = set_format(key)
    if key == 'VROT':
        if shortened:
            tmp = [[],[]]
            tmp[0] = np.hstack([[0.], profile[0]])
            tmp[1] = np.hstack([[0.], profile[1]])
            profile =  np.array(tmp,dtype=float)
        else:
            profile[:,0] = 0.

    if not fix_sbr_call:
        profile =fix_profile(Configuration, key, profile, Tirific_Template,inner_fix=inner_fixed,debug=debug)


    if debug:
        print_log(f'''SMOOTH_PROFILE: profile after smoothing.
{'':8s}{profile}
''',Configuration['OUTPUTLOG'])
    if not no_apply:
        weights = get_ring_weights(Configuration,Tirific_Template,debug=debug)
        errors = get_error(Configuration,original_profile,profile,key,weights = weights, min_error=min_error,debug=debug)
        if key not in ['VROT']:
            # Check whether it should be flat
            profile,errors =modify_flat(Configuration,profile,original_profile,errors,key,inner_fix=inner_fixed,debug=debug)

        Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in profile[0,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:{format}}' for x in profile[1,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in errors[0,:int(Configuration['NO_RINGS'])]])}")
        Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in errors[1,:int(Configuration['NO_RINGS'])]])}")
        #if key in ['INCL'] and np.mean( profile[:,int(Configuration['NO_RINGS']/2.):int(Configuration['NO_RINGS'])]) < 40.:
        #    fix_vrot_for_incl_change(Configuration,Tirific_Template,original_profile,profile,debug=debug)

        if debug:
            print_log(f'''SMOOTH_PROFILE: This has gone to the template
{'':8s}{key} = {Tirific_Template[key]}
{'':8s}{key}_2 ={Tirific_Template[f"{key}_2"]}
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration['OUTPUTLOG'])

    return profile
smooth_profile.__doc__ =f'''
 NAME:
    smooth_profile

 PURPOSE:
    Read a profile from the tirific template and smooth it using a savgol smoothing.

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Template
    key = Tirific value to be smoothed

 OPTIONAL INPUTS:
    debug = False

    min_error = 0.
    The minimum eror that should be used

    profile_in = None
    if provided this will be smoothed instead of the key being extracted from the Template

    no_apply =False
    If true the smoothed profile will only be returned but not applied to the template

    fix_sbr_call = True
    If true the smooth_profile comes from the fix_sbr routine, this is the only routine allowed to smooth the SBR profile

 KEYWORD PARAMETERS:

 OUTPUTS:
    profile = the smoothed profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 EXAMPLE:
'''

def update_disk_angles(Configuration,Tirific_Template,debug = False):
    extension = ['','_2']
    for ext in extension:
        PA = np.array(get_from_template(Configuration,Tirific_Template,[f'PA{ext}'],debug=debug),dtype=float)[0]
        inc = np.array(get_from_template(Configuration,Tirific_Template,[f'INCL{ext}'],debug=debug),dtype=float)[0]
        if debug:
            print_log(f'''UPDATE_DISK_ANGLES: abtained  this from the template
{'':8s} inc{ext} = {inc}
{'':8s} PA{ext} = {PA}
''', Configuration['OUTPUTLOG'])
        angle_adjust=np.array(np.tan((PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180.)*180./np.pi,dtype = float)
        if ext == '_2':
            angle_adjust[:] +=180.
        if debug:
            print_log(f'''UPDATE_DISK_ANGLES: adusting AZ1P{ext} with these angles
{'':8s}{angle_adjust}
''', Configuration['OUTPUTLOG'])
        Tirific_Template.insert(f'AZ1W{ext}',f'AZ1P{ext}',f"{' '.join([f'{x:.2f}' for x in angle_adjust])}")
update_disk_angles.__doc__ =f'''
 NAME:
    update_disk_angles

 PURPOSE:
    Update the AZ1W and AZ1P parameters to match the warp

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Updated template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def write_new_to_template(Configuration, filename,Tirific_Template, Variables = ['VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2'], debug = False):
    with open(filename, 'r') as tmp:
        # Separate the keyword names
        for line in tmp.readlines():
            key = str(line.split('=')[0].strip().upper())
            if key in Variables:
                Tirific_Template[key] = str(line.split('=')[1].strip())
    # If we have written the INCL we need to update the limit modifier
    if 'INCL' in Variables or 'INCL_2' in Variables:
        Inclination = np.array([(float(x)+float(y))/2. for x,y in zip(Tirific_Template['INCL'].split(),Tirific_Template['INCL_2'].split())],dtype=float)
        set_limit_modifier(Configuration,Inclination,debug= debug)
write_new_to_template.__doc__ =f'''
 NAME:
    write_new_to_template

 PURPOSE:
    Write a def file into the template

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    filename = the name of the def file
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    Variables = ['VROT', 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2']
    the parameters to be updated from the file

 OUTPUTS:
    The template is updated

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
