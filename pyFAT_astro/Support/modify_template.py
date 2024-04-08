# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to Modify the Tirific_Template



import copy

from pyFAT_astro.Support.fat_errors import InitializeError,CfluxError,\
    FunctionCallError,BadConfigurationError,FittingError,FaintSourceError
from pyFAT_astro.Support.log_functions import print_log

import pyFAT_astro.Support.support_functions as sf

import numpy as np
import warnings
from math import isclose
from scipy.optimize import curve_fit, OptimizeWarning
from scipy.signal import savgol_filter
from scipy.interpolate import CubicSpline

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Ellipse

def apply_new_size(Configuration, new_size ):
    ''' Check whether we want to apply our new size '''
    # we want to have at least a quarter beam difference to apply the size.
    apply_size = [True,True]
   
    for i in [0,1]:
        if Configuration['FIX_SIZE'][i]:
            apply_size[i] = False
            print_log(f'''APPLY_NEW_SIZE: The size on side {i} is Fixed. Not changing
''', Configuration,case=['verbose']) 
        elif abs(Configuration['SIZE_IN_BEAMS'][i]-new_size[i]) < \
            Configuration['RING_SIZE']/2.:
            apply_size[i] = False
            print_log(f'''APPLY_NEW_SIZE: The size on side {i} is equal to that of the previous fit. Not changing
''', Configuration,case=['verbose'])    
        else: 
            for sizes in Configuration['OLD_SIZE'][:-1]:
                if abs(sizes[i] - new_size[i]) < Configuration['RING_SIZE']/2.:
                    print_log(f'''APPLY_NEW_SIZE: The side {i} has been fitted before, fixing it to the new size.
''', Configuration,case=['verbose'])
                    Configuration['FIX_SIZE'][i] = True        
    return apply_size
apply_new_size.__doc__ =f'''
 NAME:
    apply new size

 PURPOSE:
    determine whether the newly calculated size should be written to the configuration
    by comparing to previous sizes

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard configuration Dictionary
    new_size = the new sizes to check

 OPTIONAL INPUTS:

 OUTPUTS:
    boolean

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def arc_tan_sdis_function(axis,center,length,amplitude,mean):
    # to prevent instant turnover
    c = axis[-1]*0.1
    #and the turnover has to be beyon 20*of the inner part
    c2 = sf.set_limits(axis[-1]*0.2,axis[2],axis[int(len(axis)/1.5)])
    return -1*np.arctan((axis-(c2+abs(center)))/(c+abs(length)))*abs(amplitude)/np.pi + abs(mean)
arc_tan_sdis_function.__doc__ =f'''
 NAME:
    arc_tan_function

 PURPOSE:
    arc tangent function for fitting the SDIS this one can only be declining

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

def arc_tan_function(axis,center,length,amplitude,mean):
    # to prevent instant turnover
    c = axis[-1]*0.2
    #and the turnover has to be beyon 20*of the inner part
    c2 = sf.set_limits(axis[-1]*0.2,axis[2],axis[int(len(axis)/1.5)])
    return np.arctan((axis-(c2+abs(center)))/(c+abs(length)))*amplitude/np.pi + abs(mean)
arc_tan_function.__doc__ =f'''
 NAME:
    arc_tan_function

 PURPOSE:
    arc tangent function that can be declining or increain

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

def calc_new_size(Configuration,Tirific_Template,radii,sbr_in,sbr_ring_limits ):
    print_log(f'''CALC_NEW_SIZE: Starting to calculate the new size
''',Configuration,case=['debug_start'])

    sbr=copy.deepcopy(sbr_in)

    for i in [0,1]:
        if sbr[i,-1] > 5*  sbr[i,-2]:
              sbr[i,-1]= 0.8*sbr[i,-2]

    sm_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',
                            min_error= sbr_ring_limits,no_apply = True,
                            fix_sbr_call = True, profile_in = sbr )
    old_size = copy.deepcopy(Configuration['SIZE_IN_BEAMS'])
    print_log(f'''CALC_NEW_SIZE: The old size = {old_size}
''',Configuration,case=['debug_add'])
    new_size = copy.deepcopy(Configuration['SIZE_IN_BEAMS'])
    difference_with_limit = np.array(sbr-sbr_ring_limits,dtype=float)
    print_log(f'''CALC_NEW_SIZE: We find the following limiting SBR 
{'':8s} sbr = {sbr}
{'':8s} sbr_limits = {sbr_ring_limits}
{'':8s} difference_with_limit = {difference_with_limit}
''',Configuration,case=['debug_add'])
    smooth_diff = smooth_profile(Configuration,{'EMPTY':None}, 'SMOOTH_DIFFERENCE' ,\
        profile_in=difference_with_limit, min_error=[0.,0.],no_fix=True,no_apply=True)
    print_log(f'''CALC_NEW_SIZE: We find the following smoothed difference profile:
{'':8s} smooth_diff = {smooth_diff}
''',Configuration,case=['debug_add'])
    for i in [0,1]:
        #If all last three values are above the limit we need to extend
        if np.all(smooth_diff[i,-3:] > 0.):
            #We fit a Gaussian to the SBR and see where it drops below the last SBR Limit
            try:
                vals = sf.fit_gaussian(Configuration,radii,sm_sbr[i,:],errors=sbr_ring_limits[i,:])
                extend_radii = np.linspace(0.,Configuration['MAX_SIZE_IN_BEAMS']*Configuration['BEAM'][0],1000)
                Gauss_profile = sf.gaussian_function(extend_radii,*vals)
                Gauss_diff = Gauss_profile-sbr_ring_limits[i,-1]
                this_size=sf.set_limits(extend_radii[np.where(Gauss_diff < 0.)[0][0]]/Configuration['BEAM'][0]+0.5, Configuration['MIN_SIZE_IN_BEAMS'], Configuration['MAX_SIZE_IN_BEAMS'])
            except FittingError:
                print_log(f'''CALC_NEW_SIZE: We failed to fit a Gaussian
''',Configuration,case=['debug_add'])
                #If we cannot fit a gaussian we extend by a beam
                this_size=sf.set_limits(radii[-1]/Configuration['BEAM'][0]+1.0, Configuration['MIN_SIZE_IN_BEAMS'], Configuration['MAX_SIZE_IN_BEAMS'])
            except IndexError:
                print_log(f'''CALC_NEW_SIZE: We failed to find where the GAussian drop below 0.
''',Configuration,case=['debug_add'])
                #If we do not find any values in the gaussian below the limits we extend 2 beams
                this_size=sf.set_limits(radii[-1]/Configuration['BEAM'][0]+2.0, Configuration['MIN_SIZE_IN_BEAMS'], Configuration['MAX_SIZE_IN_BEAMS'])

            if this_size > new_size[i]+0.5:
                new_size[i]= copy.deepcopy(this_size)
        else:
            #we need to find where we cross the zero line
            #but have to make sure it is the first from the end
            loc = 0
            counter = len(smooth_diff[i,:])-1
            while loc == 0:
                if smooth_diff[i,counter-1] < 0. and counter > 1:
                    counter -= 1
                else:
                    if counter < 1:
                        loc = 1
                    else:
                        loc = counter
            x1 = float(radii[loc-1])
            x2 = float(radii[loc])
            y1 = float(smooth_diff[i,loc-1])
            y2 = float(smooth_diff[i,loc])
            this_size = sf.set_limits((x1+y1*(x2-x1)/(y1-y2))/Configuration['BEAM'][0],Configuration['MIN_SIZE_IN_BEAMS'], Configuration['MAX_SIZE_IN_BEAMS'])


            print_log(f'''CALC_NEW_SIZE: The profile for {i} drops below zero between {x1} and {x2}.
{'':8s} with y1 = {y1} and y2 = {y2}
{'':8s} the new size is {this_size}
''',Configuration,case= ['debug_add'])
            if this_size < new_size[i]-0.5:
                new_size[i]= copy.deepcopy(this_size)
        if Configuration['FIX_RING_SIZE']:
            real_size = new_size[i]-1./5.
            rings = int(real_size/Configuration['RING_SIZE'])
            new_size[i] = rings * Configuration['RING_SIZE']+1./5.
        if new_size[i] < Configuration['TOO_SMALL_GALAXY']:
             new_size[i] = Configuration['TOO_SMALL_GALAXY']

    print_log(f'''CALC_NEW_SIZE: We have the new_size {new_size} compared to what we had before {old_size}
''',Configuration,case = ['debug_add'])

    return new_size
calc_new_size.__doc__ =f'''
 NAME:
    calc_new_size

 PURPOSE:
    Based on the current SBR profiles and sbr_limits calc ulate the size of the galaxy

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    radii = The current radii of the model
    sbr = The current sbr of the model (Both sides)
    sbr_ring_limits = The current sbr limits

 OPTIONAL INPUTS:


 OUTPUTS:
    new_size
    a 2 element array indicating the new size for both sides of the model
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified
NOTE:
    The new size is neveer applied to Configuration
'''

def calculate_change_boundary(Configuration,multiple,theta_zero,phi_zero,
                              pa_boun = None, incl_boun = None ):
    change_boundary = []
    if pa_boun is None:
       pa_boun = [Configuration['PA_CURRENT_BOUNDARY'][1][0],
                 Configuration['PA_CURRENT_BOUNDARY'][1][1]] 
    if incl_boun is None:
       incl_boun = [Configuration['INCL_CURRENT_BOUNDARY'][1][0],
                 Configuration['INCL_CURRENT_BOUNDARY'][1][1]] 
    
    poss_theta_bound= []
    poss_phi_bound = []
    #we need theta and phi for all options of PA and INCL
    limit_combined = []
    for i in [0,1]:
        for j in [0,1]:
            limit_combined.append([pa_boun[i],incl_boun[j]])

    for lim in limit_combined:
        theta_tmp,phi_tmp,multiple_tmp = sf.calculate_am_vector(Configuration,[lim[0]],[lim[1]])
        # if the multiples do not correspond in the vector can be infinite
        if multiple_tmp[0] != multiple:
            poss_theta_bound.append(float('NaN'))
            poss_phi_bound.append(float('NaN'))
        else:
            poss_theta_bound.append(float(theta_tmp))
            poss_phi_bound.append(float(phi_tmp))
    minimum = -1*np.sqrt(np.min(np.array([x-theta_zero for x in poss_theta_bound],dtype=float))**2\
                +np.min(np.array([x-phi_zero for x in poss_phi_bound],dtype=float))**2)
    maximum = np.sqrt(np.max(np.array([x-theta_zero for x in poss_theta_bound],dtype=float))**2\
                +np.max(np.array([x-phi_zero for x in poss_phi_bound],dtype=float))**2)
    if np.isnan(minimum):
        minumum = 0.
    if np.isnan(maximum):
        maximum = 0.
    
    return [minimum,maximum]

calculate_change_boundary.__doc__ =f'''
 NAME:
    calculate_change_boundary

 PURPOSE:
    when regularising the PA and incl we do this through the change in the am vector
    but it should still maintain within the boundary limit this calculates the boundary
    limits for the change vector

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Configuration
    pa_incl =  array with PA, PA_2,INCL,INCL_2
    theta_zero = the centra theta value
    phi_zero = the central phi value


 OPTIONAL INPUTS:


 OUTPUTS:
    boundaries that the for the combined theta and phi change cannot cross

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def check_angles(Configuration,Tirific_Template):
    changed_angles=False
    incl_in = sf.load_tirific(Configuration,Tirific_Template,Variables= ['INCL','INCL_2'])
    pa_in = sf.load_tirific(Configuration,Tirific_Template,Variables= ['PA','PA_2'])
    incl = copy.deepcopy(incl_in)
    pa = copy.deepcopy(pa_in)
    rad =[float(x) for x in Tirific_Template['RADI'].split()]

    print_log(f'''CHECK_ANGLES: We start with
{'':8s} radius = {rad}
{'':8s} PA appr = {pa[0]}
{'':8s} PA rec = {pa[1]}
{'':8s} INCL appr = {incl[0]}
{'':8s} INCL rec = {incl[1]}
{'':8s} Changed_Angle = {changed_angles}
''', Configuration,case=['debug_start'])

    for side in [0,1]:
        incl_too_large = (np.array(incl[side],dtype=float) > 90.)
        if incl_too_large.any():

            print_log(f'''CHECK_ANGLES: Several INCL values were too large
''', Configuration,case = ['verbose'])
            incl[side] = [180.-x if y else x for x,y in zip(incl[side],incl_too_large)]
            #changed_angles = True

        pa_too_large = (np.array(pa[side],dtype=float) > 360.)
        if pa_too_large.any():
            if np.mean(np.array(pa[side],dtype=float)[~pa_too_large]) < 180.:
                pa[side] = [x-360. if y else x for x,y in zip(pa[side],pa_too_large)]
                #changed_angles = True
                print_log(f'''CHECK_ANGLES: Several PA values were too large
''', Configuration,case=['verbose'])

        pa_too_small = (np.array(pa[side],dtype=float) < 0.)
        if pa_too_small.any():
            if np.mean(np.array(pa[side],dtype=float)[~pa_too_small]) > 180.:
                pa[side] = [x+360. if y else x for x,y in zip(pa[side],pa_too_small)]
                #changed_angles = True

                print_log(f'''CHECK_ANGLES: Several PA values were too small
''', Configuration,case=['verbose'])
        pa_tmp,incl_tmp,changed_angles = sf.check_angular_momentum_vector(Configuration,\
                                            rad,pa[side],incl[side],modified= changed_angles,\
                                            side=side)
        #pa_tmp = sf.max_profile_change(Configuration,rad,pa[side],'PA',slope = Configuration['WARP_SLOPE'][side])
        pa[side] = pa_tmp
        #incl_tmp = sf.max_profile_change(Configuration,rad,incl[side],'INCL',slope = Configuration['WARP_SLOPE'][side])
        incl[side] = incl_tmp

        #exit()

    #Ensure that we have not made PA differences
    if pa[0][0] != pa[1][0]:
        if abs( pa[0][0]-pa[1][0]) > Configuration['MIN_ERROR']['PA'][0]:
            print_log(f'''CHECK_ANGLES: The central PA differs by too much:
{'':8s} differece = { abs( pa[0][0]-pa[1][0])}, max allowed = {Configuration['MIN_ERROR']['PA'][0]}
    ''', Configuration,case=['verbose'])
            # We need to fix them if we have.
            mean_pa = np.mean([pa[0][:Configuration['INNER_FIX'][0]],pa[1][:Configuration['INNER_FIX'][1]]])
            pa[0][:Configuration['INNER_FIX'][0]+1] = mean_pa
            pa[1][:Configuration['INNER_FIX'][1]+1] = mean_pa
                 
            changed_angles = True
        else:
            print_log(f'''CHECK_ANGLES: We correct some computational changes in the PA
    ''', Configuration,case=['verbose'])
        if (pa[0][0] > 360. and pa[1][0] < 360.) or \
            (pa[0][0] < 0. and pa[1][0] > 0.):
                pa[0][0:3] = pa[1][0:3]
        elif(pa[0][0] < 360. and pa[1][0] > 360.) or \
            (pa[0][0] > 0. and pa[1][0] < 0.):
                pa[1][0:3] = pa[0][0:3]


    Tirific_Template['INCL'] = f"{' '.join(f'{x:.2e}' for x in incl[0])}"
    Tirific_Template['INCL_2'] = f"{' '.join(f'{x:.2e}' for x in incl[1])}"
    Tirific_Template['PA'] = f"{' '.join(f'{x:.2e}' for x in pa[0])}"
    Tirific_Template['PA_2'] = f"{' '.join(f'{x:.2e}' for x in pa[1])}"
    print_log(f'''CHECK_ANGLES: We wrote to the template
{'':8s} PA appr = {Tirific_Template['PA']}
{'':8s} PA rec = {Tirific_Template['PA_2']}
{'':8s} INCL appr = {Tirific_Template['INCL']}
{'':8s} INCL rec = {Tirific_Template['INCL_2']}
''', Configuration,case=['verbose'])

    return changed_angles


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


 OUTPUTS:
    Tirific_Template is modified
    output is boolean indicating changes were made or not.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      load_tirific, join

 NOTE:
'''


def check_flat(Configuration,profile,error,key, last_reliable_ring = -1,\
                inner_fix = 4):
    last_reliable_ring +=2
    if last_reliable_ring == -1 or last_reliable_ring > len(profile)-1:
        last_reliable_ring = len(profile)

    print_log(f'''CHECK_FLAT: checking flatness
{'':8s}profile = {profile}
{'':8s}error = {error}
{'':8s}last_reliable_ring = {last_reliable_ring}
''',Configuration, case= ['debug_start'])

    inner = np.mean(profile[:inner_fix])
   
    if inner_fix+1 < last_reliable_ring:
        mean_error = np.mean(error[inner_fix+1:last_reliable_ring])
        outer = np.mean(profile[inner_fix+1:last_reliable_ring])
        outer_std = np.std(profile[inner_fix+1:last_reliable_ring])
    else:
        outer = inner
        outer_std = np.mean(error)
        mean_error = np.mean(error)

    if key in ['SDIS']:
        mean_error = mean_error*0.5
        outer_std = 2.*mean_error
    if abs(outer-inner) < mean_error or outer_std  < mean_error:
        print_log(f'''CHECK_FLAT: If  {abs(outer-inner)} less than {mean_error} or
{'':8s} the outer variation {outer_std} less than the median error {mean_error} we break and set flat.
''',Configuration, case=['debug_add'])
        return True
    if key in ['INCL','INCL_2']:
        if outer < 40. or inner < 40.:
            print_log(f'''CHECK_FLAT: the outer profile is low inclination {outer} hence we set this to flat.
''',Configuration, case=['debug_add'])
            return True

    for e,x,y in zip(error[1:last_reliable_ring],profile[1:last_reliable_ring],profile[0:last_reliable_ring]):
        print_log(f'''CHECK_FLAT: x = {x}, y = {y}, e = {e}
''',Configuration, case=['debug_add'])
        if not x-e/2. < y < x+e/2.:
            print_log(f'''CHECK_FLAT: This taco is bend
''',Configuration, case=['debug_add'])
            return False

    print_log(f'''CHECK_FLAT: All values were within the error of eachother
''',Configuration, case=['debug_add'])
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


    last_reliable_ring = -1
    the last ring in the profile that can be trusted

 OUTPUTS:
       True if no variation is found false if variation is found

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      zip(),np.std,np.mean,abs,print_log

 NOTE:
'''

def check_size(Configuration,Tirific_Template, fit_type = 'Undefined', \
        stage = 'initial', Fits_Files= 'No Files' ,current_run='Not Initialized',
        no_apply=False):
    print_log(f'''CHECK_SIZE: Starting check_size with the following parameters:
{'':8s}Rings = {Configuration['NO_RINGS']}, Size in Beams = {Configuration['SIZE_IN_BEAMS']}
{'':8s}The Rings isze we use = {Configuration['RING_SIZE']}
''',Configuration,case=['debug_start','verbose'])


    sbr_ring_limits = sf.sbr_limits(Configuration,Tirific_Template)
    #get the sbr profiles

    sbr = sf.load_tirific(Configuration,Tirific_Template,Variables=['SBR','SBR_2'],\
        array=True)
    radii = sf.load_tirific(Configuration,Tirific_Template,Variables=['RADI'],\
        array=True)
    if len(sbr[0]) != len(sbr_ring_limits):
        #Check what the appropriate size should be

        print_log(f'''CHECK_SIZE: Equalizing the sizes
''',Configuration,case= ['debug_add'])
        if len(sbr[0]) != Configuration['NO_RINGS']:

            print_log(f'''CHECK_SIZE: Interpolating SBR
''',Configuration,case= ['debug_add'])
            old_radii = np.array(sf.load_tirific(Configuration,Tirific_Template, ['RADI']),dtype = float)
            for i in [0,1]:
                sbr[i] = np.interp(np.array(radii,dtype=float),np.array(old_radii[0],dtype=float),np.array(sbr[i],dtype=float))
    print_log(f'''CHECK_SIZE: These are the ring SBRs and limits we will use:
{'':8s}SBR = {sbr}
{'':8s}limits = {sbr_ring_limits}
{'':8s}No. Rings  = {Configuration['NO_RINGS']}
''',Configuration,case= ['debug_add'])
    sbr_ring_limits = np.array([sbr_ring_limits,sbr_ring_limits],dtype=float)
    size_in_beams = calc_new_size(Configuration,Tirific_Template,radii,sbr,sbr_ring_limits)
    if not no_apply:
        apply_size = apply_new_size(Configuration,size_in_beams)
    else:
        apply_size = [False,False]
    print_log(f'''CHECK_SIZE: These have been fitted before {Configuration['OLD_SIZE']}
{'':8s} This is new_size {size_in_beams}, This is what currently is in {Configuration['SIZE_IN_BEAMS']}
''',Configuration,case= ['debug_add']) 
    trig = False
    tmp_old_size=list(copy.deepcopy(Configuration['SIZE_IN_BEAMS']))
    for i in [0,1]:
        if apply_size[i]:
            trig = True
            Configuration['SIZE_IN_BEAMS'][i] = copy.deepcopy(size_in_beams[i])         
    if trig:
        Configuration['OLD_SIZE'].append(tmp_old_size)
        ring_size, number_of_rings = sf.set_ring_size(Configuration)
        print_log(f'''CHECK_SIZE: Applied the size of {Configuration['SIZE_IN_BEAMS']}, ring size {ring_size} resulting in {number_of_rings} rings
''',Configuration,case=['verbose'])


    # Set the unreliable parameters
    #the lower the inclination the sooner the RC becomes unreliable
    limit_factor = sf.set_limits(2.5*np.mean(Configuration['LIMIT_MODIFIER']) ,2.,4.)
    print_log(f'''CHECK_SIZE: Using a limit factor to calculate where the RC is reliable of  {limit_factor}
''',Configuration,case= ['debug_add'])
    Configuration['RC_UNRELIABLE'] = get_number_of_rings(Configuration,sbr,limit_factor*sbr_ring_limits)-1
    if Configuration['RC_UNRELIABLE'] == Configuration['NO_RINGS']:
        Configuration['RC_UNRELIABLE'] -= 1
    

    for i in [0,1]:
        corr_val = np.where(sbr[i,2:] > sbr_ring_limits[i,2:]*3.)[0]+2
        if corr_val.size > 0:
            Configuration['LAST_RELIABLE_RINGS'][i] = corr_val[-1]+1
        else:
            Configuration['LAST_RELIABLE_RINGS'][i] = Configuration['NO_RINGS']
    print_log(f'''CHECK_SIZE: We set these as the last reliable rings {Configuration['LAST_RELIABLE_RINGS']}
{'':8s}The RC is deemed unreliable from ring {Configuration['RC_UNRELIABLE']} on.
''',Configuration,case= ['debug_add'])


    if any(apply_size):
        # Do not move this from here else other routines such as sbr_limits are messed up
        set_new_size(Configuration,Tirific_Template,fit_type= fit_type\
            ,current_run = current_run)
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,fit_type=fit_type)
        return False
    elif no_apply:
        return False
    else:
        return True
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
def cut_low_rings(Configuration, Tirific_Template,Fits_Files, fit_type = 'Unknown', current_run=None):
    print_log(f'''CUT_LOW_RINGS: Checking for too low rings in the smoothed profile
''',Configuration,case= ['debug_start'])
    sbr_ring_limits = sf.sbr_limits(Configuration,Tirific_Template)
    #get the sbr profiles
    sbr_profiles = sf.load_tirific(Configuration,Tirific_Template,Variables=['SBR','SBR_2'],\
        array=True)
    radii = sf.load_tirific(Configuration,Tirific_Template,Variables=['RADI'],\
        array=True)
    new_size=copy.deepcopy(Configuration['SIZE_IN_BEAMS'])
    for i in [0,1]:
        for x in range(len(sbr_profiles[i,:])-1,3,-1):
            if sbr_profiles[i,x] > sbr_ring_limits[x]:
                if radii[x]/Configuration['BEAM'][0] < Configuration['SIZE_IN_BEAMS'][i]:
                    new_size[i] = radii[x]/Configuration['BEAM'][0]
                break 
    print_log(f'''CUT_LOW_RINGS: This is or new model size {new_size}
{'':8s} This was before {Configuration['SIZE_IN_BEAMS']}
''',Configuration,case= ['debug_add'])        
    apply_size = apply_new_size(Configuration,new_size)
    if any(apply_size):
        for i in [0,1]:
            if apply_size[i]:
                Configuration['SIZE_IN_BEAMS'][i] = copy.deepcopy(new_size[i])  
               
        #Need to run this to make sure Configuration['NO_RINGS'] gets updated correctly
        ring_size, number_of_rings = sf.set_ring_size(Configuration)
        set_new_size(Configuration,Tirific_Template,fit_type= fit_type\
            ,current_run = current_run)
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,fit_type=fit_type)
    print_log(f'''CUT_LOW_RINGS: At the end we get
{'':8s} {Tirific_Template}
''',Configuration,case= ['debug_add']) 
    
    




cut_low_rings.__doc__  =f'''
 NAME:
    cut_low_rings

 PURPOSE:
    if after the first smooth the outer rings are below the sbr limits they are removed.
    but only if the smoothing is reran (No_rings < 5)

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:


 OUTPUTS:
    No output the modified sbrs are writte nto the TiriFiC Template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Very simple function only to be used after smoothing
'''

def fit_arc(Configuration,radii,sm_profile,error, function_to_fit,key ):
    c2 = sf.set_limits(radii[-1]*0.2,radii[2],radii[int(len(radii)/1.5)])
    est_center = radii[-1]/2.-c2
    est_length = radii[-1]*0.1
    if key in ['SDIS']:
        est_amp = abs(np.max(sm_profile)-np.min(sm_profile))
    else:
        est_amp = np.mean(sm_profile[3:])-np.mean(sm_profile[:4])
    est_mean = np.mean(sm_profile)

    if not error.any():
        error = np.full(len(sm_profile),1.)
        absolute_sigma = False
    else:
        absolute_sigma = True
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        succes = False
        maxfev= int(100*(len(radii)))

        while not succes:

            print_log(f'''FIT_ARC: Starting the curve fit with {maxfev}
''',Configuration,case= ['debug_start'])
            try:
                arc_par,arc_cov  =  curve_fit(function_to_fit, radii, sm_profile,p0=[est_center,est_length,est_amp,est_mean]\
                                            ,sigma=error,absolute_sigma=absolute_sigma,maxfev=maxfev)
                new_profile = function_to_fit(radii,*arc_par)
                new_profile[:3] = np.mean(new_profile[:3])
                succes = True
            except OptimizeWarning:
                maxfev =  2000*(len(radii))
            except RuntimeError as e:
                split_error = str(e)
                if 'Optimal parameters not found: Number of calls to function has reached maxfev' in \
                    split_error:
                    maxfev += 100*int(len(radii))
                    print_log(f'''FIT_ARC: We failed to find an optimal fit due to the maximum number of evaluations. increasing maxfev to {maxfev}
''',Configuration,case= ['debug_add'])
                else:
                    print_log(f'''FIT_ARC: Fit arc crashed with an unknow error:
{'':8s}{split_error}
''',Configuration)
                    raise RuntimeError(split_error)
            if maxfev >  1000*(len(radii)):
                print_log(f'''FIT_ARC: We failed to find an optimal fit to dispersion, returning the smoothed profile.
''',Configuration,case = ['verbose'])
                succes = True
                new_profile = copy.deepcopy(sm_profile)

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


 OUTPUTS:
    the fitted profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fit_polynomial(Configuration,radii,profile,error, key, \
                   Tirific_Template,inner_fix = 4,zero_point=None, \
                   boundary_limits = None, allowed_order= None,
                   return_order= False ):
    if boundary_limits is None:
        boundary_limits = [0,0.]
    if allowed_order is None:
        allowed_order=[None,None]

    print_log(f'''FIT_POLYNOMIAL: starting to fit the polynomial with the following input:
{'':8s} key = {key}
{'':8s} radii = {radii}
{'':8s} profile = {profile}
{'':8s} error = {error}
''',Configuration, case = ['debug_start'])
    only_inner = False
    if key in ['PA','INCL']:
        fixed = inner_fix
        only_inner = True
    elif key in ['Z0','THETA_FACTOR','PHI_FACTOR','CHANGE_ANGLE']:
        fixed = inner_fix
        only_inner = True
        error[:fixed] = error[:fixed]/10.
    elif key in ['VROT']:
        #fixed =len(radii)-Configuration['OUTER_SLOPE_START']
        fixed =sf.set_limits(len(radii)-np.min(Configuration['LAST_RELIABLE_RINGS']),1,len(radii))
        error[np.min(Configuration['LAST_RELIABLE_RINGS']):] = Configuration['CHANNEL_WIDTH']
        error[0] = Configuration['CHANNEL_WIDTH']
        #error[1] = error[1]*3.
    else:
        fixed = 1

    if len(radii) > 10.:
        start_order = int(len(radii)/5.)
    else:
        start_order = 0



    st_fit = int(0)
    if key in ['VROT']:
        bound_fit = int(1)
        if np.mean(profile[1:3]) > 120.:
            if zero_point == None:
                zero_point = np.mean(profile[1:3])
            st_fit = int(1)

        #This needs another -1 because the 0 and 1/5. ring are more or less 1 ring
        max_order = sf.set_limits(len(radii)-fixed-2,3,8)
        #The rotation curve varies a lot so the lower limit should be as high as possible
        #But at least 3 less than max order and maximally 4
        if len(radii)-fixed-2 <= 6:
            lower_limit=sf.set_limits(3,3,max_order-2)
        elif len(radii)-fixed-2 <= 10:
            lower_limit=sf.set_limits(4,3,max_order-2)
        else:
            lower_limit=sf.set_limits(5,3,max_order-2)
        start_order = sf.set_limits(start_order,lower_limit,max_order)

    else:
        bound_fit = int(0)
        if key in ['PA']:
            max_order = sf.set_limits(len(radii)-2,4,8)
            start_order = 3  
        elif key in ['INCL']:
            max_order = sf.set_limits(len(radii)-2,3,5)
            start_order = 1               
        elif key in ['Z0']:
            max_order = sf.set_limits(len(radii)-fixed,3,5)
        elif key in ['SBR']:
            max_order = sf.set_limits(len(radii)-2,4,8)
            start_order = 3
        else:
            max_order = sf.set_limits(len(radii)-1,3,7)

    if allowed_order[0]:
        if start_order < allowed_order[0]:
            start_order=  allowed_order[0]
    if start_order >= max_order:
        max_order = start_order+1
    if allowed_order[1]:
        if max_order > allowed_order[1]:
            max_order=  allowed_order[1]

    print_log(f'''FIT_POLYNOMIAL: For {key} we start at {start_order} because we have {len(radii)} rings of which {fixed} are fixed
{'':8s} this gives us a maximum order of {max_order}
''',Configuration,case = ['debug_add'])
    found = False
    count = 0
    failed=False
    while not found:
        reduced_chi = []
        order = range(start_order,max_order+1)

        print_log(f'''FIT_POLYNOMIAL: We will fit the following radii.
    {'':8s}{radii[st_fit:]}
    {'':8s} and the following profile:
    {'':8s}{profile[st_fit:]}
    {'':8s} weights = {1./error[st_fit:]}
    {'':8s} we are doing this for the {count} times.
    ''',Configuration,case = ['debug_add'])

        #make sure there are no 0. in the errors
        zero_locations = np.where(error[st_fit:] == 0.)[0]
        if zero_locations.size > 0.:
            error[zero_locations+st_fit] = 1./np.nanmax(1./error)

        for ord in order:
            fit_prof = np.poly1d(np.polyfit(radii[st_fit:],profile[st_fit:],ord,w=1./error[st_fit:]))

            if st_fit > 0.:
                fit_profile = np.concatenate(([zero_point],[e for e in fit_prof(radii[st_fit:])]))
            else:
                fit_profile = fit_prof(radii)
            #fit_profile = fit_prof(radii)

            if key not in ['SBR','THETA_FACTOR','PHI_FACTOR']:
                fit_profile = fix_profile(Configuration, key, fit_profile, \
                    Tirific_Template=Tirific_Template,inner_fix=inner_fix,\
                    singular = True,only_inner =only_inner)
            else:
                for i in range(len(fit_profile)-5,len(fit_profile)):
                    if fit_profile[i-1] < fit_profile[i]:
                        fit_profile[i]=fit_profile[i-1]*0.9
            
            red_chi = np.sum((profile[st_fit:]-fit_profile[st_fit:])**2/error[st_fit:]**2)/(len(radii[st_fit:])-ord)
            #We penailze profiles that go outside the boundaries
            print_log(f'''FIT_POLYNOMIAL: From the final profile we get chi**2 = {red_chi}
{'':8s} {key} = {', '.join([str(x) for x in fit_profile[st_fit:]])}
''',Configuration,case = ['debug_add'])
            if np.sum(np.absolute(np.array(boundary_limits,dtype=float))) != 0.:
                    with warnings.catch_warnings():
                        warnings.filterwarnings("error",category=RuntimeWarning)
                        
                        try:                       
                            diff = np.sum(np.array([abs(x-sf.set_limits(x,\
                                                    boundary_limits[0],\
                                                    boundary_limits[1]))/abs(x) \
                                            for x in  fit_profile[bound_fit:]],dtype = float))
                        except RuntimeWarning:
                            if key in ['THETA_FACTOR','PHI_FACTOR']:
                                diff = np.where(abs(fit_profile[bound_fit:]) > 1.)[0]
                                if len(diff) > 0:
                                    diff = 1.5
                                else:
                                    diff = 0.
                            elif key in ['CHANGE_ANGLE']:
                                diff=0.
                            else:    
                                print_log(f'''FIT_POLYNOMIAL: We are throwing an error as there should not be 0's in the diff,
{'':8s} boundary_limits =  {boundary_limits}
{'':8s} fit_profile  = {fit_profile[bound_fit:]}
{'':8s} For the parameter {key}
''',Configuration,case = ['debug_add'])
                                raise RuntimeError(f"The profile in fit_polynomial should not have zeros.") 

                       
                    if diff > 1. and not np.isnan(diff):
                        print_log(f'''FIT_POLYNOMIAL: because of crossing {boundary_limits}
{'':8s} we penalize the orginal chi**2 {red_chi} of order {ord} with {diff}
{'':8s} the fitted profile is {fit_profile}
''',Configuration,case = ['debug_add'])
                        red_chi = red_chi*(diff)
            if red_chi < 1.:
                red_chi = float('NaN')
            reduced_chi.append(red_chi)
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",message="All-NaN axis encountered"\
                            ,category=RuntimeWarning)
            if not np.isnan(np.nanmin(reduced_chi)):
                found = True
            else:
                error = error/2.
                count += 1
        if count > 100 or np.sum(error) == 0.:
            print_log(f'''FIT_POLYNOMIAL: we failed to fit a polynomial
''',Configuration,case = ['debug_add'])
            found = True
            failed =True     
        #if key in ['VROT'] and Configuration['NO_RINGS'] < 2.5*max_order:
        #    reduced_chi[-1] = reduced_chi[-1]*(ord/Configuration['NO_RINGS'])**2.5
    if not failed:
        print_log(f'''FIT_POLYNOMIAL: We have fitted these:
{'':8s} order = {[x for x in order]}
{'':8s} reduced chi = {reduced_chi}
''',Configuration,case = ['debug_add'])
        reduced_chi = np.array(reduced_chi,dtype = float)
        final_order = order[np.where(np.nanmin(reduced_chi ) == reduced_chi )[0][0]]
    else:
        print_log(f'''FIT_POLYNOMIAL: we failed to fit a polynomial
{'':8s} order = 0
''',Configuration,case = ['debug_add'])
        final_order=0.

    print_log(f'''FIT_POLYNOMIAL: We have regularised {key} with a polynomial of order {final_order}.
''',Configuration,case=['verbose'])
    fit_profile = np.poly1d(np.polyfit(radii[st_fit:],profile[st_fit:],final_order,w=1./error[st_fit:]))
    if st_fit > 0.:
        new_profile = np.concatenate(([zero_point],[e for e in fit_profile(radii[st_fit:])]))
    else:
        new_profile = fit_profile(radii)
    #if key in ['VROT'] and profile[1] < profile[2]:
    #    new_profile[1] = profile[1]
    if key not in ['SBR','THETA_FACTOR','PHI_FACTOR']:
        new_profile = fix_profile(Configuration, key, new_profile, \
            Tirific_Template=Tirific_Template,inner_fix=inner_fix,\
            singular = True,only_inner=only_inner)
    else:
        for i in range(len(fit_profile)-5,len(fit_profile)):
            if fit_profile[i-1] < fit_profile[i]:
                fit_profile[i]=fit_profile[i-1]*0.9
    if return_order:
        return new_profile,final_order
    else:
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
    error = errors on the profile
    key = parameter that is being regularised
    Tirific_Template = standard tirific template

 OPTIONAL INPUTS:
    zero_point = The inner point for the RC 

 OUTPUTS:
    polynomial fitted profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fix_outer_rotation(Configuration,profile):

    print_log(f'''FIX_OUTER_ROTATION: adjust last rings of VROT profile:
{'':8s}{profile}
{'':8s} from ring {Configuration['RC_UNRELIABLE']} we do not trust the rings.
''',Configuration,case = ['debug_start'])
    profile = np.array(profile,dtype=float)
    # if the outer parts are less then 5 channels there is something wrong And we just take a flat curve from max
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",message="Mean of empty slice."\
                            ,category=RuntimeWarning) 
        warnings.filterwarnings("ignore",message="invalid value encountered in scalar divide"\
                            ,category=RuntimeWarning)
        if np.mean(profile[Configuration['RC_UNRELIABLE']:]) < 5.*Configuration['CHANNEL_WIDTH']:
            Configuration['RC_UNRELIABLE'] = int(np.where(np.max(profile) == profile)[0][0])+1
            if Configuration['RC_UNRELIABLE'] < Configuration['NO_RINGS']-1:
                profile[Configuration['RC_UNRELIABLE']:] = profile[Configuration['RC_UNRELIABLE']-1]

                print_log(f'''FIX_OUTER_ROTATION: we adjusted the unreliable part.
    {'':8s}{profile}
    {'':8s} from ring {Configuration['RC_UNRELIABLE']} we do not trust the rings.
        ''',Configuration, case= ['debug_add'])

    #inner_slope = int(round(sf.set_limits(NUR*(4.-Configuration['LIMIT_MODIFIER'][0])/4.,round(NUR/2.),NUR-2)))
    if Configuration['RC_UNRELIABLE'] < Configuration['NO_RINGS']-1 and np.mean(profile[1:3]) > 180.:
        profile[Configuration['RC_UNRELIABLE']:] = profile[Configuration['RC_UNRELIABLE']-1]

    for i in range(int(Configuration['NO_RINGS']*3./4),Configuration['NO_RINGS']-1):
        if profile[i+1] > profile[i]*1.3:
            profile[i+1] = profile[i]*1.3

    print_log(f'''FIX_OUTER_ROTATION: this is corrected profile:
{profile}
''',Configuration, case= ['debug_add'])

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


 OUTPUTS:
    profile = the corrected profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Declining rotation curves are dealt with in no_declining_vrot
'''

def fix_profile(Configuration, key, profile, Tirific_Template = None,\
                 inner_fix = 4,singular = False,only_inner = False ):
  
    if isinstance(inner_fix,int):
        if singular:
            inner_fix = [inner_fix]
        else:
            inner_fix = [inner_fix,inner_fix]
       
    print_log(f'''FIX_PROFILE: Starting to fix {key} with the input values:
{'':8s}{profile}
''', Configuration, case= ['debug_add'])
    if key  == 'VROT' and Tirific_Template is None:
        print_log(f'''FIX_PROFILE: To fix VROT  we require a Tirific Template 
''', Configuration,case = ['main','screen'])
        raise FunctionCallError("FIX_PROFILE:  To fix VROT  we require a Tirific Template.")

    if key == 'SBR':
        print_log(f'''FIX_PROFILE: To fix sbr profiles use FIX_SBR.
''', Configuration,case = ['main','screen'])
        raise FunctionCallError("FIX_PROFILE: To fix sbr profiles use check SBR.")
    if singular:
        indexes = [0]
        profile = np.array([profile,profile])
        if key not in ['THETA_FACTOR','PHI_FACTOR','CHANGE_ANGLE']:
            inner_mean = np.nanmean([profile[0,:inner_fix[0]]])
        else:
            inner_mean = 0.
    else:
        indexes = [0,1]
        if np.sum(inner_fix) != 0. and key not in ['THETA_FACTOR','PHI_FACTOR','CHANGE_ANGLE'] :
            inner_mean = np.nanmean(np.concatenate((profile[0,:inner_fix[0]],profile[1,:inner_fix[1]])))
        else:
            inner_mean= 0.

    profile = np.array(profile,dtype=float)
    if key in ['VROT']:
        indexes = [0]
        inner_mean = 0.
        profile[0,0] = 0.
    else:
        for i in indexes:
            profile[i,:inner_fix[i]] = inner_mean
        print_log(f'''FIX_PROFILE: the  {inner_fix} inner rings are fixed for the profile:
{'':8s} profile = {profile[i,:]}
''', Configuration,case=['debug_add'])
    print_log(f'''FIX_PROFILE: the  inner mean is {inner_mean}.
''', Configuration,case=['debug_add'])

    for i in indexes:
        print_log(f'''FIX_PROFILE: From ring {Configuration['LAST_RELIABLE_RINGS'][i]} on we do not trust these rings.
''', Configuration,case=['debug_add'])
        if Configuration['LAST_RELIABLE_RINGS'][i] < len(profile[i,:]) and not only_inner:
            #profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:] = profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:]*0.25+profile[i,Configuration['LAST_RELIABLE_RINGS'][i]-1]*0.75
            profile[i,Configuration['LAST_RELIABLE_RINGS'][i]:] = profile[i,Configuration['LAST_RELIABLE_RINGS'][i]-1]
        print_log(f'''FIX_PROFILE:After fixing the last reliable rings
{'':8s} profile = {profile[i,:]}
''', Configuration,case=['debug_add'])
        if key == 'VROT':
            profile[i] =fix_outer_rotation(Configuration,profile[i])
        elif inner_fix[i] != 0:
            xrange = sf.set_limits((int(round(len(profile[0])-5.)/4.)),1,4)
        # need to make sure this connects smoothly
            for x in range(0,xrange):
                if inner_fix[i]+x < len(profile[i,:]):
                    profile[i,inner_fix[i]+x] = 1/(x+4./xrange)*inner_mean+ (1-1/(x+4./xrange))*profile[i,inner_fix[i]+x]
            #if key in ['PA','INCL']:
            #    profile[i,:] = sf.max_profile_change(Configuration,rad,profile[i,:],key)
        print_log(f'''FIX_PROFILE:After smoothed transition
{'':8s} profile = {profile[i,:]}
''', Configuration,case=['debug_add'])
        #profile[:,:Configuration['INNER_FIX']] = np.nanmean(profile[:,:Configuration['INNER_FIX']])
        if key in ['SDIS']:
            inner_max = np.nanmax(profile[i,:int(len(profile[i,:])/2.)])
            mean = np.nanmean(profile[i,:])
            if inner_mean < mean or inner_mean < np.nanmean(profile[i,-3:]):
                ind = np.where(profile[i,:] == inner_max)[0]
                if ind.size > 1:
                    ind = int(ind[0])
                else:
                    ind = int(ind)
                profile[i,:ind] = inner_max
            profile[i] =np.hstack([[profile[i,0]],[y if y <= x else x*0.95 for x,y in zip(profile[i,:],profile[i,1:])]])
    if key == 'VROT':
        tmp  = no_declining_vrot(Configuration,Tirific_Template,profile = profile[0])
        profile[0] = tmp
        profile[1] = tmp
    if singular:
        profile = profile[0]


    print_log(f'''FIX_PROFILE: The final profile for {key} is:
{'':8s}{profile}
''', Configuration,case=['debug_add'])
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

def fix_sbr(Configuration,Tirific_Template, smooth = False, initial = False ):
    print_log(f'''FIX_SBR: Starting a SBR check.
''',Configuration,case=['debug_start'])

    # get the cutoff limits
    cutoff_limits = sf.sbr_limits(Configuration,Tirific_Template)
    cutoff_limits = np.array([cutoff_limits,cutoff_limits],dtype=float)
    print_log(f'''FIX_SBR: Using these cutoff limits.
{'':8s}cutoff_limits = {cutoff_limits}
''',Configuration,case=['debug_add'])
    # Then get the profile from the template
    radii = sf.load_tirific(Configuration,Tirific_Template,Variables=['RADI'],array=True)
    sbr = sf.load_tirific(Configuration,Tirific_Template,Variables=['SBR','SBR_2'], array=True)
    original_sbr=copy.deepcopy(sbr)
    print_log(f'''FIX_SBR: Before modify.
{'':8s}sbr from template = {sbr}
''',Configuration,case=['debug_add'])
    # First make a correction on the inner 2 values
    sbr = inner_sbr_fix(Configuration,sbr,cutoff_limits)
    # Let's use a smoothed profile for the fittings
    sm_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',
                            min_error= cutoff_limits,no_apply = True,
                            fix_sbr_call = True,profile_in = sbr )

    print_log(f'''FIX_SBR: Before modify.
{'':8s}sbr  = {sbr}
{'':8s}sm_sbr  = {sm_sbr}
''',Configuration,case=['debug_add'])

    # We interpolate negative values as well as values below the limits inner part with a cubi
    errors = get_error(Configuration,sbr,sm_sbr,'SBR',min_error=cutoff_limits)
    print_log(f'''FIX_SBR: retrieved errors.
{'':8s}errors  = {errors}
''',Configuration,case=['debug_add'])
    no_zero_sbr = copy.deepcopy(sm_sbr)
    zeros = np.where(no_zero_sbr < 1e-9)
    no_zero_sbr[zeros] = 0.1*cutoff_limits[zeros]
    error_weights = cutoff_limits/no_zero_sbr*1./np.nanmin(cutoff_limits[:,2:]/no_zero_sbr[:,2:])
    error_weights[:,0] = 3.
    for i in [0,1]:
        error_weights[i] = [x if x <100. else 100. for x in error_weights[i]]
    errors =errors*error_weights
    print_log(f'''FIX_SBR: weighed errors errors.
{'':8s}sm_sbr = {sm_sbr}
{'':8s}cutoff_limits  = {cutoff_limits}
{'':8s}errors  = {errors}
{'':8s}weights  = {error_weights}
''',Configuration,case=['debug_add'])


    store_gaussian = []
    for i in [0,1]:
        corr_val = np.where(sbr[i,2:] > cutoff_limits[i,2:])[0]+2
        # If we have enough safe values in the profile we attempt to fit it
        if corr_val.size > 3.:
            fit_sbr = sm_sbr[i,corr_val]
            print_log(f'''FIX_SBR: The values used for fitting are {fit_sbr}.
''',Configuration,case=['debug_add'])
            try:
                if 'SBR' in Configuration['FIXED_PARAMETERS'][0] \
                    or initial:
                    vals = sf.fit_gaussian(Configuration,radii[corr_val],\
                                fit_sbr,errors=errors[i,corr_val])
                    gaussian = sf.gaussian_function(radii,*vals)
                else:
                    gaussian = fit_polynomial(Configuration,radii,sbr[i,:],\
                        errors[i,:],'SBR', Tirific_Template)
                # if the peak of this gaussian is in the inner two points replace it with the smoothed profile
                if np.any(np.where(np.max(gaussian) == gaussian)[0] < 2):
                    print_log(f'''FIX_SBR: We are trying to replace the inner gaussian.
{'':8s} gaussian = {gaussian[[0,1]]}
{'':8s} sm_sbr = {sm_sbr[i,[0,1]]}
''',Configuration,case=['debug_add'])
                    gaussian[[0,1]] = sm_sbr[i,[0,1]]
            except FittingError:
                # If we fail we try a CubicSpline interpolation
                try:
                    tmp = CubicSpline(radii[corr_val],fit_sbr,extrapolate = True,bc_type ='natural')
                    gaussian = tmp(radii)
                except RuntimeError:
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
        sbr[np.where(sbr<cutoff_limits)] = 2.*cutoff_limits[np.where(sbr<cutoff_limits)]
        #if our warp_slope is very small we want to increase the SBR significantly
        for i in [0,1]:
            if Configuration['WARP_SLOPE'][i] <  0.75*(Configuration['NO_RINGS']) and Configuration['WARP_SLOPE'][i] != None:
                sbr[i,Configuration['WARP_SLOPE'][i]:] = 0.5*sbr[i,Configuration['WARP_SLOPE'][i]-1]+2.5*sbr[i,Configuration['WARP_SLOPE'][i]:]
  

      

    # and ensure that both sides the inner two rings are the same
    sbr[:,[0,1]] = np.mean(sbr[:,[0,1]])
    print_log(f'''FIX_SBR: after ensuring both ring are the same
{'':8s}sbr = {sbr[i, :]}
''',Configuration,case=['debug_add'])
    #And that they are less then 3/4 of ring 3 unless ring 3 is blanked
    ring_three = np.min(sbr[[0,1],2])
    if ring_three < 1e-15:
        ring_three = sbr[0,1]*4

    if sbr[0,1] > 3/4.* ring_three:
        sbr[:,[0,1]] = 3/4.* ring_three

    print_log(f'''FIX_SBR: after ring three
{'':8s}sbr = {sbr[i, :]}
''',Configuration,case=['debug_add'])
    if smooth and len(sbr[0]) > 6:
        for j in [0,1]:
            fit_profile = sbr[j]
            for i in range(len(fit_profile)-5,len(fit_profile)):
                if fit_profile[i-1] < fit_profile[i]:
                    fit_profile[i]=fit_profile[i-1]*0.9
            sbr[j] = fit_profile

    # And finally we need to make sure all values are positive and reasonable
    # And that values that are not fitted are set to 0.   
    not_too_faint = False             
    for i in [0,1]:
        #In setting rings to fit the Tirific rings start at one but python is 0 based
        last_ring_to_fit = np.where(radii <= \
            Configuration['SIZE_IN_BEAMS'][i]*Configuration['BEAM'][0])\
                [0][-1]   #+1            
        for n in range(len(sbr[i,:])):
            if n > last_ring_to_fit: 
                sbr[i,n] = 0. 
            elif sbr[i,n] < cutoff_limits[i,n]:
                sbr[i,n] = cutoff_limits[i,n]*1.5
            else:
                not_too_faint = True
        if sbr[i,last_ring_to_fit-1] < sbr[i,last_ring_to_fit]:
                last = sbr[i,last_ring_to_fit]
                second = sbr[i,last_ring_to_fit-1]
                sbr[i,last_ring_to_fit-1] = last
                sbr[i,last_ring_to_fit] = second         
    if smooth and not not_too_faint:
        print_log(f'''FIX_SBR: All SBR rings are below the cutoff limits
{'':8s}original sbr = {original_sbr}
{'':8s}cutoff_limits  = {cutoff_limits}
''',Configuration,case=['debug_add'])
        raise FaintSourceError(f"After correcting the SBRs all values we're below the cutoff limit, your source is too faint.")
           


    print_log(f'''FIX_SBR: After modify.
{'':8s}{sbr}
''',Configuration,case=['debug_add'])
    Tirific_Template['SBR'] = f"{' '.join([f'{x:.2e}' for x in sbr[0]])}"
    Tirific_Template['SBR_2'] = f"{' '.join([f'{x:.2e}' for x in sbr[1]])}"
    print_log(f'''FIX_SBR: We checked the surface brightness profiles.
''',Configuration,case=['verbose'])
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

def fix_vrot_for_incl_change(Configuration,Tirific_Template, \
                            incl_original,incl_modified ):
    vrot = sf.load_tirific(Configuration,Tirific_Template,\
            Variables=["VROT","VROT_2"],array=True)
    new_vrot = copy.deepcopy(vrot)
    for i in [0,1]:
        vobs = np.array([x*np.sin(np.radians(y)) for x,y in zip(vrot[i],incl_original[i])],dtype=float)
        new_vrot[i] = np.array([x/np.sin(np.radians(y)) for x,y in zip(vobs,incl_modified[i])],dtype=float)
    vrot = np.array([np.mean([x,y]) for x,y in zip(new_vrot[0],new_vrot[1])],dtype=float)
    format = sf.set_format("VROT")
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


 OUTPUTS:
    The template is corrected

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def flatten_the_curve(Configuration,Tirific_Template):
    to_flatten = ['INCL','Z0','PA','SDIS']
    for key in to_flatten:
        profile = sf.load_tirific(Configuration,Tirific_Template,Variables=\
            [key,f'{key}_2'],array=True)
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

    fit_type = 'Undefined'

 OUTPUTS:
    the flattened profiles for PA, INLC, SDIS and Z0 are written to the template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_error(Configuration,profile,sm_profile,key,min_error = 0.,\
    singular = False,weights= None,apply_max_error = False ):
    if weights is None:
        weights= [1.]
    print_log(f'''GET_ERROR: starting;
''',Configuration, case=['debug_start'])
    if singular:
        if len(weights) == 1:
            weights = np.full((len(profile)),weights[0])
        profile = np.array([profile],dtype=float)
        sm_profile = np.array([sm_profile],dtype=float)
        error =[[]]
        sides = [0]
     
    else:
        error = [[],[]] 
        profile = np.array(profile,dtype=float)
        if len(weights) == 1:
            weights = np.full((len(profile[0]),len(profile[1])),1.)
        sides =[0,1]
    min_error = get_min_error(Configuration,min_error_in=min_error,profile=profile,parameter=key)
    if singular:
        if min_error.shape[0] != 1:
            min_error = np.array([min_error],dtype=float)
    print_log(f'''GET_ERROR: Using these values;
{'':8s}original profile = {profile}
{'':8s}new profile = {sm_profile}
{'':8s}weights = {weights}
{'':8s}singular = {singular}
{'':8s}min_error = {min_error}
{'':8s}max_error = {apply_max_error}
''',Configuration, case=['debug_add'])
    
    for i in sides:
        error[i] = abs(profile[i]-sm_profile[i])/2.
        error[i]= error[i]/weights[i]
       
        error[i] = [np.max([y,x]) for x,y in zip(error[i],min_error[i])]
        
        if apply_max_error:
            if len(Configuration['MAX_ERROR'][key]) == len(error[i]):
                error[i] = [np.nanmin([x,y]) for x,y in zip(error[i],Configuration['MAX_ERROR'][key])]
            else:
                error[i] = [np.nanmin([x,float(Configuration['MAX_ERROR'][key][0])]) for x in error[i]]
        if key in ['PA','INCL','Z0','ARBITRARY']:
            error[i][:Configuration['INNER_FIX'][i]+1] = [np.min(min_error) for x in error[i][:Configuration['INNER_FIX'][i]+1]]
    if singular:
        error = np.array(error[0],dtype=float)
    else:
        error = np.array(error,dtype=float)
    
    low = np.where(error < min_error)
    if np.sum(low) > 0.:
        error[low] =min_error[low]

  
    print_log(f'''GET_ERROR: error =
{'':8s}{error}
''',Configuration,case=['debug_add'])
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
def get_min_error(Configuration,profile = None,min_error_in = None, parameter = None):
    if profile is None:
        raise FunctionCallError(f'You have to provide a profile to match')
    if min_error_in is None:
        if parameter is None:
            raise FunctionCallError(f'You have to provide a parameter (parameter=) or input min_error (min_error_in=)')
        else:
            if parameter in Configuration['MIN_ERR']:
                min_error = np.full(profile.shape,Configuration['MIN_ERR'][parameter][0])
            else:
                min_error = 0.
    else:
        min_error = copy.deepcopy(min_error_in)
    
    #If they have the same shape we stop here
   
    
    # else we first check it is equivalant
    if sf.isiterable(min_error):
        min_error= np.array(min_error,dtype=float)
        if min_error.shape == profile.shape:
            return min_error
        print_log(f'''GET_MIN_ERROR:  These are the input shapes:
{'':8s} min_errror =  {min_error.shape} 
{'':8s} profile =  {profile.shape}
''',Configuration,case=['debug_add'])
        try:
            if min_error.shape[0] == profile.shape[1]:
                min_error = np.array([min_error,min_error],dtype=float)
            elif min_error.shape[0] == 1:
                min_error = np.full(profile.shape,min_error[0])
            else:
                tmp_err1,profile1 = sf.make_equal_length(min_error[0],profile[0]) 
                tmp_err2,profile2 = sf.make_equal_length(min_error[1],profile[1])
                min_error = np.array([tmp_err1,tmp_err2],dtype=float)
        
        except IndexError:
            print_log(f'''GET_MIN_ERROR: Assuming singular
''',Configuration,case=['debug_add'])
            if min_error.shape[0] == 1:
                min_error = np.full(profile.shape,min_error[0])
            else:
                tmp_err1,profile1 = sf.make_equal_length(min_error[0],profile[0]) 
                min_error = np.array([tmp_err1],dtype=float)
    else:
        print_log(f'''GET_MIN_ERROR:  min_error is not iterable
''',Configuration,case=['debug_add'])
        min_error = np.full(profile.shape,min_error)
    if min_error.shape != profile.shape:
        print_log(f'''GET_MIN_ERROR:  These are the input shapes:
{'':8s} min_errror =  {min_error.shape} 
{'':8s} profile =  {profile.shape}
''',Configuration,case=['debug_add'])
        raise FunctionCallError('GET_MIN_ERROR: the min error and profile shape are not matching')
    return min_error          

def get_number_of_rings(Configuration,sbr,sbr_ring_limits ):
    '''Determine whether the amount of rings is good for the limits or not'''
    new_rings = Configuration['NO_RINGS']
    difference_with_limit = np.array(sbr-sbr_ring_limits,dtype=float)
    if np.all(difference_with_limit[:,-1] < 0.):
        print_log(f'''GET_NUMBER_OF_RINGS: both last rings are below the limit
''',Configuration,case=['debug_start'])
        for i in range(len(difference_with_limit[0,:])-1,int(new_rings/2.),-1):
            print_log(f'''GET_NUMBER_OF_RINGS: Checking ring {i}
''',Configuration,case=['debug_add'])
            if np.all(difference_with_limit[:,i] < 0.):
                #check that 1 any of the lesser rings are bright enough
                if np.any(sbr[:,i-1] > 1.5 *sbr_ring_limits[:,i-1]):
                    new_rings = i+1
                    print_log(f'''GET_NUMBER_OF_RINGS: we find that the previous rings are bright enough, rings = {new_rings}
''',Configuration,case=['debug_add'])
                    break
                else:
                    print_log(f'''GET_NUMBER_OF_RINGS: the previous rings are not bright enough so we reduce 1, old_rings = {new_rings}, new_rings = {i}
''',Configuration,case=['debug_add'])
                    new_rings = i
            else:
                #if not both values are below than this is the extend we want
                new_rings = i+1
                print_log(f'''GET_NUMBER_OF_RINGS: Not both rings warrant cutting, rings = {new_rings}
''',Configuration,case=['debug_add'])
                break
    else:
        print_log(f'''GET_NUMBER_OF_RINGS: Not both last rings are below the limit
''',Configuration,case=['debug_add'])
        # if they are not we first check wether both second to last rings are
        if ((difference_with_limit[0,-2] < 0.) and (sbr[0,-1] < 2*sbr_ring_limits[0,-1]) and (difference_with_limit[1,-1] < 0.)) or \
            ((difference_with_limit[1,-2] < 0.) and (sbr[1,-1] < 2*sbr_ring_limits[1,-1]) and (difference_with_limit[0,-1] < 0.)) or\
            ((difference_with_limit[0,-2] < 0.) and (difference_with_limit[1,-2] < 0.) and (sbr[0,-1] < 3*sbr_ring_limits[0,-1]) and (sbr[1,-1] < 3*sbr_ring_limits[1,-1])):
            new_rings -= 1
            print_log(f'''GET_NUMBER_OF_RINGS: A second ring is too faint, rings = {new_rings}
''',Configuration,case=['debug_add'])
        elif np.all(difference_with_limit[:,-2] < 0.) and np.all(sbr[:,-1] < 5*sbr_ring_limits[:,-1]):
            new_rings -= 1
            print_log(f'''GET_NUMBER_OF_RINGS: Both second rings are too faint, rings = {new_rings}
''',Configuration,case=['debug_add'])
        else:
            print_log(f'''GET_NUMBER_OF_RINGS: The second rings are too bright and do not allow for cutting.
''',Configuration,case=['debug_add'])
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


 OUTPUTS:
    new_rings = the required number of rings

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_ring_weights(Configuration,Tirific_Template):
    print_log(f'''GET_RING_WEIGTHS: Getting the importance of the rings in terms of SBR.
''',Configuration,case=['debug_start'])
    sbr = sf.load_tirific(Configuration,Tirific_Template,Variables=["SBR",f"SBR_2"],\
                array=True )
    
    print_log(f'''GET_RING_WEIGTHS: retrieve this sbr.
{'':8s} sbr = {sbr}
''',Configuration,case=['debug_add'])
    cut_off_limits = sf.sbr_limits(Configuration,Tirific_Template )
    print_log(f'''GET_RING_WEIGTHS: retrieved these cut_off_limits.
{'':8s} col = {cut_off_limits}
''',Configuration,case=['debug_add'])
    sm_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',
                            min_error= [cut_off_limits,cut_off_limits],no_apply = True,
                            fix_sbr_call = True, profile_in = sbr )
    print_log(f'''GET_RING_WEIGTHS: retrieve this sbr.
{'':8s} sbr = {sbr}
{'':8s} using sm_sbr {sm_sbr}
''',Configuration,case=['debug_add'])

    weights= [[],[]]
    for i in [0,1]:
        weights[i] = [sf.set_limits(x/y,0.1,10.) for x,y in zip(sm_sbr[i],cut_off_limits)]
        weights[i] = weights[i]/np.nanmax(weights[i])
        weights[i][0:2] = np.nanmin(weights[i])
        weights[i] = [sf.set_limits(x,0.1,1.) for x in weights[i]]
    print_log(f'''GET_RING_WEIGTHS: Obtained the following weights.
{'':8s}{weights}
''',Configuration,case=['debug_add'])
    return np.array(weights,dtype = float)

get_ring_weights.__doc__=f'''
 NAME:
    get_ring_weights

 PURPOSE:
    Get the importance of the rings based on how much the SBR lies above the noise limit for each ring

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard Tirific template containg the SBR profiles

 OPTIONAL INPUTS:


 OUTPUTS:
    numpy array with the weight normalized to the the maximum.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    Weight 1 is most important, weight 0. least important.
    Errors should be divided by these weights to reflect the importance
'''
def get_warp_slope(Configuration,Tirific_Template):
    print_log(f'''GET_WARP_SLOPE: We have {Tirific_Template['NUR']} rings in the template. and this should be {Configuration['NO_RINGS']}
''', Configuration,case=['debug_start'])
    sbr_ring_limits = sf.sbr_limits(Configuration,Tirific_Template)
    #get the sbr profiles
    sbr = sf.load_tirific(Configuration,Tirific_Template,Variables=['SBR','SBR_2'],\
                array=True)
    print_log(f'''GET_WARP_SLOPE: We have {len(sbr_ring_limits)} rings in our limits.
{'':8s}GET_WARP_SLOPE: And we have {len(sbr[0])} rings in our profiles.
''', Configuration,case=['debug_add'])
    warp_slope = [int(Tirific_Template['NUR']),int(Tirific_Template['NUR'])]
    sbr_ring_limits = 1.5*np.array([sbr_ring_limits,sbr_ring_limits],dtype=float)
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
            final = Configuration['LAST_RELIABLE_RINGS'][i]-1
            for parameter in ['INCL',"PA",'SDIS','Z0']:
                if parameter not in Configuration['FIXED_PARAMETERS'][0]:
                    Configuration['FIXED_PARAMETERS'][0].append(parameter)
        else:
            final = len(slope)

        if final > Configuration['LAST_RELIABLE_RINGS'][i]-1:
            final = Configuration['LAST_RELIABLE_RINGS'][i]-1
        # we have to vary some rings else tirific will break with a segmentation fault
        if final < 2:
            final = 2

        warp_slope[i] = int(final)
    print_log(f'''GET_WARP_SLOPE: We find a slope of {warp_slope}.
''', Configuration,case=['debug_add'])
    Configuration['WARP_SLOPE'] = warp_slope
    incl = sf.load_tirific(Configuration,Tirific_Template,Variables=['INCL','INCL_2'],\
            array=True)
    if np.mean(incl[:,:int(Configuration['NO_RINGS']/2.)]) < 35. :
        if 'INCL' not in Configuration['FIXED_PARAMETERS'][0]:
            Configuration['FIXED_PARAMETERS'][0].append('INCL')
    else:
        if 'INCL' in Configuration['FIXED_PARAMETERS'][0] and 'INCL' not in Configuration['FIXED_PARAMETERS'][1]:
            Configuration['FIXED_PARAMETERS'][0].remove('INCL')
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


 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def inner_sbr_fix(Configuration,sbr,cutoff_limits ):
    print_log(f'''INNER_SBR_FIX: Checking the SBR inner points for runaway values
{'':8s} sbr in  = {sbr}
''',Configuration,case=['debug_start'])

    if np.all(sbr[:,0] > 2*sbr[:,2]) or np.all(sbr[:,1] > 2*sbr[:,2]):
        print_log(f'''INNER_SBR_FIX: We need to correct
{'':8s} sbr 0  = {sbr[:,0]} sbr 1  = {sbr[:,1]} sbr 2  = {sbr[:,2]}
''',Configuration,case=['debug_add'])
        if np.mean(sbr[:,2]) > cutoff_limits[0,2]:
            sbr[:,[0,1]] = np.mean(sbr[:,2])
            print_log(f'''INNER_SBR_FIX: We need to correct with mean
{'':8s} mean = {np.mean(sbr[:,2])}
''',Configuration,case=['debug_add'])
        else:
            print_log(f'''INNER_SBR_FIX: We need to correct with cut_off_limits
{'':8s} limit = {1.5*cutoff_limits[0,2]}
''',Configuration,case=['debug_add'])
            sbr[:,[0,1,2]] = 1.5*cutoff_limits[0,2]
    if np.any(sbr[:,0] > sbr[:,1]):
        print_log(f'''INNER_SBR_FIX: We correct 0 point
{'':8s} mean 1 = {np.mean(sbr[:,1])}
''',Configuration,case=['debug_add'])
        sbr[:,0] = np.mean(sbr[:,1])

    for i in [0,1]:
        if np.any(sbr[:,i] < cutoff_limits[:,2]):
            print_log(f'''INNER_SBR_FIX: correcting ring {i}
''',Configuration,case=['debug_add'])
            sbr[:,i] = 1.5*cutoff_limits[0,2]

    print_log(f'''INNER_SBR_FIX: the fixed sbr {sbr}
''',Configuration,case=['debug_add'])

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


 OUTPUTS:
    sbr = profile with the modified inner points if required

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def modify_flat(Configuration,profile,original_profile,errors,key,inner_fix = None):
    if inner_fix is None:
        inner_fix = [4,4]
    print_log(f'''MODIFY_FLAT: These {key} profiles are checked to be flat.
{'':8s} profile = {profile}
{'':8s} original_profile = {original_profile}
{'':8s} errors = {errors}
''',Configuration,case=['debug_start'])

    flatness = []


    for side in [0,1]:
         flatness.append(check_flat(Configuration,profile[side],errors[side],key,inner_fix=inner_fix[side]\
                                    ,last_reliable_ring= Configuration['LAST_RELIABLE_RINGS'][side]))
    print_log(f'''MODIFY_FLAT: Side 0 is flat = {flatness[0]}
{'':8s} Side 1 is flat = {flatness[1]}
''',Configuration,case=['debug_add'])

    if all(flatness):
        if key in ['PA','INCL']:
            profile[:] = profile[0,0]
        else:
            profile[:] = np.nanmedian(original_profile[:,:round(len(original_profile)/2.)])
        errors = get_error(Configuration,original_profile,profile,key,apply_max_error = True,min_error =np.nanmin(errors))
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
                        errors[side] = get_error(Configuration,\
                            original_profile[side],profile[side],key,\
                            apply_max_error = True,\
                            min_error =np.nanmin(errors[side]),singular = True )
                profile[:,0:3] = flat_val
                profile[:,4] = (flat_val+profile[:,4])/2.
            else:
                profile[:] = np.nanmedian(original_profile[:,:round(len(original_profile)/2.)])
                errors = get_error(Configuration,original_profile,profile,key,\
                        apply_max_error = True,min_error =np.nanmin(errors))

    print_log(f'''MODIFY_FLAT: Returning:
{'':8s} profile = {profile}
{'':8s} errors = {errors}
''',Configuration,case=['debug_add'])
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

    fit_type = 'Undefined'

 OUTPUTS:
    the final profile and errors

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def no_declining_vrot(Configuration, Tirific_Template, profile = None):
    print_log(f'''NO_DECLINING_VROT: make RC flat from highest point on.
{'':8s}NO_DECLINING_VROT: But only for low value RCs
''',Configuration, case=['debug_start'])
    no_input = False
    if profile is None:
        no_input = True
        profile = sf.load_tirific(Configuration,Tirific_Template,Variables=['VROT'],\
                    array=True)
    RCval = np.mean(profile[2:])
    RCmax = np.where(profile == np.max(profile))[0]
    if len(RCmax) > 1:
        RCmax = RCmax[0]
    print_log(f'''NO_DECLINING_VROT: We find the maximum at ring {RCmax}
{'':8s}NO_DECLINING_VROT: And a mean value of {RCval}.
''',Configuration,case=['debug_add'])
    Configuration['OUTER_SLOPE_START'] = Configuration['NO_RINGS']-1
    if RCmax < len(profile)/2. or RCval > 180.:
        print_log(f'''NO_DECLINING_VROT: We shan't adapt the RC
''',Configuration,case=['debug_add'])
    else:
        for i in range(int(len(profile)/2.),len(profile)-1):
            if profile[i+1] < profile[i]:
                profile[i:] =profile[i]
                print_log(f'''NO_DECLINING_VROT: Flattening from ring {i} on.)
    ''',Configuration,case=['debug_add'])
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
    format = sf.set_format('VROT')
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


def regularise_profile(Configuration,Tirific_Template, key,min_error= None, \
                            no_apply =False,warp=False):
    if min_error is None:
        if key in Configuration['MIN_ERROR']:
            min_error= Configuration['MIN_ERROR'][key]
        else: 
            min_error = [0.]
    if key in ['PA','INCL'] and not warp:
        raise FunctionCallError('The warp is regularised in regularise_warp. ')
    # We start by getting an estimate for the errors
    min_error=np.array(min_error,dtype=float)
    weights = get_ring_weights(Configuration,Tirific_Template)
    #For a lot of things it is important to use a unmodified input profile but we have to avoid the declining vrot
    if key == 'VROT':
        no_declining_vrot(Configuration, Tirific_Template)
    profile = sf.load_tirific(Configuration,Tirific_Template, [key,f"{key}_2"],array=True)
    diff = abs(np.sum(profile[0]-profile[1]))#Check that we have two profiles

    

    #First if we have an RC we flatten the curve
    print_log(f'''REGULARISE_PROFILE: profile of {key} before regularistion
{'':8s}{profile[0]}
{'':8s}{profile[1]}
{'':8s}The minimal error is
{'':8s}{min_error}
''',Configuration,case=['debug_start'])
    #if the profiles is fixed we do not have to regularise and we return a flatted profile
    if key in Configuration['FIXED_PARAMETERS'][0]:
        for i in [0,1]:
            profile[i][:]=np.mean(profile[i])  
        set_errors(Configuration,Tirific_Template,key,min_error=min_error)
        return profile

    # get a smoothed profiles
    sm_profile = smooth_profile(Configuration,Tirific_Template, key \
                                ,min_error=min_error,no_apply=True)

   
    error = get_error(Configuration,profile,sm_profile,key,\
                        min_error=min_error,weights = weights)
    if key in ['PA']:
        for i in [0,1]:
            error[i,:] = [np.mean([np.sqrt(x),min_error[0]]) if \
                          np.sqrt(x) > min_error[0] else x for x in error[i,:]]

  
    if key in ['SDIS','VROT']:
        diff = 0.
    if diff > 1e-8:
        print_log(f'''REGULARISE_PROFILE: Treating both sides independently.
''',Configuration,case=['debug_add'])
        sides = [0,1]
    else:
        print_log(f'''REGULARISE_PROFILE: Found symmetric profiles.
''',Configuration,case=['debug_add'])
        sides = [0]
        error[0] = np.array([np.mean([x,y]) for x,y in \
                                zip(error[0],error[1])],dtype=float)

    radii =sf.set_rings(Configuration)
    for i in sides:
        if key in ['SDIS']:
            function_to_fit =arc_tan_sdis_function
            fit_profile = fit_arc(Configuration,radii,sm_profile[i],\
                                    error[i],function_to_fit,key)
        else:

            fit_profile = fit_polynomial(Configuration,radii,\
                profile[i],error[i],key, Tirific_Template,\
                inner_fix = Configuration['INNER_FIX'][i],\
                boundary_limits= Configuration[f"{key}_CURRENT_BOUNDARY"][i+1])
        profile[i] = fit_profile

    if diff < 1e-8:
        profile[1] = profile[0]

    original = np.array(sf.load_tirific(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype=float)
    error = get_error(Configuration,original,profile,key,weights=weights,apply_max_error = True,min_error=min_error)
    print_log(f'''REGULARISE_PROFILE: This the fitted profile without corrections:
{'':8s}{profile}
''',Configuration,case=['debug_add'])
#then we want to fit the profiles with a polynomial
    if key not in ['SBR','VROT','SDIS','PA','INCL']:
        #We should not fix the profile again as the fitted profile is fixed should be good
        profile,error = modify_flat(Configuration, profile, original, error,key,inner_fix= Configuration['INNER_FIX'])
    format = sf.set_format(key)

    if not no_apply:
        Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in profile[0,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:{format}}' for x in profile[1,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in error[0,:int(Configuration['NO_RINGS'])]])}")
        Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in error[1,:int(Configuration['NO_RINGS'])]])}")
        #if key in ['INCL'] and np.mean( profile[:,int(Configuration['NO_RINGS']/2.):int(Configuration['NO_RINGS'])]) < 40.:
        #    fix_vrot_for_incl_change(Configuration,Tirific_Template,original,profile)

        print_log(f'''REGULARISE_PROFILE: And this has gone to the template.
{'':8s}{key} = {Tirific_Template[key]}
{'':8s}{key}_2 ={Tirific_Template[f"{key}_2"]}
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration,case=['debug_add'])
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


def regularise_warp(Configuration,Tirific_Template, min_error = None, \
                        no_apply =False,smooth_only=False):
    ''' 
    Initially we did this by fitting a polynomial to the change_angle and
    The factors but that sometimes led to very bad results. The we tried
    to merely use smoothing on the factors but that leaves saw tooth
    if we smooth and fit Theta and Phi individually this can lead to weird results
    We cannot smooth the factors either they have to fitted with a polynomial
    The fixed part  should be equal to fix+1 
    Well it seems that the best option is to regularise the the profiles and 
    only use an angular momentum check during the fitting.
    but as we are running it every iteration we are not exactly sure what happens 
    in the fitting
    '''
     # We start by getting an estimate for the errors
    if min_error is None:
        if 'PA' in Configuration['MIN_ERROR'] and 'INCL' in \
            Configuration['MIN_ERROR']:
            min_error = [Configuration['MIN_ERROR']['PA'],\
                         Configuration['MIN_ERROR']['INCL']]
        else:
            min_error = [0.,0.]

    min_error=np.array(min_error,dtype=float)
    weights = get_ring_weights(Configuration,Tirific_Template)
   
    #if the PA and INCL are fixed we stop here
    if 'PA' in Configuration['FIXED_PARAMETERS'][0] and 'INCL' in \
          Configuration['FIXED_PARAMETERS'][0]:
        profile = np.array(sf.load_tirific(Configuration,Tirific_Template, [f"PA",f"PA_2",f"INCL",f"INCL_2"]),dtype=float)
        for i in [0,1,2,3]:
            profile[i][:]=np.mean(profile[i])
        for i,key in enumerate(['PA','INCL']):  
            set_errors(Configuration,Tirific_Template,key,min_error=min_error[i])
        return profile


   
    if smooth_only:
        pa_both_sides = smooth_profile(Configuration,Tirific_Template, 'PA' ,\
                    min_error=min_error[0]/3.,no_fix=False,no_apply=True)
    
        incl_both_sides = smooth_profile(Configuration,Tirific_Template, 'INCL' ,\
                    min_error=min_error[1]/3.,no_fix=True,no_apply=True)

    else:
        pa_both_sides = regularise_profile(Configuration,Tirific_Template,\
                        'PA',warp=True,no_apply=True) 
    
        incl_both_sides = regularise_profile(Configuration,Tirific_Template,\
                        'INCL',warp=True,no_apply=True) 
    
    profile= np.concatenate((pa_both_sides,incl_both_sides))  
    diff = abs(np.sum(profile[0]-profile[1]))+abs(np.sum(profile[2]-profile[3]))  #Check that we have two profiles
    if diff < 1e-8:
        print_log(f'''REGULARISE_WARP: Found symmetric profiles.
''',Configuration,case=['debug_add'])
        sides = [0]  
    else:
        print_log(f'''REGULARISE_WARP: Treating both sides independently.
''',Configuration,case=['debug_add'])      
        sides = [0,1]
    radii =sf.set_rings(Configuration)    
    for i in sides:
        #put the pa and inc in a single profile
        one_side_profile = [profile[i],profile[i+2]] 
        print_log(f'''REGULARISE_WARP: profile of the warp for side {i} before regularistion
{'':8s}PA = {one_side_profile[0]}
{'':8s}INCL = {one_side_profile[1]}
{'':8s}The minimal error is
{'':8s}PA min Error = {min_error[0]}, INCL min Error = {min_error[1]}
''',Configuration,case=['debug_start'])
        #we need to fix them in the inner part before transforming
        pa = fix_profile(Configuration, 'PA', \
            one_side_profile[0], inner_fix = Configuration['INNER_FIX'][i],\
            singular = True,only_inner = True )
        inclination = fix_profile(Configuration, 'INCL', \
            one_side_profile[1], inner_fix = Configuration['INNER_FIX'][i],\
            singular = True,only_inner = True )
        #first fit PA and INCL indpendently
        #pa,inclination,changed_angles = sf.check_angular_momentum_vector(Configuration,\
        #    radii,one_side_profile[0],one_side_profile[1],side=i,regularise=True)

        
        if smooth_only:
        # get Phi and Theta for both sides
            Theta,Phi,multiple= sf.calculate_am_vector(Configuration,
                                    one_side_profile[0],one_side_profile[1])
            
            
            high_theta,high_phi,high_multiple = sf.calculate_am_vector(\
                Configuration,one_side_profile[0]+min_error[0],\
                one_side_profile[1]+min_error[1])
            if not np.array_equiv(high_multiple,multiple):
                high_theta,high_phi,high_multiple = sf.calculate_am_vector(\
                Configuration,one_side_profile[0]-min_error[0],\
                one_side_profile[1]-min_error[1])

            min_change_error = np.sqrt(np.min(np.array([abs(x-y)\
                for x,y in zip(Theta,high_theta)],dtype=float))**2+\
                np.min(np.array([abs(x-y) for x,y in zip(Phi,high_phi)],dtype=float))**2)

            print_log(f'''REGULARISE_WARP: These converted arrays for side {i}
{'':8s}Theta = {Theta}
{'':8s}Phi = {Phi}
{'':8s}The minimal error is
{'':8s}min change Error = {min_change_error}
''',Configuration,case=['debug_add'])
            multiple_zero = multiple[0]
            
            #Calculate the change angle
            change_angle = sf.calculate_change_angle(Configuration,Theta,Phi)
                #if this side is not warping we have to skip
            if np.sum(change_angle['CHANGE_ANGLE']) == 0.:
                continue
            change_boundary = calculate_change_boundary(Configuration,multiple_zero,\
                        change_angle['THETA']['ZERO'],change_angle['PHI']['ZERO'],\
                        pa_boun =Configuration['PA_CURRENT_BOUNDARY'][i+1],\
                        incl_boun =Configuration['INCL_CURRENT_BOUNDARY'][i+1])
            sm_change_angle = smooth_profile(Configuration,Tirific_Template, 'CHANGE_ANGLE' \
                ,profile_in=change_angle['CHANGE_ANGLE'],\
                min_error=min_change_error,no_apply=True,no_fix =True,singular=True)
            

            sm_theta_factor,sm_phi_factor = smooth_factor_profile(Configuration, \
                change_angle,radii,inner_fix=Configuration['INNER_FIX'][i],\
                Tirific_Template=Tirific_Template)
          
           
            error_change_angle = get_error(Configuration,change_angle['CHANGE_ANGLE'],\
                sm_change_angle,'CHANGE_ANGLE',min_error=min_change_error,\
                weights = weights,singular=True)
            
            print_log(f'''REGULARISE_WARP: For side {i} we  regularise the following profile.
        {'':8s} change_angle = {change_angle['CHANGE_ANGLE'][i]}
        {'':8s} sm_change_angle = {sm_change_angle[i]}
        ''',Configuration,case=['debug_add'])
            #if smooth_only:
            new_change_angle = sm_change_angle
            #    new_theta_factor = sm_theta_factor
            #    new_phi_factor = sm_phi_factor
            #else:
                
            #    new_change_angle,fit_order = fit_polynomial(Configuration,radii,\
            #        change_angle['CHANGE_ANGLE'],error_change_angle,'CHANGE_ANGLE', Tirific_Template,\
            #        inner_fix = Configuration['INNER_FIX'][i],\
            #        allowed_order= [3,6],boundary_limits= change_boundary, return_order=True)
                
            print_log(f'''REGULARISE_WARP:
        {'':8s} new_change_angle = {new_change_angle}
        {'':8s} new_theta_factor = {sm_theta_factor}
        {'':8s} new_phi_factor  = {sm_phi_factor}
        ''',Configuration,case=['debug_add'])
            new_change_dict = copy.deepcopy(change_angle)
            new_change_dict['CHANGE_ANGLE'] = new_change_angle
            new_change_dict['THETA']['FACTOR'] = sm_theta_factor
            new_change_dict['PHI']['FACTOR'] = sm_phi_factor
            #tyr= True
            #if tyr:
            #    plot_output(radii,change_angle,new_change_dict,i)
                

                                
            print_log(f'''REGULARISE_WARP: For side {i} we get
        {'':8s} Original Theta = {Theta}
        {'':8s} Original Phi = {Phi}
        ''',Configuration,case=['debug_add']) 
                    
            Theta,Phi = sf.revert_change_angle(Configuration,new_change_dict)
            
            print_log(f'''REGULARISE_WARP: For side {i} we get
        {'':8s} Theta = {Theta}
        {'':8s} Phi = {Phi}
        ''',Configuration,case=['debug_add'])
            
            pa,inclination=sf.calculate_am_vector(Configuration,Theta,\
                                Phi,multiple = multiple, invert=True)
    
        pa[:Configuration['INNER_FIX'][i]] = np.mean(pa[:Configuration['INNER_FIX'][i]])
        inclination[:Configuration['INNER_FIX'][i]] = np.mean(inclination[:Configuration['INNER_FIX'][i]])
        profile[i]=pa
        profile[i+2]=inclination
        #if i == 0:
        #    print('Decidedly odd')
        #    print(pa,inclination)
            
           
    
    if np.sum(sides) == 0.:
        profile[1] = profile[0]
        profile[3] = profile[2]
    else:
        #print(profile,Configuration['INNER_FIX'])
        fixed = np.mean(np.concatenate((profile[0,:Configuration['INNER_FIX'][0]],
                          profile[1,:Configuration['INNER_FIX'][1]])))
        profile[0,:Configuration['INNER_FIX'][0]]= fixed
        profile[1,:Configuration['INNER_FIX'][1]]= fixed
        fixed = np.mean(np.concatenate((profile[2,:Configuration['INNER_FIX'][0]],
                          profile[3,:Configuration['INNER_FIX'][1]])))
        profile[2,:Configuration['INNER_FIX'][0]]= fixed
        profile[3,:Configuration['INNER_FIX'][1]]= fixed
       
    error= copy.deepcopy(profile)
    error[:] = 0
   
    for i,key in enumerate(['PA','INCL']):
        original = np.array(sf.load_tirific(Configuration,Tirific_Template,\
                                             [key,f"{key}_2"]),dtype=float)
        sm_profile = smooth_profile(Configuration,Tirific_Template, key ,\
                    profile_in=profile[2*i:2*i+2], min_error=min_error[i]/3.,\
                    no_fix=True,no_apply=True)
        sm_error =  get_error(Configuration,original,sm_profile,key,\
                weights=weights,apply_max_error = True,\
                min_error=min_error[i]/3.)
        
        # As there is too much in flux we do not want to flatten the warp in the first iteration
        #if Configuration['ITERATIONS'] > 1:       
        #    sm_profile,sm_error = modify_flat(Configuration, profile[2*i:2*i+2],\
        #             original, sm_error,key,inner_fix= Configuration['INNER_FIX'])
    
        profile[2*i:2*i+2]=sm_profile
        error[2*i:2*i+2] = get_error(Configuration,original,profile[2*i:2*i+2],\
                    key,weights=weights,apply_max_error = True,\
                    min_error=min_error[i])
        format = sf.set_format(key)
        if not no_apply:
            Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in profile[2*i,:int(Configuration['NO_RINGS'])]])}"
            Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:{format}}' for x in profile[2*i+1,:int(Configuration['NO_RINGS'])]])}"
            Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in error[2*i,:int(Configuration['NO_RINGS'])]])}")
            Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in error[2*i+1,:int(Configuration['NO_RINGS'])]])}")

            print_log(f'''REGULARISE_WARP: This has gone to the template.
{'':8s}{key} = {Tirific_Template[key]}
{'':8s}{key}_2 ={Tirific_Template[f"{key}_2"]}
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration,case=['debug_add'])
    
    return profile

regularise_warp.__doc__ =f'''
 NAME:
    regularise_warp

 PURPOSE:
    Regularise a the PA and INCL by smoothing the AM vector

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    key = parameter to fix

 OPTIONAL INPUTS:


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
def plot_output(radius, input_parameters, output_parameters,i):
        colors = ['r','b','g','k']
        plt.figure(89,figsize=(16,16),dpi=300,facecolor = 'w', edgecolor = 'k')
        plt.subplot(2,1,1)
        count=0
                    

            
        plt.plot(radius,input_parameters['CHANGE_ANGLE'],c=colors[count],label='Change Angle' )
        plt.plot(radius,output_parameters['CHANGE_ANGLE'],c=colors[count],linestyle='--',lw=3 )

        plt.subplot(2,1,2)
        plt.plot(radius,input_parameters['THETA']['FACTOR'],c=colors[count],label='THETA_FACTOR' )
        plt.plot(radius,output_parameters['THETA']['FACTOR'],c=colors[count],linestyle='--',lw=3 )
        count =1
        plt.plot(radius,input_parameters['PHI']['FACTOR'],c=colors[count],label='PHI_FACTOR' )
        plt.plot(radius,output_parameters['PHI']['FACTOR'],c=colors[count],linestyle='--',lw=3 )


        plt.subplot(2,1,1)    
        plt.legend()
        plt.subplot(2,1,2)    
        plt.legend()
        plt.savefig(f"/home/peter/FAT_Main/Component_Test/Angular_Momentum/Comparison_Angle_{i}.png", bbox_inches='tight')
        plt.close()

def set_boundary_limits(Configuration,Tirific_Template,key,values = None, \
                        tolerance = 0.01, fixed = False):
    if values is None:
        values = [0.,0.]
    print_log(f'''SET_BOUNDARY_LIMITS: checking limits for {key},
{'':8s} current Boundaries = {Configuration[f"{key}_CURRENT_BOUNDARY"]}
{'':8s} values = {values}
''',Configuration,case=['debug_start'])
    profile = np.array(sf.load_tirific(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype = float)
    current_boundaries = Configuration[f"{key}_CURRENT_BOUNDARY"]

    if key == 'VROT' and np.sum(values) != 0.:
        current_boundaries =[[values[0]-values[1]*5.,values[0]+values[1]*5.] for x in range(3)]
    if np.sum(current_boundaries) == 0. and np.sum(values) != 0.:
        current_boundaries =[[values[0]-values[1],values[0]+values[1]] for x in range(3)]
        print_log(f'''SET_BOUNDARY_LIMITS: set boundaries from values
{'':8s} current Boundaries = {current_boundaries}
''',Configuration,case=['debug_add'])
    elif np.sum(current_boundaries) == 0. and np.sum(values) == 0.:
        raise FunctionCallError(f'SET_BOUNDARY_LIMITS: if boundaries are not set in the Configuration call set_boundary_limits with values')
    else:
        print_log(f'''SET_BOUNDARY_LIMITS: Checking
{'':8s} current Boundaries = {current_boundaries}
{'':8s} Applying to the profiles {profile}
''',Configuration,case=['debug_add'])

        if fixed:
            range_to_check = [0]
        else:
            range_to_check = [0,1,2]
        for i in range_to_check:

            buffer = np.array(float(current_boundaries[i][1]-current_boundaries[i][0])*2,dtype=float)\
                    *np.array([tolerance,0.25],dtype=float)
            print_log(f'''SET_BOUNDARY_LIMITS: Using a buffer of {buffer[0]} and a change of {buffer[1]}.
''',Configuration,case=['debug_add'])
            if i == 0:
                profile_part = profile[0,:int(np.min(Configuration['INNER_FIX']))+1]
                infix = np.min(Configuration['INNER_FIX'])
            else:
                profile_part = profile[i-1,Configuration['INNER_FIX'][i-1]:]
                infix = Configuration['INNER_FIX'][i-1]
                #Configuration['LAST_RELIABLE_RINGS'][i-1]]
            print_log(f'''SET_BOUNDARY_LIMITS: Checking {profile_part}.
{'':8s} because for this part the inner fix = {infix}
    ''',Configuration,case=['debug_add'])
            #check the upper bounderies
            on_boundary = np.where(profile_part > float(current_boundaries[i][1])-buffer[0])[0]
            print_log(f'''SET_BOUNDARY_LIMITS: Found the following on the upper {on_boundary}.
    ''',Configuration,case=['debug_add'])
            if len(on_boundary) > 0:
                if on_boundary[0] != len(profile[0])-1:
                    current_boundaries[i][1] = current_boundaries[i][1] + buffer[1]
            #check the lower boundaries.
            on_boundary = np.where(profile_part < float(current_boundaries[i][0])+buffer[0])[0]
            print_log(f'''SET_BOUNDARY_LIMITS: Found the following on the lower {on_boundary}.
    ''',Configuration,case=['debug_add'])
            if len(on_boundary) > 0:
                if on_boundary[0] != len(profile[0])-1:
                    current_boundaries[i][0] = current_boundaries[i][0] - buffer[1]
    low = [x[0] for x in current_boundaries]
    high= [x[1] for x in current_boundaries]
    if key == 'Z0':
        inc = np.array(sf.load_tirific(Configuration,Tirific_Template, ['INCL']),dtype = float)
         # Linear increase of the maximum from 30-70 inc
        min,max = sf.convertskyangle(Configuration,[0.2,sf.set_limits(float(inc[0])/20.-1,0.5,2.5)], physical = True)
        high = [sf.set_limits(x,min,max) for x in high]
        min,max = sf.convertskyangle(Configuration,[0.05,0.2], physical = True)
        low = [sf.set_limits(x,min,max) for x in low]
    elif key == 'INCL':
        # for inclination the other boundaries are bracketed by the input boundary
        high = [x if x > 50. else 50. for x in high]
        low = [x if x < 60. else 60. for x in low]
    elif key == 'PA':
        # for inclination the other boundaries are bracketed by the input boundary
        high = [x if x > 170. else 170. for x in high]
        low = [x if x < 190. else 190. for x in low]
        if abs(high[1]-high[2]) > 30:
            high[1:] = [np.max(high[1:]) for x in high[1:]]
        if abs(low[1]-low[2]) > 30:
            low[1:] = [np.min(low[1:]) for x in low[1:]]

    print_log(f'''SET_BOUNDARY_LIMITS:We use these low = {low} and these high {high} to set {Configuration[f"{key}_CURRENT_BOUNDARY"]}.
''',Configuration,case=['debug_add'])
    sf.set_boundaries(Configuration,key,low,high)
    print_log(f'''SET_BOUNDARY_LIMITS: We have adjusted the boundaries to  {Configuration[f"{key}_CURRENT_BOUNDARY"]}.
''',Configuration,case=['debug_add'])

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


    tolerance = 0.01,
    The minimum difference with the boundary

    fixed = False
    If true it is assumed there is only one boundary value to check else three


 OUTPUTS:
    The new Boundaries, they are also updated in Configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_cflux(Configuration,Tirific_Template):

    if any(np.isnan(Configuration['NO_POINTSOURCES'])):
        print_log(f'''SET_CFLUX: We detected an infinite number of model point sources.
{"":8s}SET_CFLUX: This must be an error. Exiting the fitting.
''',Configuration,case=['main','screen'])
        raise CfluxError('The model had infinite point sources')
    if np.max(Configuration['SIZE_IN_BEAMS']) < 15:
        factor = 1.
    else:
        factor=(np.max(Configuration['SIZE_IN_BEAMS'])/15.)**1.5
    triggered = 0
    if not 0.5e6 < Configuration['NO_POINTSOURCES'][0] < 2.2e6:
        new_cflux = sf.set_limits(float(Tirific_Template['CFLUX'])*Configuration['NO_POINTSOURCES'][0]/(factor*1e6),1e-7,5e-3)
        print_log(f'''SET_CFLUX: CFLUX is adapted from {Tirific_Template['CFLUX']} to {new_cflux:.2e}
''',Configuration,case=['verbose'])
        Tirific_Template['CFLUX'] = f"{new_cflux:.2e}"
        triggered = 1
    if not 0.5e6 < Configuration['NO_POINTSOURCES'][1] < 2.2e6:
        new_cflux = sf.set_limits(float(Tirific_Template['CFLUX_2'])*Configuration['NO_POINTSOURCES'][1]/(factor*1e6),1e-7,5e-3)
        print_log(f'''SET_CFLUX: CFLUX_2 is adapted from {Tirific_Template['CFLUX_2']} to {new_cflux:.2e}
''',Configuration,case=['verbose'])
        Tirific_Template['CFLUX_2'] = f"{new_cflux:.2e}"
        triggered = 1
    if not triggered:
        print_log(f'''SET_CFLUX: CFLUXES are within the required limits.
''',Configuration,case=['verbose'])
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


 OUTPUTS:
    Updated CFLUX values in template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_errors(Configuration,Tirific_Template,key,min_error = 0.):
    error = np.full((2,int(Configuration['NO_RINGS'])),min_error)
    format=sf.set_format(key)
    Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in error[0,:int(Configuration['NO_RINGS'])]])}")
    Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in error[1,:int(Configuration['NO_RINGS'])]])}")

    print_log(f'''SET_ERRORS: This has gone to the template.
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration,case=['debug_start'])
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


    min_error =0 .
    the error to be used for all rings

 OUTPUTS:
    Updated Tirific Template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_fitting_parameters(Configuration, Tirific_Template, \
        parameters_to_adjust= None, modifiers = None, stage = 'initial',\
         initial_estimates = None):
    if parameters_to_adjust is None:
        parameters_to_adjust = ['NO_ADJUSTMENT']
    if modifiers is None:
        modifiers = ['EMPTY']
    if initial_estimates is None:
        initial_estimates = ['EMPTY']
    print_log(f'''SET_FITTING_PARAMETERS: We are starting with these modifiers.
{'':8s} {modifiers}
''',Configuration,case=['debug_start'])
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
    if Configuration['INSTALLATION_CHECK']:
        Tirific_Template['LOOPS'] = "1"
    elif Configuration['ITERATIONS'] < 3:
        Tirific_Template['LOOPS'] = "10"
    else:
        Tirific_Template['LOOPS'] = f"{Configuration['LOOPS']}"
    fitting_settings = {}
    fitting_keys = ['VARY','VARINDX','MODERATE','DELEND','DELSTART','MINDELTA','PARMAX','PARMIN']

    if 'INCL' not in initial_estimates:
        try:
            profile = np.array([np.mean([x,y]) for x,y in \
                zip(sf.load_tirific(Configuration,Tirific_Template, ['INCL']),\
                sf.load_tirific(Configuration,Tirific_Template, [f"INCL_2"]) )],dtype=float)  
            diff = abs(np.max(profile)-np.min(profile))/10.
        except ValueError:
            profile = sf.load_tirific(Configuration,Tirific_Template, ['INCL'],array=True)
            diff = abs(np.max(profile)-np.min(profile))/10.
        initial_estimates['INCL'] = [profile[0],sf.set_limits(diff,1.,5.)/np.sin(np.radians(profile[0]))]

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
            if Configuration['ITERATIONS'] == 0:
                parameters_to_adjust = ['VSYS','XPOS','YPOS','SBR','VROT','SDIS','INCL','PA','Z0']
            else:
                if initial_estimates['INCL'][0] < 30.:
                    parameters_to_adjust = ['SBR','PA','SDIS','INCL','VROT']
                elif initial_estimates['INCL'][0] < 50.:
                    parameters_to_adjust = ['SBR','PA','SDIS','VROT','INCL']
                elif initial_estimates['INCL'][0] > 75.:
                    parameters_to_adjust = ['VROT','SBR','PA','INCL','SDIS','Z0']
                else:
                    parameters_to_adjust = ['SBR','VROT','PA','INCL','SDIS','Z0']

                if Configuration['CENTRAL_CONVERGENCE']:
                    parameters_to_adjust = parameters_to_adjust+['VSYS','XPOS','YPOS']
                else:
                    parameters_to_adjust = ['VSYS','XPOS','YPOS']+parameters_to_adjust
        else:
            print_log(f'''SET_FITTING_PARAMETERS: No default adjustment for unknown stage.
''',Configuration)
            raise InitializeError('No default adjustment for unknown stage. ')

    for center in Configuration['CENTRAL_FIX']:
        if center in parameters_to_adjust:
            parameters_to_adjust.remove(center)

    #VSYS should always be in the initial estimates as sbr setting uses it
    for key in ['VSYS']:
        if key not in initial_estimates:
            profile = np.array(sf.load_tirific(Configuration,Tirific_Template, ['VSYS']),dtype=float)
            initial_estimates['VSYS'] = [profile[0],Configuration['CHANNEL_WIDTH']/2.]

    for key in parameters_to_adjust:
        if key not in initial_estimates:
            profile = np.array([np.mean([x,y]) for x,y in \
                                zip(sf.load_tirific(Configuration,Tirific_Template, [key]),\
                                sf.load_tirific(Configuration,Tirific_Template, [f"{key}_2"]))],dtype=float)
            diff = abs(np.max(profile)-np.min(profile))/10.
            if key == 'PA':
                initial_estimates['PA'] = [profile[0],sf.set_limits(diff,0.5,10)]
            elif key == 'VROT':
                initial_estimates['VROT'] = [np.nanmax(profile),sf.set_limits(np.nanstd(profile[1:]),np.nanmax(profile)-np.nanmin(profile[1:]),np.nanmax(profile))]
            elif key == 'XPOS':
                initial_estimates['XPOS'] = [profile[0],0.1*Configuration['BEAM_IN_PIXELS'][1]*Configuration['PIXEL_SIZE']]
            elif key == 'YPOS':
                initial_estimates['YPOS'] = [profile[0],0.1*Configuration['BEAM_IN_PIXELS'][1]*Configuration['PIXEL_SIZE']]
            elif key == 'SDIS':
                initial_estimates['SDIS'] = [np.mean(profile),Configuration['CHANNEL_WIDTH']]
            elif key == 'Z0':
                initial_estimates['Z0'] = sf.convertskyangle(Configuration,[0.2,0.05], physical = True)

        if key not in modifiers:
            print_log(f'''SET_FITTING_PARAMETERS: Adding {key} to the modifiers
''',Configuration,case=['debug_add'])
            if key == 'Z0': modifiers['Z0'] = [0.5,0.5,2.]
            elif key in ['XPOS','YPOS']:
                if Configuration['CENTRAL_CONVERGENCE']:
                    modifiers[key] = [0.5,0.5,1.]
                else:
                    modifiers[key] = [1.,1.,2.]
            elif key == 'VSYS':
                if Configuration['CENTRAL_CONVERGENCE']:
                    modifiers[key] = [0.5,0.5,1.]
                else:
                    modifiers[key] = [2.,0.5,0.5]
            elif stage in  ['initial','run_cc','after_cc','after_ec','after_os','final_os']:
                if key == 'INCL': modifiers['INCL'] = [1.,1.,1.]
                elif key == 'PA': modifiers['PA'] = [3.,1.,1.]
                elif key == 'SDIS': modifiers['SDIS'] =  [1.,1.,2.]
            else:
                if key == 'INCL': modifiers['INCL'] = [2.0,0.5,0.5]
                elif key == 'PA': modifiers['PA'] =[1.0,1.0,2.0]
                elif key == 'SDIS': modifiers['SDIS'] =  [1.,1.,0.5]

    for key in modifiers:
        print_log(f'''SET_FITTING_PARAMETERS: This {key} is in modifiers
{'':8s} With these values {modifiers[key]}
''',Configuration,case=['debug_add'])
    ###############-- Modfier adaptations for large galaxies ####################
    if Configuration['OUTER_RINGS_DOUBLED']:
        for par_to_change in ['XPOS','YPOS']:
            if par_to_change in modifiers:
                modifiers[par_to_change]  = np.array(modifiers[par_to_change],dtype=float)*\
                                            np.array([4.,3.,5.],dtype=float)
    if  Configuration['VEL_SMOOTH_EXTENDED']:
        if 'VSYS' in modifiers:
            modifiers['VSYS']  = np.array(modifiers['VSYS'],dtype=float)*\
                                        np.array([2.,2.,3.],dtype=float)
    ###############-- Modfier adaptations based on the inclination ####################
    if stage not in ['final_os']:
        if initial_estimates['INCL'][0] < 30.:
            if 'Z0' in modifiers:
                modifiers['Z0'][0:1] = np.array(modifiers['Z0'][0:1],dtype=float)*(0.2/0.5)
                modifiers['Z0'][2] = float(modifiers['Z0'][2])*1.5
            if 'INCL' in modifiers:
                modifiers['INCL'][0:1] =np.array(modifiers['INCL'][0:1],dtype=float)*(0.1/1.0)
                modifiers['INCL'][2] = float(modifiers['INCL'][2])*5.
            if 'SDIS' in modifiers:
                modifiers['SDIS'][0:1] = np.array(modifiers['SDIS'][0:1],dtype=float)*1.5
                modifiers['SDIS'][2] = float(modifiers['SDIS'][2])*0.5
            if 'VSYS' in modifiers:
                modifiers['VSYS'][0] = modifiers['VSYS'][0]*1.5
            print_log(f'''SET_FITTING_PARAMETERS: These are the  modifiers after correcting < 30.
{'':8s} {modifiers}
''',Configuration,case=['debug_add'])
        elif initial_estimates['INCL'][0] < 50.:
            if 'Z0' in modifiers:
                modifiers['Z0'][0:1] = np.array(modifiers['Z0'][0:1],dtype=float)*(0.4/0.5)
                modifiers['Z0'][2] = float(modifiers['Z0'][2])*1.2
            if 'INCL' in modifiers:
                modifiers['INCL'][0:1] =np.array( modifiers['INCL'][0:1],dtype=float)*(0.5/1.0)
                modifiers['INCL'][2] = float(modifiers['INCL'][2])*1.5
            print_log(f'''SET_FITTING_PARAMETERS: These are the  modifiers after correcting < 50.
{'':8s} {modifiers}
''',Configuration,case=['debug_add'])
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
            print_log(f'''SET_FITTING_PARAMETERS: These are the modifiers after correcting > 75.
{'':8s} {modifiers}
''',Configuration,case=['debug_add'])
        else:
            pass
    else:
        if initial_estimates['INCL'][0] < 40.:
            if 'INCL' in modifiers:
                modifiers['INCL'] = np.array(modifiers['INCL'],dtype=float)*(0.2)
    if stage in ['initial','initialize_os']:
        no_boun_check = ['SBR']
    else:
        no_boun_check =  ['SBR', 'INCL', 'PA']
    for key in parameters_to_adjust:
        if key in Configuration['FIXED_PARAMETERS'][0]:
            fixed=True
        else:
            fixed =False

        if key not in no_boun_check:
            bracket_values = initial_estimates[key]
            if key == 'PA':
                if initial_estimates['INCL'][0] > 35:
                    upper_limit = 20.
                else:
                    upper_limit= 90.
                bracket_values[1] = sf.set_limits(bracket_values[1]*10., 10.,upper_limit)
            if key == 'INCL':
                if bracket_values[0] > 40.:
                    bracket_values[1] = sf.set_limits(bracket_values[1]*2., 5.,15.)
            if key == 'VROT':
                bracket_values[1] = sf.set_limits(bracket_values[1], 10., 75.)
            # PA and INCL are set in check_angles after the initial
            set_boundary_limits(Configuration,Tirific_Template,key, values=bracket_values , tolerance = 0.01\
                            ,fixed = fixed)
        if key == 'VROT':
            fitting_settings['VROT'] = set_vrot_fitting(Configuration,stage = stage, rotation = initial_estimates['VROT'] )
        elif key == 'SBR':
            fitting_settings['SBR'] = set_sbr_fitting(Configuration,Tirific_Template,stage = stage)
        else:
            flat_slope = False
            symmetric = False
            if key in ['PA','INCL','Z0']:
                slope = Configuration['WARP_SLOPE']
                profile = np.array(sf.load_tirific(Configuration,Tirific_Template, [key,f"{key}_2"]),dtype=float)
                for i in [0,1]:
                    if slope[i] != None:
                        with warnings.catch_warnings():
                            warnings.simplefilter("error")
                            try:
                                if slope[i] < Configuration['NO_RINGS'] and \
                                    (np.mean(profile[i,:3]) > np.mean(profile[i,3:slope[i]]) < np.mean(profile[i,slope[i]:]) or \
                                    np.mean(profile[i,:3]) < np.mean(profile[i,3:slope[i]]) > np.mean(profile[i,slope[i]:])):
                                    flat_slope = True
                            except RuntimeWarning:
                                flat_slope = False


                if stage in ['initialize_os']:
                    inner = int(sf.set_limits(Configuration['NO_RINGS']*1./3., 3,Configuration['NO_RINGS']-2 ))
                else:
                    inner =  Configuration['INNER_FIX']
            elif key in ['SDIS']:
                flat_slope = True
                symmetric = True
                inner = int(sf.set_limits(Configuration['NO_RINGS']*0.33,4,Configuration['NO_RINGS']*0.5))
                slope = [int(Configuration['NO_RINGS']*0.66),int(Configuration['NO_RINGS']*0.66)]
            else:
                inner = 4
                slope = [None ,None ]

            if key in ['SDIS','INCL']:
                fact = sf.set_limits(float(initial_estimates['INCL'][0])/20.-2.5,1.,2.)
                try:
                    inner[0] = int(sf.set_limits(inner[0]*fact,4,Configuration['NO_RINGS']/2.))
                    inner[1] = int(sf.set_limits(inner[1]*fact,4,Configuration['NO_RINGS']/2.))
                except TypeError:
                    inner = int(sf.set_limits(inner*fact,4,Configuration['NO_RINGS']/2.))
            # set the moderate value
            moderate = set_generic_moderate(Configuration,key)
            #If we do a single disk then symmetric is always true
            if Configuration['NUMBER_OF_DISKS'] == 1:
                symmetric = True

            fitting_settings[key] =  set_generic_fitting(Configuration,key,\
                stage = stage, basic_variation = initial_estimates[key][1],\
                slope= slope, flat_slope = flat_slope,fixed =fixed, \
                flat_inner = inner, step_modifier = modifiers[key], \
                symmetric = symmetric,moderate = moderate)

    # Reset the fitting parameters
    for fit_key in fitting_keys:
        Tirific_Template[fit_key]= ''
    #write the new parameters
    for key in parameters_to_adjust:
        if key in fitting_settings:
            for fit_key in fitting_keys:
                if  fit_key in fitting_settings[key]:
                    format = sf.set_format(key)

                    if fit_key in ['DELSTART','DELEND','MINDELTA']:
                    #if fit_key in ['DELEND','MINDELTA']:
                        # These should never be 0.
                        if fit_key == 'MINDELTA':
                            print_log(f'''SET_FITTING_PARAMETERS: The MINDELTA for {key} should never be below {Configuration['MIN_ERROR'][key][0]*0.1}.
''',Configuration,case=['debug_add'])

                        for i,x in enumerate(fitting_settings[key][fit_key]):
                            while float(f'{fitting_settings[key][fit_key][i]:{format}}') == 0.:
                                if float(fitting_settings[key][fit_key][i]) == 0.:
                                    fitting_settings[key][fit_key][i] += 0.01
                                else:
                                    fitting_settings[key][fit_key][i] *= 2.
                            if fit_key in ['MINDELTA','DELEND'] and key not in ['SBR']:
                                if fitting_settings[key][fit_key][i] < Configuration['MIN_ERROR'][key][0]*0.05:
                                    fitting_settings[key][fit_key][i] = Configuration['MIN_ERROR'][key][0]*0.05
                            if not Configuration['ACCEPTED_TIRIFIC'] and fit_key == 'MINDELTA':
                                fitting_settings[key][fit_key][i] = fitting_settings[key][fit_key][i]*sf.set_limits(Configuration['LOOPS']-9,1,4)







                    if fit_key == 'VARY':
                        if len(Tirific_Template[fit_key]) == 0:
                            Tirific_Template[fit_key] = ', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])
                        else:
                            Tirific_Template[fit_key] = f"{Tirific_Template[fit_key]}, {', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])}"
                    else:
                        if fit_key == 'VARINDX':
                            format = '<10s'
                        else:
                            format = sf.set_format(key)
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



def set_generic_fitting(Configuration, key , stage = 'initial', \
        basic_variation = 5.,slope = None, flat_slope = False, symmetric = False,\
                        fixed = True, moderate = 3, step_modifier = None,\
                        flat_inner = 3):
    if slope is None:
        slope = [None, None]
    if step_modifier is None:
        step_modifier = [1.,1.,1.]

    print_log(f'''SET_GENERIC_FITTING: We are processing {key}.
''', Configuration, case=['debug_start'])
    if isinstance(flat_inner,int):
        flat_inner = [flat_inner,flat_inner]
    NUR = Configuration['NO_RINGS']
    input= {}
    print_log(f'''SET_GENERIC_FITTING: flat is {fixed}
''', Configuration,case=['debug_add'])
    if (stage in ['after_os','final_os','after_cc','after_ec','parameterized']) or fixed:
        print_log(f'''SET_GENERIC_FITTING: Fitting all as 1.
''', Configuration,case=['debug_add'])
        input['VARY'] =  np.array([f"{key} 1:{NUR} {key}_2 1:{NUR}"],dtype=str)
        input['PARMAX'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][0][1]],dtype=float)
        input['PARMIN'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][0][0]],dtype=float)
        input['MODERATE'] = np.array([moderate],dtype=int) #How many steps from del start to del end
        input['DELSTART'] = np.array([basic_variation*step_modifier[0]],dtype=float) # Starting step
        input['DELEND'] = np.array([0.1*basic_variation*step_modifier[1]],dtype=float) #Ending step
        input['MINDELTA'] = np.array([0.05*basic_variation*step_modifier[2]],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
    else:
        if not symmetric:
            print_log(f'''SET_GENERIC_FITTING: implementing a varying non-symmetric profile.
{'':8s} step_modifier = {step_modifier}
{'':8s} variations = {basic_variation}
{'':8s} limits = {Configuration[f'{key}_CURRENT_BOUNDARY']}
''', Configuration,case=['debug_add'])
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
            input['PARMAX'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][1][1],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][2][1],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][0][1]],dtype=float)

            input['PARMIN'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][1][0],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][2][0],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][0][0]],dtype=float)
            input['MODERATE'] =np.array([moderate,moderate,moderate],dtype=float) #How many steps from del start to del end
            input['DELSTART'] =np.array([2.,2.,0.5],dtype=float)*step_modifier[0]*basic_variation# Starting step
            input['DELEND'] = np.array([0.1,0.1,0.05],dtype=float)*step_modifier[1]*basic_variation#Ending step
            input['MINDELTA'] = np.array([0.1,0.1,0.075],dtype=float)*step_modifier[2]*basic_variation #saturation criterum when /SIZE SIZE should be 10 troughout the code
        else:
            print_log(f'''SET_GENERIC_FITTING: implementing a varying symmetric profile
{'':8s} step_modifier = {step_modifier}
{'':8s} basic_variation = {basic_variation}
{'':8s} limits = {Configuration[f'{key}_CURRENT_BOUNDARY']}
''', Configuration,case=['debug_add'])
            flat_inner = int(np.min(flat_inner))
            if flat_inner+1 >= NUR:
                input['VARY'] =  np.concatenate((np.array([f"!{key} {NUR} {key}_2 {NUR}"],dtype=str),\
                                                 np.array([f"{key} 1:{NUR-1} {key}_2 1:{NUR-1}"],dtype=str)))
            else:
                input['VARY'] =  np.concatenate((np.array([f"!{key} {NUR}:{flat_inner+1} {key}_2 {NUR}:{flat_inner+1}"],dtype=str),\
                                                 np.array([f"{key} 1:{flat_inner} {key}_2 1:{flat_inner}"],dtype=str)))
            input['PARMAX'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][1][1],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][0][1]],dtype=float)
            input['PARMIN'] = np.array([Configuration[f'{key}_CURRENT_BOUNDARY'][1][0],\
                                        Configuration[f'{key}_CURRENT_BOUNDARY'][0][0]],dtype=float)
            input['MODERATE'] =np.array([moderate,moderate],dtype=float) #How many steps from del start to del end
            input['DELSTART'] =np.array([2.,0.5],dtype=float)*step_modifier[0]*basic_variation # Starting step
            input['DELEND'] = np.array([0.1,0.05],dtype=float)*step_modifier[1]*basic_variation #Ending step
            input['MINDELTA'] = np.array([0.1,0.075],dtype=float)*step_modifier[2]*basic_variation #saturation criterum when /SIZE SIZE should be 10 troughout the code
        # then we need to set the warp slope

        forvarindex = ''
        if Configuration['NO_RINGS'] > 5:
            implement_keys = [key, f"{key}_2"]
            for i,cur_key in enumerate(implement_keys):
                if slope[i] < NUR and slope[i] != None:
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


    stage = 'initial'
    STAGE OF THE FITTING WE ARE IN

    basic_variation = 5.
        the basic variation parameters which is the base for delstart
        delstart =  basic_variation  without modifiers
        delend and mindelt = 0.1*basic_variation without modifiers

    slope = [None, None]
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
def set_generic_moderate(Configuration,key):
    ''' set the moderate value for a parameter
    This does not include VROT and SBR '''
    if not Configuration['CENTRAL_CONVERGENCE'] and key in ['VSYS','XPOS','YPOS']:
        moderate = 3
    else:
        moderate = 3
    return moderate


#Function
def set_model_parameters(Configuration, Tirific_Template,Model_Values, \
        stage = 'initial'):
    parameters_to_set = ['RADI','VROT_profile','Z0','SBR_profile','INCL','PA','XPOS','YPOS','VSYS','SDIS']


    check_parameters = []
    if 'VSYS' in Model_Values:
        vsys = Model_Values['VSYS']
    else:
        vsys=100.
    scramble = np.zeros(len(parameters_to_set))

    for key in parameters_to_set:
        if key in Model_Values:
            if key in ['VROT_profile','SBR_profile']:
                key_to_set = key.split('_')[0]
            else:
                key_to_set = key
            # if 2 long we have a value and error
            format = sf.set_format(key_to_set)
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
                rad = sf.set_rings(Configuration)
                Tirific_Template['RADI']= f"{' '.join([f'{x:.2f}' for x in rad])}"
                Tirific_Template['NUR']=str(len(rad))
                check_parameters.append('RADI')
            elif key == 'Z0':
                check_parameters.append('Z0')
                if Model_Values['INCL'][0] > 80:
                    Tirific_Template['Z0'] = f"{np.max([sf.convertskyangle(Configuration,0.2,physical= True),Configuration['BEAM'][0]/4.]):.3f}"

                else:
                    Tirific_Template['Z0'] = f"{sf.convertskyangle(Configuration,0.2,physical= True):.3f}"


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
    no_declining_vrot(Configuration, Tirific_Template)
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


    stage = 'initial'
    stage of the fitting process

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
#function to check that all parameters in template have the proper length.
def set_new_size(Configuration,Tirific_Template, fit_type = 'Undefined',
                    current_run='Not Initialized', Variables = None):

    if Variables is None:
        Variables =['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS',\
            'VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2',\
            'SDIS_2', 'AZ1P', 'AZ1W' ,'AZ1P_2','AZ1W_2', 'RADI']

    print_log(f'''SET_NEW_SIZE: Starting the adaptation of the size in the template
''',Configuration,case=['debug_start'])
    for key in Variables:
        if not key in Tirific_Template:
            Variables.remove(key)
    parameters = sf.load_tirific(Configuration,Tirific_Template,Variables)
    # do a check on whether we need to modify all
    radii = sf.set_rings(Configuration)
    if not Configuration['NEW_RING_SIZE']:
        print_log(f'''SET_NEW_SIZE: The rings size is stable and we are not interpolating
''',Configuration,case=['debug_add'])
        interpolate = False
    else:
        print_log(f'''SET_NEW_SIZE: The rings size is updated and we are interpolating the old values.
''',Configuration,case= ['verbose'])
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
        print_log(f'''SET_NEW_SIZE: We are processing {key}
{'':8s}SET_NEW_SIZE: We have a parameter of length {len(parameters[i])}.
{'':8s}SET_NEW_SIZE: Our current number of rings in the model is {Configuration['NO_RINGS']}.
{'':8s}SET_NEW_SIZE: {parameters[i]}
''',Configuration,case=['debug_add'])
        if key == 'RADI':
            Tirific_Template[key] = f" {' '.join([f'{x:.2f}' for x in radii])}"
        else:
            if interpolate:
                radii_int=np.array(radii,dtype=float)
                old_radii_int = np.array(old_radii,dtype=float)
                parameters_int = copy.deepcopy(parameters[i])
                print_log(f'''SET_NEW_SIZE: We are interpolating par = {parameters_int} old radii={old_radii_int} new radii={radii_int}
    ''',Configuration,case=['debug_add'])
                if len(parameters_int) > len(old_radii_int):
                    print_log(f'''SET_NEW_SIZE: The parameters have more values than the radii. Cutting the end.
    ''',Configuration,case=['debug_add'])
                    parameters_int = parameters_int[:len(old_radii)-1]
                elif len(parameters_int) < len(old_radii_int):
                    print_log(f'''SET_NEW_SIZE: The parameters have less values than the radii. Adding the last value until match.
    ''',Configuration,case=['debug_add'])
                    while len(parameters_int) < len(old_radii_int):
                        parameters_int.append(parameters_int[-1])

                parameters[i] = list(np.interp(radii_int,old_radii_int,np.array(parameters_int,dtype=float)))
                del radii_int
                del old_radii_int
                del parameters_int
            format = sf.set_format(key)

          
            if len(parameters[i]) > Configuration['NO_RINGS']-1:
                # if we are cutting a ring it is likely the outer ring have done weird stuff so we flatten the curve
                Tirific_Template[key] =f" {' '.join([f'{x:{format}}' for x in parameters[i][:int(Configuration['NO_RINGS'])]])}"

            elif len(parameters[i]) <= Configuration['NO_RINGS']-1:
                # We simply add the last value
                counter = len(parameters[i])
                while counter !=  Configuration['NO_RINGS']:
                    Tirific_Template[key] = Tirific_Template[key] + f" {parameters[i][-1]:{format}}"
                    counter += 1
                   
                       
        print_log(f'''SET_NEW_SIZE: We wrote the following line {Tirific_Template[key]}
''',Configuration,case=['debug_add'])

    #Replace the old ring numbers in VARY and VARINDX
    #old_rings = sf.calc_rings(Configuration,size_in_beams = int(round(np.max([float(Configuration['OLD_SIZE'][-1][0]),float(Configuration['OLD_SIZE'][-1][1])]))), ring_size  = 0.)
    old_rings = int(copy.deepcopy(Tirific_Template['NUR']))
    current_rings = sf.calc_rings(Configuration)
    #This could lead to replacing a value smaller than the other side
    Tirific_Template['VARY'] = Tirific_Template['VARY'].replace(f"{old_rings}",f"{current_rings}")
    Tirific_Template['VARINDX'] = Tirific_Template['VARINDX'].replace(f"{old_rings-1}",f"{current_rings-1}")
    Tirific_Template['VARINDX'] = Tirific_Template['VARINDX'].replace(f"{old_rings}",f"{current_rings}")
    Tirific_Template['NUR'] = f"{current_rings}"
    # if we cut we want to flatten things
    #if current_rings < old_rings:
    #    flatten_the_curve(Configuration,Tirific_Template)
    # Check whether the galaxy has become to small for variations
    if np.sum(Configuration['SIZE_IN_BEAMS']) > Configuration['MINIMUM_WARP_SIZE']:
        Configuration['FIXED_PARAMETERS'][0] = Configuration['FIXED_PARAMETERS'][1]
    else:
        for parameter in ['INCL','Z0','SDIS','PA']:
            if parameter not in Configuration['FIXED_PARAMETERS'][0]:
                Configuration['FIXED_PARAMETERS'][0].append(parameter)
        flatten_the_curve(Configuration,Tirific_Template)
    print_log(f'''The following parameters will be fixed'
{'':8s} {Configuration['FIXED_PARAMETERS'][0]}
''',Configuration,case=['debug_add'])

    #update the limit_modifier
    sf.set_limit_modifier(Configuration,Tirific_Template)
    # Maybe increase the amount of loops in smaller galaxies
    #set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,fit_type=fit_type)
    # if we change the radii we need to restart tirific
    sf.finish_current_run(Configuration,current_run)
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


    fit_type = 'Undefined'
    type of fitting


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


def set_overall_parameters(Configuration, Fits_Files,Tirific_Template,\
                fit_type='Undefined', flux = None):
    print_log(f'''SET_OVERALL_PARAMETERS: startin to set the overall parameters.
''',Configuration,case=['debug_add'])
    if Configuration['OPTIMIZED']:
        Tirific_Template['INSET'] = f"{Fits_Files['OPTIMIZED_CUBE']}"
    else:
        Tirific_Template['INSET'] = f"{Fits_Files['TIR_RUN_CUBE']}"

    if Configuration['NO_RINGS'] < 20:
        Tirific_Template['INIMODE'] = '1'
    elif Configuration['NO_RINGS'] < 30:
        Tirific_Template['INIMODE'] = '2'
    else:
        Tirific_Template['INIMODE'] = '3'

    if Configuration['NO_RINGS'] > 8.:

        Tirific_Template['INTY'] = 0
        #Tirific_Template['INDINTY'] = 0
    else:
        #This should not be used with 8 <  rings.
        Tirific_Template['INTY'] = 2
        #preferably we'd use the akima spline but there is an issue with that where the final ring does not get modified
        #Tirific_Template['INDINTY'] = 2
    Tirific_Template['NUR'] = f"{Configuration['NO_RINGS']}"
    #current_cwd = os.getcwd()
    short_log = Configuration['LOG_DIRECTORY'].replace(Configuration['FITTING_DIR'],'')
    Tirific_Template['RESTARTNAME'] = f"{short_log}restart_{fit_type}.txt"
    #this could be fancier
    Tirific_Template['NCORES'] = Configuration['PER_GALAXY_NCPU']

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

    if Configuration['CHANNEL_DEPENDENCY'].lower() == 'hanning':
        instrumental_vres = (Configuration['CHANNEL_WIDTH']*2)/(2.*np.sqrt(2.*np.log(2)))
    elif Configuration['CHANNEL_DEPENDENCY'].lower() == 'sinusoidal':
        instrumental_vres = (Configuration['CHANNEL_WIDTH']*1.2)/(2.*np.sqrt(2.*np.log(2)))
    elif Configuration['CHANNEL_DEPENDENCY'].lower() == 'independent':
        instrumental_vres = 0.
    else:
        raise BadConfigurationError('Something went wrong in the Configuration setup')

    Tirific_Template['CONDISP'] = f"{instrumental_vres:.2f}"
    if flux:
        Tirific_Template['CFLUX'] = f"{sf.set_limits(flux/1.5e7,1e-6,1e-3):.2e}"
        Tirific_Template['CFLUX_2'] = f"{sf.set_limits(flux/1.5e7,1e-6,1e-3):.2e}"
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

def set_sbr_fitting(Configuration,Tirific_Template, stage = 'no_stage'):
    print_log(f'''SET_SBR_FITTING: We are setting the SBR limits.
{'':8s} No_Rings = {Configuration['NO_RINGS']}
''',Configuration,case=['debug_start'])
    sbr_input = {}
    inner_ring = 2
    sbr_profile= sf.load_tirific(Configuration,Tirific_Template,Variables= ['SBR','SBR_2'],array=True)
    radii= sf.load_tirific(Configuration,Tirific_Template, Variables=['RADI'],array=True)
    last_ring_to_fit = []
    for i in [0,1]:
        last_ring_to_fit.append(np.where(radii <= \
            Configuration['SIZE_IN_BEAMS'][i]*Configuration['BEAM'][0])[0][-1])
        
  

    if stage in ['initial','run_cc','initialize_ec','run_ec','initialize_os','run_os']:
        sbr_ring_limits = sf.sbr_limits(Configuration,Tirific_Template)
        if stage in ['run_ec','run_os']:
            sbr_ring_limits[-4:]=[x/5 for x in sbr_ring_limits[-4:]]
        print_log(f'''SET_SBR_FITTING: Using these SBR limits.
{'':8s} limits = {sbr_ring_limits}
{'':8s} no of limits = {len(sbr_ring_limits)}
''',Configuration,case=['debug_add'])
        if stage in ['initial','run_cc']:
            max_size = 4
        elif stage in ['initialize_ec','run_ec','initialize_os','run_os']:
            max_size = 2.
        #Make sure the SBR profile is not dropping below the minimum we are setting

        fact= [2.,2.]
        format = sf.set_format('SBR')
        pmax = np.full((2,len(sbr_profile[0,:])),1.)
        pmin = np.full((2,len(sbr_profile[0,:])),0.)
        for i in [0,1]:
            if stage in ['initialize_os']:
                fact[i] = 0.75
            elif Configuration['SIZE_IN_BEAMS'][i] < max_size:
                fact[i]=2.5

            for x in range(len(radii)-1,inner_ring-1,-1):
                if radii[x] <= Configuration['SIZE_IN_BEAMS'][i]*Configuration['BEAM'][0]:
                    sbr_profile[i,x] = sf.set_limits(sbr_profile[i,x],sbr_ring_limits[x]/fact[i]*2.,1.)
                else:
                    sbr_profile[i,x] = 0.
            sbr_profile[i,:inner_ring] = [sf.set_limits(x,np.min(sbr_ring_limits),1.) for x in sbr_profile[i,:inner_ring]]
            if i == 0:
                ext=''
            else:
                ext='_2'
            Tirific_Template[f"SBR{ext}"]= f"{' '.join([f'{x:{format}}' for x in sbr_profile[i]])}"

        sbr_smoothed_profile = smooth_profile(Configuration,Tirific_Template,'SBR',
                                min_error= [sbr_ring_limits,sbr_ring_limits],no_apply = True,
                                fix_sbr_call = True,profile_in = sbr_profile )
        sbr_av_smoothed = [(x+y)/2. for x,y in zip(sbr_smoothed_profile[0],sbr_smoothed_profile[1])]
        print_log(f'''SET_SBR_FITTING:
{'':8s} This is the mean SNR {Configuration['SNR']}
beamarea = {Configuration['BEAM_AREA']}, channelwidth = {Configuration['CHANNEL_WIDTH']}, noise = {Configuration['NOISE']}
''',Configuration,case=['debug_add'])

        #if x < 4:
        limits_for_max = True
        if Configuration['SNR'] > 10:
            limits_for_max = False
            print_log(f'''SET_SBR_FITTING:
{'':8s} The SNR is so high that we will not use the ring limits for parmax
''',Configuration,case=['debug_add'])

        for i in [0,1]:
            mean_signal = Configuration['SNR']*Configuration['NOISE']/\
                    (Configuration['BEAM_AREA']*Configuration['SIZE_IN_BEAMS'][i]/2.)\
                    *Configuration['CHANNEL_WIDTH'] #in Jy/arcsec *km/s
            print_log(f'''SET_SBR_FITTING:
{'':8s} mean_signal = {mean_signal}
''',Configuration,case=['debug_add'])
            for x in range(len(radii)-1,inner_ring-1,-1):

                #    min_max = sf.set_limits(sbr_ring_limits[x]*30.,1e-3,0.9 )
                #else:
                #    min_max = sbr_ring_limits[x]*30.

                if  sbr_profile[i,x] > 0.:
                    if limits_for_max:
                        pmax[i,x] = sf.set_limits(sbr_av_smoothed[x-1]*20.,sbr_ring_limits[x]*30.,1.)
                    else:
                        pmax[i,x] = sf.set_limits(sbr_av_smoothed[x-1]*20.,\
                            mean_signal*100./(radii[x])**0.5,1.)
                else:
                    if limits_for_max:
                        pmax[i,x] = sf.set_limits(sbr_av_smoothed[x-1]*5.\
                        ,sbr_ring_limits[x]*10.,sbr_ring_limits[x]*30.)
                    else:
                        pmax[i,x] = sf.set_limits(sbr_av_smoothed[x-1]*5.\
                        , mean_signal*50./(radii[x])**0.5, mean_signal*200./(radii[x])**0.5)
                pmin[i,x] = sbr_ring_limits[x]/fact[i]


     
       
        if np.mean(Configuration['SIZE_IN_BEAMS']) < max_size or Configuration['NUMBER_OF_DISKS'] == 1:
            extend = np.max(last_ring_to_fit)
               
            sbr_input['VARY'] =  np.array([f"SBR {x+1} SBR_2 {x+1}" for x \
                in range(extend,inner_ring-1,-1)],dtype=str)


            sbr_input['PARMAX'] = np.array([sf.set_limits(sbr_av_smoothed[x-1]\
                *10.,np.mean(pmax[:,x]),1.) for x in range(len(radii)-1,\
                    inner_ring-1,-1)],dtype=float)
            #if stage in ['initial','run_cc']:
            #    sbr_input['PARMIN'] = np.array([sbr_ring_limits[x]/2. if x <= (3./4.)*len(radii) else 0 for x in range(len(radii)-1,inner_ring-1,-1)])
            #elif stage in ['initialize_ec','run_ec']:

            sbr_input['PARMIN'] = np.array([sbr_ring_limits[x]/np.max(fact) for x in range(extend,inner_ring-1,-1)],dtype=float)
            sbr_input['MODERATE'] = np.array([5 for x in range(extend,inner_ring-1,-1)],dtype=float) #How many steps from del start to del end
            sbr_input['DELSTART'] = np.array([1e-4 for x in range(extend,inner_ring-1,-1)],dtype=float) # Starting step
            sbr_input['DELEND'] = np.array([2.5e-6 for x in range(extend,inner_ring-1,-1)],dtype=float) #Ending step
            sbr_input['MINDELTA'] = np.array([5e-6 for x in range(extend,inner_ring-1,-1)],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
        else:
            input_variables = ['VARY','PARMAX','PARMIN','MODERATE','DELSTART','DELEND','MINDELTA']
            for key in input_variables:
                sbr_input[key] = []

            ext = ''
            for x in  range(np.max(last_ring_to_fit),inner_ring-1,-1):
                for i in [0,1]:
                    if x <= last_ring_to_fit[i]:
                        if i == 0:
                            ext = ''
                        else:
                            ext='_2'
                        sbr_input['VARY'].append(f'SBR{ext} {x+1}') 
                        sbr_input['PARMAX'].append(pmax[i,x]) 
                        sbr_input['PARMIN'].append(pmin[i,x]) 
                        sbr_input['MODERATE'].append(5.)
                        sbr_input['DELSTART'].append(sbr_smoothed_profile[i,x]/10.) 
                        sbr_input['DELEND'].append(sbr_ring_limits[x]/20.) 
                        sbr_input['MINDELTA'].append(sbr_ring_limits[x]/20.)
            for key in input_variables:
                sbr_input[key] =np.array(sbr_input[key]) 
            '''           
            sbr_input['VARY'] =  np.array([[f"SBR {x+1}",f"SBR_2 {x+1}"] for x in range(extend,inner_ring-1,-1)],dtype=str).reshape((len(radii)-inner_ring)*2)
            #pmax = np.array([np.full(2,sf.set_limits(sbr_av_smoothed[x-1]*20.,sbr_ring_limits[x]*30.,1.)) for x in range(len(radii)-1,inner_ring-1,-1)], dtype=float)
            sbr_input['PARMAX'] = np.array([[pmax[0,x],pmax[1,x]] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            #sbr_input['PARMAX'] = np.array([[sf.set_limits(sbr_smoothed_profile[0,x-1]*20.,sbr_ring_limits[x]*30.,1.),sf.set_limits(sbr_smoothed_profile[1,x-1]*20.,sbr_ring_limits[x]*30.,1.)] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
        #if stage in ['initial','run_cc']:
            #    sbr_input['PARMIN'] = np.array([[sbr_ring_limits[x]/2.,sbr_ring_limits[x]/2.] if x <= (3./4.)*len(radii) else [0.,0.] for x in range(len(radii)-1,inner_ring-1,-1)]).reshape((len(radii)-inner_ring)*2)
            #elif stage in ['initialize_ec','run_ec']:
            sbr_input['PARMIN']  = np.array([[pmin[0,x],pmin[1,x]] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            #sbr_input['PARMIN'] = np.array([[sbr_ring_limits[x]/fact[0],sbr_ring_limits[x]/fact[1]] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            sbr_input['MODERATE'] = np.array([[5,5] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2) #How many steps from del start to del end
            sbr_input['DELSTART'] = np.array([[1e-4,1e-4] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2) # Starting step
            sbr_input['DELEND'] = np.array([[sbr_ring_limits[x]/20.,sbr_ring_limits[x]/20.] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            sbr_input['MINDELTA'] = np.array([[sbr_ring_limits[x]/20.,sbr_ring_limits[x]/20.] for x in range(len(radii)-1,inner_ring-1,-1)],dtype=float).reshape((len(radii)-inner_ring)*2)
            '''
        sbr_input['VARY'] = np.concatenate((sbr_input['VARY'],\
            [f"SBR {' '.join([str(int(x)) for x in range(1,inner_ring+1)])} SBR_2 {' '.join([str(int(x)) for x in range(1,inner_ring+1)])}"]\
            ),axis=0)
        if limits_for_max:
            sbr_input['PARMAX'] = np.concatenate((sbr_input['PARMAX'],\
                [sf.set_limits(np.mean([sbr_smoothed_profile[0,2:4],\
                sbr_smoothed_profile[1,2:4]])*4.,sbr_ring_limits[2],\
                np.max(sbr_profile))]))
        else:
            sbr_input['PARMAX'] = np.concatenate((sbr_input['PARMAX'],\
                [sf.set_limits(np.mean([sbr_smoothed_profile[0,2:4],\
                sbr_smoothed_profile[1,2:4]])*4.,mean_signal*100,1.)]))

        if Configuration['CENTRAL_CONVERGENCE']:
            sbr_input['PARMIN'] = np.concatenate((sbr_input['PARMIN'],[np.min(sbr_ring_limits)]))
        else:
            sbr_input['PARMIN'] = np.concatenate((sbr_input['PARMIN'],[sbr_ring_limits[inner_ring]*1.5]))
        sbr_input['MODERATE'] = np.concatenate((sbr_input['MODERATE'],[5]))
        sbr_input['DELSTART'] = np.concatenate((sbr_input['DELSTART'],[1e-5]))
        sbr_input['DELEND'] = np.concatenate((sbr_input['DELEND'],[1e-6]))
        sbr_input['MINDELTA'] = np.concatenate((sbr_input['MINDELTA'],[2e-6]))
    elif stage in ['after_cc','after_ec','after_os'] :
        #Used in Fit_Smoothed_Check
        # Initialize the parameters
        par_to_fill = ['PARMAX','PARMIN','MODERATE','DELSTART',\
                'DELEND','MINDELTA']
        value =       [np.max(sbr_profile[0,2:])*3., 0., 5. , 1e-5, 1e-6,2e-6]
        # initialize
        for var in par_to_fill:
            sbr_input[var] = []

        if Configuration['NUMBER_OF_DISKS'] == 2:
            sbr_input['VARY'] = [f"SBR {inner_ring+1}:{last_ring_to_fit[0]+1}, SBR_2 {inner_ring+1}:{last_ring_to_fit[1]+1}"]
        else:
            sbr_input['VARY'] = [f"SBR {inner_ring+1}:{last_ring_to_fit[0]+1} SBR_2 {inner_ring+1}:{last_ring_to_fit[1]+1}"]

        for i in range(Configuration['NUMBER_OF_DISKS']):
            for i,var in enumerate(par_to_fill):
                sbr_input[var].append(value[i])
        for var in par_to_fill:
            sbr_input[var] = np.array(sbr_input[var])

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

def set_vrot_fitting(Configuration, stage = 'initial', rotation = None):
    if rotation is None:
        rotation = [100,5.]
    NUR = Configuration['NO_RINGS']
    vrot_input = {}

    print_log(f'''SET_VROT_FITTING: We are setting the VROT limits.
{'':8s} No_Rings = {Configuration['NO_RINGS']}
{'':8s} Limits = {Configuration['VROT_CURRENT_BOUNDARY'][0]}
''',Configuration,case=['debug_start'])
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
    vrot_input['PARMAX'] = np.array([Configuration['VROT_CURRENT_BOUNDARY'][0][1]],dtype=float)
    vrot_input['PARMIN'] = np.array([Configuration['VROT_CURRENT_BOUNDARY'][0][0]],dtype=float)
    vrot_input['MODERATE'] = np.array([5],dtype=float) #How many steps from del start to del end
    vrot_input['DELSTART'] = np.array([2.*Configuration['CHANNEL_WIDTH']*modifier[0]],dtype=float) # Starting step
    #These were lower in the original fat
    vrot_input['DELEND'] = np.array([0.1*Configuration['CHANNEL_WIDTH']*modifier[1]],dtype=float) #Ending step
    vrot_input['MINDELTA'] = np.array([0.1*Configuration['CHANNEL_WIDTH']*modifier[2]],dtype=float) #saturation criterum when /SIZE SIZE should be 10 troughout the code
    #if there is not values in the center we connect the inner ring to the next ring
    forvarindex = ''
    if Configuration['NO_RINGS'] > 5:
        if Configuration['EXCLUDE_CENTRAL']: # Let's remove this for now as from ESO 92_G021 and serveral masssive tests it seems a bad idea or rotation[0] > 150.:
            forvarindex = 'VROT 2 VROT_2 2 '
        if Configuration['OUTER_SLOPE_START'] == NUR:
            if Configuration['NO_RINGS'] > 5:
                inner_slope =  int(round(sf.set_limits(Configuration['RC_UNRELIABLE'],round(NUR/2.),NUR-1)))
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

def smooth_factor_profile(Configuration, change_angle,radius,inner_fix= None,\
                          Tirific_Template= None):
    print_log(f'''SMOOTH_PROFILE: Starting to smooth the phi and theta factor profile.
Original Phi = {change_angle['PHI']['FACTOR']}
Original Theta = {change_angle['THETA']['FACTOR']}             
''',Configuration,case= ['debug_start'])
    if inner_fix  is None:
        inner_fix = Configuration['INNER_FIX'][0]
    theta_factor  = change_angle['THETA']['FACTOR'][inner_fix:]   
    phi_factor =  change_angle['PHI']['FACTOR'][inner_fix:]
    
    #Smoothing doesn't work
    theta_factor = np.concatenate((np.full(inner_fix,theta_factor[0]),theta_factor))
    phi_factor = np.concatenate((np.full(inner_fix,phi_factor[0]),phi_factor))
    if np.mean( theta_factor) !=  theta_factor[0]:
        theta_factor = fit_polynomial(Configuration,radius,\
                    theta_factor,np.full(len(theta_factor),0.1),'THETA_FACTOR', Tirific_Template,\
                    inner_fix = inner_fix,allowed_order= [1,6],
                    boundary_limits= [-1,1])
    if np.mean(phi_factor) !=  phi_factor[0]:  
        phi_factor = fit_polynomial(Configuration,radius,\
                phi_factor,np.full(len(phi_factor),0.1),'PHI_FACTOR', Tirific_Template,\
                inner_fix = inner_fix,allowed_order= [1,6],
                boundary_limits= [-1,1])
    
    #theta_factor = np.concatenate((np.zeros(inner_fix),theta_factor))
    #reapply the signs
    #with warnings.catch_warnings():
    #    warnings.filterwarnings("ignore",message="invalid value encountered in divide"\
     #                       ,category=RuntimeWarning)
    #    theta_factor = theta_factor*change_angle['THETA']['FACTOR']\
    #        /abs(change_angle['THETA']['FACTOR'])
    #theta_factor[np.isnan(theta_factor)] = 0.

    #if len (phi_factor) > 3:
    #    phi_factor= savgol_filter(phi_factor, 3, 1)
    #phi_factor = np.concatenate((np.zeros(inner_fix),phi_factor))
    #with warnings.catch_warnings():
    #    warnings.filterwarnings("ignore",message="invalid value encountered in divide"\
    #                        ,category=RuntimeWarning)
    #    phi_factor = phi_factor*change_angle['PHI']['FACTOR']\
    #        /abs(change_angle['PHI']['FACTOR'])
    #phi_factor[np.isnan(phi_factor)] = 0.

    #Then the factor together need to 1.
    theta_factor,phi_factor=ensure_total_dist(Configuration,\
        theta_factor,phi_factor,inner_fix=inner_fix)    
    return theta_factor,phi_factor

def ensure_total_dist(Configuration,theta_factor,phi_factor,inner_fix=3):
    print_log(f'''ensure_total_dist: Starting to level the factors to add to 1 in square
Original Phi = {phi_factor}
Original Theta = {theta_factor}             
''',Configuration,case= ['debug_start'])

    for i in range(len(theta_factor)):
        #if i < inner_fix:
        #    if phi_factor[i]+theta_factor[i] != 0:
        #        raise RuntimeError('The fixed part of the factors is not 0')
        #else:
                new_comb =  phi_factor[i]**2+theta_factor[i]**2
                diff = np.sqrt(abs(1.-new_comb)/2.)
                if phi_factor[i] == 0.:
                    tmp_phi = abs(1-abs(theta_factor[i]))
                else:
                    tmp_phi = phi_factor[i]
                if theta_factor[i] == 0.:
                    tmp_theta = abs(1-abs(phi_factor[i]))
                else:
                    tmp_theta = theta_factor[i]
                count = 0
                while not isclose(new_comb,1.,abs_tol=1e-3):
                    if new_comb >   1.:
                        diff = -1.*abs(diff/2.)
                    else:
                        diff = abs(diff/2.)
                    tmp_phi =  tmp_phi+diff*phi_factor[i]/abs(phi_factor[i])
                    tmp_theta =  tmp_theta+diff*theta_factor[i]/abs(theta_factor[i])
                    new_comb =  tmp_phi**2+tmp_theta**2
                    if abs(diff) < 1e-7:
                        tmp_phi = phi_factor[i]
                        tmp_theta = theta_factor[i] 
                        diff = float(np.random.default_rng().random())
                        new_comb =  tmp_phi**2+tmp_theta**2 
                    count += 1
                    if count > 1000:
                        break   
                if np.isnan(tmp_phi):
                    tmp_phi=0.
                if np.isnan(tmp_theta):
                    tmp_theta=0.
                phi_factor[i]=tmp_phi
                theta_factor[i]=tmp_theta
    return theta_factor,phi_factor

def smooth_profile(Configuration,Tirific_Template,key,min_error = 0.  \
                    ,profile_in = None, no_apply =False,no_fix = False,\
                    fix_sbr_call=False,singular=False,inner_fix=None):
    print_log(f'''SMOOTH_PROFILE: Starting to smooth the {key} profile.
''',Configuration,case= ['debug_start'])
    if key == 'SBR' and not fix_sbr_call:
        error_message = f'''SMOOTH_PROFILE: Do not use smooth_profile for the SBR, SBR is regularised in fix_sbr'''
        print_log(error_message,Configuration,\
            case= ['main','screen'])
        raise FunctionCallError(error_message)
    if fix_sbr_call:
        no_fix=True
    if key in 'ARBITRARY' and (profile_in is None or no_apply is not True):
        error_message = f'''SMOOTH_PROFILE: We cannot read or write an ARBITRARY profile from the template.'''
        print_log(error_message,Configuration,\
            case= ['main','screen'])
        raise FunctionCallError(error_message)


    if profile_in is None:
        profile = np.array(sf.load_tirific(Configuration,Tirific_Template,[key,f"{key}_2"]),dtype = float)
    else:
        profile= copy.deepcopy(profile_in)

    original_profile = copy.deepcopy(profile)
    min_error = np.array(min_error,dtype=float)
    if min_error.size == 1:
        min_error = np.full(profile.size,min_error)
    if inner_fix is None:
        if key in ['INCL','Z0', 'PA','ARBITRARY','CHANGE_ANGLE','THETA_FACTOR','PHI_FACTOR']:
            inner_fixed = Configuration['INNER_FIX']
        elif key in ['SDIS']:
            inner_fixed = [4,4]
        else:
            inner_fixed = [0,0]
    else:
        inner_fixed = inner_fix
    #he sbr profile is already fixed before getting to the smoothing
    if not no_fix:
        profile =fix_profile(Configuration, key, profile,\
                    Tirific_Template=Tirific_Template,\
                    inner_fix = inner_fixed,singular=singular)
       

    if key == 'VROT':
        #if profile[0,1] > profile[0,2] or np.mean(profile[1:3]) > 120.:
        if np.mean([profile[0][1:3],profile[1][1:3]]) > 120.:
            shortened =True
            profile = np.delete(profile, 0, axis = 1)
        else:
            shortened = False
    print_log(f'''SMOOTH_PROFILE: retrieved profile.
{'':8s}{profile}
''',Configuration,case= ['debug_add'])
    # savgol filters do not work for small array
    
    for i in [0,1]:
        if singular:
            if i ==1:
              continue
            else:
                singular_profile= copy.deepcopy(profile)  
        else:
            singular_profile= copy.deepcopy(profile[i]) 
        if key in ['THETA_FACTOR','PHI_FACTOR']:
            sign = [x/abs(x) if x !=0. else 1. for x in singular_profile]
            singular_profile = abs(singular_profile)
        if len(singular_profile) <= 8:
            #In this case we are using a cubic spline so no smoothing required
            pass
        elif len(singular_profile) < 15:
            singular_profile= savgol_filter(singular_profile, 3, 1)
        elif len(singular_profile) < 20:
            singular_profile= savgol_filter(singular_profile, 5, 2)
        elif len(singular_profile) < 25:
            singular_profile = savgol_filter(singular_profile, 7, 3)
        else:
            singular_profile = savgol_filter(singular_profile, 9, 4)
        if fix_sbr_call:
            below_zero = np.where(singular_profile[:int(len(profile)/2.)] < 0.)[0]
            if below_zero.size > 0.:
                singular_profile[below_zero] = min_error[i,below_zero]
            below_zero = np.where(singular_profile[int(len(profile)/2.):] < 0.)[0]+1
            if below_zero.size > 0.:
                singular_profile[below_zero] = singular_profile[below_zero-1]/2.
        if key in ['THETA_FACTOR','PHI_FACTOR']:
            singular_profile =  singular_profile*sign
        if singular:
            profile=singular_profile
        else:
            profile[i]=singular_profile

    # Fix the settings
    format = sf.set_format(key)
    if key == 'VROT':
        if shortened:
            tmp = [[],[]]
            tmp[0] = np.hstack([[0.], profile[0]])
            tmp[1] = np.hstack([[0.], profile[1]])
            profile =  np.array(tmp,dtype=float)
        else:
            profile[:,0] = 0.

    if not no_fix:
        profile =fix_profile(Configuration, key, profile,\
            Tirific_Template=Tirific_Template,\
            inner_fix=inner_fixed,singular=singular)


    print_log(f'''SMOOTH_PROFILE: profile after smoothing.
{'':8s}{profile}
''',Configuration,case= ['debug_add'])
    if not no_apply:
        weights = get_ring_weights(Configuration,Tirific_Template)
        errors = get_error(Configuration,original_profile,profile,key,weights = weights, min_error=min_error)
        if key not in ['VROT']:
            # Check whether it should be flat
            profile,errors =modify_flat(Configuration,profile,original_profile,errors,key,inner_fix=inner_fixed)

        Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in profile[0,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template[f"{key}_2"]= f"{' '.join([f'{x:{format}}' for x in profile[1,:int(Configuration['NO_RINGS'])]])}"
        Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x:{format}}' for x in errors[0,:int(Configuration['NO_RINGS'])]])}")
        Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x:{format}}' for x in errors[1,:int(Configuration['NO_RINGS'])]])}")
        #if key in ['INCL'] and np.mean( profile[:,int(Configuration['NO_RINGS']/2.):int(Configuration['NO_RINGS'])]) < 40.:
        #    fix_vrot_for_incl_change(Configuration,Tirific_Template,original_profile,profile)

        print_log(f'''SMOOTH_PROFILE: This has gone to the template
{'':8s}{key} = {Tirific_Template[key]}
{'':8s}{key}_2 ={Tirific_Template[f"{key}_2"]}
{'':8s}# {key}_ERR ={Tirific_Template[f"# {key}_ERR"]}
{'':8s}# {key}_2_ERR ={Tirific_Template[f"# {key}_2_ERR"]}
''',Configuration,case= ['debug_add'])

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

def update_disk_angles(Configuration,Tirific_Template):
    extension = ['','_2']
    for ext in extension:
        PA = np.array(sf.load_tirific(Configuration,Tirific_Template,[f'PA{ext}']),dtype=float)
        inc = np.array(sf.load_tirific(Configuration,Tirific_Template,[f'INCL{ext}']),dtype=float)
        print_log(f'''UPDATE_DISK_ANGLES: obtained  this from the template
{'':8s} inc{ext} = {inc}
{'':8s} PA{ext} = {PA}
''', Configuration,case= ['debug_start'])
        angle_adjust=np.array(np.tan((PA[0]-PA)*np.cos(inc*np.pi/180.)*np.pi/180.)*180./np.pi,dtype = float)
        if ext == '_2':
            angle_adjust[:] +=180.
        print_log(f'''UPDATE_DISK_ANGLES: adusting AZ1P{ext} with these angles
{'':8s}{angle_adjust}
''', Configuration,case= ['debug_add'])
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


 OUTPUTS:
    Updated template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def write_center(Configuration,Tirific_Template,center ):
    extension = ['','_2']
    variables = ['XPOS','YPOS','VSYS']
    for i,var in enumerate(variables):
        values= [center[i]]*int(Tirific_Template['NUR'])
        format= sf.set_format(var)
        for ext in extension:
            Tirific_Template[f"{var}{ext}"] = f"{' '.join([f'{x:{format}}' for x in values])}"

write_center.__doc__ =f'''
 NAME:
    write_center

 PURPOSE:
    Update the center in the template

 CATEGORY:
    modify_template

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template
    center = [xpos,ypos,vsys]

 OPTIONAL INPUTS:


 OUTPUTS:
    Updated template

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def write_new_to_template(Configuration, filename,Tirific_Template,\
        Variables = None):
    if Variables is None:
        Variables = ['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS'\
                    ,'VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','XPOS_2','YPOS_2'\
                    ,'VSYS_2','SDIS_2']
    '''Write a def file into the template'''
    with open(filename, 'r') as tmp:
        # Separate the keyword names
        for line in tmp.readlines():
            key = str(line.split('=')[0].strip().upper())
            if key in Variables:
                Tirific_Template[key] = str(line.split('=')[1].strip())
    # If we have written the INCL we need to update the limit modifier
    if 'INCL' in Variables or 'INCL_2' in Variables:
        Inclination = np.array([(float(x)+float(y))/2. for x,y in zip(Tirific_Template['INCL'].split(),Tirific_Template['INCL_2'].split())],dtype=float)
        sf.set_limit_modifier(Configuration,Tirific_Template)
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
