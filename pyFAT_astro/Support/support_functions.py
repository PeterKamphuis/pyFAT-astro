# -*- coding: future_fstrings -*-

from pyFAT_astro.Support.fat_errors import SupportRunError,SmallSourceError,\
                                              FileNotFoundError,TirificKillError,\
                                              InputError,ProgramError,DefFileError,\
                                              BadHeaderError,FittingError,TirificOutputError,\
                                              SofiaMissingError
from pyFAT_astro.Support.log_functions import print_log,get_usage_statistics,write_config
from pyFAT_astro import Templates as templates
from collections import OrderedDict #used in Proper_Dictionary
from numpy.linalg import LinAlgError
from scipy.optimize import curve_fit, OptimizeWarning
from scipy import ndimage
from scipy.signal import savgol_filter
from astropy.wcs import WCS
from astropy.io import fits

try:
    from importlib.resources import files as import_pack_files
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    # For Py<3.9 files is not available
    from importlib_resources import files as import_pack_files
try:
    from importlib.resources import open_text as pack_open_txt
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    from importlib_resources import open_text as pack_open_txt
import matplotlib.pyplot as plt
import os
import psutil as psu
import signal
import time
import traceback
import numpy as np
import copy
import warnings
import re
import subprocess
import omegaconf
from datetime import datetime

# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def __setitem__(self, key, value):
        if key not in self:
            # If it is a new item we only allow it if it is not Configuration or Original_Cube or if we are in setup_configuration
            try:
                function,variable,empty = traceback.format_stack()[-2].split('\n')
            except ValueError: 
                function,variable = traceback.format_stack()[-2].split('\n')
            function = function.split()[-1].strip()
            variable = variable.split('[')[0].strip()
            if variable == 'Original_Configuration' or variable == 'Configuration':
                if function != 'setup_configuration':
                    raise ProgramError("FAT does not allow additional values to the Configuration outside the setup_configuration in support_functions.")
        OrderedDict.__setitem__(self,key, value)
    #    "what habbens now")
    def insert(self, existing_key, new_key, key_value):
        time.sleep(0.05)
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    f"----!!!!!!!! YOUR {new_key} was appended at the end as you provided the non-existing {existing_key} to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")

Proper_Dictionary.__doc__=f'''
A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
'''

def calc_rings(Configuration,size_in_beams = 0., ring_size  = 0.):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    if size_in_beams == 0.:
        size_in_beams = np.max(Configuration['SIZE_IN_BEAMS'])
    print_log(f'''CALC_RINGS: Calculating the number of rings in the model.
{'':8s} size in beams = {size_in_beams}
{'':8s} ring_size = {ring_size}
{'':8s} the maximum amount of rings = {Configuration['MAX_SIZE_IN_BEAMS']}
''',Configuration,case= ['debug_start'])

    est_rings = round((size_in_beams-1./5.)/(ring_size)+2.)
    #if est_rings > 20 and Configuration['MAX_SIZE_IN_BEAMS'] > 25:
    if est_rings > 20:
        Configuration['OUTER_RINGS_DOUBLED'] = True
        # 0.5 should always be rounding up else we get problemsin the reverse
        no_rings = 2+10+np.floor((est_rings-12)/2.+0.5)
    else:
        Configuration['OUTER_RINGS_DOUBLED'] = False
        no_rings = est_rings
    print_log(f'''CALC_RINGS: We started with ring size = {ring_size} and the size in beams {size_in_beams}.
{'':8s}The model will have {no_rings} rings and the rings are {'doubled' if Configuration['OUTER_RINGS_DOUBLED'] else 'not doubled'}.
''',Configuration,case= ['debug_add'])
    return int(no_rings)

calc_rings.__doc__ =f'''
 NAME:
    calc_rings

 PURPOSE:
    Calculate the actual number of rings in the model from ring size and the size in beams:

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    size_in_beams = 0.
    Radius of the model in beams. Default from Configuration

    ring_size  = 0.
    size of the rings in the model. Default from Configuration

 OUTPUTS:
    the number of rings in the model

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def calculate_number_processes(Configuration):
    cores_set = False
    no_galaxies = set_limits(Configuration['CATALOGUE_END_ID']-Configuration['CATALOGUE_START_ID'],1,float('inf'))
    while not cores_set:
        no_process= int(np.floor(Configuration['NCPU']/ Configuration['PER_GALAXY_NCPU']))
        while no_process > no_galaxies:
            no_process -= 1
            cores_set = True
        extra_cores=  Configuration['NCPU']-  Configuration['PER_GALAXY_NCPU']*no_process
        if not cores_set and \
            extra_cores >= Configuration['PER_GALAXY_NCPU']*no_process/2. and \
            Configuration['PER_GALAXY_NCPU'] > 3.:
            Configuration['PER_GALAXY_NCPU'] -= 1
        else:
            cores_set = True


    while no_process < extra_cores:
        extra_cores -= no_process
        Configuration['PER_GALAXY_NCPU'] += 1

    input_cores = copy.deepcopy(Configuration['PER_GALAXY_NCPU'])
    while no_process <= extra_cores:
        Configuration['PER_GALAXY_NCPU'] += 1
        extra_cores -= no_process
    if Configuration['PER_GALAXY_NCPU'] > 2.*input_cores:
        Configuration['PER_GALAXY_NCPU'] =  2.*input_cores
    sofia_processes = Configuration['NCPU']
    if sofia_processes > no_galaxies:
        sofia_processes = no_galaxies
    print(f"We use {no_process} processes for {no_galaxies} galaxies.")
    print(f"With {Configuration['PER_GALAXY_NCPU']} cpus per galaxy.")
    print(f"For SoFiA and initial estimates we use {sofia_processes} processes.")
    #The updated Configuration should be modified automatically
    return no_process,sofia_processes

calculate_number_processes.__doc__ =f'''
 NAME:
    calculate_number_processes

 PURPOSE:
    Calculate the optimal number of processs and create Configurations for each of them.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard Original FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    a set of Configurations for each process.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def complex_am_invert(Configuration,Theta_in,Phi_in):
    print_log(f'''COMPLEX_AM_INVERT:  as PA or INCL in {Configuration['FIXED_PARAMETERS']} we are doing a complex invert.
''',Configuration,case= ['debug_start'])

    Theta_out = copy.deepcopy(Theta_in)
    Phi_out = copy.deepcopy(Phi_in)
    Theta_out[0]=Theta_in[0]
    Phi_out[0]=Phi_in[0]
    if 'INCL' in Configuration['FIXED_PARAMETERS'] and 'PA' in Configuration['FIXED_PARAMETERS']:
        Theta_out[:]=Theta_in[0]
        Phi_out[:]=Phi_in[0]
        print_log(f'''COMPLEX_AM_INVERT:  both PA and INCL are fixed returning the centra Phi and Theta.
''',Configuration,case= ['debug_add'])
    else:
        #Create the PA and Inclination plane
        Theta = np.linspace(0,np.pi,1000)
        Phi = np.linspace(0,np.pi,1000)


        if 'INCL' in Configuration['FIXED_PARAMETERS']:
            Inclination_plane =   np.zeros((len(Theta),len(Phi)))
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",message="divide by zero encountered in true_divide"\
                                    ,category=RuntimeWarning)
                warnings.filterwarnings("ignore",message="divide by zero encountered in divide"\
                                    ,category=RuntimeWarning)
                for i in range(len(Theta)):
                    Inclination_plane[:,i]= 90-np.arctan(1./(np.cos(Theta[i])*np.tan(Phi[:])))*(360./(2*np.pi))
            Inclination_plane[np.where(Inclination_plane > 90.)] = 180 - Inclination_plane[np.where(Inclination_plane > 90.)]
            Inclination_plane[np.where(Inclination_plane < 0.)] = -1 * Inclination_plane[np.where(Inclination_plane < 0.)]
            current_inclination =  90-np.arctan(1./(np.cos(Theta_in[:])*np.tan(Phi_in[:])))*(360./(2*np.pi))
            current_inclination[np.where(current_inclination > 90.)] = 180 - current_inclination[np.where(current_inclination > 90.)]
            current_inclination[np.where(current_inclination < 0.)] = -1 * current_inclination[np.where(current_inclination < 0.)]
            current_diff = np.array([abs(x-current_inclination[0]) for x in current_inclination],dtype=float)
            not_flat = np.where(current_diff > 0.25)[0]
            if not_flat.size != 0:
                options_are = np.where(np.logical_and(current_inclination[0]-0.25 < Inclination_plane, Inclination_plane < current_inclination[0]+0.25))
                options_Phi = Phi[options_are[0][:]]
                options_Theta = Theta[options_are[1][:]]
                for i in range(not_flat.size):
                    diff = [abs(x-Phi_in[not_flat[i]])+abs(y-Theta_in[not_flat[i]]) for x,y in zip(options_Phi,options_Theta)]
                    location= np.where(np.min(diff) == diff)[0]
                    Phi_out[not_flat[i]]=options_Phi[location[0]]
                    Theta_out[not_flat[i]]=options_Theta[location[0]]
        if 'PA' in Configuration['FIXED_PARAMETERS']:
            PA_plane =   np.zeros((len(Theta),len(Phi)))
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore",message="divide by zero encountered in true_divide"\
                                    ,category=RuntimeWarning)
                warnings.filterwarnings("ignore",message="divide by zero encountered in divide"\
                                    ,category=RuntimeWarning)
                for i in range(len(Theta)):
                    PA_plane[:,i] = abs(np.arctan(np.sin(Theta[i])*np.tan(Phi[:]))*(360./(2*np.pi)))

            current_pa = abs(np.arctan(np.sin(Theta_in[:])*np.tan(Phi_in[:]))*(360./(2*np.pi)))
            current_diff = np.array([abs(x-current_pa[0]) for x in current_pa],dtype=float)
            not_flat = np.where(current_diff > 0.25)[0]
            if not_flat.size != 0:
                print(f'Starting the search for a constant PA')
                options_are = np.where(np.logical_and(current_pa[0]-0.25 < PA_plane, PA_plane < current_pa[0]+0.25))
                options_Phi = Phi[options_are[0][:]]
                options_Theta = Theta[options_are[1][:]]
                for i in range(not_flat.size):
                    diff = [abs(x-Phi_in[not_flat[i]])+abs(y-Theta_in[not_flat[i]]) for x,y in zip(options_Phi,options_Theta)]
                    location= np.where(np.min(diff) == diff)[0]
                    Phi_out[not_flat[i]]=options_Phi[location[0]]
                    Theta_out[not_flat[i]]=options_Theta[location[0]]


    return Theta_out,Phi_out


def calculate_am_vector(Configuration,PA_in,Inclination_in, multiple = None,
        invert=False):
    if multiple is None and invert:
        print_log(f'''CALCULATE_AM_VECTOR: inverting phi and theta without the multiple will result in all PAs being in the first quadrant.
''',Configuration,case= ['debug_start'])
        multiple=np.zeros(len(PA_in))

    if invert:
        Theta=np.array(PA_in,dtype=float)
        Phi=np.array(Inclination_in,dtype=float)
        if 'INCL' in Configuration['FIXED_PARAMETERS'] or 'PA' in Configuration['FIXED_PARAMETERS']:
            Theta,Phi = complex_am_invert(Configuration,Theta,Phi)

        Inclination=90-np.arctan(1./(np.cos(Theta)*np.tan(Phi)))*(360./(2*np.pi))
    # return inclination boundary adjustements
        Inclination[np.where(Inclination > 90.)] = 180 - Inclination[np.where(Inclination > 90.)]
        Inclination[np.where(Inclination < 0.)] = -1 * Inclination[np.where(Inclination < 0.)]
        if (Phi[0] < np.pi/2.) and (Phi[-1] > np.pi/2) and (Inclination[0] < 5.):
            PA= np.arctan(np.sin(Theta)*np.tan(Phi-np.pi/2.))*(360./(2*np.pi))\
                        +multiple*90.+np.arctan(np.sin(Theta[0])*np.tan(Phi[0]))*(360./(2*np.pi))
        else:
            PA= abs(np.arctan(np.sin(Theta)*np.tan(Phi))*(360./(2*np.pi)))+multiple*90.
        # And make sure the inner 3 rings are the same
        Inclination[0:3] = np.mean(Inclination[0:3])
        PA[0:3] = np.mean(PA[0:3])
        return PA,Inclination
    else:
        # For this the PA has to be between 0-90
        PA=np.array(PA_in,dtype=float)
        Inclination=np.array(Inclination_in,dtype=float)
        multiple=np.floor(PA/90.)
        PA= PA-multiple*90.
        #if any(PA > 90.):

        Inclination= 90.-Inclination
        PA[np.where(PA == 0.)] = 0.001
        Inclination[np.where(Inclination == 0.)] = 0.001
        Theta=np.arctan(abs(np.tan(Inclination*(np.pi/180.)))*abs(np.tan(PA*(np.pi/180.))))
        Phi = np.arctan(abs(np.tan(PA*(np.pi/180.)))/np.sin(Theta))

        return Theta,Phi,multiple
calculate_am_vector.__doc__ =f'''
 NAME:
    calculate_am_vector

 PURPOSE:
    Calculate the optimal number of processs and create Configurations for each of them.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard Original FAT configuration

 OPTIONAL INPUTS:

 OUTPUTS:
    a set of Configurations for each process.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

'This is for singular arrays not both sides'
def calculate_change_angle(Configuration,Theta,Phi):
    theta_zero = Theta[0]
    phi_zero = Phi[0]
       
    theta_change= np.array([float(x-theta_zero) for x in Theta]\
                           ,dtype=float)
    
    phi_change=  np.array([float(x-phi_zero) for x in Phi]\
                           ,dtype=float)
   
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore",message="invalid value encountered in true_divide"\
                            ,category=RuntimeWarning)
        #pYthon being really incredibly dumb again as changed the run time warning 
        warnings.filterwarnings("ignore",message="invalid value encountered in divide"\
                            ,category=RuntimeWarning)
        theta_factor = np.sqrt(theta_change**2/(theta_change**2+phi_change**2))\
                    *(theta_change)/abs(theta_change)
        
        theta_factor[np.where(np.array(theta_change) == 0.)] = 0.
        phi_factor = np.sqrt(phi_change**2/(theta_change**2+phi_change**2))*(phi_change)/abs(phi_change)
        phi_factor[np.where(np.array(phi_change) == 0.)] = 0.

    in_zero = np.where(np.array(theta_change+phi_change) == 0.)
    phi_factor[in_zero]=0.
    theta_factor[in_zero]=0.

    change_angle = np.sqrt(theta_change**2+phi_change**2)
    change_angle[in_zero] =0.

    return {'CHANGE_ANGLE': change_angle, 'PHI': {'ZERO': phi_zero, 'FACTOR':phi_factor},\
                                     'THETA':{'ZERO': theta_zero, 'FACTOR':theta_factor}}

def revert_change_angle(Configuration,change_angle):
    new_theta = change_angle['THETA']['ZERO']+change_angle['CHANGE_ANGLE']\
            *change_angle['THETA']['FACTOR']
    new_phi =  change_angle['PHI']['ZERO']+change_angle['CHANGE_ANGLE']\
            *change_angle['PHI']['FACTOR']
    return new_theta,new_phi

def check_angular_momentum_vector(Configuration,radius_in,pa_in,inclination_in,\
                                    modified= False,side=0):

    inclination = np.array(inclination_in)
    pa = np.array(pa_in)
    inclination_in = np.array(inclination_in)
    pa_in = np.array(pa_in)
    radius = np.array(radius_in,dtype=float)
    radkpc = convertskyangle(Configuration,radius)
    #radkpc = np.array([convertskyangle(Configuration,float(x)) for x in radius],dtype=float)
    #phi is dependent on theta but the max change should be the same
    max_shift = np.arctan(np.tan(Configuration['MAX_CHANGE']['INCL']*(np.pi/180.))\
                    *np.tan(Configuration['MAX_CHANGE']['PA']*(np.pi/180.)))
  
    print_log(f'''CHECK_ANGULAR_MOMENTUM_VECTOR: The maximum allowed shift  = {max_shift}.
''',Configuration,case= ['debug_start'])
    succes = False
    counter = 0.
    while not succes:
        #print_log(f'''CHECK_ANGULAR_MOMENTUM_VECTOR: Calculating Phi and Theta from these PA and inclination
#PA = {pa}
#Inclination = {inclination}
#''',Configuration,case= ['debug_add'])
        Theta,Phi,quadrant = calculate_am_vector(Configuration,pa,inclination )
        change_angle = calculate_change_angle(Configuration,Theta,Phi) 
        new_change_angle = max_profile_change(Configuration,radkpc,change_angle['CHANGE_ANGLE'],'CHANGE_ANGLE',\
            slope = Configuration['WARP_SLOPE'][side],max_change=max_shift, kpc_radius=True, quadrant=quadrant )
        
        change_angle['CHANGE_ANGLE'] = new_change_angle
        new_theta,new_phi = revert_change_angle(Configuration,change_angle)
       
     
        diff_phi = np.array(np.where(np.array([abs(x-y) for x,y in zip(Phi,new_phi) ],dtype=float) > 1e-6))
        diff_theta = np.array(np.where(np.array([abs(x-y) for x,y in zip(Theta,new_theta)],dtype=float) > 1e-6))

        if diff_phi.size != 0. or \
            diff_theta.size != 0.:
            #print_log(f'''CHECK_ANGULAR_MOMENTUM_VECTOR
#{'':8s} Phi = {Phi}, new_phi = {new_phi}
#{'':8s} Phi diff = {np.array([abs(x-y) for x,y in zip(Phi,new_phi) ],dtype=float)}{'':8s}
#{'':8s} Theta = {Theta}, new_theta = {new_theta}
#{'':8s} Theta diff = {np.array([abs(x-y) for x,y in zip(Theta,new_theta) ],dtype=float)}
#''',Configuration,case= ['debug_add'])


            pa,inclination = calculate_am_vector(Configuration,new_theta,new_phi,\
                                multiple=quadrant,invert=True )
            # Because when the quadrant changes the PA and inclination are not checked against maximum variations 
            # Hence this stille needs to be done
            pa,inclination = check_quadrant_change(Configuration,pa,inclination,quadrant)
            #for i in range(len(pa)):
            #    print(f'PA = {pa[i]}, inclination = {inclination[i]}, quadrant = {quadrant[i]} ')
#            print_log(f'''CHECK_ANGULAR_MOMENTUM_VECTOR: new PA and inclination
#PA = {pa}
#Inclination = {inclination}
#''',Configuration,case= ['debug_add'])
            counter += 1
            if counter > 1000.:
                #It is not really a succes but we have to exit the loop and one point.
                succes = True
        else:
            succes = True
        #exit()

    print_log(f'''CHECK_ANGULAR_MOMENTUM_VECTOR: Before checking we have modified = {modified}
''',Configuration,case= ['debug_add'])
    #if len(np.where(np.array([abs(x-y) for x,y in zip(pa_in,pa) ],dtype=float) != 0.)) != 0. or \
    #    len(np.where(np.array([abs(x-y) for x,y in zip(inclination_in,inclination) ],dtype=float) != 0.)) != 0.:

    if not np.isclose(pa_in,pa, atol=float(Configuration['MIN_ERROR']['PA'][0])*3.).any() or \
        not np.isclose(inclination_in,inclination, atol=float(Configuration['MIN_ERROR']['INCL'][0])*3.).any():
        modified= True
    return pa,inclination, modified
check_angular_momentum_vector.__doc__ =f'''
 NAME:
    check_angular_momentum_vector

 PURPOSE:
    Check the output warp by ensuring that the angular momentum is varying smoothly

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard Original FAT configuration
    radius = radius of the model
    pa = PA of the model
    inclination = inclination of the model

 OPTIONAL INPUTS:
    angle_check = False
    indicates whether we have modified the angles

 OUTPUTS:
    pa_new = the modfied PA
    incl_new= the modified incl
    modified = boolean indicating whether the PA and incl ave beet modfied.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: PA and inclination are singular sides
'''

def check_quadrant_change(Configuration,pa,inclination,quadrant):
    diff_quadrant = [x -quadrant[0] for x in quadrant]
    if np.sum(diff_quadrant) != 0.:
        for i,x in enumerate(diff_quadrant):
            if x != 0. and (1 < i <len(diff_quadrant)-1):
                pa[i] = np.mean([pa[i-1],pa[i+1]]) 
                inclination[i] = np.mean([inclination[i-1],inclination[i+1]])
            elif i == len(diff_quadrant)-1:
                pa[i] = pa[i-1]
                inclination[i] = inclination[i-1]
            else:
                pass
    return pa,inclination
check_angular_momentum_vector.__doc__ =f'''
 NAME:
    check_angular_momentum_vector

 PURPOSE:
    Check the output warp by ensuring that the pa and inclination vary smoothly at the rings where the 
    angular momentum is switching quadrants

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard Original FAT configuration
    pa = PA of the model
    inclination = inclination of the model
    quadrants = the array that specifies the quadrants into which the am vector belongs

 OPTIONAL INPUTS:
   
 OUTPUTS:
    pa_new = the modfied PA
    incl_new= the modified incl
   
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: PA and inclination are singular sides
'''             
        
def check_sofia(Configuration,Fits_Files):
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits','_cat.txt']
    for file in files:
        if os.path.exists(f'{Configuration["FITTING_DIR"]}Sofia_Output/{Configuration["BASE_NAME"]}{file}'):
            pass
        else:
            raise InputError()

check_sofia.__doc__ =f'''
 NAME:
    check_sofia

 PURPOSE:
    Check that the sofia files are in place when no sofia stages is used

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:


 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_tiltogram(Configuration, tiltogram,inner_min=3):
    for i in [0,1]:
        theta_inner = np.array([np.mean(tiltogram[i,0:x+1,0:x+1]) for x in range(tiltogram.shape[1])],dtype=float)
        theta_mutual = np.array([np.mean(tiltogram[i,x:,0:x+1]) for x in range(tiltogram.shape[1])],dtype=float)
        #theta_outer = np.array([np.mean(tiltogram[i,x+1:,x+1:]) for x in range(tiltogram.shape[1])],dtype=float)
        #And then we want to apply the rules
        # (i) the difference between theta_inner and theta_mutual is larger than the differences observed at other radii
        # (ii) theta_inner < 5 deg
        # (iii) thetamut > 15 deg
        diff = np.array(abs(theta_inner-theta_mutual),dtype = float)
        print_log(f'''CHECK_TILTOGRAM: Checking the tiltogram in side {i}.
{'':8s} Diff = {diff}
''',Configuration,case= ['debug_start'])
        rings_found = False
        while not rings_found:
            ring_location = np.where(np.nanmax(diff) == diff)[0]
            if ring_location.size > 1:
                ring_location = ring_location[0]
            if theta_inner[ring_location] < 5. and theta_mutual[ring_location] > 15.:
                Configuration['INNER_FIX'][i] = int(set_limits(ring_location-1,inner_min,Configuration['NO_RINGS']*3./4.-1))
                rings_found = True
            else:
                diff[ring_location] = 0.
            if np.nansum(diff) == 0.:
                Configuration['INNER_FIX'][i] = int(Configuration['NO_RINGS']*3./4.-1)
                rings_found = True

check_tiltogram.__doc__ =f'''
NAME:
    check_tiltogram

PURPOSE:
   Set the inner_fix values based on the tiltograms

CATEGORY:
   support_functions

INPUTS:
   Configuration = standard FAT Configuration
   tiltogram = the array containing the tiltogram for both sides

OPTIONAL INPUTS:

   inner_min = minimum set of inner rings that should be fixed

OUTPUTS:
   Updates Configuration['INNER_FIX']

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:

'''
# clean the header
def clean_header(Configuration,hdr_in,two_dim=False,mask_file=False):
    hdr = copy.deepcopy(hdr_in)
    print_log(f'''CLEAN_HEADER: Starting to clean the header.
''',Configuration,case= ['debug_start'])
    keywords = ['CDELT','CUNIT','CRPIX','CRVAL','CTYPE','CROTA']
    for key in keywords:
        try:
            del hdr[f'{key}4']
        except KeyError:
            pass
    if not two_dim:
        hdr['NAXIS'] = 3
        if not 'CUNIT3' in hdr:
            if abs(hdr['CDELT3']) > 500:
                hdr['CUNIT3'] = 'm/s'
            else:
                hdr['CUNIT3'] = 'km/s'
            print_log(f'''CLEAN_HEADER:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
Your header did not have a unit for the third axis, that is bad policy.
{"":8s} We have set it to {hdr['CUNIT3']}. Please ensure that this is correct.'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration)
        if hdr['CUNIT3'].upper() == 'HZ' or hdr['CTYPE3'].upper() == 'FREQ':
            print_log('CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.', Configuration, case=['main','screen'])
            raise BadHeaderError('The Cube has frequency as a velocity axis this is not supported')

        #vel_types = ['VELO-HEL','VELO-LSR','VELO', 'VELOCITY','VOPT','VRAD']
        vel_types = ['VELO','VOPT','VRAD']
        if hdr['CTYPE3'].upper() not in vel_types:
            if hdr['CTYPE3'].split('-')[0].upper() in ['RA','DEC']:
                print_log(f'''CLEAN_HEADER: Your zaxis is a spatial axis not a velocity axis.
{"":8s}CLEAN_HEADER: Please arrange your cube logically
''',Configuration,case= ['main','screen'])
                raise BadHeaderError("The Cube's third axis is not a velocity axis")
            if hdr['CTYPE3'].split('-')[0].upper() == 'VELO':
                frame = hdr['CTYPE3'].split('-')[1].upper()
                hdr['CTYPE3'] = 'VELO'
                if frame in ['LSR','LSRK']:
                    hdr['SPECSYS3']='LSRK'
                elif frame in ['HEL','HELO']:
                    hdr['SPECSYS3']='HELIOCEN'
            else:
                print_log(f'''CLEAN_HEADER: Your velocity projection is not standard (CTYPE3 = {hdr["CTYPE3"]}). 
Please note that FAT can only deal with velocity cubes, not frequency.
''',Configuration, case=['main','screen'])
                raise BadHeaderError(f'The Cube has a non-standard velocity projection (CTYPE3 = {hdr["CTYPE3"]})')
        #Tirific will only run with VELO so if VRAD of VOPT then we need to run with VELO and fix it in the end
        Configuration['HDR_VELOCITY'] = hdr['CTYPE3']
        if  hdr['CTYPE3'].upper() in ['VRAD','VOPT']:
            print_log(f'''CLEAN_HEADER: you are use a velocity type of {hdr['CTYPE3'].upper()} however tirific takes issue with this.
As such we are replacing it with VELO before every tirific run. 
!!!!!!!!!!!! ---- We are not converting anything so your fitting and all products are are of type {hdr['CTYPE3'].upper()}-------!!!!!!!!                                              
''',Configuration)
            #hdr['CTYPE3'] = 'VELO'

        if 'SPECSYS3' not in hdr:
            if 'SPECSYS' in hdr:
                hdr['SPECSYS3'] = hdr['SPECSYS']
            else:
                print_log(f'''CLEAN_HEADER: There is no reference frame defined in your header. This is ok for FAT but WCS might complain.
''',Configuration)

        if hdr['CUNIT3'].lower() == 'km/s':
            print_log( f'''CLEAN_HEADER: The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'
''',Configuration)
            hdr['CUNIT3'] = 'm/s'
            hdr['CDELT3'] = hdr['CDELT3']*1000.
            hdr['CRVAL3'] = hdr['CRVAL3']*1000.

        #because astropy is truly stupid

        if hdr['CUNIT3'] == 'M/S':
            hdr['CUNIT3'] = 'm/s'

    if not 'EPOCH' in hdr:
        if 'EQUINOX' in hdr:
            print_log(f'''CLEAN_HEADER: Your cube has no EPOCH keyword but we found EQUINOX.
{"":8s}We have set EPOCH to {hdr['EQUINOX']}
''',Configuration)
            hdr['EPOCH'] = hdr['EQUINOX']
            del hdr['EQUINOX']
        else:
            print_log(f'''CLEAN_HEADER: Your cube has no EPOCH keyword
{"":8s}CLEAN_HEADER: We assumed J2000
''',Configuration)
            hdr['EPOCH'] = 2000.


    #FAT cannot deal with the PC nonsense
    if f'PC01_01' in hdr:
        del hdr['PC0*_0*']
    if f'PC1_1' in hdr:
        del hdr['PC*_*']


    # Check for the beam
    if not 'BMAJ' in hdr and not mask_file:
        if 'BMMAJ' in hdr:
            hdr['BMAJ']= hdr['BMMAJ']/3600.
        else:
            found = False
            if 'HISTORY' in hdr:
                for line in hdr['HISTORY']:
                    tmp = [x.strip().upper() for x in line.split()]
                    if 'BMAJ=' in tmp:
                        hdr['BMAJ'] = float(tmp[tmp.index('BMAJ=') + 1])
                        found = True
                    if 'BMIN=' in tmp:
                        hdr['BMIN'] = float(tmp[tmp.index('BMIN=') + 1])
                    if 'BPA=' in tmp:
                        hdr['BPA'] = float(tmp[tmp.index('BPA=') + 1])
                    if found:
                        break
            if not found:
                print_log(f'''CLEAN_HEADER: WE CANNOT FIND THE MAJOR AXIS FWHM IN THE HEADER
''',Configuration,case= ['main','screen'])
                raise BadHeaderError("The Cube has no major axis FWHM in the header.")
    if not 'CTYPE1' in hdr or not 'CTYPE2' in hdr:
        print_log(f'''CLEAN_HEADER: Your spatial axes have no ctype. this can lead to errors.
''',Configuration,case= ['main','screen'])
        raise BadHeaderError("The Cube header has no ctypes.")

    if hdr['CTYPE1'].split('-')[0].upper() in ['DEC']:
        print_log(f'''CLEAN_HEADER: !!!!!!!!!!Your declination is in the first axis. !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
''',Configuration,case= ['main','screen'])
        raise BadHeaderError("Your spatial axes are reversed")
    if hdr['CTYPE2'].split('-')[0].upper() in ['RA']:
        print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your right ascension is on the second axis. !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
''',Configuration,case= ['main','screen'])
        raise BadHeaderError("Your spatial axes are reversed")
    if hdr['CRVAL1'] < 0.:
        print_log(f'''CLEAN_HEADER: your RA crval is negative, this can lead to errors. Adding 360. deg
''',Configuration)
        hdr['CRVAL1'] = hdr['CRVAL1']+360.
    if not 'BMIN' in hdr and not mask_file:
        if 'BMMIN' in hdr:
            hdr['BMIN']= hdr['BMMIN']/3600.
        else:
            print_log(f'''CLEAN_HEADER: We cannot find the minor axis FWHM. Assuming a circular beam.
''',Configuration)
            hdr['BMIN'] = hdr['BMAJ']
    if not 'BPA' in hdr and not mask_file:
            print_log(f'''CLEAN_HEADER: We cannot find the Beam PA assuming it to be 0
''',Configuration)
            hdr['BPA'] = 0.

    try:
        if len(hdr['HISTORY']) > 10:
            del hdr['HISTORY']
            if Configuration:
                print_log( f'''CLEAN_HEADER: Your cube has a significant history attached we are removing it for easier interpretation.
''',Configuration)
    except KeyError:
        pass
    #The blank value only applies to integer value data
    try:
        if int(hdr['BITPIX']) != 16:
            del hdr['BLANK']
    except KeyError:
        pass

    try:
        if abs(hdr['BMAJ']/hdr['CDELT1']) < 2:
            print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your cube has less than two pixels per beam major axis.!!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!           This will lead to bad results.              !!!!!!!!!!!!!!!!'
''',Configuration)

        if abs(hdr['BMAJ']/hdr['CDELT1']) > hdr['NAXIS1']:
            print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your cube is smaller than the beam major axis. !!!!!!!!!!!!!!!!!
    {"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
    ''',Configuration,case= ['main','screen'])
            raise BadHeaderError("Your cube is too small for your beam")
    except KeyError:
        pass
    #Make sure the header is WCS compatible
    if not mask_file:
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                check_wcs = WCS(hdr)
        except:
            raise BadHeaderError("The header of your cube contains keywords that confuse astropy. Please make your header astropy WCS compatible")

    return hdr
clean_header.__doc__ =f'''
 NAME:
    clean_header

 PURPOSE:
    Clean up the cube header and make sure it has all the right
    variables that we require in the process of fitting

 CATEGORY:
    supprot_functions

 INPUTS:
   Configuration = Standard FAT configuration
    hdr = header to be cleaned

 OPTIONAL INPUTS:


 OUTPUTS:
    the updated header

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def columndensity(Configuration,levels,systemic = 100.,beam=None,\
        channel_width=-1.,column= False,arcsquare=False,solar_mass_input =False\
        ,solar_mass_output=False):
    if beam is None:
        if not arcsquare:
            beam = Configuration['BEAM'][:2]
    if channel_width == -1:
        if arcsquare:
            channel_width = 1.
        else:
            channel_width = Configuration['CHANNEL_WIDTH']
    print_log(f'''COLUMNDENSITY: Starting conversion from the following input.
{'':8s}Levels = {levels}
{'':8s}Beam = {beam}
{'':8s}channel_width = {channel_width}
''',Configuration,case= ['debug_start'])
    beam=np.array(beam,dtype=float)
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg
    print_log(f'''COLUMNDENSITY: We have the following input for calculating the columns.
{'':8s}COLUMNDENSITY: level = {levels}, channel_width = {channel_width}, beam = {beam}, systemic = {systemic})
''',Configuration,case= ['debug_add'])

    f = f0 * (1 - (systemic / c)) #Systemic frequency
    if arcsquare:
        #Should we have the (f0/f)**2 factor here????
        HIconv = 605.7383 * 1.823E18 * (2. *np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass_input:
                levels=np.array([x*solarmass/(mHI*pc**2) for x in levels],dtype=float)
            #levels=levels/(HIconv*channel_width)

            levels = levels/(HIconv*channel_width)
        else:

            levels = HIconv*levels*channel_width
            if solar_mass_output:
                levels=levels*mHI/solarmass*pc*pc
    else:
        if beam.size <2:
            beam= [beam,beam]
        b=beam[0]*beam[1]
        if column:
            if solar_mass_input:
                levels=levels*solarmass/(mHI*pc**2)
            TK = levels/(1.823e18*channel_width)
            levels = TK/(((605.7383)/(b))*(f0/f)**2)
        else:
            TK=((605.7383)/(b))*(f0/f)**2*levels
            levels = TK*(1.823e18*channel_width)
    if not column and solar_mass_input:
        levels = levels*mHI*pc**2/solarmass
    return levels

columndensity.__doc__ =f'''
 NAME:
    columndensity

 PURPOSE:
    Convert the various surface brightnesses to other units

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    levels = the values to convert

 OPTIONAL INPUTS:


    systemic = 100.
    the systemic velocity of the source

    beam  = [-1.,-1.]
    the FWHM of the beam in arcsec, if unset taken from Configuration

    channelwidth = -1. width of a channel in km/s
    channelwidth of the observation if unset taken from Configuration

    column = false
    if True input is columndensities else in mJy

    arcsquare=False
    If true then  input is assumed to be in Jy/arcsec^2.
    If the input is in Jy/arcsec^2*km/s then channelwidth must be 1.
    This is assumed when channelwidth is left unset

    solar_mass_input =False
    If true input is assumed to be in M_solar/pc^2

    solar_mass_output=False
    If true output is provided in M_solar/pc^2

 OUTPUTS:
    The converted values

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Cosmological column densities are taken from Meyer et al 2017
'''

        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(Configuration,RAin,DECin,invert=False, colon=False):
    print_log(f'''CONVERTRADEC: Starting conversion from the following input.
    {'':8s}RA = {RAin}
    {'':8s}DEC = {DECin}
''',Configuration,case= ['debug_start'])
    RA = copy.deepcopy(RAin)
    DEC = copy.deepcopy(DECin)
    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA= [RA]
            DEC =[DEC]
        for i in range(len(RA)):
            xpos=RA
            ypos=DEC
            xposh=int(np.floor((xpos[i]/360.)*24.))
            xposm=int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh=int(np.floor(np.absolute(ypos[i]*1.)))
            yposm=int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs=(((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign=ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i]="{}:{}:{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}:{}:{:2.2f}".format(yposh,yposm,yposs)
            else:
                RA[i]="{}h{}m{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}d{}m{:2.2f}".format(yposh,yposm,yposs)
            if sign < 0.: DEC[i]='-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA,str):
            RA=[RA]
            DEC=[DEC]

        xpos=RA
        ypos=DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+",xpos[i])
            RA[i]=(float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+",ypos[i])
            if float(tmp[0]) != 0.:
                DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))
            else:
                DEC[i] = float(np.absolute(float(tmp[0])) + ((float(tmp[1]) + (float(tmp[2]) / 60.)) / 60.))
                if tmp[0][0] == '-':
                    DEC[i] = float(DEC[i])*-1.
        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])
        else:
            RA =np.array(RA,dtype=float)
            DEC = np.array(DEC,dtype=float)
    return RA,DEC

convertRADEC.__doc__ =f'''
 NAME:
    convertRADEC

 PURPOSE:
    convert the RA and DEC in degre to a string with the hour angle

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    RAin = RA to be converted
    DECin = DEC to be converted

 OPTIONAL INPUTS:


    invert=False
    if true input is hour angle string to be converted to degree

    colon=False
    hour angle separotor is : instead of hms

 OUTPUTS:
    converted RA, DEC as string list (hour angles) or numpy float array (degree)

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def convert_type(array, type = 'float'):
    if type =='int':
        return  [int(x) for x in array]
    elif type =='str':
        return  [str(x) for x in array]
    else:
        return  [float(x) for x in array]

convert_type.__doc__ =f'''
 NAME:
    convert_type

 PURPOSE:
    Convert a list of variables from one type to another and be able to have them unpack into single varaiable again.

 CATEGORY:
    support_functions

 INPUTS:
    array = array to be converted

 OPTIONAL INPUTS:
    type = 'float'
    type to convert to

 OUTPUTS:
    converted list with requested types

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# function for converting kpc to arcsec and vice versa
def convertskyangle(Configuration, angle, distance=-1., unit='arcsec', \
        distance_unit='Mpc', physical=False):
    if distance == -1.:
        distance = Configuration['DISTANCE']
    print_log(f'''CONVERTSKYANGLE: Starting conversion from the following input.
    {'':8s}Angle = {angle}
    {'':8s}Distance = {distance}
''',Configuration,case= ['debug_start'])

    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle,dtype=float)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print_log(f'''CONVERTSKYANGLE: {distance_unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: please use Mpc, kpc or pc.
''',Configuration,case= ['main','screen'])
        raise SupportRunError('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.')
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print_log(f'''CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: arcsec, arcmin or degree.
''',Configuration,case= ['main','screen'])
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')


        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle * (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle / (10 ** 3)
        else:
            print_log(f'''CONVERTSKYANGLE: {unit} is an unknown unit to convertskyangle.
{'':8s}CONVERTSKYANGLE: please use Mpc, kpc or pc.
''',Configuration,case= ['main','screen'])
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')

        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc

convertskyangle.__doc__ =f'''
 NAME:
    convertskyangle

 PURPOSE:
    convert an angle on the sky to a distance in kpc or vice versa

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    angle = the angles or lengths to be converted

 OPTIONAL INPUTS:


    distance=1.
    Distance to the galaxy for the conversion

    unit='arcsec'
    Unit of the angle or length options are arcsec (default),arcmin, degree, pc, kpc(default) and Mpc

    distance_unit='Mpc'
    Unit of the distance options are pc, kpc and Mpc

    physical=False
    if true the input is a length converted to an angle

 OUTPUTS:
    converted value or values

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def copy_homemade_sofia(Configuration,no_cat=False):
    create_directory('Sofia_Output',Configuration['FITTING_DIR'])
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits']
    if not no_cat:
        files.append('_cat.txt')
    for file in files:
        source = get_system_string(f"{Configuration['SOFIA_DIR']}{Configuration['SOFIA_BASENAME']}{file}")
        target = get_system_string(f"{Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']}{file}")
        os.system(f'''rm -f {target}''')
        os.system(f'''cp {source} {target}''')
        if not os.path.exists(f"{Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']+file}"):
            if file in ['_mask.fits','_cat.txt']:
                print_log(f'''COPY_HOMEMADE_SOFIA: Something went wrong copying the file  {Configuration['SOFIA_DIR']}{Configuration['SOFIA_BASENAME']+file}
{'':8s} to the file {Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']+file}.
{'':8s}We are aborting this fit as we cannot make it without Sofia input.
''',Configuration,case= ['main','screen'])
                raise SofiaMissingError("Either the sofia mask or catalogue is missing. We can not run without it")
            else:
                pass

        elif file != '_cat.txt':
            #make sure the header is decent
            two_dim = True
            mask_file=False

            if file == '_mask.fits':
                two_dim = False
                mask_file = True

            initial = fits.open(f"{Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']}{file}")
            hdr_in = initial[0].header
            data=initial[0].data
            initial.close()
            hdr = clean_header(Configuration,hdr_in,two_dim=two_dim,mask_file=mask_file )
            update = False
            for key in hdr:
                try:
                    if hdr[key] != hdr_in[key]:
                        update =True
                except KeyError:
                    update = True
            if file == '_mom0.fits':
                if hdr['BUNIT'].strip().lower() == 'jy/beam*m/s':
                    update = True
                    hdr['BUNIT'] = 'Jy/beam*km/s'
                    data = data/1000.
            if file in ['_mom1.fits','_mom2.fits'] :
                if hdr['BUNIT'].strip().lower() == 'm/s':
                    update = True
                    hdr['BUNIT'] = 'km/s'
                    data = data/1000.

            if update:
                 fits.writeto(f"{Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']}{file}",data,hdr,overwrite=True)

copy_homemade_sofia.__doc__ =f'''
 NAME:
    copy_homemade_sofia

 PURPOSE:
    Copy user provided Sofia files to the FAT specified directory such that the original are a) kept in place, b) Never modified.

 CATEGORY:
    clean_functions

 INPUTS:
    Configuration = Standard FAT configuration


 OPTIONAL INPUTS:

    check = False
    if set to true FAT merely checks for the existings of the sofia files where they expect them
 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def create_directory(directory,base_directory):
    split_directory = [x for x in directory.split('/') if x]
    split_directory_clean = [x for x in directory.split('/') if x]
    split_base = [x for x in base_directory.split('/') if x]
    #First remove the base from the directory but only if the first directories are the same
    if split_directory[0] == split_base[0]:
        for dirs,dirs2 in zip(split_base,split_directory):
            if dirs == dirs2:
                split_directory_clean.remove(dirs2)
            else:
                if dirs != split_base[-1]:
                    raise InputError(f"You are not arranging the directory input properly ({directory},{base_directory}).")
    for new_dir in split_directory_clean:
        if not os.path.isdir(f"{base_directory}/{new_dir}"):
            os.mkdir(f"{base_directory}/{new_dir}")
        base_directory = f"{base_directory}/{new_dir}"
create_directory.__doc__ =f'''
 NAME:
    create_directory

 PURPOSE:
    create a directory recursively if it does not exists and strip leading directories when the same fro the base directory and directory to create

 CATEGORY:
    support_functions

 INPUTS:
    directory = string with directory to be created
    base_directory = string with directory that exists and from where to start the check from

 OPTIONAL INPUTS:


 OUTPUTS:

 OPTIONAL OUTPUTS:
    The requested directory is created but only if it does not yet exist

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def deproject(Configuration,map,angle, center = 0., invert = False):
    axis = np.array(range(len(map[:,0])),dtype=float)-center
    if invert:
        newaxis = axis/np.cos(np.radians(angle))
    else:
        newaxis = axis*np.cos(np.radians(angle))
    for x in range(len(map[0,:])):
        profile = copy.deepcopy(map[:,x])
        new_profile = np.interp(np.array(newaxis,dtype=float),np.array(axis,dtype=float),np.array(profile,dtype=float))
        map[:,x] = new_profile
    return map
deproject.__doc__ =f'''
 NAME:
    deproject

 PURPOSE:
    Deproject an image assuming circular orienation and an inclination given by the angle

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map = the image to deproject
    angle = the inclination angle

 OPTIONAL INPUTS:


    center = 0.
    height on the yaxis to rotate around

    invert = False
    If true the image is inclined by tangle instead of deprojected

 OUTPUTS:
    the deprojected map

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def ensure_list(variable):
    '''Make sure that variable is a list'''
    if not isinstance(variable,list):
        if not isiterable(variable):
            variable=[variable]
        else:
            variable=[x for x in variable]
    return variable
ensure_list.__doc__ =f'''
 NAME:
    ensure_list

 PURPOSE:
    Make sure that variable is a list

 CATEGORY:
    support_functions

 INPUTS:
    variable = variable to check and transform

 OPTIONAL INPUTS:

 OUTPUTS:
    variable = variable in list form

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def find_program(name,search):
    '''check whether a program is available for use.'''
    found = False
    while not found:
        try:
            run = subprocess.Popen([name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            run.stdout.close()
            run.stderr.close()
            os.kill(run.pid, signal.SIGKILL)
            found = True
        except OSError:
            name = input(f'''You have indicated to use {name} for using {search} but it cannot be found.
Please provide the correct name : ''')
    return name
find_program.__doc__ =f'''
 NAME:
    find_program

 PURPOSE:
    check whether a program is available for use.

 CATEGORY:
    support_functions

 INPUTS:
    name = command name of the program to run
    search = Program we are looking for
 OPTIONAL INPUTS:

 OUTPUTS:
    the correct command for running the program

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def finish_current_run(Configuration,current_run):
    print_log(f'''FINISH_CURRENT_RUN: Is Tirific Running? {Configuration['TIRIFIC_RUNNING']}.
''',Configuration)
    if Configuration['TIRIFIC_RUNNING']:
        try:
            current_run.stdout.close()
            current_run.stderr.close()
        except:
            print_log(f'''FINISH_CURRENT_RUN: We failed to close the pipe to the current run even though there should be one.
''',Configuration)
            pass
        try:
            os.kill(Configuration['TIRIFIC_PID'], signal.SIGKILL)
            print_log(f'''FINISH_CURRENT_RUN: We killed PID = {Configuration['TIRIFIC_PID']}.
''',Configuration)
        except:
            try:
                current_run.kill()
                print_log(f'''FINISH_CURRENT_RUN: We killed the current run although we failed on the PID = {Configuration['TIRIFIC_PID']}.
''',Configuration)
            except AttributeError:
                print_log(f'''FINISH_CURRENT_RUN: We failed to kill the current run with PID {Configuration['TIRIFIC_PID']} even though we have tirific running
''',Configuration,case=['main','screen'])
                raise TirificKillError('FINISH_CURRENT_RUN: Despite having an initialized tirific we could not kill it. This should not happen.')
        Configuration['TIRIFIC_RUNNING'] = False
        Configuration['TIRIFIC_PID'] = 'Not Initialized'
    else:
        print_log(f'''FINISH_CURRENT_RUN: No run is initialized.
''',Configuration)
finish_current_run.__doc__ =f'''
 NAME:
    finish_current_run
 PURPOSE:
    make sure that the initiated tirific is cleaned when pyFAT stops
 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    current_run = subprocess structure for the current tirific run

 OPTIONAL INPUTS:


 OUTPUTS:
    kills tirific if initialized or raises an error when it fails to do so while tirific is running

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fit_center_ellipse(Configuration,map,inclination= -1, pa = -1):
    max = np.max(map)
    ellipses = np.linspace(max/2.,max,num=10)
    parameters =[]
    for el_out,el_in in zip(ellipses,ellipses[1:]):
        map_to_fit = np.zeros(map.shape)
        map_to_fit[map < el_in ] = 1.
        map_to_fit[map < el_out ] = 0.
        y,x = np.where(map_to_fit == 1.)
        cx,cy,a,b, orient = fit_ellipse(Configuration,x,y)
        if np.isnan(cx) or np.isnan(cx):
            continue
        else:
            if a < b:
                tmp = copy.deepcopy(a)
                a = copy.deepcopy(b)
                b = copy.deepcopy(tmp)
                orient = orient-90.


            ratio = float(set_limits(b/a,0.2,1.0))
            inc = np.degrees(np.arccos(np.sqrt((ratio**2-0.2**2)/0.96)))
            parameters.append([cx,cy,inc , orient,len(y)])

    x_in = [x[0] for x  in parameters]
    y_in = [x[1] for x  in parameters]
    if inclination == -1. or pa == -1.:
        weight = [x[-1]  for x in parameters ]
    else:
        weight = [x[-1]/((abs(x[2]-inclination)+abs(x[3]-inclination)))  for x in parameters ]
    x_center = np.sum([x*y for x,y in zip(x_in,weight)])/np.sum(weight)
    y_center = np.sum([x*y for x,y in zip(y_in,weight)])/np.sum(weight)

    return [x_center,y_center]
fit_center_ellipse.__doc__= '''
NAME:
   fit_center_ellipse
PURPOSE:
   determine the center through ellipse fits

CATEGORY:
   support_functions

INPUTS:
  Configuration = stadard configurstion
  map = moment 0 map
  inclination= current inclination
  pa = current pa

OPTIONAL INPUTS:



OUTPUTS:
   center
   the fitted center in pixels

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def fit_ellipse(Configuration,x,y):
    x=x[:,None]
    y=y[:,None]
    diag = np.hstack([x*x,x*y,y*y,x,y,np.ones(x.shape)])
    dmult = np.dot(diag.T,diag)
    constant=np.zeros([6,6])
    constant[0,2]=constant[2,0]=2
    constant[1,1]=-1
    try:
        eigen,vec=np.linalg.eig(np.dot(np.linalg.inv(dmult),constant))
    except LinAlgError:
        return float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')
    n= np.argmax(eigen)
    # if we have more then 3 arguments then we have complex numbers which is pointless for our center
    if n > 3:
        return float('NaN'),float('NaN'),float('NaN'),float('NaN'),float('NaN')
    a=vec[:,n]
        #-------------------Fit ellipse-------------------
    b,c,d,f,g,a=a[1]/2., a[2], a[3]/2., a[4]/2., a[5], a[0]
    bb=b**2
    num=b**2-a*c
    new_x=(c*d-b*f)/num
    new_y=(a*f-b*d)/num

    angle=0.5*np.arctan(2*b/(a-c))*180/np.pi+90.
    up = 2*(a*f**2+c*d**2+g*bb-2*b*d*f-a*c*g)
    down1=(bb-a*c)*( (c-a)*np.sqrt(1+4*bb/((a-c)*(a-c)))-(c+a))
    down2=(bb-a*c)*( (a-c)*np.sqrt(1+4*bb/((a-c)*(a-c)))-(c+a))
    a=np.sqrt(abs(up/down1))
    b=np.sqrt(abs(up/down2))

    return new_x,new_y,a,b,angle

fit_ellipse.__doc__= '''
NAME:
   fit_ellipse
PURPOSE:
   fit ellipse to a bunch of coordinates

CATEGORY:
   support_functions

INPUTS:
  Configuration = stadard configurstion
  x = x -coordinates to evaluated
  y = y -coordinates to evaluated

OPTIONAL INPUTS:



OUTPUTS:
   new_x = new x center
   new_y = new y center
   a =  major axis
   b = minor axis
   angle= orientaion of ellipse


OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''

def fit_sine(Configuration,x,y):
    print_log(f'''FIT_SINE: Starting to fit a sin.
{'':8s}x = {x}
{'':8s}y = {y}
{'':8s} x size = {x.size} y size = {y.size}
''', Configuration,case=['debug_start'])
    # Make sure we have numpy arrays
    x= np.array(x,dtype=float)
    y= np.array(y,dtype=float)
    est_peak = np.nanmax(y)-np.mean(y)
    peak_location = np.where(y == np.nanmax(y))[0]
    if peak_location.size > 1:
        peak_location = int(peak_location[0])
    min_location =  np.where(y == np.nanmin(y))[0]
    if min_location.size > 1:
        min_location = int(min_location[0])
    est_width = float(abs(x[peak_location]-x[min_location])/(2.*np.pi))
    est_amp = np.mean(y)
    peak_location = np.where(y == np.nanmax(y))[0]
    if peak_location.size > 1:
        peak_location = int(peak_location[0])

    est_center = float(x[peak_location])

    est_width=est_width*(2.*np.pi/180.)
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            sin_parameters, sin_covariance = curve_fit(sine, x[~np.isnan(y)], y[~np.isnan(y)],p0=[est_peak,est_center,est_width,est_amp])
        except RuntimeError:
            sin_parameters = [float('NaN') for z in range(4)]
        except OptimizeWarning:
            print_log(f'''FIT_SINE: the covariance could not be estimated. Using initial estimates
{'':8s}peak est = {est_peak}
{'':8s}center est = {est_center}
{'':8s}width est = {est_width}
{'':8s}amp est = {est_amp}
''', Configuration,case=['debug_add'])
            sin_parameters = [0.,0.,0.,0.]

    if not 0.4 < sin_parameters[2] <0.6:
        ratios = y
        sin_parameters = [est_peak,est_center,est_width,est_amp]
    else:
        ratios = sine(x,*sin_parameters)
    return ratios,sin_parameters

fit_sine.__doc__= '''
NAME:
   fit_sine
PURPOSE:
   Fit a modified sine function to a profile, with initial estimates

CATEGORY:
   support_functions

INPUTS:
   x = x-axis of profile
   y = y-axis of profile
   Configuration = Standard FAT configuration

OPTIONAL INPUTS:



OUTPUTS:
   ratios
   the fitted sin profile or when the  width is too small or to wide the original profiles

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def fit_gaussian(Configuration,x,y, covariance = False,errors = None):
    print_log(f'''FIT_GAUSSIAN: Starting to fit a Gaussian.
{'':8s}x = {x}
{'':8s}y = {y}
''', Configuration,case=['debug_start'])
    # Make sure we have numpy arrays
    x= np.array(x,dtype=float)
    y= np.array(y,dtype=float)
    # First get some initial estimates
    est_peak = np.nanmax(y)
    if not errors.any():
        errors = np.full(len(y),1.)
        absolute_sigma = False
    else:
        absolute_sigma = True
    peak_location = np.where(y == est_peak)[0]
    if peak_location.size > 1:
        peak_location = peak_location[0]
    est_center = float(x[peak_location])
    est_sigma = np.nansum(y*(x-est_center)**2)/np.nansum(y)

    with warnings.catch_warnings():
        warnings.simplefilter("error")
        succes = False
        maxfev= int(100*(len(x)))

        while not succes:
            print_log(f'''FIT_GAUSSIAN: Starting the curve fit with {maxfev}
''',Configuration,case=['debug_add'])
            try:
                gauss_parameters, gauss_covariance = curve_fit(gaussian_function, \
                        x, y,p0=[est_peak,est_center,est_sigma],sigma= errors,\
                        absolute_sigma= absolute_sigma,maxfev=maxfev)
                succes = True
            except OptimizeWarning:
                maxfev =  2000*(len(x))
            except RuntimeError as e:
                split_error = str(e)
                if 'Optimal parameters not found: Number of calls to function has reached maxfev' in \
                    split_error:
                    maxfev += 100*int(len(x))
                    print_log(f'''FIT_GAUSSIAN: We failed to find an optimal fit due to the maximum number of evaluations. increasing maxfev to {maxfev}
''',Configuration,case=['debug_add'])
                else:
                    print_log(f'''FIT_GAUSSIAN: failed due to the following error {split_error}
''',Configuration,case=['main','screen'])
                    raise RuntimeError(split_error)
            if maxfev >  1000*(len(x)):
                raise FittingError("FIT_GAUSSIAN: failed to find decent gaussian parameters")

    #gauss_parameters, gauss_covariance = curve_fit(gaussian_function, x, y,p0=[est_peak,est_center,est_sigma],sigma= errors,absolute_sigma= absolute_sigma)
    if covariance:
        return gauss_parameters, gauss_covariance
    else:
        return gauss_parameters

fit_gaussian.__doc__ =f'''
 NAME:
    fit_gaussian
 PURPOSE:
    Fit a gaussian to a profile, with initial estimates
 CATEGORY:
    supprt_functions

 INPUTS:
    x = x-axis of profile
    y = y-axis of profile
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    covariance = false
    return to covariance matrix of the fit or not



 OUTPUTS:
    gauss_parameters
    the parameters describing the fitted Gaussian

 OPTIONAL OUTPUTS:
    gauss_covariance
    The co-variance matrix of the fit

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def gaussian_function(axis,peak,center,sigma):
    return peak*np.exp(-(axis-center)**2/(2*sigma**2))

gaussian_function.__doc__ =f'''
 NAME:
    gaussian_function
 PURPOSE:
    Describe a Gaussian function

 CATEGORY:
    support_functions

 INPUTS:
    axis = the points where to evaluate the gaussian
    peak = amplitude of the peak of the Gaussian
    center = location of the peak on axis
    sigma = dispersion of the gaussian

 OPTIONAL INPUTS:

 KEYWORD PARAMETERS:

 OUTPUTS:
    The gaussian function

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 EXAMPLE:

 NOTE:
'''

def get_fit_groups(Configuration,Tirific_Template):
    parameter_groups = []
    block = []
    rings = []
    groups = Tirific_Template['VARY'].split(',')
    variation_type = []
    variation = []
    cut_off_limits = sbr_limits(Configuration,Tirific_Template )
    sbr_standard = np.mean(cut_off_limits) * 5. 
    paramater_standard_variation = {'PA': [10.,'a'],
                                   'INCL': [10.,'a'],
                                   'VROT': [Configuration['CHANNEL_WIDTH']*3.,'a'],
                                   'VSYS': [Configuration['CHANNEL_WIDTH'],'a'],
                                   'XPOS': [Configuration['BEAM'][0]/3600.,'a'],
                                   'YPOS': [Configuration['BEAM'][0]/3600.,'a'],
                                   'SBR': [sbr_standard,'a'],
                                   'Z0': [Configuration['BEAM'][0]*2.,'a'],
                                   'SDIS': [Configuration['CHANNEL_WIDTH']*1.5,'a'],
                                   }
    print_log(f'''GET_FIT_GROUPS: We have found the following unformatted groups from VARY:
{'':8s}{groups}
''',Configuration,case=['debug_add'])
    for group in groups:
        print_log(f'''GET_FIT_GROUPS: We are processing {group}
''',Configuration,case=['debug_add'])
        parts = group.split()
        if parts[0][0] == '!':
            block.append(False)
            parts[0] = parts[0][1:]
        else:
            block.append(True)
        found_rings = []
        preliminary_group = []
        for part in parts:
            if part[0].isnumeric():
                if ':' in part:
                    in_rings = [int(x) for x in part.split(':')]
                    if in_rings.sort():
                        in_rings.sort()
                    in_rings[1] += 1
                    current = [int(x) for x in range(*in_rings)]
                else:
                    current = [int(part)]
                if len(found_rings) == 0:
                    found_rings = current
                else:
                    if np.array_equal(np.array(found_rings),np.array(found_rings)):
                        pass
                    else:
                        raise DefFileError("The VARY settings in this deffile are not acceptable you have different ring for one block.")
            else:
                preliminary_group.append(f'{part}')
        if len(found_rings) == 1:
            block[-1] = False
        parameter_groups.append(preliminary_group)
        rings.append(found_rings)
        block_str = 'as a block.'
        if not block[-1]:
            block_str = 'individually.'
        print_log(f'''GET_FIT_GROUPS: We determined the group {parameter_groups[-1]}
{'':8s}with rings {rings[-1]}
{'':8s} which are varied {block_str}
''',Configuration,case=['debug_add'])
        # Then get current errors that are present
        per_par_variation  = []
        current_par = ''
        for par in parameter_groups[-1]:
            if current_par == '':
                current_par = par
                if current_par[-1] == '2':
                    current_par=current_par[:-2]
            par = [f'# {par}_ERR']
            all_errors = load_tirific(Configuration,Tirific_Template, par,\
                array=True )
            current_rings = np.array(rings[-1],dtype=int)-1
            if current_rings.size == 1:
                current_rings = int(current_rings)
            if len(all_errors) ==0. and current_par != 'SBR':
                if block[-1]:
                    per_par_variation.append(paramater_standard_variation[current_par][0]/3.)
                else:
                    per_par_variation.append(paramater_standard_variation[current_par][0])
            elif current_par == 'SBR':
                print_log(f'''GET_FIT_GROUPS: For SBR we are setting
{'':8s}the cut_off_limits {cut_off_limits}
{'':8s} in ring {current_rings}  == {cut_off_limits[current_rings]}
''',Configuration,case=['debug_add'])
                if block[-1]:
                    per_par_variation.append(np.mean(cut_off_limits[current_rings])*3.)
                else:
                    per_par_variation.append(np.max(cut_off_limits[current_rings])*3.)
            else:

                if block[-1]:
                    per_par_variation.append(np.mean(all_errors[current_rings])*3.)
                else:
                    per_par_variation.append(np.max(all_errors[current_rings]))
        if block[-1]:
            tmp_variation = np.mean(per_par_variation)*3.
        else:
            tmp_variation = np.max(per_par_variation)
        if tmp_variation > paramater_standard_variation[current_par][0]:
            variation.append(tmp_variation)
            variation_type.append('a')
        else:
            variation.append(paramater_standard_variation[current_par][0])
            variation_type.append(paramater_standard_variation[current_par][1])


    return parameter_groups,rings,block,variation,variation_type

get_fit_groups.__doc__ =f'''
 NAME:
    get_fit_groups
 PURPOSE:
    get the groups that are fitting, whether they are a block or not and their expected errors

 CATEGORY:
    support_functions

 INPUTS:
    Tirific_Template =  the def file to get errors.
 OPTIONAL INPUTS:

 KEYWORD PARAMETERS:

 OUTPUTS:
     tirshaker settings
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 EXAMPLE:

 NOTE:
'''

def get_inclination_pa(Configuration, Image, center, cutoff = 0.):
    map = copy.deepcopy(Image[0].data)
    for i in [0,1]:
        print_log(f'''GET_INCLINATION_PA: Doing iteration {i} in the estimates
''',Configuration,case=['debug_start'])
        # now we need to get profiles under many angles let's say 100
        #extract the profiles under a set of angles
        angles = np.linspace(0, 360, 180)
        ratios, maj_extent = obtain_ratios(Configuration,map, center, angles,noise = cutoff)
        if np.any(np.isnan(ratios)):
            return [float('NaN'),float('NaN')],  [float('NaN'),float('NaN')],float('NaN')
        sin_ratios,sin_parameters = fit_sine(Configuration,angles,ratios)
        if np.any(np.isnan(sin_parameters)):
            return [float('NaN'),float('NaN')],  [float('NaN'),float('NaN')],float('NaN')
        
        ratios=sin_ratios

        max_index = np.where(ratios == np.nanmax(ratios))[0]
        if max_index.size > 1:
            max_index =int(max_index[0])
        min_index = np.where(ratios == np.nanmin(ratios))[0]
        if min_index.size > 1:
            min_index =int(min_index[0])
        #From here we know it is 1 value so lets set it to an integer    
        max_index = set_limits(int(max_index),2,177)
        min_index = set_limits(int(min_index),2,177)

        if min_index > 135 and max_index < 45:
            min_index=min_index-90
        if max_index > 135 and min_index < 45:
            max_index=max_index-90
        if i == 0:
            print_log(f'''GET_INCLINATION_PA: We initially find these indeces min {min_index } {angles[min_index]} max {max_index} {angles[max_index]}.
''',Configuration,case=['debug_add'])
        else:
            print_log(f'''GET_INCLINATION_PA: From the cleaned map we find these indeces min {min_index }  {angles[min_index]} max {max_index} {angles[max_index]}.
''',Configuration,case=['debug_add'])
        #get a 10% bracket

        tenp_max_index = np.where(ratios > np.nanmax(ratios)*0.9)[0]
        tenp_min_index = np.where(ratios < np.nanmin(ratios)*1.1)[0]

        if tenp_max_index.size <= 1 and 2 <= max_index <=177 :
            tenp_max_index= [max_index-2,max_index+2]

        if tenp_min_index.size <= 1:
            tenp_min_index= [min_index-2,min_index+2]
       
        if angles[min_index] > angles[max_index]:
            pa = float(np.nanmean(np.array([angles[min_index]-90,angles[max_index]],dtype=float)))
        else:
            pa = float(np.nanmean(np.array([angles[min_index]+90,angles[max_index]],dtype=float)))
        if 180. < pa:
            pa = pa -180

        if i == 0:
            print_log(f'''GET_INCLINATION_PA: The initial  pa = {pa}.
''',Configuration,case=['debug_add'])
        else:
            print_log(f'''GET_INCLINATION_PA: The pa from the cleaned map  = {pa}.
''',Configuration,case=['debug_add'])
        pa_error = set_limits(np.nanmean([abs(float(angles[int(tenp_min_index[0])])-float(angles[min_index])),\
                            abs(float(angles[int(tenp_min_index[-1])])-float(angles[min_index])),\
                            abs(float(angles[int(tenp_max_index[0])])-float(angles[max_index])), \
                            abs(float(angles[int(tenp_min_index[-1])])-float(angles[max_index]))]), \
                            0.5,15.)
        ratios[ratios < 0.204] = 0.204
        ratios[1./ratios < 0.204] = 1./0.204
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            inclination = np.nanmean([np.degrees(np.arccos(np.sqrt((float(ratios[min_index])**2-0.2**2)/0.96))) \
                              ,np.degrees(np.arccos(np.sqrt(((1./float(ratios[max_index]))**2-0.2**2)/0.96))) ])
          

        if i == 0:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                inclination_error = set_limits(np.nanmean([abs(np.degrees(np.arccos(np.sqrt((ratios[min_index]**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt((ratios[tenp_min_index[0]]**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt((ratios[min_index]**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt((ratios[tenp_min_index[-1]]**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt(((1./ratios[max_index])**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt(((1./ratios[tenp_max_index[0]])**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt(((1./ratios[max_index])**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt(((1./ratios[tenp_max_index[0]])**2-0.2**2)/0.96))))]), \
                                     1.5,10.)
            if not np.isfinite(inclination_error):
                inclination_error = 90./inclination

        if ratios[max_index]-ratios[min_index] < 0.4:
            inclination = float(inclination-(inclination*0.01/(ratios[max_index]-ratios[min_index])))
            if i == 0:
                inclination_error = float(inclination_error*0.4/(ratios[max_index]-ratios[min_index]))
        if maj_extent*3600./Configuration['BEAM'][0] < 4:
            inclination = float(inclination+(inclination/10.*np.sqrt(4./(maj_extent/Configuration['BEAM'][0]*3600.))))
            if i == 0:
                inclination_error = float(inclination_error*4./(maj_extent*3600./Configuration['BEAM'][0]))

        if i == 0:
            print_log(f'''GET_INCLINATION_PA: The initial inclination = {inclination} +- {inclination_error}.
''',Configuration,case=['debug_add'])
        else:
            print_log(f'''GET_INCLINATION_PA: From the cleaned map we find inclination = {inclination}+- {inclination_error}.
''',Configuration,case=['debug_add'])
        # this leads to trouble for small sources due to uncertain PA and inclination estimates
        # at low inclinations the error is too uncertain to have a proper deprojection    
        inclination_error = inclination_error/np.sin(np.radians(inclination)) 
        if inclination < 70. and maj_extent*3600./Configuration['BEAM'][0] > 4 and inclination_error < 10.:
            Image_clean = remove_inhomogeneities(Configuration,Image,inclination=inclination, pa = pa,iteration=i, center = center,WCS_center = False )
            map = Image_clean[0].data
            Image_clean.close()
        else:
            break
    
    return [inclination,inclination_error], [pa,pa_error],maj_extent

get_inclination_pa.__doc__ =f'''
 NAME:
    get_inclination_pa
 PURPOSE:
    Get estimates for the PA, inclination and radius from a intensity map.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Image = the intensity map to be evaluated, this should be an astropy structure
    center = center of the galaxy in pixels

 OPTIONAL INPUTS:


    cutoff = 0.
    Limit to trust the map to, below this values are not evaluated

 OUTPUTS:
    inclination,inclination_error = inclination and error
    pa,pa_error =  pa and errror
    maj_extent = the estimated radius of the galaxy

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_new_center(Configuration, map_in, inclination=-1,pa=-1, noise= 0.):

    new_center = fit_center_ellipse(Configuration,map_in,inclination=inclination,pa=pa )


    return new_center,True,False


get_new_center.__doc__ =f'''
 NAME:
    get_new-center
 PURPOSE:
    Get the best fitting center for the galaxy.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map_in = the intensity map to be evaluated, this should be an astropy structure
    inclination = the current inclination for weighing purposes
    pa =  the current inclination for weighing purposes

 OPTIONAL INPUTS:



 OUTPUTS:
    center = the newly determined center
    center_stable = whether the center has shifted or now.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


# Function to get the amount of inner rings to fix
def get_inner_fix(Configuration,Tirific_Template ):
    print_log(f'''GET_INNER_FIX: Attempting to get the inner rings to be fixed.
''',Configuration,case=['debug_start'])
    sbr_av = np.array([(float(x)+float(y))/2. for x,y  in zip(Tirific_Template['SBR'].split(),Tirific_Template['SBR_2'].split())],dtype = float)
    column_levels = columndensity(Configuration,sbr_av, arcsquare = True )
    column_levels[0]= 1e21
    tmp = np.where(column_levels > 1e20)[0]
    if Configuration['OUTER_RINGS_DOUBLED']:
        inner_min =set_limits(Configuration['NO_RINGS']/3.,5.,11)
    else:
        inner_min =set_limits(Configuration['NO_RINGS']/5.,4.,Configuration['NO_RINGS']*0.7)
    inner_min = int(set_limits(np.floor(tmp[-1]/1.5-1), inner_min, Configuration['NO_RINGS']*0.9))


    tiltogram = make_tiltogram(Configuration,Tirific_Template )
    check_tiltogram(Configuration,tiltogram,inner_min=inner_min )
get_inner_fix.__doc__ =f'''
 NAME:
    get_inner_fix

 PURPOSE:
    Obtain the number of rings that should be fixed to a single value in the inner parts
    All ring > 1e20 column density should be fixed upt to a maximum of 90% and a minimu of 4 rings.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard Tirific template containing the SBR profiles

 OPTIONAL INPUTS:


 OUTPUTS:
    Update INNER_FIX in the configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    !!!!!!!!!!!!!This appears to currently not be working well.
'''

def get_kinematical_center(Configuration,map,angle,center = None ):
    if center is None:
        center = [len(map[0,:])/2.,len(map[:,0])/2.]
    if angle > 180.:
        angle= angle-180.

    if angle < 90:
        angle = angle+90.
    else:
        angle= angle-90.
    buffer = int(round(np.mean(Configuration['BEAM_IN_PIXELS'][:2])/3.))
    min_axis_prof, min_axis, min_res = get_profile(Configuration,map,angle,center= center )
    found_vsys = [np.nanmean(min_axis_prof),0,0]
    found_diff=  abs(np.nanmin(min_axis_prof)-np.nanmax(min_axis_prof))+ np.nansum([abs(x-found_vsys[0]) for x in min_axis_prof])
    for x in range(-buffer,buffer):
        for y in range(-buffer,buffer):
                var_center = [int(round(center[0]+x)),int(round(center[1]+y))]
                min_axis_prof, min_axis, min_res = get_profile(Configuration,map,angle,center= var_center )
                var_diff=  abs(np.nanmin(min_axis_prof)-np.nanmax(min_axis_prof))+np.nansum([abs(x-np.nanmean(min_axis_prof)) for x in min_axis_prof])
                if var_diff < found_diff:
                    found_diff = copy.deepcopy(var_diff)
                    found_vsys = [np.nanmean(min_axis_prof),x,y]
    return found_vsys

get_kinematical_center.__doc__=f'''
 NAME:
    get_kinematical_center

 PURPOSE:
    Determine at which position a profile has the least variation

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map = the 2D array with the map
    angle = the angle of the major axis, with error

 OPTIONAL INPUTS:

    center= [0.,0.]
    default is the cenetr of the map given in pixels

 OUTPUTS:
    profile, the corresponding axis in pixel size, resolution of a step on the axis

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def get_profile(Configuration,map,angle,center = None):

    if center is None:
        center = [len(map[0,:])/2.,len(map[:,0])/2.]
    x1,x2,y1,y2 = obtain_border_pix(Configuration,angle,center )
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    resolution = np.sqrt((x2-x1)**2+(y2-y1)**2)/1000.
    #maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(pa[0])))+abs(abs(y2-y1)/1000.*np.cos(np.radians(pa[0])))
    profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    axis =  np.linspace(0,1000*resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(angle)))+abs(abs(center[1])*np.cos(np.radians(angle))))
    return profile,axis,resolution
get_profile.__doc__=f'''
 NAME:
    get_profile

 PURPOSE:
    extract a profile under an arbitry angle from a map.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map = the 2D array with the map
    pa = the angle

 OPTIONAL INPUTS:

    center= [0.,0.]
    default is the cenetr of the map given in pixels

 OUTPUTS:
    profile, the corresponding axis in pixel size, resolution of a step on the axis

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def get_system_string(string):
    '''Escape any spaces in string with backlash'''
    if len(string.split()) > 1:
        string = "\ ".join(string.split())
    return string
get_system_string.__doc__=f'''
 NAME:
    get_system_string

 PURPOSE:
    Escape any spaces in string with backlash

 CATEGORY:
    support_functions

 INPUTS:
    string = is input string with spaces

 OPTIONAL INPUTS:


 OUTPUTS:
    string = string with escaped spaces

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def get_tirific_output_names(Configuration,work_dir,deffile,shaker=False):
    print_log(f'''GET_TIRIFIC_OUTPUT_NAMES: Starting the extraction of the output names in {deffile} tot be ran in {work_dir})
''',Configuration,case=['debug_start'])

    fitsfile,deffile = load_tirific(Configuration,f'{work_dir}/{deffile}' ,\
        Variables=['OUTSET','TIRDEF'])
    output_deffile = f'{work_dir}/{deffile[0]}'
    #Tirshakerdoes not produce the fits outset
    if not shaker:   
        output_fits = f'{work_dir}/{fitsfile[0]}'
    else:
        output_fits = f'{work_dir}/{deffile[0]}'
    return output_fits,output_deffile
    
get_tirific_output_names.__doc__=f'''
 NAME:
    get_tirific_output_names

 PURPOSE:
    get the  name of the output fits file and def file as defined in the running directory + deffile

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    work_dir = The directory where the def file is ran as tirific input can be relative
    deffile = the tirific input deffile


 OPTIONAL INPUTS:


 OUTPUTS:
    output_fits = the full output name
    output_deffile = the full output def file name

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_vel_pa(Configuration,velocity_field,center= [0.,0.] ):
    print_log(f'''GET_VEL_PA: This is the center we use {center}
''',Configuration,case=['debug_start'])
    if center == None:
        center = [x/2. for x in velocity_field.shape]
    else:
         # because python is stupid the ceneter should be reversed to match the map
        center.reverse()
    sigma = [3.,3.]
    sm_velocity_field = ndimage.gaussian_filter(velocity_field, sigma=(sigma[1], sigma[0]), order=0)
    while len(sm_velocity_field[~np.isnan(sm_velocity_field)]) < 10 and sigma[0] > 0.5:
        sigma = [x-0.5 for x in sigma]
        sm_velocity_field = ndimage.gaussian_filter(velocity_field, sigma=(sigma[1], sigma[0]), order=0)
    if len(sm_velocity_field[~np.isnan(sm_velocity_field)]) < 10:
        sm_velocity_field = copy.deepcopy(velocity_field)

    max_pos = np.where(np.nanmax(sm_velocity_field) == sm_velocity_field)
    #Python is a super weird language so make a decent list of np output

    max_pos = [float(max_pos[0]),float(max_pos[1])]
    min_pos = np.where(np.nanmin(sm_velocity_field) == sm_velocity_field)
    min_pos = [float(min_pos[0]),float(min_pos[1])]

    print_log(f'''GET_VEL_PA: This is the location of the maximum {max_pos} and minimum {min_pos}
''',Configuration,case=['debug_add'])
    try:
        pa_from_max = np.arctan((center[1]-max_pos[1])/(center[0]-max_pos[0]))
    except ZeroDivisionError:
        if center[1]-max_pos[1] >= 0.:
            pa_from_max = np.radians(90.)
        else:
            pa_from_max = np.radians(-90.)
    try:
        pa_from_min = np.arctan((center[1]-min_pos[1])/(center[0]-min_pos[0]))
    except ZeroDivisionError:
        if center[1]-min_pos[1] >= 0.:
            pa_from_min = np.radians(90.)
        else:
            pa_from_min = np.radians(-90.)
    try:
        pa_max_to_min = np.arctan((max_pos[1]-min_pos[1])/(max_pos[0]-min_pos[0]))
    except ZeroDivisionError:
        if max_pos[1]-min_pos[1] >= 0.:
            pa_max_to_min = np.radians(90.)
        else:
            pa_max_to_min = np.radians(-90.)

    for i in [0,1,2]:
        if i == 0:
            pa = pa_from_max
            pos1 = center
            pos2 = max_pos
        elif i == 1:
            pa_from_max = pa
            pa = pa_from_min
            pos1 = center
            pos2 = max_pos
        else:
            pa_from_min = pa
            pa = pa_max_to_min
            pos1 = min_pos
            pos2 = max_pos
        if pos1[1]-pos2[1] == 0:
            if pos1[0]-pos2[0] < 0.:
                pa = 0.0175
            else:
                pa = np.radians(180.)
        elif pos1[1]-pos2[1] < 0:
            if pa < 0.:
                pa = abs(pa)+np.radians(180.)
            else:
                pa = np.radians(360.) - pa
        elif pos1[1]-pos2[1] > 0:
            if pa < 0.:
                pa = abs(pa)
            else:
                pa = np.radians(180.) - pa
    print_log(f'''GET_VEL_PA: This is the PA
{'':8s} pa = {pa} pa_from_max = {pa_from_max} pa_from_min = {pa_from_min}
''',Configuration,case=['debug_add'])
    if np.degrees(abs(pa-pa_from_max)) > 170. or   np.degrees(abs(pa-pa_from_min)) > 170:
        pa = pa
    else:
        pa = np.nanmean([pa,pa_from_max,pa_from_min])
    center.reverse()
    print_log(f'''GET_VEL_PA: This is the PA we extract from the velpa {np.degrees(pa)}
''',Configuration,case=['debug_add'])
    return np.degrees([pa, np.nanstd([pa,pa_from_max,pa_from_min])])
get_vel_pa.__doc__ =f'''
 NAME:
    get_vel_pa

 PURPOSE:
    Routine to obtain the PA from a moment 1 map of a galaxy

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    velocity_field = the VF

 OPTIONAL INPUTS:

    center = center of the galaxy

 OUTPUTS:
    mean and standard deviation of the pa from minimum value to center, maximum to center and minimum to maximum

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def isiterable(variable):
    '''Check whether variable is iterable'''
    #First check it is not a string as those are iterable
    if isinstance(variable,str):
        return False
    try:
        iter(variable)
    except TypeError:
        return False

    return True
isiterable.__doc__ =f'''
 NAME:
    isiterable

 PURPOSE:
    Check whether variable is iterable

 CATEGORY:
    support_functions

 INPUTS:
    variable = variable to check

 OPTIONAL INPUTS:

 OUTPUTS:
    True if iterable False if not

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
#The functions load_tirifiic,load_template and get_from template were extremely similar
#This replaces all with a single function
def load_tirific(Configuration,def_input,Variables = None,array = False,\
        ensure_rings = False,brightness_check=False ):
    #Cause python is the dumbest and mutable objects in the FAT_defaults
    # such as lists transfer
    if Variables is None:
        Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI',\
                     'VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS',\
                     'SDIS','VROT_2',  'Z0_2','SBR_2','INCL_2','PA_2','XPOS_2',\
                     'YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']
    sbr_added =[False,False]
    if brightness_check:
        if 'SBR' not in Variables:
            Variables.append('SBR')
            sbr_added[0] =True
        if 'SBR_2' not in Variables:
            Variables.append('SBR_2')
            sbr_added[1] =True
    # if the input is a string we first load the template
    if isinstance(def_input,str):
        def_input = tirific_template(filename = def_input )

    out = []
    for key in Variables:
        try:
            out.append([float(x) for x  in def_input[key].split()])
        except KeyError:
            out.append([])
        except ValueError:
            out.append([x for x  in def_input[key].split()])

    #Because lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand

    if array:
        tmp = out
        #We can ensure that the output has the same number of values as there are rings
        if ensure_rings:
            length=int(def_input['NUR'])
        else:
            #or just take the longest input as the size
            length = max(map(len,out))
        #let's just order this in variable, values such that it unpacks properly into a list of variables
        out = np.zeros((len(Variables),length),dtype=float)
        for i,variable in enumerate(tmp):
            if len(variable) > 0.:
                out[i,0:len(variable)] = variable[0:len(variable)]
    if brightness_check:
        for variable in Variables:
            tmp = variable.split('_')
            ext=''
            if tmp[-1] == '2':
                ext = '_2'
            elif tmp[-1] == 'ERR':
                if tmp[-2] == 2:
                   ext = '_2'  
           
            cut = np.where(out[Variables.index(f'SBR{ext}'),:] < 1e-15)
            if np.sum(cut) > 0.:
                out[Variables.index(variable),cut] = float('NaN')
        

    if len(Variables) == 1:
        out= out[0]
    print_log(f'''LOAD_TIRIFIC: We extracted the following profiles from the Template.
{'':8s}Requested Variables = {Variables}
{'':8s}Extracted = {out}
''',Configuration,case=['debug_start'])
    #Beware that lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    # However if you make a np. array from it make sure that you specify float  or have lists of the same length else you get an array of lists which behave just as dumb
    return out
load_tirific.__doc__ =f'''
 NAME:
    load_tirific

 PURPOSE:
    Load values from variables set in the tirific files

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    def_input = Path to the tirific def file or a FAT tirific template dictionary

 OPTIONAL INPUTS:
    Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']


    array = False
        Specify that the output should be an numpy array with all varables having the same length

    ensure_rings =false
        Specify that the output array should have the length of the NUR parameter in the def file
 OUTPUTS:
    outputarray array with all the values of the parameters requested

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    This replaces load_tirific,load_template and get_from_template.
'''


def make_equal_length(list1_in,list2_in):
    '''ensure that 2 lists have the same length'''
    #python is a silly language
    list1=list1_in.copy()
    list2=list2_in.copy()

    #We can only do this with lists so first we ensure lists
    list1 = ensure_list(list1)
    list2 = ensure_list(list2)

    #We are matching the lists on the highest level only so if
    if isiterable(list1[0]) or isiterable(list2[0]):
        raise TypeError('make_equal_length only works on 1D arrays or scalars')

    while len(list1) < len(list2):
        list1.append(list1[-1])
    while len(list1) > len(list2):
        list2.append(list2[-1])
    return list1, list2
make_equal_length.__doc__ =f'''
 NAME:
    make_equal_length

 PURPOSE:
    ensure that 2 lists have the same length

 CATEGORY:
    support_functions

 INPUTS:
    list1, list2

 OPTIONAL INPUTS:


 OUTPUTS:
    the lists with the same length.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    Lists are mutable objects so if the input are lists the return could give the impression that
    the original lists remain unmodified. To ensure proper behaviour and not python magic the input
    is copied at the start.
'''


def make_tiltogram(Configuration,Tirific_Template):
    print_log(f'''MAKE_TILTOGRAM: Starting tiltogram.
''',Configuration,case=['debug_start'])
    pa_incl = load_tirific(Configuration,Tirific_Template,\
            Variables=['PA','PA_2','INCL','INCL_2'],array=True)

    add = [[],[]]
    Theta = [[],[]]
    phi = [[],[]]
    x = [[],[]]
    y = [[],[]]
    z = [[],[]]
    tiltogram = [[],[]]
    for i in [0,1]:
        add[i] = [0. if x < 90 else 90. if  90 <= x < 180. else 180. if 180<= x < 270 else 270. for x in pa_incl[i]]
        pa_incl[i] = pa_incl[i]-add[i]
        pa_incl[i] = np.array([x if x!= 0. else 0.00001 for x in pa_incl[i]],dtype=float)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            Theta[i] = np.arctan(np.tan(np.radians(pa_incl[i+2]))*np.tan(np.radians(pa_incl[i])))
            phi[i] = np.arctan(np.tan(np.radians(pa_incl[i]))/np.sin(Theta[i]))
            x[i]=np.sin(Theta[i])*np.cos(phi[i])
            y[i]=np.sin(Theta[i])*np.sin(phi[i])
            z[i]=np.cos(Theta[i])
            print_log(f'''MAKE_TILTOGRAM: For the cartesian coordinates we find in side {i}
{'':8s} x = {x[i]}
{'':8s} y = {y[i]}
{'':8s} z = {z[i]}
''',Configuration,case=['debug_add'])

            tiltogram[i] = np.degrees(np.arccos(np.multiply.outer(x[i], x[i]).ravel().reshape(len(x[i]),len(x[i])) + \
                                         np.multiply.outer(y[i], y[i]).ravel().reshape(len(x[i]),len(x[i])) + \
                                         np.multiply.outer(z[i], z[i]).ravel().reshape(len(x[i]),len(x[i]))))
    tiltogram = np.array(tiltogram,dtype = float)
    tiltogram[np.isnan(tiltogram)]= 0.

    return tiltogram

make_tiltogram.__doc__ =f'''
 NAME:
     make_tiltogram

 PURPOSE:
    Make the tiltogram of the current Template

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = standard FAT Configuration
    Tirific_Template = Standard Tirific template

 OPTIONAL INPUTS:


 OUTPUTS:
    a multidimensionals array containig the tiltograms for both sides

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:

'''



def max_profile_change(Configuration,radius,profile,key,max_change=None,\
                        slope = None, kpc_radius=False , quadrant = None):
    if not kpc_radius:
        radkpc = convertskyangle(Configuration,radius)
    else:
        radkpc = copy.deepcopy(radius)
    new_profile = copy.deepcopy(profile)
  
    sm_profile = savgol_filter(new_profile, 3, 1)
    if key in ['ARBITRARY','CHANGE_ANGLE'] and not max_change:
        raise InputError(f'For {key} checks you have to set the max_change.')
    if max_change == None:
        max_change=Configuration['MAX_CHANGE'][key]
    diff_rad =  [float(y-x) for x,y in zip(radkpc,radkpc[1:])]
    diff_profile = [float(y-x) for x,y in zip(profile,profile[1:])]
    diff_sm_profile = [float(y-x) for x,y in zip(sm_profile,sm_profile[1:])]
    if quadrant is None:
        diff_q = np.zeros(len(diff_profile))
    else:
        diff_q = [float(y-x) for x,y in zip(quadrant,quadrant[1:])]

    print_log(f'''MAX_CHANGE_PROFILE: The profile {key} starts with.
{'':8s} {key} = {new_profile} and max change = {max_change}
{'':8s} smoothed {key}  = {sm_profile}
{'':8s} diff per ring = {diff_profile}
{'':8s} smoothed dif per ring = {diff_sm_profile}
{'':8s} slope = {slope}
{'':8s} diff/kpc = {[x/y for x,y in zip(diff_profile,diff_rad)]}
''', Configuration,case=['debug_start'])
    #print(f''' For the PA and the inlcination we get a max shift per ring
    #PA = {[Configuration['MAX_CHANGE']['PA']*y for x,y in zip(diff_profile,diff_rad)]}   
    #INCL =  {[Configuration['MAX_CHANGE']['INCL']*y for x,y in zip(diff_profile,diff_rad)]}            
    #          ''')
    #    exit()
    for i,diff in enumerate(diff_profile):
        if abs(diff)/diff_rad[i] > max_change and diff_q[i] == 0.:
            print_log(f'''MAX_CHANGE_PROFILE: for ring {i} with {new_profile[i+1]} we find the difference {diff/diff_rad[i]}
''', Configuration,case=['debug_add'])
            #if it is the last point we simply limit it
            if i == len(diff_profile)-1:
                new_profile[i+1] = new_profile[i]+ diff/abs(diff)*(max_change*0.5*diff_rad[i])
                #print_log(f'''MAX_CHANGE_PROFILE: This is the last ring putting {new_profile[i+1] }.
#''', Configuration,case=['debug_add'])
            elif i+1 == slope:
                # if we have the start of the slope on the max_change then put it to value of the previous ring
                new_profile[i+1] = new_profile[i]
                print_log(f'''MAX_CHANGE_PROFILE: This is the start of the slope putting {new_profile[i+1] }.
''', Configuration,case=['debug_add'])
            elif abs(diff_sm_profile[i])/diff_rad[i] < max_change*0.5:
                new_profile[i+1] = sm_profile[i+1]
                print_log(f'''MAX_CHANGE_PROFILE: The smoothed profile is acceptable using {new_profile[i+1] }.
''', Configuration,case=['debug_add'])
                sm_profile = savgol_filter(new_profile, 3, 1)
                diff_sm_profile = [float(y-x) for x,y in zip(sm_profile,sm_profile[1:])]
            elif diff_profile[i+1] == 0:
                new_profile[i+1] = profile[i]+ diff/abs(diff)*(max_change*0.9*diff_rad[i]) 
                print_log(f'''MAX_CHANGE_PROFILE: The next ring is the same using {new_profile[i+1] }.
''', Configuration,case=['debug_add'])
                #If the change does not reverse we simply limit
            elif diff/abs(diff) == diff_profile[i+1]/abs(diff_profile[i+1]):
                new_profile[i+1] = new_profile[i]+ diff/abs(diff)*(max_change*0.9*diff_rad[i])
                print_log(f'''MAX_CHANGE_PROFILE: We are simply applying the maximum change using {new_profile[i+1] }.
''', Configuration,case=['debug_add'])
            else:
                if abs(diff) > abs(diff_profile[i+1]):
                    gapped_diff = radkpc[i+2] - radkpc[i]
                    if abs(new_profile[i]-new_profile[i+2])/gapped_diff < max_change:
                        new_profile[i+1] = np.mean([new_profile[i],new_profile[i+2]])
                    else:
                        new_profile[i+1] = new_profile[i]
                    print_log(f'''MAX_CHANGE_PROFILE: Bridged the gap using {new_profile[i+1] }.
''', Configuration,case=['debug_add'])
                else:

                    new_profile[i+1] = new_profile[i]+ diff/abs(diff)*(max_change*0.9*diff_rad[i])

            if i < len(diff_profile)-1:
                diff_profile[i+1] = float(new_profile[i+2]-new_profile[i+1])

    print_log(f'''MAX_CHANGE_PROFILE: The returned profile is:
{'':8s}{key} = {new_profile}
''', Configuration,case=['debug_add'])
    #for i,diff in enumerate([float(y-x) for x,y in zip(new_profile,new_profile[1:])]):
    #    print(f'The change {diff} allowed {max_change*diff_rad[i]}')
    #
    #exit()
    return new_profile
max_profile_change.__doc__ =f'''
 NAME:
     max_profile_change

 PURPOSE:
    Check that the profile is not change more than the maximum per kpc and correct if it does

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = standard FAT Configuration
    radius = radius corresponding to the profile
    profile = profile to examine
    key = the parameter to indicate the type of profile

 OPTIONAL INPUTS:


 OUTPUTS:
    new_profile
    the modified profile

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    Currently this function is only applied to the angular momentum vector

'''


def obtain_border_pix(Configuration,angle,center):
    rotate = False
    # only setup for 0-180 but 180.-360 is the same but -180
    if angle > 180.:
        angle -= 180.
        rotate = True

    if angle < 90.:
        x1 = center[0]-(Configuration['NAXES'][1]-center[1])*np.tan(np.radians(angle))
        x2 = center[0]+(center[1])*np.tan(np.radians(angle))
        if x1 < 0:
            x1 = 0
            y1 = center[1]+(center[0])*np.tan(np.radians(90-angle))
        else:
            y1 = Configuration['NAXES'][1]
        if x2 > Configuration['NAXES'][0]:
            x2 = Configuration['NAXES'][0]
            y2 = center[1]-(center[0])*np.tan(np.radians(90-angle))
        else:
            y2 = 0
    elif angle == 90:
        x1 = 0 ; y1 = center[1] ; x2 = Configuration['NAXES'][0]; y2 = center[1]
    else:
        x1 = center[0]-(center[1])*np.tan(np.radians(180.-angle))
        x2 = center[0]+(Configuration['NAXES'][1]-center[1])*np.tan(np.radians(180-angle))
        if x1 < 0:
            x1 = 0
            y1 = center[1]-(center[0])*np.tan(np.radians(angle-90))
        else:
            y1 = 0
        if x2 > Configuration['NAXES'][0]:
            x2 = Configuration['NAXES'][0]
            y2 = center[1]+(center[0])*np.tan(np.radians(angle-90))
        else:
            y2 = Configuration['NAXES'][1]
    # if the orginal angle was > 180 we need to give the line 180 deg rotation
    x = [x1,x2]
    y = [y1,y2]
    if rotate:
        x.reverse()
        y.reverse()
    return (*x,*y)

obtain_border_pix.__doc__ =f'''
 NAME:
    obtain_border_pix
 PURPOSE:
    Get the pixel locations of where a line across a map exits the map

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = standard FAT Configuration
    hdr = header of the map
    angle = the angle of the line running through the map
    center = center of the line running through the map

 OPTIONAL INPUTS:


 OUTPUTS:
    x,y
    pixel locations of how the line runs through the map

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    np.tan,np.radians,.reverse()

 NOTE:
'''

def obtain_ratios(Configuration, map, center, angles, noise = 0. ):
    ratios = []
    max_extent = 0.
    for angle in angles:
        #major axis
        maj_profile,maj_axis,maj_resolution = get_profile(Configuration,map,angle,center=center )
        tmp = np.where(maj_profile > noise)[0]
        if tmp.shape[0] == 0.:
             width_maj = 0.
        else:
            width_maj = (tmp[-1]-tmp[0])*maj_resolution
            if width_maj**2 > Configuration['BEAM_IN_PIXELS'][0]**2:
                width_maj = np.sqrt(width_maj**2 - Configuration['BEAM_IN_PIXELS'][0]**2)
            else:
                if width_maj < Configuration['BEAM_IN_PIXELS'][0]:
                    width_maj = Configuration['BEAM_IN_PIXELS'][0]

        if width_maj > max_extent:
            max_extent = width_maj
        #minor axis
        if angle < 90:
            min_profile,min_axis,min_resolution = get_profile(Configuration,map,angle+90,center=center )
        else:
            min_profile,min_axis,min_resolution = get_profile(Configuration,map,angle-90,center=center )
        tmp = np.where(min_profile > noise)[0]
        if tmp.shape[0] == 0.:
             width_min = 0.
        else:
            width_min = (tmp[-1]-tmp[0])*min_resolution
            if width_min**2 > Configuration['BEAM_IN_PIXELS'][0]**2:
                width_min = np.sqrt(width_min**2 - Configuration['BEAM_IN_PIXELS'][0]**2)
            else:
                if width_min < Configuration['BEAM_IN_PIXELS'][1]:
                    width_min = Configuration['BEAM_IN_PIXELS'][1]
            if width_min > max_extent:
                max_extent = width_min
        if width_min != 0. and width_maj != 0.:
            ratios.append(width_maj/width_min)
        else:
            ratios.append(float('NaN'))
    #as the extend is at 25% let's take 2 time the sigma of that
    #max_extent = (max_extent/(2.*np.sqrt(2*np.log(2))))*2.
    max_extent = max_extent/2.
    return np.array(ratios,dtype=float), max_extent*Configuration['PIXEL_SIZE']

obtain_ratios.__doc__ = '''
 NAME:
       obtain_ratios

 PURPOSE:
       Obtain the ratio between an axis along the specified angle as well as one rotated 90 degrees to determine the pa and inclination.
       Additionally keep track of the maximum width.

 CATEGORY:
       support_functions

 INPUTS:
       Configuration
       map = moment 0 map
       hdr = hdr of the map
       center = estimated center
       angles = angles to look for the ratios

 OPTIONAL INPUTS:
       noise = 0.
       noise level in the map



 OUTPUTS:
         ratios =  ratios corresponding to the input angles
         max_extent = the maximum extend of any profile analysed (in degree)

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      np.mean, ndimage.map_coordinates, np.linspace, np.vstack, np.cos, np.radians, np.sin

 NOTE:
'''

def remove_inhomogeneities(Configuration,fits_map_in,inclination=30., pa = 90.\
        ,center = None, WCS_center = True, iteration= 0 ):
    if center is None:
        center = [0.,0.]
    fits_map = copy.deepcopy(fits_map_in)
    print_log(f'''REMOVE_INHOMOGENEITIES: These are the values we get as input
{'':8s}Inclination = {inclination}
{'':8s}pa = {pa}
{'':8s}center = {center}, WCS = {WCS_center}
{'':8s}map shape = {np.shape(fits_map[0].data)}
''',Configuration,case=['debug_start'])
    map = fits_map[0].data
    # first rotate the pa
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        map_wcs = WCS(fits_map[0].header)
    # convert the boundaries to real coordinates
    if WCS_center:
        x, y = map_wcs.wcs_world2pix(*center, 1)
    else:
        x = center[0]
        y = center[1]
    rot_map = rotateImage(Configuration,map,pa-90,[x,y] )
    # deproject
    dep_map = deproject(Configuration,copy.deepcopy(rot_map),inclination,center = y )

    if Configuration['DEBUG']:
        fits.writeto(f"{Configuration['FITTING_DIR']}rot_map_{int(iteration)}.fits",rot_map,fits_map[0].header,overwrite = True)
        fits.writeto(f"{Configuration['FITTING_DIR']}dep_map_{int(iteration)}.fits",dep_map,fits_map[0].header,overwrite = True)

    angles = np.linspace(5.,360.,71)
    minimum_map = copy.deepcopy(dep_map)
    for angle in angles:
        rot_dep_map =  rotateImage(Configuration,copy.deepcopy(dep_map),angle,[x,y] )

        #tmp = np.where(rot_dep_map < minimum_map)[0]
        minimum_map[rot_dep_map < minimum_map] =rot_dep_map[rot_dep_map < minimum_map]
    clean_map = rotateImage(Configuration,deproject(Configuration,copy.deepcopy(minimum_map),inclination,center = y,invert= True ),-1*(pa-90),[x,y] )

    if Configuration['DEBUG']:
        fits.writeto(f"{Configuration['FITTING_DIR']}minimum_map_{int(iteration)}.fits",minimum_map,fits_map[0].header,overwrite = True)
        fits.writeto(f"{Configuration['FITTING_DIR']}clean_map_{int(iteration)}.fits",clean_map,fits_map[0].header,overwrite = True)
    fits_map[0].data = clean_map
    return fits_map

remove_inhomogeneities.__doc__ =f'''
 NAME:
remove_inhomogeneities

 PURPOSE:
    Remove the inhomogeneities from a galaxy by deprojecting the map and rotating the map in steps of 5 deg over 360 degrees.
    Only keep the minimum values from every angle to ensure obtaining the underlying disk only not the overdensities.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    fits_map = intensity map of the galaxy as an astropy structure, i.e fits_map[0].header and fits_map[0].data

 OPTIONAL INPUTS:


    inclination=30.
    inclination of the galaxy

    pa = 90.
    PA of the galaxy

    center = [0.,0.]
    center of the galaxy, used for rotation and deprojection

    WCS_center = True
    Boolean, if True the center is assumed to be in RA, DEC (degrees)
    If false the center is in pixels.

    iteration= 0.
    counter for outputting the intermediate maps (See optional outputs)

 OUTPUTS:
    The cleaned minimum value map

 OPTIONAL OUTPUTS:
    If debug is True all intemediate step map are written to fits file in the fitting directory.
    With rot_map_{{iteration}}.fits = the input map rotated such that the major axis as defined by pa is along the x-axis
    dep_map_{{iteration}}.fits = the deprojected original map
    minimum_map_{{iteration}}.fits = the deprojected map with only minimum values from each angle
    clean_map_{{iteration}}.fits = the final map reprojected and rotated back to the original angle

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def rename_fit_products(Configuration,stage = 'initial', fit_type='Undefined'):
    extensions = ['def','log','ps','fits']
    for filetype in extensions:
        source = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}")
        if filetype == 'log':
            if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                target = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev.{filetype}")
                os.system(f"cp {source} {target}")
        else:

            if filetype == 'def':
                if stage in ['run_ec','run_os','run_cc']:
                    Loopnr = f"Iteration_{Configuration['ITERATIONS']-1}"
                elif stage in ['after_cc','after_ec','after_os'] :
                    Loopnr = f"Iteration_{Configuration['ITERATIONS']}"
                elif fit_type == Configuration['USED_FITTING'] and stage in ['final_os']:
                    Loopnr = f"Smoothed_1"
                else:
                    Loopnr = 'Output_Before_'+stage
                if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                    target = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_{Loopnr}.{filetype}")
                    os.system(f"mv {source} {target}")

            elif os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                target = get_system_string(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev.{filetype}")
                os.system(f"mv {source} {target}")

rename_fit_products.__doc__ =f'''
 NAME:
    rename_fit_products
 PURPOSE:
    rename the tirific products from the previous stage.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    stage = 'initial'
    stage of the type of fitting

    fit_type='Undefined'
    Type of fitting

 OUTPUTS:
    renamed files
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#function to rotate an image without losing info
def rotateImage(Configuration,image, angle, pivot):
    padX = [int(image.shape[1] - pivot[0]), int(pivot[0])]
    padY = [int(image.shape[0] - pivot[1]), int(pivot[1])]
    imgP = np.pad(image, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, angle, axes=(1, 0), reshape=False)
    return imgR[padY[0]: -padY[1], padX[0]: -padX[1]]

rotateImage.__doc__ =f'''
 NAME:
    rotateImage

 PURPOSE:
    Rotate an image around a specified center

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    image = Image to rotate
    angle =  the angle to rotate the image by
    pivot = the center around which to rotate

 OPTIONAL INPUTS:


 OUTPUTS:
    The rotated image is returned

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def run_tirific(Configuration, current_run, stage = 'initial',\
        fit_type = 'Undefined',deffile='Undefined',max_ini_time = 600):
    print_log(f'''RUN_TIRIFIC: For the galaxy {Configuration['ID']} in directory {Configuration['SUB_DIR']} we are starting a RUN_TIRIFIC.
We are in in stage {stage} and fit_type {fit_type} and have done {Configuration['ITERATIONS']} loops
''',Configuration,case=['debug_start'])
    if deffile == 'Undefined':
        deffile=f"{fit_type}_In.def"

    # First move the previous fits
    rename_fit_products(Configuration,fit_type = fit_type, stage=stage )
    #If we are debug let's write the config file
    if Configuration['DEBUG']:
        if stage == 'run_os':
            write_config(
                f'{Configuration["LOG_DIRECTORY"]}CFG_Before_Iteration_{Configuration["ITERATIONS"]}.txt', Configuration)
    if fit_type == 'Error_Shaker':
        work_dir = os.getcwd()
        restart_file = f"restart_Error_Shaker.txt"
        shaker = True
    else:
        restart_file = f"{Configuration['LOG_DIRECTORY']}restart_{fit_type}.txt"
        work_dir = Configuration['FITTING_DIR']
        shaker = False
    #Get the output fits file and def file defined in workdir+ deffile
    output_fits,output_deffile = get_tirific_output_names(Configuration,work_dir,deffile\
                                ,shaker=shaker )
    # Then if already running change restart file
    if Configuration['TIRIFIC_RUNNING']:
        print_log(f'''RUN_TIRIFIC: We are using an initialized tirific in {Configuration['FITTING_DIR']} with the file {deffile}
''',Configuration)

        with open(restart_file,'a') as file:
            file.write("Restarting from previous run \n")
    else:
        print_log(f'''RUN_TIRIFIC: We are inizializing a new TiRiFiC in {Configuration['FITTING_DIR']} with the file {deffile}
''',Configuration)
        with open(restart_file,'w') as file:
            file.write("Initialized a new run \n")
        current_run = subprocess.Popen([Configuration['TIRIFIC'],f"DEFFILE={deffile}","ACTION= 1"],\
                                       stdout = subprocess.PIPE, stderr = subprocess.PIPE,\
                                       cwd=work_dir,universal_newlines = True)


        Configuration['TIRIFIC_RUNNING'] = True
        Configuration['TIRIFIC_PID'] = current_run.pid
    currentloop =1
    max_loop = 0
    counter = 0

    current_process= psu.Process(current_run.pid)

    initialized = datetime.now()
    if Configuration['TIMING']:
        time.sleep(0.1)
        with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
            file.write(f"# TIRIFIC: Initializing Tirific at stage = {fit_type}, Loop = {Configuration['ITERATIONS']} {datetime.now()} \n")
            CPU,mem = get_usage_statistics(Configuration,current_process)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
    else:
        time.sleep(0.1)

    if Configuration['VERBOSE_SCREEN']:
        print(f"\r{'':8s}RUN_TIRIFIC: 0 % Completed", end =" ",flush = True)
    triggered = False
    for tir_out_line in current_run.stdout:
        tmp = re.split(r"[/: ]+",tir_out_line.strip())
        counter += 1
        if (counter % 50) == 0:
            if Configuration['TIMING']:
                with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
                    if tmp[0] == 'L' and not triggered:
                        if tmp[1] == '1':
                            file.write(f"# TIRIFIC: Started the actual fitting {datetime.now()} \n")
                    CPU,mem = get_usage_statistics(Configuration,current_process)
                    file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
        if tmp[0] == 'L':
            if not triggered:
                triggered = True
            if int(tmp[1]) != currentloop and Configuration['VERBOSE_SCREEN']:
                print(f"\r{'':8s}RUN_TIRIFIC: {set_limits(float(tmp[1])-1.,0.,float(max_loop))/float(max_loop)*100.:.1f} % Completed", end =" ",flush = True)
            currentloop  = int(tmp[1])
            if max_loop == 0:
                max_loop = int(tmp[2])
            Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)

        if tmp[0].strip() == 'Finished':
            break
        if tmp[0].strip() == 'Abort':
            break
        if not triggered:
            #Check that the initialization doesn't take to long
            check = datetime.now()
            diff = (check-initialized).total_seconds()
            if diff > max_ini_time:
                print_log(f'''RUN_TIRIFIC: After {diff/60.} min we could not find the expected output from the tirific run. 
running in the directory = {work_dir} 
and the file deffile = {deffile}                         
''',Configuration,case=['main','screen'])
                raise TirificOutputError(f'''{diff/60.} minutes after initialization the fitting has still not started.
We were running {deffile} and failed to find the output {output_fits} or {output_deffile}.
''')
    if Configuration['VERBOSE_SCREEN']:
        print(f'\n')
    if Configuration['TIMING']:
        with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt",'a') as file:
            file.write(f"# TIRIFIC: Finished this run {datetime.now()} \n")
            CPU,mem = get_usage_statistics(Configuration,current_process)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb for TiRiFiC \n")
    if Configuration['VERBOSE_SCREEN']:
        print(f"{'':8s}RUN_TIRIFIC: Finished the current tirific run.")

    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    time.sleep(1.0)
    wait_counter = 0
    if fit_type != 'Error_Shaker':
        while not os.path.exists(output_fits) and wait_counter < 100.:
            time.sleep(0.3)
            wait_counter += 1
            if wait_counter/10. == int(wait_counter/10.):
                print(f"\r Waiting for {output_fits}. \n", end = "", flush = True)
                print_log(f'''RUN_TIRIFIC: we have waited {0.3*wait_counter} seconds for the output of tirific but it is not there yet.
''',Configuration)
        if not os.path.exists(output_fits) \
            or not  os.path.exists(output_deffile):
            print_log(f'''RUN_TIRIFIC: After 30 seconds we could not find the expected output from the tirific run. We are raising an error for this galaxy.
''',Configuration,case=['main','screen'])
            raise TirificOutputError(f'''The tirific subprocess did not produce the correct output, most likely it crashed.
We were running {deffile} and failed to find the output {output_fits} or {output_deffile}.
''')
    #If we are running tirific with an incorrect velocity type we should addapt the velocity type of the output
    if Configuration['HDR_VELOCITY'] != 'VELO':
        cube = fits.open(output_fits)
        cube[0].header['CTYPE3'] = Configuration['HDR_VELOCITY']
        fits.writeto(output_fits,cube[0].data,cube[0].header,overwrite=True)
    if currentloop != max_loop:
        return 1,current_run
    else:
        return 0,current_run
run_tirific.__doc__ =f'''
 NAME:
    run_tirific

 PURPOSE:
    Check whether we have an initialized tirific if not initialize and run else restart the initialized run.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    current_run = subprocess structure of tirific

 OPTIONAL INPUTS:


    fit_type = 'Undefined'
    type of fitting

    stage = 'initial'
    stage of the fitting process

    max_ini_time = 600
    maximum time it can take for tirific to initialize 
    Higher ini times take longer

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def set_boundaries(Configuration,key,lower,upper,input=False):
    print_log(f'''SET_BOUNDARIES: Starting to implement the boundaries for {key} with lower = {lower} and upper = {upper}.
''',Configuration,case=['debug_start'])

    check=False
    key = key.upper()
    if np.sum(Configuration[f'{key}_INPUT_BOUNDARY']) != 0.:
        check = True
    if input:
        boundary = f'{key}_INPUT_BOUNDARY'
    else:
        boundary = f'{key}_CURRENT_BOUNDARY'
    try:
        _ = (e for e in lower)
    except TypeError:
        lower= [lower]
    try:
        _ = (e for e in upper)
    except TypeError:
        upper= [upper]

    if len(lower) == 1:
        lower = lower*3
    if len(upper) == 1:
        upper = upper*3
    for i in range(len(Configuration[boundary])):
        if check:
            if lower[i] < Configuration[f'{key}_INPUT_BOUNDARY'][i][0]:
                lower[i] = Configuration[f'{key}_INPUT_BOUNDARY'][i][0]
            if upper[i] > Configuration[f'{key}_INPUT_BOUNDARY'][i][1]:
                upper[i] = Configuration[f'{key}_INPUT_BOUNDARY'][i][1]
        Configuration[boundary][i][0] = float(lower[i])
        Configuration[boundary][i][1] = float(upper[i])
set_boundaries.__doc__ =f'''
 NAME:
    set_boundaries

 PURPOSE:
    set the bopundaries for a current parameter

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    key = name of the parameter
    lower = the lower boundary
    upper = upper boundary

 OPTIONAL INPUTS:

    input = False
    if set to true the user input parameters are updated. this should only be used
    when reading the cube to set the input VSYS, XPOS and YPOS limits.

 OUTPUTS:
    updated boundary limits. When user input is set these can never be loewer/higher than those limits.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
#Set the indivdual configuration such that the loop can run without the full catalogue
def set_individual_configuration(current_galaxy_index,Full_Catalogue,\
        Original_Configuration):
    Configuration = copy.deepcopy(Original_Configuration)
    Configuration['ID'] = Full_Catalogue['ID'][current_galaxy_index]
    if Full_Catalogue['DISTANCE'][current_galaxy_index] != -1.:
        Configuration['DISTANCE'] = Full_Catalogue['DISTANCE'][current_galaxy_index]
    Configuration['SUB_DIR'] = Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]
    Configuration['BASE_NAME'] = Full_Catalogue['CUBENAME'][current_galaxy_index]+'_FAT'
    if not Configuration['SOFIA_BASENAME']:
        if 'BASENAME' in Full_Catalogue['ENTRIES']:
            Configuration['SOFIA_BASENAME'] = Full_Catalogue['BASENAME'][current_galaxy_index]
        else:
            Configuration['SOFIA_BASENAME'] = Configuration['BASE_NAME']
    #Add our fitting directory to the Configuration
    #Maindir always ends in slash already
    if Full_Catalogue['DIRECTORYNAME'][current_galaxy_index] == './':
        Configuration['FITTING_DIR'] = f"{Configuration['MAIN_DIRECTORY']}"
    else:
        Configuration['FITTING_DIR'] = f"{Configuration['MAIN_DIRECTORY']}{Full_Catalogue['DIRECTORYNAME'][current_galaxy_index]}/"
    if 'sofia_catalogue' in Configuration['FITTING_STAGES']:
        Configuration['INPUT_CUBE']= f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}_FAT.fits"
    else:        
        Configuration['INPUT_CUBE']= f"{Full_Catalogue['CUBENAME'][current_galaxy_index]}.fits"
    return(Configuration)
set_individual_configuration.__doc__ =f'''
 NAME:
    set_individual_configuration

 PURPOSE:
    Fill the parameters of the configuration file that come from the catalogue

 CATEGORY:
    support_functions

 INPUTS:
    current_galaxy_index = index in the full catalogue

    Full_Catalogue = the full catalogue

    Original_Configuration = The original Configuration that applies to each galaxy


 OPTIONAL INPUTS:

 OUTPUTS:
    A configuration for an individual galaxy

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:Configf
'''

#Function to read FAT configuration file into a dictionary
def setup_configuration(cfg):
    if cfg.installation_check:
        cfg.fitting.fixed_parameters=['INCL','PA','SDIS']
        cfg.advanced.max_iterations= 1
        cfg.advanced.loops= 1
        cfg.advanced.shaker_iterations = 2
        cfg.input.main_directory= f'{cfg.input.main_directory}/FAT_Installation_Check/'
        cfg.output.output_quantity = 0
        cfg.fitting.fitting_stages = ['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC','Tirshaker']
        #cfg.fitting.fitting_stages = ['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC']
        cfg.cube_name = 'NGC_2903.fits'
        #cfg.output.debug = False
        cfg.output.timing = True
        test_files = ['NGC_2903.fits','ModelInput.def']
        if not os.path.isdir(cfg.input.main_directory):
            os.mkdir(cfg.input.main_directory)
        else:
            for file in test_files:
                try:
                    os.remove(f'{cfg.input.main_directory}/{file}')
                except FileNotFoundError:
                    pass

        my_resources = import_pack_files('pyFAT_astro.Installation_Check')
        for file in test_files:
            data = (my_resources / file).read_bytes()
            with open(cfg.input.main_directory+file,'w+b') as tmp:
                tmp.write(data)
    if cfg.cube_name:
        cfg.input.catalogue = None
        cfg.advanced.multiprocessing = False
        cfg.fitting.catalogue_start_id= '-1'
        cfg.fitting.catalogue_end_id= '-1'

    Configuration = Proper_Dictionary({})

    for key in cfg._content:
        input_key = getattr(cfg,key)

        if isinstance(input_key,omegaconf.dictconfig.DictConfig):
            for sub_key in input_key._content:
                if str(key) == 'output' and str(sub_key) == 'catalogue':
                    Configuration['OUTPUT_CATALOGUE'] =  getattr(input_key,sub_key)
                elif str(key) == 'fitting' and str(sub_key) == 'fixed_parameters':
                    value = getattr(input_key,sub_key)
                    for req_key in ['XPOS','YPOS','VSYS']:
                        if req_key not in value:
                            value.append(req_key)
                    Configuration[str(sub_key).upper()] =  [value,value]
                elif isinstance(getattr(input_key,sub_key),omegaconf.listconfig.ListConfig):
                    Configuration[str(sub_key).upper()] = list(getattr(input_key,sub_key))
                    for i,oyt in enumerate(Configuration[str(sub_key).upper()]):
                        if isinstance(oyt,omegaconf.listconfig.ListConfig):
                            Configuration[str(sub_key).upper()][i] = list(oyt)
                else:
                    Configuration[str(sub_key).upper()] =  getattr(input_key,sub_key)
        else:
            Configuration[str(key).upper()] = input_key


    # None cfg additions, that is additions that should be reset for every galaxy

    boolean_keys = ['OPTIMIZED', # Are we fitting an optimized cube
                'TIRIFIC_RUNNING', # Is there a tirific initialized
                'OUTER_RINGS_DOUBLED', #Do the outer rings have twice the size of the inner rings
                'NEW_RING_SIZE',  #Have we update the size of the ring while not yet updating the template
                'VEL_SMOOTH_EXTENDED', # Is the velocity smoothing extended ????
                'EXCLUDE_CENTRAL', # Do we exclude the central part of the fitting due to blanks/an absorption source
                'ACCEPTED',
                'FIX_RING_SIZE', #Do we allow FAT to vary the rings size (False) or not (True)
                'ACCEPTED_TIRIFIC', #Did Tirific run the full loops (False) or not (True)
                'SOFIA_RAN', #Check if we have ran Sofia
                'NO_RADEC',
                'CENTRAL_CONVERGENCE' #Have we found the center
                ]
#
    for key in boolean_keys:
        Configuration[key] = False
    if Configuration['RING_SIZE'] < 0.:
        Configuration['FIX_RING_SIZE'] = True
        Configuration['RING_SIZE'] = abs(Configuration['RING_SIZE'])
    Configuration['FIX_SIZE'] = [False,False]# If we have before fitted a size we want to fix to this size to avoid looping
    other_keys={'ID': 'Unset', # ID of the galaxy in the catalogue , set from the catalogue at start of loop
               'SUB_DIR': 'Unset', # Name of the directory in which galaxy resides, set from the catalogue at start of loop
               'FITTING_DIR': 'Unset', # Full path of the directory in which the fitting takes place, set at start of loop
               'BASE_NAME': 'Unset', #Basename for FAT products, typically {input_cube}_FAT, set at start of loop
               'LOG_DIR': 'Unset', #Directory to put log files from run, set at start of loop
               'INPUT_CUBE': 'Unset', # Name of the input cube as listed in full catalogue
               'STOP_INDIVIDUAL_ERRORS': ['SmallSourceError','BadSourceError'\
                                        ,'FaintSourceError','BadHeaderError',\
                                        'BadCubeError','BadMaskError',\
                                        'BadCatalogueError'],
               'MAXIMUM_DIRECTORY_LENGTH': len('Directory Name'),
               'FULL_TIME': [datetime.now(),'Not completed'],
               'INITIALIZATION_TIME': ['Not started','Not completed'],
               'PREPARATION_TIME': ['Not started','Not completed'],
               'SOFIA_TIME': ['Not started','Not completed'],
               'INITIAL_GUESSES_TIME': ['Not started','Not completed'],
               'FIT_TIME': ['Not started','Not completed'],
               'TIRSHAKER_TIME': ['Not started','Not completed'],
               'OUTPUTLOG': None,
               'RUN_COUNTER': 0,
               'CENTRAL_CONVERGENCE_COUNTER': 1.,
               'RP_COUNTER': 1,
               'ITERATIONS': 0,
               'CURRENT_STAGE': 'initial', #Current stage of the fitting process, set at switiching stages
               'USED_FITTING': None,
               'TIRIFIC_PID': 'Not Initialized', #Process ID of tirific that is running
               'FAT_PID': os.getpid(), #Process ID of FAT that is running
               'FAT_PSUPROCESS': 'cant copy',
               'FINAL_COMMENT': "This fitting stopped with an unregistered exit.",

               'MAX_SIZE_IN_BEAMS': 30, # The galaxy is not allowed to extend beyond this number of beams in radius, set in check_source
               'MIN_SIZE_IN_BEAMS': 0., # Minimum allowed radius in number of beams of the galaxy, set in check_source
               'SIZE_IN_BEAMS': np.full(2,0.,dtype=float),  #The radius of the galaxy in number of beams, adapted after running Sofia
               'NO_RINGS': 0., # The number of rings in the fit,
               'LAST_RELIABLE_RINGS': [0.,0.], # Location of the rings where the SBR drops below the cutoff limits, adapted after every run. Should only be modified in check_size
               'LIMIT_MODIFIER': [1.], #Modifier for the cutoff limits based on the inclination , adapted after every run.
               'OLD_SIZE': [[0.,0.]], # List to keep track of the size in beams.

               'NO_POINTSOURCES': 0. , # Number of point sources, set in run_tirific

               'INNER_FIX': [4,4], #Number of rings that are fixed in the inner part for the INCL and PA, , adapted after every run in get_inner_fix in support_functions and for both sides
               'CENTRAL_FIX': [], #Parameters of the center to fix. This set in check_central_convergence
               'WARP_SLOPE': [None, None], #Ring numbers from which outwards the warping should be fitted as a slope, set in get_warp_slope in modify_template
               'OUTER_SLOPE_START': 1, # Ring number from where the RC is fitted as a slope
               'RC_UNRELIABLE': 1, # Ring number from where the RC values are set flat. Should only be set in check_size

               'NOISE': 0. , #Noise of the input cube in Jy/beam, set in read_cube
               'SNR': 0, # the mean SNR in the moment 0 map sofia mask these are set in guess_orientation
               'BEAM_IN_PIXELS': [0.,0.,0.], #FWHM BMAJ, BMIN in pixels and total number of pixels in beam area, set in read_cube in read_functions
               'BEAM': [0.,0.,0.], #  FWHM BMAJ, BMIN in arcsec and BPA, set in main
               'BEAM_AREA': 0., #BEAM_AREA in arcsec set in main
               'HDR_VELOCITY': 'VELO', #Track the velocity frame as tirific will only run with VELO
               'TIR_RUN_CUBE': None, # if the type of velocity axis is not VELO we need to run on a cube with that key word adjusted.
               'NAXES': [0.,0.,0.], #  Size of the cube in pixels x,y,z arranged like sane people not python, set in main
               'MAX_ERROR': {}, #The maximum allowed errors for the parameters, set in main derived from cube
               'MIN_ERROR': {}, #The minumum allowed errors for the parameters, initially set in check_source but can be modified through out INCL,PA,SDSIS,Z0 errors change when the parameters is fixed or release
               'MAX_CHANGE': {'INCL': 6.25, 'PA': 25.}, #The maximum change in a parameter in unit/kpc
               'CHANNEL_WIDTH': 0., #Width of the channel in the cube in km/s, set in main derived from cube
               'PIXEL_SIZE': 0., #'Size of the pixels in degree'
               }
    #####!!!!!!!!!!!!!!!!!!!!!!!!!!!! The use of min_error is not full implemented yet!!!!!!!!!!!!!!!!!!!!!!!!!
    for key in other_keys:
        Configuration[key] = other_keys[key]

    #If the input boundaries are unset we will set some reasonable limits


# The parameters that need boundary limits are set here
    boundary_limit_keys = ['PA','INCL', 'SDIS', 'Z0','VSYS','XPOS','YPOS','VROT','SBR']
    for key in boundary_limit_keys:
        if np.sum(Configuration[f"{key}_INPUT_BOUNDARY"]) == 0.:
            if key == 'INCL':
                Configuration[f"{key}_INPUT_BOUNDARY"] = [[5.,90.] for x in range(3)]
            elif key == 'PA':
                Configuration[f"{key}_INPUT_BOUNDARY"] = [[-10.,370.] for x in range(3)]
            else:
                #These are set when we have read the cube as they depend on it or the distance (Z0)
                pass

        Configuration[f"{key}_CURRENT_BOUNDARY"] = [[0.,0.],[0.,0.],[0.,0.]]

    #Make the input idiot safe
    if Configuration['MAIN_DIRECTORY'][-1] != '/':
        Configuration['MAIN_DIRECTORY'] = f"{Configuration['MAIN_DIRECTORY']}/"

    while not os.path.isdir(Configuration['MAIN_DIRECTORY']):
        Configuration['MAIN_DIRECTORY'] = input(f'''
Your main fitting directory ({Configuration['MAIN_DIRECTORY']}) does not exist.
Please provide the correct directory.
:  ''')


    if Configuration['CATALOGUE']:
        if not os.path.exists(Configuration['CATALOGUE']) and len(Configuration['CATALOGUE'].split('/')) == 1:
            Configuration['CATALOGUE']= f'''{Configuration['MAIN_DIRECTORY']}{Configuration['CATALOGUE']}'''

        while not os.path.exists(Configuration['CATALOGUE']):
            Configuration['CATALOGUE'] = input(f'''
Your input catalogue ({Configuration['CATALOGUE']}) does not exist.
Please provide the correct path and file name.
: ''')
    #Make sure there is only one Fit_ stage

    # Make sure all selected stages exist
    possible_stages = ['Fit_Tirific_OSC','Create_FAT_Cube',\
        'Run_Sofia','Restart_Fitting','External_Sofia','Sofia_Catalogue',\
        'Tirshaker','Fit_Make_Your_Own']
    possible_stages_l = [x.lower() for x in possible_stages]
    approved_stages = []
    fit_count = 0
    for stage in Configuration['FITTING_STAGES']:
        while stage.lower() not in possible_stages_l:
            stage = input(f'''
The stage {stage} is not supported by FAT.
Please pick one of the following {', '.join(possible_stages)}.
: ''')

        if 'fit_' in stage.lower():
            fit_count += 1
            if fit_count == 1:
                Configuration['USED_FITTING'] = possible_stages[possible_stages_l.index(stage.lower())]


        if fit_count > 1 and 'fit_' in stage.lower():
            print(f''' FAT only supports one single fitting stage. You have already selected one and hence {stage} will be ignored.
    ''')
        else:
            approved_stages.append(stage.lower())
    if Configuration['DEBUG']:
        Configuration['DEBUG_FUNCTION'] = [x.upper() for x in \
                                        Configuration['DEBUG_FUNCTION']]

    Configuration['FITTING_STAGES'] =approved_stages
    if 'sofia_catalogue' in Configuration['FITTING_STAGES'] and 'external_sofia' in Configuration['FITTING_STAGES']:
        Configuration['FITTING_STAGES'].remove('external_sofia')
    if ('sofia_catalogue' in Configuration['FITTING_STAGES'] or 'external_sofia' in Configuration['FITTING_STAGES']) and 'run_sofia' in Configuration['FITTING_STAGES']:
        Configuration['FITTING_STAGES'].remove('run_sofia')
    #if 'sofia_catalogue' in Configuration['FITTING_STAGES'] and 'create_fat_cube' in Configuration['FITTING_STAGES']:
    #    Configuration['FITTING_STAGES'].remove('create_fat_cube')
    if 'sofia_catalogue' in Configuration['FITTING_STAGES']:
        if Configuration['SOFIA_DIR'][-1] != '/':
            Configuration['SOFIA_DIR'] = f"{Configuration['SOFIA_DIR']}/"

  


    if 'run_sofia' in Configuration['FITTING_STAGES'] and 'external_sofia' in Configuration['FITTING_STAGES']:
        raise InputError(f"You are both providing existing sofia input and ask for sofia to be ran. This won't work exiting.")
    if 'fit_tirific_osc' in Configuration['FITTING_STAGES'] and \
        Configuration['NUMBER_OF_DISKS'] > 2:
        raise InputError(f"Fit_Tirific_OSC only allows for 2 or 1 disk. This won't work exiting.")

    #Check that the channel_dependency is acceptable
    while Configuration['CHANNEL_DEPENDENCY'].lower() not in ['sinusoidal','independent','hanning']:
        Configuration['CHANNEL_DEPENDENCY'] = input(f'''
The channel dependency  {'CHANNEL_DEPENDENCY'} is not supported by FAT.
Please pick one of the following {', '.join(['sinusoidal','independent','hanning'])}.
: ''')


    #The output catalogue only needs to be in a valid directory as we create it
    if Configuration['OUTPUT_CATALOGUE']:
        output_catalogue_dir = Configuration['OUTPUT_CATALOGUE'].split('/')
        if len(output_catalogue_dir) > 1:
            check_dir = '/'.join(output_catalogue_dir[:-1])
            while not os.path.isdir(check_dir):
                check_dir= input(f'''
                        The directory for your output catalogue ({Configuration['OUTPUT_CATALOGUE']}) does not exist.
                        Please provide the correct directory name.
                        ''')
                Configuration['OUTPUT_CATALOGUE'] = f"{check_dir}/{output_catalogue_dir[-1]}"
        elif len(output_catalogue_dir) == 1:
            Configuration['OUTPUT_CATALOGUE'] = f"{Configuration['MAIN_DIRECTORY']}{Configuration['OUTPUT_CATALOGUE']}"

    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['OUTPUT_QUANTITY'] == 5:
        Configuration['OUTPUT_QUANTITY'] = 4
    if Configuration['SHAKER_ITERATIONS'] < 2:
        Configuration['SHAKER_ITERATIONS'] = 2
    if Configuration['NCPU'] == 1:
        if Configuration['MULTIPROCESSING']:
             print(f'''You are only requesting 1 CPU to be used, we turn off multiprocessing.
''')
        Configuration['MULTIPROCESSING'] = False
    

    return Configuration
setup_configuration.__doc__ =f'''
 NAME:
    setup_configuration

 PURPOSE:
    Read the FAT config file and write into a dictionary

 CATEGORY:
    support_functions

 INPUTS:
    cfg = OmegaConf input object

 OPTIONAL INPUTS:

    fit_type = 'Undefined'

 OUTPUTS:
    Configuration = dictionary with the config file input

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Only this function can add entries to Original_Configuration or Configuration
'''



def sbr_limits(Configuration, Tirific_Template ):
    radii = load_tirific(Configuration,Tirific_Template,Variables=['RADI'],array=True)
    if len(radii) == 0.:
        radii = set_rings(Configuration)
        Tirific_Template['RADI']= ' '.join([str(x) for x in radii])
    print_log(f'''SBR_LIMITS: Got {len(radii)} rings at the start
''',Configuration, case=['debug_start'])
    level = Configuration['NOISE']*1000
    noise_in_column = columndensity(Configuration,level,systemic = \
        float(load_tirific(Configuration,Tirific_Template,Variables=['VSYS'] )[0]))
    J2007col=9.61097e+19
    #J2007scl= 2. #arcsec in a sech^2 layer
    ratio=(noise_in_column/J2007col)**0.5
    #*np.sqrt(J2007scl/scaleheight)
    beamsolid=(np.pi*Configuration['BEAM'][0]*Configuration['BEAM'][1])/(4.*np.log(2.))
    ringarea= [0 if radii[0] == 0 else np.pi*((radii[0]+radii[1])/2.)**2]
    ringarea = np.hstack((ringarea,
                         [np.pi*(((y+z)/2.)**2-((y+x)/2.)**2) for x,y,z in zip(radii,radii[1:],radii[2:])],
                         [np.pi*((radii[-1]+0.5*(radii[-1]-radii[-2]))**2-((radii[-1]+radii[-2])/2.)**2)]
                         ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sbr_ring_limits=9e-4*(ringarea/beamsolid)**(-0.82)*ratio
    if ringarea[0] == 0.:
         sbr_ring_limits[0]=np.nanmin(sbr_ring_limits)
         sbr_ring_limits[1]=sbr_ring_limits[2]/2.
    print_log(f'''SBR_LIMITS: Applying the modifiier:
{'':8s}{Configuration['LIMIT_MODIFIER']}
''',Configuration,case=['debug_add'])
    if len(Configuration['LIMIT_MODIFIER']) == 1:
        sbr_ring_limits= sbr_ring_limits*float(Configuration['LIMIT_MODIFIER'][0])
    else:
        mod_list = list(Configuration['LIMIT_MODIFIER'])
        while len(mod_list) < len(sbr_ring_limits):
            mod_list.append(Configuration['LIMIT_MODIFIER'][-1])
        Configuration['LIMIT_MODIFIER'] = np.array(mod_list,dtype=float)
        sbr_ring_limits=[x*y for x,y in zip(sbr_ring_limits,Configuration['LIMIT_MODIFIER'])]
    print_log(f'''SBR_LIMITS: Retrieved these radii and limits:
{'':8s}{radii}
{'':8s}{sbr_ring_limits}
''',Configuration,case=['debug_add'])
    return sbr_ring_limits

sbr_limits.__doc__ =f'''
 NAME:
    sbr_limits

 PURPOSE:
    Calculate the sbr limits below which rings can not be trusted.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard Tirific template


 OPTIONAL INPUTS:


 OUTPUTS:
    radii = the radii of the rings in arcsec
    sbr_ring_limits = the limits of reliability for each ring based on the noise in the cube in Jy/arcsec^2 x km/s

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_format(key):
    if key in ['SBR', 'SBR_2']:
        format = '.5e'
    elif key in ['XPOS', 'YPOS','XPOS_2', 'YPOS_2']:
        format = '.7e'
    else:
        format = '.2f'
    return format

set_format.__doc__ =f'''
 NAME:
    set_format

 PURPOSE:
    Get the format code for specific tirific parameter

 CATEGORY:
    support_functions

 INPUTS:
    key = Tirific parameter to get code for

 OPTIONAL INPUTS:

 OUTPUTS:
    The format code

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#simple function keep track of how to modify the edge limits
def set_limit_modifier(Configuration,Tirific_Template):
    '''Write to the Configuration a value or list of values to modify sbr dependent values.'''
    Profiles = load_tirific(Configuration,Tirific_Template,
        Variables=['INCL','INCL_2','Z0','Z0_2'],array=True )
    Inclination = [np.mean([x,y]) for x,y in zip(Profiles[0],Profiles[1])]
    Z0_av = [np.mean([x,y]) for x,y in zip(Profiles[2],Profiles[3])]
    kpc_radius=convertskyangle(Configuration,load_tirific(Configuration,\
        Tirific_Template,['RADI'] ))
    modifier_list = []
    # This should be per kpc
    if len(Inclination) > 1:
        if np.abs(Inclination[-1]-Inclination[-2])/(kpc_radius[-1]-kpc_radius[-2]) > Configuration['MAX_CHANGE']['INCL']:
            Inclination[-1] = Inclination[-2]
    # Correction because the noise applies to non-face on rings while the SBR is face on,correction is normalized to the average inclination
    # The square root is because of the square root between the noise in cube and the noise in J2007
    # The lower limit corresponds to a inclination 80 above which the cos becomes too steep
    for i,inc in enumerate(Inclination):
        if inc > 90.:
            inc  = 180.-inc
        if inc < 0.:
            inc= abs(inc)
        if i < 10.:
            modifier_list.append(set_limits(np.sqrt(np.cos(np.radians(inc))/np.cos(np.radians(60.))),0.75,2.))
        else:
            modifier_list.append(set_limits(np.sqrt(np.cos(np.radians(inc))/np.cos(np.radians(60.)))*i/10.,0.8,1.75))

    if Configuration['OUTER_RINGS_DOUBLED']:
        if len(modifier_list) > 10:
            modifier_list[10:]= np.sqrt(modifier_list[10:])
    #We also need a correction based on Z0
    #Get the avergae scaleheight for each ring
    Z0_av,modifier_list = make_equal_length(Z0_av,modifier_list)
    Z0_kpc = convertskyangle(Configuration,Z0_av)

    if not isiterable(Z0_kpc):
        Z0_kpc = [Z0_kpc]
    #Scale the limits with the deviation away from 0.2 kpc as this is more or less the unmodified scale height

    modifier_list=[set_limits(x*(1.125-0.625*y)*set_limits(((Configuration['RING_SIZE']*Configuration['BEAM'][0])/45.)**0.25,0.75,1.25),0.5,3.) for x,y in zip(modifier_list,Z0_kpc)]
    #The 1.1 is emprically determined and might need modification. Lets put it in advanced configuration
    Configuration['LIMIT_MODIFIER'] = np.array(modifier_list,dtype=float)*Configuration['LIMIT_MODIFIER_FACTOR']
    print_log(f'''SET_LIMIT_MODIFIER: Based on a Z0 in kpc {Z0_kpc}
{'':8s} and Inclination = {Inclination}
{'':8s} We updated the LIMIT_MODIFIER to {Configuration['LIMIT_MODIFIER']}.
''', Configuration,case=['debug_start'])


set_limit_modifier.__doc__ =f'''
 NAME:
    set_limit_modifier

 PURPOSE:
    Write to the Configuration a value or list of values to modify sbr dependent values. Foremost the sbr_limits

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard Tirific template

 OPTIONAL INPUTS:


 OUTPUTS:
    'LIMIT_MODIFIER' is updated in the Configuration.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_limits(value,minv,maxv):
    if value < minv:
        return minv
    elif value > maxv:
        return maxv
    else:
        return value

set_limits.__doc__ =f'''
 NAME:
    set_limits
 PURPOSE:
    Make sure Value is between min and max else set to min when smaller or max when larger.
 CATEGORY:
    support_functions

 INPUTS:
    value = value to evaluate
    minv = minimum acceptable value
    maxv = maximum allowed value

 OPTIONAL INPUTS:


 OUTPUTS:
    the limited Value

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
def set_ring_size(Configuration,requested_ring_size = 0., size_in_beams = 0., \
        check_set_rings = False):
    double_size = False
    if size_in_beams == 0.:
        size_in_beams =  np.min(Configuration['SIZE_IN_BEAMS'])
        double_size=True
    if requested_ring_size == 0.:
         requested_ring_size =  Configuration['RING_SIZE']

    print_log(f'''SET_RING_SIZE: Starting with the following parameters.
{'':8s}minimum size in beams = {size_in_beams}
{'':8s}requested ring size = {requested_ring_size}
''', Configuration,case=['debug_start'])
    if not Configuration['FIX_RING_SIZE']:
        no_rings = 0.
        ring_size = 100.
        while (requested_ring_size > 0.5 and no_rings < Configuration['MINIMUM_RINGS']) or no_rings == 0.:

            est_rings = set_limits(round((size_in_beams-1./5.)/(requested_ring_size)),8.,float('NaN'))

            ring_size = set_limits((size_in_beams-1./5.)/est_rings,0.5,float('NaN') )

            no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams )

            if no_rings < Configuration['MINIMUM_RINGS'] and requested_ring_size > 0.5:
                previous_ringsize = copy.deepcopy(requested_ring_size)
                requested_ring_size = set_limits(requested_ring_size/1.25,0.5,float('NaN') )
                print_log(f'''SET_RING_SIZE: Because we had less than {Configuration['MINIMUM_RINGS']} rings we have reduced the requested ring size from {previous_ringsize} to {requested_ring_size}
    ''',Configuration)

        if double_size:
            # we had two values in the size_in beams we get the number of rings from the maximum
            #As we used the minimum to calate the ring size we need to recalculated the number of rings
            no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=np.max(Configuration['SIZE_IN_BEAMS']) )
            size_in_beams = copy.deepcopy(Configuration['SIZE_IN_BEAMS'])
    else:
        ring_size = requested_ring_size
        no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=np.max(Configuration['SIZE_IN_BEAMS']) )

    print_log(f'''SET_RING_SIZE: After checking the size we get
{'':8s}size in beams = {size_in_beams}, ring size = {ring_size} and the number of ring = {no_rings}
''', Configuration,case=['debug_add'])

    if not check_set_rings:
        if no_rings <  Configuration['MINIMUM_RINGS']:
            print_log(f'''SET_RING_SIZE: With a ring size of {ring_size} we still only find {no_rings}.
{"":8s}SET_RING_SIZE: This is not enough for a fit.
''',Configuration,case=['main','screen'])
            raise SmallSourceError('This source is too small to reliably fit.')
        Configuration['NO_RINGS'] = int(no_rings)
        Configuration['RING_SIZE'] = ring_size
    return ring_size,int(no_rings)

set_ring_size.__doc__ =f'''
 NAME:
    set_ring_size

 PURPOSE:
    Calculate the size and number of rings and the galaxy size in beams and update them in the Configuration dictionary

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    size_in_beams = 0.
    The size of the galaxy in beams across the major axis radius, if unset it is taken from the Configuration

    check_set_rings = False
    Set this parameter to True to not apply the calaculated ring size and number and size of the model to the
    Configuration but to return the parameters as  size_in_beams,ring_size,no_rings

 OUTPUTS:
    Will raise an error when the number of rings goes below the minimum

 OPTIONAL OUTPUTS:
    size_in_beams,ring_size,no_rings

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_rings(Configuration,ring_size = 0. ):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    print_log(f'''SET_RINGS: Starting with the following parameters.
{'':8s}RING_SIZE = {ring_size}
''', Configuration, case=['debug_start'])
    no_rings = calc_rings(Configuration )
    print_log(f'''SET_RINGS: We find {no_rings} rings.
''', Configuration, case=['debug_add'])
    if Configuration['OUTER_RINGS_DOUBLED']:
        print_log(f'''SET_RINGS: This is a large galaxy. Therefore we use twice the ring size in the outer parts.
''',Configuration,case=['verbose'])
        radii = [0.,1./5.*Configuration['BEAM'][0]]
        radii = np.hstack((radii,(np.linspace(Configuration['BEAM'][0]*ring_size,Configuration['BEAM'][0]*10.*ring_size, \
                                        10)+1./5*Configuration['BEAM'][0])))
        radii = np.hstack((radii,(np.linspace(Configuration['BEAM'][0]*12.*ring_size, \
                                              Configuration['BEAM'][0]*ring_size*(10+(no_rings-12.)*2.), \
                                              no_rings-12)) \
                                              +1./5*Configuration['BEAM'][0]))

    else:
        radii = [0.,1./5.*Configuration['BEAM'][0]]
        radii = np.hstack((radii,(np.linspace(Configuration['BEAM'][0]*ring_size,Configuration['BEAM'][0]*ring_size*(no_rings-2.), \
                                        no_rings-2)+1./5.*Configuration['BEAM'][0])))
    if Configuration['DEBUG']:
        if Configuration['OUTER_RINGS_DOUBLED']:
            req_outer_ring = 2.*Configuration['BEAM'][0]*ring_size
        else:
            req_outer_ring = Configuration['BEAM'][0]*ring_size
        print_log(f'''SET_RINGS: Got the following radii.
{'':8s}{radii}
{'':8s}We should have {Configuration['NO_RINGS']} or we incorrectly updated and should have {no_rings}
{'':8s}We have {len(radii)} rings.
{'':8s}The last ring should be around {Configuration['BEAM'][0]*np.max(Configuration['SIZE_IN_BEAMS'])}
{'':8s}The rings should be size {Configuration['BEAM'][0]*ring_size} and outer rings {req_outer_ring}
{'':8s}They are {radii[2]-radii[1]} and {radii[-1]-radii[-2]}
''', Configuration, case=['debug_add'])
    return np.array(radii,dtype = float)

set_rings.__doc__ =f'''
 NAME:
    set_rings
 PURPOSE:
    Calculate the radii of all rings required for the model

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    ring_size = 0.
    Width of the rings

 OUTPUTS:
    numpy array with the central locations of the rings in arcsec.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def sine(x,amp,center,f,mean):
    return amp*np.sin(np.radians(x/f+abs(center)))+mean

sine.__doc__ =f'''
 NAME:
    sine

 PURPOSE:
    A sin function that can be modified

 CATEGORY:
    support_functions

 INPUTS:
    x = xaxis
    amp = amplitude size
    center = the location of the first highest point of the sin
    f = width of the size profile, should correspond to 2pi on the xaxis
    mean = the average value of the sin, basically the offset from 0.

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def sofia_output_exists(Configuration,Fits_Files):
    print_log(f'''SOFIA_OUTPUT_EXISTS: Starting check
''', Configuration, case=['debug_start'])
    req_files= ['MOMENT1','MOMENT0','MASK']
    for file in req_files:
        if os.path.exists(f"{Configuration['FITTING_DIR']}{Fits_Files[file]}"):
            continue
        else:
            log_statement=f'''SOFIA_OUTPUT_EXISTS: The file {Configuration['FITTING_DIR']}{Fits_Files[file]} is not found.'''
            print_log(log_statement, Configuration,case=['main','screen'])
            raise FileNotFoundError(log_statement)

    if not os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'):
        log_statement=f'''SOFIA_OUTPUT_EXISTS: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'} is not found.'''
        print_log(log_statement, Configuration,case=['main','screen'])
        raise FileNotFoundError(log_statement)

sofia_output_exists.__doc__ =f'''
 NAME:

 PURPOSE:
    Simple function to make sure all sofia output is present as expected
 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:


 OUTPUTS:
    Raises a FileNotFoundError if the files are not found.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def tirific_template(filename = ''):
    if filename == '':
        with pack_open_txt(templates, 'template.def') as tmp:
            template = tmp.readlines()
    elif filename == 'Installation_Check':
        from pyFAT_astro import Installation_Check as IC
        with pack_open_txt(IC, 'ModelInput.def') as tmp:
            template = tmp.readlines()
    else:
        with open(filename, 'r') as tmp:
            template = tmp.readlines()
    result = Proper_Dictionary()
    counter = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip().upper())
        if key == '':
            result[f'EMPTY{counter}'] = line
            counter += 1
        else:
            result[key] = str(line.split('=')[1].strip())
   
    return result
tirific_template.__doc__ ='''
 NAME:
    tirific_template

 PURPOSE:
    Read a tirific def file into a dictionary to use as a template.
    The parameter ill be the dictionary key with the values stored in that key

 CATEGORY:
    read_functions

 INPUTS:
    filename = Name of the def file

 OPTIONAL INPUTS:
    filename = ''
    Name of the def file, if unset the def file in Templates is used



 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      split, strip, open

 NOTE:
'''

