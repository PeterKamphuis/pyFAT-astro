# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in FAT to read input files

from pyFAT_astro.Support.support_functions import Proper_Dictionary,print_log,convertRADEC,set_limits, remove_inhomogeneities, \
                                obtain_border_pix, get_inclination_pa,get_vel_pa,columndensity,get_profile, get_kinematical_center,\
                                create_directory,copy_homemade_sofia,clean_header
from pyFAT_astro.Support.fits_functions import check_mask,clean_header

from astropy.io import fits
from astropy.wcs import WCS
from scipy import ndimage
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Ellipse
    import matplotlib.axes as maxes
import os
import time
import copy
import numpy as np
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res

import shutil
import traceback


from pyFAT_astro import Templates as templates
#Function to read a FAT input Catalogue
class BadCatalogueError(Exception):
    pass
class BadConfigurationError(Exception):
    pass
class NoConfigFile(Exception):
    pass

def catalogue(filename, debug = False):
    Catalogue = Proper_Dictionary({})
    with open(filename,'r') as tmpfile:
        #Define the exsiting catalogue input()
        input_columns = [x.strip().upper() for x in tmpfile.readline().split('|')]
        Catalogue['ENTRIES'] = ['ENTRIES']
        Catalogue['ENTRIES'].extend(input_columns)
        for key in input_columns:
            Catalogue[key] = []

        for line in tmpfile.readlines():
            input = [x.strip() for x  in line.split('|')]
            for i,key in enumerate(input_columns):
                if key == 'DISTANCE':
                    Catalogue[key].append(float(input[i]))
                else:
                    Catalogue[key].append(input[i])
    if 'NUMBER' in Catalogue['ENTRIES']:
        Catalogue['NUMBER'] = np.array(Catalogue['NUMBER'],dtype=int)

    return Catalogue
catalogue.__doc__ =f'''
 NAME:
    catalogue

 PURPOSE:
    Read the FAT input catalogue and write into the a dictionary

 CATEGORY:
    read_functions

 INPUTS:
    filename = name of the catalogue to read

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Catalogue = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_edge_limits(xmin,xmax,ymin,ymax,zmin,zmax,Configuration,debug=False ,beam_edge = 0.5, vel_edge = 0.5 ):
    diff = np.array([xmin,abs(xmax - Configuration['NAXES'][0]),
                     ymin,abs(ymax - Configuration['NAXES'][1])],dtype = float)
    if debug:
        print_log(f'''CHECK_EDGE_LIMIT: We find these differences and this edge size {beam_edge*Configuration['BEAM_IN_PIXELS'][0]}
{'':8s} diff  = {diff}
''',Configuration['OUTPUTLOG'],debug= True)
    if np.where(diff < beam_edge*Configuration['BEAM_IN_PIXELS'][0])[0].size:
        return True
    diff = np.array([zmin,abs(zmax-Configuration['NAXES'][2])],dtype=float)
    if debug:
        print_log(f'''CHECK_EDGE_LIMIT: And for velocity edge =  {vel_edge}
{'':8s} diff  = {diff}
''',Configuration['OUTPUTLOG'])
    if np.where(diff < vel_edge)[0].size:
        print(f"On the edge")
        return True
    else:
        print(f"Off the edge")
        return False
check_edge_limits.__doc__ =f'''
 NAME:
    check_edge_limits

 PURPOSE:
    Check whether a Sofia source is properly separated from the edge of the cube.
    In order to read

 CATEGORY:
    read_functions

 INPUTS:
    xmin,xmax,ymin,ymax,zmin,zmax,header,Configuration
    these are the parameters that describe the mask + header + Configuration

 OPTIONAL INPUTS:
    debug=False

    beam_edge = 0.5
    minimum tolerance for spatial edges

    vel_edge = 0.5
    minimum tolerance for velocity

 OUTPUTS:
    Boolean which is True when the source is closer to the edge than beam_edge and vel_edge False other wise

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    abs, np.where, np.array, print_log,

 NOTE:
'''

def extract_vrot(Configuration,map ,angle,center, debug= False):
    if debug:
        print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
{'':8s} PA= {angle}
{'':8s} center= {center}
''',Configuration['OUTPUTLOG'], debug = True)
    maj_profile,maj_axis,maj_resolution = get_profile(Configuration,map,angle,center=center,debug=debug)
    # We should base extracting the RC on where the profile is negative and positive to avoid mistakes in the ceneter coming through
    neg_index = np.where(maj_profile < 0.)[0]
    pos_index = np.where(maj_profile > 0.)[0]

    if debug:
        print_log(f'''EXTRACT_VROT: The resolution on the extracted axis
{'':8s} resolution = {maj_resolution}
''',Configuration['OUTPUTLOG'])
    avg_profile = []
    neg_profile = []
    #pos_profile = []
    diff = 0.
    counter = 1.
    for i in range(np.nanmin([neg_index.size,pos_index.size])):
        if not np.all(np.isnan([maj_profile[neg_index[neg_index.size-i-1]],-1*maj_profile[pos_index[i]]])):
            avg_profile.append(np.nanmean([-1*maj_profile[neg_index[neg_index.size-i-1]],maj_profile[pos_index[i]]]))
            #neg_profile.append(-1*maj_profile[neg_index[neg_index.size-i-1]])
            #pos_profile.append(maj_profile[pos_index[i]])
        else:
            avg_profile.append(avg_profile[-1])
            #neg_profile.append(float('NaN'))
            #pos_profile.append(float('NaN'))
        # When they start declining we just keep them flat
        if i > 1:
            if avg_profile[-1] < avg_profile[-2]:
                avg_profile[-1] = avg_profile[-2]
        #correct for beam smearing in the center
        beam_back = -1*int(Configuration['BEAM_IN_PIXELS'][0]/maj_resolution)

        if Configuration['BEAM_IN_PIXELS'][0]/maj_resolution < i < -2.*beam_back:
            avg_profile[beam_back] = avg_profile[beam_back]+0.25*avg_profile[int(beam_back/2.)]+0.1*avg_profile[-1]
        elif -2.*beam_back <= i < -3.*beam_back:
            avg_profile[beam_back] = avg_profile[beam_back]+0.25/counter*avg_profile[int(beam_back/2.)]+0.1/counter*avg_profile[-1]
            counter += 1
            #neg_profile[beam_back] = neg_profile[beam_back]+0.5*neg_profile[int(beam_back/2.)]+0.1*neg_profile[-1]
            #pos_profile[beam_back] = pos_profile[beam_back]+0.5*pos_profile[int(beam_back/2.)]+0.1*pos_profile[-1]
    if debug:
        print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
''',Configuration['OUTPUTLOG'])
    ring_size_req = Configuration['BEAM_IN_PIXELS'][0]/maj_resolution
    if debug:
        print_log(f'''EXTRACT_VROT: We need a rings size of
{'':8s} ringsize= {ring_size_req}
{'':8s} because bmaj in pixels  ={Configuration['BEAM_IN_PIXELS'][0]}  and the resolution of the profile = {maj_resolution} pixels
''',Configuration['OUTPUTLOG'])
    profile = np.array(avg_profile[0::int(ring_size_req)],dtype=float)
    if debug:
        print_log(f'''EXTRACT_VROT: Constructing the final RC
{'':8s} initial RC= {profile}
{'':8s} at step width= {ring_size_req}
''',Configuration['OUTPUTLOG'])
    try:
        profile[np.isnan(profile)] = profile[~np.isnan(profile)][-1]
    except IndexError:
        profile = []

    if len(profile) < 1 and len(avg_profile) > 0.:
        profile = [np.max(avg_profile)]
    return profile
extract_vrot.__doc__ =f'''
 NAME:
    extract_vrot

 PURPOSE:
    Extract a profile from the velocity field and average it over both sides from the center.

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map = the velocity field
    angle = PA of major axis
    center = center of the galaxy in pixel coordinates

 OPTIONAL INPUTS:
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:
    profile = The RC at given ring locations averaged over the approaching and receding side

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_DHI(Configuration,Model='Finalmodel' ,debug=False):
    #Get the sbrs
    radi,sbr,sbr_2,systemic = load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{Model}/{Model}.def",Variables = ['RADI','SBR','SBR_2','VSYS'],debug=debug)
    #convert to solar_mass/pc^2
    sbr_msolar = columndensity(Configuration,sbr*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = columndensity(Configuration,sbr_2*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    # interpolate these to ~1" steps
    new_radii = np.linspace(0,radi[-1],int(radi[-1]))
    new_sbr_msolar = np.interp(new_radii,radi,sbr_msolar)
    new_sbr_2_msolar = np.interp(new_radii,radi,sbr_2_msolar)

    index_1 = np.where(new_sbr_msolar > 1.)[0]
    index_2 = np.where(new_sbr_2_msolar > 1.)[0]
    if index_1.size > 0 and index_2.size > 0:
        DHI = float(new_radii[index_1[-1]]+new_radii[index_2[-1]])
    elif index_1.size > 0:
        DHI = float(new_radii[index_1[-1]])
    elif index_2.size > 0:
        DHI = float(new_radii[index_2[-1]])
    else:
        DHI = float('NaN')
    return DHI
get_DHI.__doc__ =f'''
 NAME:
    get_DHI

 PURPOSE:
    get the DHI as determined by the SBR profiles in the fit from the Tirific Template

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

    Model = 'Finalmodel'
    location of the def file to get DHI from. it should be in the fitting dir in the {{Model}}/{{Model}}.def

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_totflux(Configuration,map_name,debug=False):
    image = fits.open(f"{Configuration['FITTING_DIR']}{map_name}")
    flux = np.nansum(image[0].data)
    error = np.sqrt((np.where(image[0].data> 0.)[0].size)/Configuration['BEAM_IN_PIXELS'][2])*Configuration['NOISE']
    image.close()
    return [float(flux/Configuration['BEAM_IN_PIXELS'][2]),error]
get_totflux.__doc__ =f'''
 NAME:
    get_totflux

 PURPOSE:
    Get the total flux from a intensity map

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map_name = name of the intensity fits file

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    total flux in the map in Jy*km/s

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# Function to get the PA and inclination from the moment 0 for initial estimates
def guess_orientation(Configuration,Fits_Files, v_sys = -1 ,center = None, debug = False):
    #open the moment 0
    if debug:
        print_log(f'''GUESS_ORIENTATION: starting extraction of initial parameters.
''',Configuration['OUTPUTLOG'], debug = True)
    Image = fits.open(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['MOMENT0']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = Image[0].data
    hdr = Image[0].header
    mom0 = copy.deepcopy(Image)
    Image.close()

    if not center:
        center = [hdr['NAXIS1']/2.-1,hdr['NAXIS2']/2.-1]
    Image = fits.open(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['CHANNEL_MAP']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    noise_map = np.sqrt(Image[0].data)*Configuration['NOISE']*Configuration['CHANNEL_WIDTH']

    SNR = np.nanmean(map[noise_map > 0.]/noise_map[noise_map > 0.])
    noise_hdr = Image[0].header
    Image.close()
    median_noise_in_map = np.nanmedian(noise_map[noise_map > 0.])
    minimum_noise_in_map = np.nanmin(noise_map[noise_map > 0.])

    map[0.5*minimum_noise_in_map > noise_map] = 0.
    mom0[0].data= map
    scale_factor = set_limits(SNR/3.*minimum_noise_in_map/median_noise_in_map, 0.05, 1.)
    if debug:
        print_log(f'''GUESS_ORIENTATION: We find SNR = {SNR} and a scale factor {scale_factor} and the noise median {median_noise_in_map}
{'':8s} minimum {minimum_noise_in_map}
''',Configuration['OUTPUTLOG'])
    beam_check=[Configuration['BEAM_IN_PIXELS'][0],Configuration['BEAM_IN_PIXELS'][0]/2.]
    inclination_av, pa_av, maj_extent_av = get_inclination_pa(Configuration, mom0, center, cutoff = scale_factor* median_noise_in_map, debug = debug)

    inclination_av = [inclination_av]
    int_weight = [2.]
    pa_av = [pa_av]
    maj_extent_av = [maj_extent_av]
    for mod in beam_check:
        for i in [[-1,-1],[-1,1],[1,-1],[1,1]]:
            center_tmp = [center[0]+mod*i[0],center[1]+mod*i[1]]
            inclination_tmp, pa_tmp, maj_extent_tmp= get_inclination_pa(Configuration, mom0, center_tmp, cutoff = scale_factor* median_noise_in_map, debug = debug)
            inclination_av.append(inclination_tmp)
            pa_av.append(pa_tmp)
            int_weight.append(mod/beam_check[0]*0.25)
            maj_extent_av.append(maj_extent_tmp)
    int_weight = np.array(int_weight)
    weight = np.array([1./x[1] for x in inclination_av],dtype= float)*int_weight
    inclination = np.array([np.nansum(np.array([x[0] for x in inclination_av],dtype=float)*weight)/np.nansum(weight),\
                            np.nansum(np.array([x[1] for x in inclination_av],dtype=float)*weight)/np.nansum(weight)],dtype=float)
    weight = np.array([1./x[1] for x in pa_av],dtype= float)
    pa = np.array([np.nansum(np.array([x[0] for x in pa_av],dtype=float)*weight)/np.nansum(weight),\
                            np.nansum(np.array([x[1] for x in pa_av],dtype=float)*weight)/np.nansum(weight)],dtype=float)

    maj_extent= np.nansum(maj_extent_av*weight)/np.nansum(weight)
    # For very small galaxies we do not want to correct the extend
    if maj_extent/Configuration['BEAM'][0]/3600. > 3.:
        maj_extent = maj_extent+(Configuration['BEAM'][0]/3600.*0.2/scale_factor)

    if debug:
        print_log(f'''GUESS_ORIENTATION: From the maps we find
{'':8s} inclination = {inclination}
{'':8s} pa = {pa}
{'':8s} size in beams = {maj_extent/(Configuration['BEAM'][0]/3600.)}
''',Configuration['OUTPUTLOG'])
            #map[3*minimum_noise_in_map > noise_map] = 0.
    # From these estimates we also get an initial SBR
    maj_profile,maj_axis,maj_resolution = get_profile(Configuration,map, pa[0],center,debug=debug)
    # let's get an intensity weighted center for the extracted profile.

    center_of_profile = np.sum(maj_profile*maj_axis)/np.sum(maj_profile)
    neg_index = np.where(maj_axis < 0.)[0]
    pos_index = np.where(maj_axis > 0.)[0]
    avg_profile = []
    neg_profile = []
    pos_profile = []
    diff = 0.
    for i in range(np.nanmin([neg_index.size,pos_index.size])):
        avg_profile.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],maj_profile[pos_index[i]]]))
        neg_profile.append(maj_profile[neg_index[neg_index.size-i-1]])
        pos_profile.append(maj_profile[pos_index[i]])
        #diff = diff+abs(avg_profile[i]-neg_profile[i])+abs(avg_profile[i]-pos_profile[i])
        diff = diff+abs(pos_profile[-1]-neg_profile[-1])*abs(np.mean([pos_profile[-1],neg_profile[-1]]))
    diff = diff/np.nanmin([neg_index.size,pos_index.size])
    if debug:
        print_log(f'''GUESS_ORIENTATION:'BMAJ in pixels, center of profile, center, difference between pos and neg
{'':8s}{Configuration['BEAM_IN_PIXELS'][0]*0.5} {center_of_profile} {center} {diff}
''',Configuration['OUTPUTLOG'],screen=True)

    # if the center of the profile is more than half a beam off from the Sofia center let's see which on provides a more symmetric profile
    if (abs(center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution) > Configuration['BEAM_IN_PIXELS'][0]*0.5 \
        or abs(center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution) > Configuration['BEAM_IN_PIXELS'][0]*0.5) and SNR > 3.:
        if debug:
                print_log(f'''GUESS_ORIENTATION: The SoFiA center and that of the SBR profile are separated by more than half a beam.
{'':8s}GUESS_ORIENTATION: Determining the more symmetric profile.
''',Configuration['OUTPUTLOG'])
        neg_index = np.where(maj_axis < center_of_profile)[0]
        pos_index = np.where(maj_axis > center_of_profile)[0]
        avg_profile_new = []
        neg_profile_new = []
        pos_profile_new = []
        diff_new =0.
        for i in range(np.nanmin([neg_index.size,pos_index.size])):
            avg_profile_new.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],maj_profile[pos_index[i]]]))
            neg_profile_new.append(maj_profile[neg_index[neg_index.size-i-1]])
            pos_profile_new.append(maj_profile[pos_index[i]])
            #diff_new = diff_new+abs(avg_profile_new[i]-neg_profile_new[i])+abs(avg_profile_new[i]-pos_profile_new[i])
            diff_new = diff_new+abs(pos_profile_new[-1]-neg_profile_new[-1])*abs(np.mean([pos_profile_new[-1],neg_profile_new[-1]]))
        diff_new = diff_new/np.nanmin([neg_index.size,pos_index.size])
        if debug:
            print_log(f'''GUESS_ORIENTATION: We have this old difference {diff} and this new difference {diff_new}
''',Configuration['OUTPUTLOG'])
        if diff_new < diff:
            if debug:
                print_log(f'''GUESS_ORIENTATION: We are updating the center from {center[0]},{center[1]} to {center[0]-center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution},{center[1]+center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution}
''',Configuration['OUTPUTLOG'])
            avg_profile = avg_profile_new
            center[0] = center[0]-center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution
            center[1] = center[1]+center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution
            maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa[0])))+abs(abs(center[1])*np.cos(np.radians(pa[0]))))


    ring_size_req = Configuration['BEAM_IN_PIXELS'][0]/maj_resolution
    SBR_initial = avg_profile[0::int(ring_size_req)]/(np.pi*Configuration['BEAM'][0]*Configuration['BEAM'][1]/(4.*np.log(2.))) # Jy*km/s
    SBR_initial =np.hstack((SBR_initial[0],SBR_initial,SBR_initial[-1]))

    SBR_initial[0:3] = SBR_initial[0:3] * (1.2 -float(inclination[0])/90.)


    #We need to know which is the approaching side and which is receding



    Image = fits.open(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['MOMENT1']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = copy.deepcopy(Image[0].data)
    hdr = Image[0].header
    if debug:
        print_log(f'''GUESS_ORIENTATION: This is the amount of values we find initially {len(map[noise_map > 0.])}
''',Configuration['OUTPUTLOG'])
    #Image.close()
    #map[3*minimum_noise_in_map > noise_map] = float('NaN')
    map[3*minimum_noise_in_map > noise_map] = float('NaN')
    if debug:
        print_log(f'''GUESS_ORIENTATION: This is the amount of values we find after blanking low SNR {len(map[~np.isnan(map)])}
''',Configuration['OUTPUTLOG'])
    if len(map[~np.isnan(map)]) < 5:
        no_values  = True
        noise_level = 2.5
        sigma = [0.5,0.5]
        while no_values:
            if debug:
                print_log(f'''GUESS_ORIENTATION: We smooth the velocity field to use a lower Noise threshold
''',Configuration['OUTPUTLOG'])
            tmp = copy.deepcopy(Image[0].data)
            tmp =  ndimage.gaussian_filter(Image[0].data, sigma=(sigma[1], sigma[0]), order=0)
            tmp[noise_level*minimum_noise_in_map > noise_map] =  float('NaN')
            if debug:
                print_log(f'''GUESS_ORIENTATION: We used the threshold {noise_level} and sigme = {sigma}
{'':8s}We find the len {len(tmp[~np.isnan(tmp)])}
''',Configuration['OUTPUTLOG'])
            if len(tmp[~np.isnan(tmp)]) < 5:
                sigma = [sigma[0]+0.5,sigma[0]+0.5]
                noise_level -= 0.5
            else:
                no_values = False
                map = tmp
            if noise_level < 1. and no_values:
                VROT_initial = [0]
                SBR_initial = [0]
                Image.close()
                return np.array(pa,dtype=float),np.array(inclination,dtype=float),SBR_initial,maj_extent,center[0],center[1],VROT_initial

    Image.close()
    noise_map = []

    # First we look for the kinematics center
    #found_vsys = get_kinematical_center(Configuration,map,float(pa[0]),center=center)
    #map_vsys = found_vsys[0]
    #if np.sum([abs(found_vsys[1]), abs(found_vsys[2])]) != 0.:
    #    center = [center[0]+found_vsys[1], center[1]+found_vsys[2]]
    #    print_log(f'''GUESS_ORIENTATION: We are updating the center with x,y = {found_vsys[1:]}
#{'':8s}to {center}
#''',Configuration['OUTPUTLOG'])
    #if debug:
    #    print_log(f'''GUESS_ORIENTATION: We found the following initial VSYS:
#{'':8s}vsys = {map_vsys}, at {center}
#''',Configuration['OUTPUTLOG'])
    try:
        vel_pa = get_vel_pa(Configuration,map,center=center,debug=debug)
    except:
        vel_pa = [float('NaN'),float('NaN')]

    maj_profile,maj_axis,maj_resolution = get_profile(Configuration,map,pa[0], center=center,debug=debug)
    loc_max = np.mean(maj_axis[np.where(maj_profile == np.nanmax(maj_profile))[0]])

    if loc_max > 0.:
        pa[0] = pa[0]+180
        print_log(f'''GUESS_ORIENTATION: We have modified the pa by 180 deg as we found the maximum velocity west of the center.
''' , Configuration['OUTPUTLOG'])
    if abs(pa[0]-vel_pa[0]) > 25. and ~np.isnan(vel_pa[0]):
        pa=vel_pa
    else:
        if ~np.isnan(vel_pa[0]):
            pa = [np.nansum([vel_pa[0]/vel_pa[1],pa[0]/pa[1]])/np.nansum([1./vel_pa[1],1./pa[1]]),2.*1./np.sqrt(np.nansum([1./vel_pa[1],1./pa[1]]))]
    #As python is utterly moronic the center goes in back wards to the map
    buffer = int(round(np.mean(Configuration['BEAM_IN_PIXELS'][:2])/2.))
    if debug:
        print_log(f'''GUESS_ORIENTATION: We start with vsys = {v_sys:.2f} km/s
''' , Configuration['OUTPUTLOG'])

    if v_sys == -1:
        map_vsys = np.nanmean(map[int(round(center[1]-buffer)):int(round(center[1]+buffer)),int(round(center[0]-buffer)):int(round(center[0]+buffer))])
    else:
        map_vsys = v_sys
    if debug:
        print_log(f'''GUESS_ORIENTATION: We subtract {map_vsys:.2f} km/s from the moment 1 map to get the VROT
''' , Configuration['OUTPUTLOG'])
    map = map  - map_vsys
    VROT_initial = extract_vrot(Configuration,map ,pa[0],center, debug= debug)
    min_RC_length= len(VROT_initial)
    if pa[1] < 10.:
        for x in [pa[0]-pa[1],pa[0]-pa[1]/2.,pa[0]+pa[1]/2.,pa[0]+pa[1]]:
            tmp  = extract_vrot(Configuration,map ,x,center, debug= debug)
            if len(tmp) > 0.:
                RC_length = len(tmp)
                if RC_length < min_RC_length:
                    min_RC_length = RC_length
                if debug:
                    print_log(f'''GUESS_ORIENTATION: We found the lengths for the angled rotation curves:
{'':8s}GUESS_ORIENTATION: RC length = {RC_length}, min_RC length = {min_RC_length}, Vrot ini = {len(VROT_initial)}, tmnp = {len(tmp)}
{'':8s}GUESS_ORIENTATION: Vrot {VROT_initial}
{'':8s}GUESS_ORIENTATION: tmp {tmp}
''',Configuration['OUTPUTLOG'])
                VROT_initial = np.vstack((VROT_initial[:min_RC_length],tmp[:min_RC_length]))
                VROT_initial = np.mean(VROT_initial,axis=0)

    map= []
    if len(VROT_initial) < 1:
        VROT_initial = 0.
    elif len(VROT_initial) == 1:
        VROT_initial = np.array([0.,VROT_initial[0]])
    else:
        VROT_initial[0] = 0

    VROT_initial = np.abs(VROT_initial/np.sin(np.radians(inclination[0])))
    if debug:
        print_log(f'''GUESS_ORIENTATION: We found the following initial rotation curve:
{'':8s}GUESS_ORIENTATION: RC = {VROT_initial}
''',Configuration['OUTPUTLOG'])
    return np.array(pa,dtype=float),np.array(inclination,dtype=float),SBR_initial,maj_extent,center[0],center[1],VROT_initial
guess_orientation.__doc__ =f'''
 NAME:
    guess_orientation

 PURPOSE:
    Obtain the initial estimates for the galaxy from the SoFiA products

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:
    debug = False

    center = None
    Will assume the center of the map provided when unset.
    In pixel coordinates
 OUTPUTS:

 OPTIONAL OUTPUTS:
    np.array(pa) = [pa, pa errors]
    np.array(inclination) = [inc, inc error]
    SBR_initial = Initial guess of the SBR profile from the moment 0 map
    maj_extent = Initial guess of the extend of the galaxy
    center[0] = RA in pixels
    center[1] = DEC in pixels
    VROT_initial = initial RC

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# load the basic info file to get the Sofia FAT Initial_Estimates
def load_basicinfo(Configuration,filename, Variables = ['RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask','Tot FLux','D_HI','Distance','HI_Mass' ,'D_HI_kpc' ], unpack = True, debug = False):
    outputarray = np.zeros((2,len(Variables)), dtype=float)
    try:
        with open(filename, 'r') as tmp:
            fileline = tmp.readlines()
    except FileNotFoundError:
        print("That Basic info file does not exist. Returning empty array.")
        if unpack:
            return (*outputarray.T,)
        else:
            return outputarray
    Var_inFile = [f.strip() for f in fileline[3].split('  ') if f != '']
    del Var_inFile[0]
    invalues = [f.strip() for f in fileline[6].split('  ') if f != '']
    for i,var in enumerate(Variables):
        if var == 'RA' or var == 'DEC':
            tmp_str = invalues[Var_inFile.index('RA')].split('+/-')
            tmp_str2 = invalues[Var_inFile.index('DEC')].split('+/-')
            RA, DEC = convertRADEC(Configuration,tmp_str[0],tmp_str2[0],invert=True,debug=debug)
            if var == 'RA':
                outputarray[:,i] = [RA,float(tmp_str[1])]
            if var == 'DEC':
                outputarray[:,i] = [DEC,float(tmp_str2[1])]
        else:
            outputarray[:,i] = invalues[Var_inFile.index(var)].split('+/-')
    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray
load_basicinfo.__doc__ =f'''
 NAME:
    load_basicinfo

 PURPOSE:
    Load the file with the basic info of the different steps to get the SoFiA initial estimates stored there.

 CATEGORY:
    read_functions

 INPUTS:
    filename = The name of the basic info file
 OPTIONAL INPUTS:
    Variables = ['RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask','Tot FLux','D_HI','Distance','HI_Mass' ,'D_HI' ]
    Variable to extract from the file

    unpack = True
    Unpack the values when returning or not

    debug = False

 OUTPUTS:
    outputarray = Array with values of all requested parameters

 OPTIONAL OUTPUTS:
       -

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Function for loading the variables of a tirific def file into a set of variables to be used
def load_template(Configuration,Template,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True  ,debug = False):
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    if debug:
        print_log(f'''LOAD_TEMPLATE: We get the following number of rings from the template {Template['NUR']}
''',Configuration['OUTPUTLOG'], debug=True)
    numrings = int(Template['NUR'])
    outputarray=np.zeros((numrings,len(Variables)),dtype=float)
    counter = 0
    for var in Variables:
        if debug:
            print_log(f'''LOAD_TEMPLATE: We are processing {var}.
{'':8s}LOAD_TEMPLATE: With the following values {Template[var]}
''',Configuration['OUTPUTLOG'])
        tmp =  np.array(Template[var].rsplit(),dtype=float)
        outputarray[0:len(tmp),counter] = tmp[0:len(tmp)]
        counter +=1

    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray
load_template.__doc__ =f'''
 NAME:
    load_template

 PURPOSE:
    Load values from variables set in the tirific templates

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Template = The tirific template to extract values from, this is a dictionary where the def file is set as  Template['parameter'] = values

 OPTIONAL INPUTS:
    Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']

    unpack = True
    if true unpack the values

    debug = False

 OUTPUTS:
    outputarray array with all the values of the parameters requested

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


#Function for loading the variables of a tirific def file into a set of variables to be used
def load_tirific(Configuration,filename,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True , debug = False ):
    if debug:
        print_log(f'''LOAD_TIRIFIC: Starting to extract the following paramaters:
{'':8s}{Variables}
''',Configuration['OUTPUTLOG'], debug = True)
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    numrings = []
    while len(numrings) < 1:
        time.sleep(0.1)
        with open(filename, 'r') as tmp:
            numrings = [int(e.split('=')[1].strip()) for e in tmp.readlines() if e.split('=')[0].strip().upper() == 'NUR']

    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
    with open(filename, 'r') as tmp:
        unarranged = tmp.readlines()
    # Separate the keyword names
    for line in unarranged:
        var_concerned = str(line.split('=')[0].strip().upper())
        if len(var_concerned) < 1:
            var_concerned = 'xxx'
        varpos = np.where(Variables == var_concerned)[0]
        if varpos.size > 0:
            tmp =  np.array(line.split('=')[1].rsplit(),dtype=float)
            if len(outputarray[:,0]) < len(tmp):
                tmp_out=outputarray
                outputarray = np.zeros((len(tmp), len(Variables)), dtype=float)
                outputarray[0:len(tmp_out),:] = tmp_out
            outputarray[0:len(tmp),int(varpos)] = tmp[0:len(tmp)]
        else:
            if var_concerned[0] == '#':
                varpos = np.where(var_concerned[2:].strip() == Variables)[0]
                if varpos.size > 0:
                    tmp = np.array(line.split('=')[1].rsplit(),dtype=float)
                    if len(outputarray[:, 0]) < len(tmp):
                        tmp_out = outputarray
                        outputarray = np.zeros((len(tmp), len(Variables)), dtype=float)
                        outputarray[0:len(tmp_out), :] = tmp_out
                    outputarray[0:len(tmp),int(varpos)] = tmp[:]
    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray
load_tirific.__doc__ =f'''
 NAME:
    load_tirific

 PURPOSE:
    Load values from variables set in the tirific files

 CATEGORY:
    read_functions

 INPUTS:
    filename = Path to the tirific def file

 OPTIONAL INPUTS:
    Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2']

    unpack = True
    if true unpack the values

    debug = False

 OUTPUTS:
    outputarray array with all the values of the parameters requested

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def read_cube(Configuration,cube,debug=False):
    cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}{cube}")
    Configuration['NOISE'] = cube_hdr['FATNOISE']
    Configuration['CHANNEL_WIDTH'] = cube_hdr['CDELT3']/1000.
    Configuration['PIXEL_SIZE'] = np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])])
    Configuration['NAXES'] = [cube_hdr['NAXIS1'],cube_hdr['NAXIS2'], cube_hdr['NAXIS3']]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        coordinate_frame = WCS(cube_hdr)
        xlow,ylow,zlow = coordinate_frame.wcs_pix2world(1,1,1., 1.)
        xhigh,yhigh,zhigh = coordinate_frame.wcs_pix2world(*Configuration['NAXES'], 1.)
        Configuration['NAXES_LIMITS'] = [np.sort([xlow,xhigh]),np.sort([ylow,yhigh]),np.sort([zlow,zhigh])/1000.]
    # We write the pixels per beam info to Configuration such that it is easily accesible
    beamarea=(np.pi*abs(cube_hdr['BMAJ']*cube_hdr['BMIN']))/(4.*np.log(2.))
    Configuration['BEAM_AREA'] = beamarea*3600.**2 # beamarea in arcsec
    Configuration['BEAM_IN_PIXELS'] = [cube_hdr['BMAJ']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                       cube_hdr['BMIN']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                       beamarea/(abs(cube_hdr['CDELT1'])*abs(cube_hdr['CDELT2']))]
    # Ad the major beam to configuration as we need it in many places
    Configuration['BEAM'] = [float(cube_hdr['BMAJ']*3600.), float(cube_hdr['BMIN']*3600.),float(cube_hdr['BPA'])]
    #Let's set some maximum errors based on the input cube
    Configuration['MAX_ERROR'] = {'VROT': [Configuration['CHANNEL_WIDTH']*5.], \
                                  'VSYS': [Configuration['CHANNEL_WIDTH']*1.5], \
                                  'SBR': [Configuration['NOISE']/Configuration['BEAM_AREA']*Configuration['CHANNEL_WIDTH']*3.],\
                                  'PA' : [15.],\
                                  'INCL': [15.],\
                                  'SDIS': [Configuration['CHANNEL_WIDTH']*2.5],\
                                  'Z0' : [Configuration['BEAM'][0]],\
                                  'XPOS': [cube_hdr['BMAJ']],\
                                  'YPOS': [cube_hdr['BMAJ']],\
    }
    Configuration['MIN_ERROR'] = {'VROT': [Configuration['CHANNEL_WIDTH']*0.5], \
                                  'VSYS': [Configuration['CHANNEL_WIDTH']*0.1], \
                                  'SBR': [Configuration['NOISE']/Configuration['BEAM_AREA']*Configuration['CHANNEL_WIDTH']*0.3],\
                                  'PA' : [1.],\
                                  'INCL': [2.],\
                                  'SDIS': [Configuration['CHANNEL_WIDTH']*0.1],\
                                  'Z0' : [Configuration['BEAM'][0]*0.1],\
                                  'XPOS': [cube_hdr['BMAJ']*0.1],\
                                  'YPOS': [cube_hdr['BMAJ']*0.1],\
                                }
# function to read the sofia catalogue
def sofia_catalogue(Configuration,Fits_Files, Variables =['id','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra',\
                    'dec','v_app','f_sum','kin_pa','w50','err_f_sum','err_x','err_y','err_z'], debug = False):
    if debug:
        print_log(f'''SOFIA_CATALOGUE: Reading the source from the catalogue.
''',Configuration['OUTPUTLOG'],debug= True)
    outlist = [[] for x in Variables]
    with open(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt') as sof_cat:
        for line in sof_cat.readlines():
            tmp =line.split()
            if line.strip() == '' or line.strip() == '#':
                pass
            elif tmp[0] == '#' and len(tmp) > 1:
                if tmp[1].strip().lower() in ['name','id']:
                    # get the present columns
                    input_columns  = [x.strip() for x in tmp[1:]]
                    #determin their location in the line
                    column_locations = []
                    for col in input_columns:
                        column_locations.append(line.find(col)+len(col))
                    # check that we found all parameters
                    for value in Variables:
                        if value.lower() in input_columns:
                            continue
                        else:
                            print_log(f'''SOFIA_CATALOGUE: We cannot find the required column for {value} in the sofia catalogue.
{"":8s}SOFIA_CATALOGUE: This can happen because a) you have tampered with the sofiainput.txt file in the Support directory,
{"":8s}SOFIA_CATALOGUE: b) you are using an updated version of SoFiA2 where the names have changed and FAT is not yet updated.'
{"":8s}SOFIA_CATALOGUE:    In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
{"":8s}SOFIA_CATALOGUE: c) You are using pre processed SoFiA output of your own and do not have all the output'
{"":8s}SOFIA_CATALOGUE:    Required output is {','.join(Variables)})
''',Configuration['OUTPUTLOG'],screen= True)
                            raise BadCatalogueError("SOFIA_CATALOGUE: The required columns could not be found in the sofia catalogue.")
            else:
                for col in Variables:
                    if input_columns.index(col) == 0:
                        start = 0
                    else:
                        start = column_locations[input_columns.index(col)-1]
                    end = column_locations[input_columns.index(col)]
                    outlist[Variables.index(col)].append(line[start:end].strip())
    #if we have a header we convert w50 or w20
    if ('w50' in Variables or 'w20' in Variables or 'f_sum' in Variables or 'err_f_sum' in Variables):
        if 'w50' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('w50')][i] = float(outlist[Variables.index('w50')][i])*Configuration['CHANNEL_WIDTH']
        if 'w20' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('w20')][i] = float(outlist[Variables.index('w20')][i])*Configuration['CHANNEL_WIDTH']
        if 'f_sum' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('f_sum')][i] = float(outlist[Variables.index('f_sum')][i])/Configuration['BEAM_IN_PIXELS'][2]*Configuration['CHANNEL_WIDTH']
        if 'err_f_sum' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('err_f_sum')][i] = float(outlist[Variables.index('err_f_sum')][i])/Configuration['BEAM_IN_PIXELS'][2]*Configuration['CHANNEL_WIDTH']

    # we want to fit a specific source
    if len(outlist[0]) > 1 and 'f_sum' in Variables:
        if debug:
            print_log(f'''SOFIA_CATALOGUE: Multiple sources were found we will try to select the correct one.
''',Configuration['OUTPUTLOG'])
        if Configuration['SOFIA_RAN']:
            found = False
            beam_edge=2.
            if Configuration['VEL_SMOOTH_EXTENDED'] or Configuration['HANNING_SMOOTHED'] :
                vel_edge = 1.
                min_vel_edge = 0.
            else:
                vel_edge = 2.
                min_vel_edge = 1.

            while not found:
                many_sources  = copy.deepcopy(outlist)
                # We want to exclude any edge sources
                for i in range(len(many_sources[0])):

                    edge = check_edge_limits(float(many_sources[Variables.index('x_min')][i]),
                                    float(many_sources[Variables.index('x_max')][i]),
                                    float(many_sources[Variables.index('y_min')][i]),
                                    float(many_sources[Variables.index('y_max')][i]),
                                    float(many_sources[Variables.index('z_min')][i]),
                                    float(many_sources[Variables.index('z_max')][i]),
                                    Configuration, debug = debug,beam_edge = beam_edge, vel_edge= vel_edge)
                    if edge:
                        many_sources[Variables.index('f_sum')][i]=0.
                if np.nansum(many_sources[Variables.index('f_sum')]) == 0.:
                    if beam_edge > 0.5:
                        beam_edge = beam_edge/2.
                    elif vel_edge > min_vel_edge:
                        vel_edge = vel_edge/2.
                        if vel_edge < 1.:
                            vel_edge= 0.
                    else:
                        # if our sources are all close to the edge we check whether there is one which is more than half of the spatial size in the channels it exists
                        for i in range(len(many_sources[0])):
                            cube= float(Configuration['NAXES'][0])*float(Configuration['NAXES'][1])* \
                                    (float(many_sources[Variables.index('z_max')][i])-float(many_sources[Variables.index('z_min')][i]))
                            source_size=(float(many_sources[Variables.index('z_max')][i])-float(many_sources[Variables.index('z_min')][i]))* \
                                        (float(many_sources[Variables.index('y_max')][i])-float(many_sources[Variables.index('y_min')][i]))* \
                                        (float(many_sources[Variables.index('x_max')][i])-float(many_sources[Variables.index('x_min')][i]))
                            #print(source_size,cube)
                            #cube= np.array([float(Configuration['NAXES'][0]),float(Configuration['NAXES'][1]), \
                            #        (float(many_sources[Variables.index('z_max')][i])-float(many_sources[Variables.index('z_min')][i]))])
                            #source_size=np.array([(float(many_sources[Variables.index('z_max')][i])-float(many_sources[Variables.index('z_min')][i])), \
                            #            (float(many_sources[Variables.index('y_max')][i])-float(many_sources[Variables.index('y_min')][i])), \
                            #            (float(many_sources[Variables.index('x_max')][i])-float(many_sources[Variables.index('x_min')][i]))])
                            #print(source_size,cube)
                            if source_size/cube > 0.5:
                                print_log(f'''SOFIA_CATALOGUE: We discarded a very large source, so we will restore is and try for that.
    !!!!!!!!!!!!!!!!!!!!!!!!! This means your original cube is in principle too small!!!!!!!!!!!!!!!!!!!!!!!!
    ''',Configuration['OUTPUTLOG'],screen=True)
                                many_sources[Variables.index('f_sum')][i]=outlist[Variables.index('f_sum')][i]
                        if np.nansum(many_sources[Variables.index('f_sum')]) == 0.:
                            print_log(f'''SOFIA_CATALOGUE:The found sources are too close to the edges of the cube. And not large enough to warrant trying them.
    {'':8s} The edge limits were {beam_edge} beams spatially and {vel_edge} channels.
    ''',Configuration['OUTPUTLOG'],screen=True)
                            raise BadCatalogueError("The found sources are too close to the edges of the cube. And not large enough to warrant trying them.")
                        else:
                            found = True
                else:
                    found = True
            # We need to check we are not throwing away a source that is infinitely brighter
            fluxes = np.array(many_sources[Variables.index('f_sum')],dtype =float)
            if np.any(fluxes == 0.):
                no_edge_fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
                if np.nanmax(no_edge_fluxes) > 10.* np.nanmax(fluxes):
                    if debug:
                        print_log(f'''SOFIA_CATALOGUE: We discarded a very bright source, let's check wether it satisfies our minimum boundaries.
    ''',Configuration['OUTPUTLOG'])
                    index = np.where(np.nanmax(no_edge_fluxes) == no_edge_fluxes)[0][0]
                    edge = check_edge_limits(float(outlist[Variables.index('x_min')][index]),
                                    float(outlist[Variables.index('x_max')][index]),
                                    float(outlist[Variables.index('y_min')][index]),
                                    float(outlist[Variables.index('y_max')][index]),
                                    float(outlist[Variables.index('z_min')][index]),
                                    float(outlist[Variables.index('z_max')][index]),
                                    Configuration,debug = debug,vel_edge=min_vel_edge)

                    if edge:
                        print_log(f'''SOFIA_CATALOGUE: The bright source is very close to limits
    ''',Configuration['OUTPUTLOG'])
                    else:
                        print_log(f'''SOFIA_CATALOGUE: The bright source is acceptable, restoring its flux
    ''',Configuration['OUTPUTLOG'])
                        many_sources  = copy.deepcopy(outlist)
                        fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
            if debug:
                print_log(f'''SOFIA_CATALOGUE: after checking edges we find these fluxes
{'':8s}{many_sources[Variables.index('f_sum')]}
''',Configuration['OUTPUTLOG'])
        else:
            #If we did not run SoFia we simply assume we want to fit the brightest source in the cube
            many_sources  = copy.deepcopy(outlist)
            fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
        outlist = []
        #We want the source with the most total flux.
        index = np.where(np.nanmax(fluxes) == fluxes)[0][0]
        print_log(f'''SOFIA_CATALOGUE: We select the {index} source of this list.
''',Configuration['OUTPUTLOG'])
        outlist = [x[index] for x in many_sources]
        # check that our mask has the selected source
    else:
        outlist = [x[0] for x in outlist]

    check_mask(Configuration,outlist[Variables.index('id')],Fits_Files,debug=debug)
    if debug:
        print_log(f'''SOFIA_CATALOGUE: we found these values
{'':8s}{outlist}
''',Configuration['OUTPUTLOG'])
    return outlist
sofia_catalogue.__doc__ =f'''
 NAME:
    sofia_catalogue

 PURPOSE:
       Read the sofia catalogue and extract several basic parameters from into a list.
       In order to read

 CATEGORY:
    read_functions

 INPUTS:
    Configuration
    Fits_Files = Fits file names

 OPTIONAL INPUTS:
    Variables =['id','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra',\
                        'dec','v_app','f_sum','kin_pa','w50','err_f_sum','err_x','err_y','err_z']
    Variables to read from the catalogue

    header = None
    header of the initial cube to translate units

    debug = False

 OUTPUTS:
    outlist = list with the requested values for the brightest source in the catalogue provided it is not on the edge of the cube.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def sofia_input_catalogue(Configuration,debug=False):
    Catalogue = Proper_Dictionary({})

    Catalogue = {'ENTRIES':  ['ENTRIES','DISTANCE','NUMBER','DIRECTORYNAME','CUBENAME'],\
                'DISTANCE': [],'DIRECTORYNAME': [],\
                'NUMBER': [], 'CUBENAME': []}
    basename = Configuration['CATALOGUE'].split('_cat.txt')[0]
    Variables =['id','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra',\
                    'dec','v_app','f_sum','kin_pa','w50','err_f_sum','err_x','err_y','err_z']
    headerlines=[]
    #Read the sofia catalogue
    with open(Configuration['CATALOGUE']) as sof_cat:
        for line in sof_cat.readlines():
            tmp =line.split()

            if line.strip() == '' or line.strip() == '#':
                headerlines.append(line)
                pass
            elif tmp[0] == '#' and len(tmp) > 1:
                headerlines.append(line)
                if tmp[1].strip().lower() in ['name']:

                    input_columns = [x.strip() for x in tmp[1:]]
                    column_locations = []
                    for col in input_columns:
                        column_locations.append(line.find(col)+len(col))

                    # check that we found all parameters
                    for value in Variables:
                        if value.lower() in input_columns:
                            continue
                        else:
                            print_log(f'''SOFIA_CATALOGUE: We cannot find the required column for {value} in the sofia catalogue.
{"":8s}SOFIA_CATALOGUE: This can happen because a) you have tampered with the sofiainput.txt file in the Support directory,
{"":8s}SOFIA_CATALOGUE: b) you are using an updated version of SoFiA2 where the names have changed and FAT is not yet updated.'
{"":8s}SOFIA_CATALOGUE:    In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
{"":8s}SOFIA_CATALOGUE: c) You are using pre processed SoFiA output of your own and do not have all the output'
{"":8s}SOFIA_CATALOGUE:    Required output is {','.join(Variables)})
''',Configuration['OUTPUTLOG'],screen= True)
                            raise BadCatalogueError("SOFIA_CATALOGUE: The required columns could not be found in the sofia catalogue.")
            else:

                outlist = ['' for x in input_columns]
                for col in input_columns:
                    if input_columns.index(col) == 0:
                        start = 0
                    else:
                        start = column_locations[input_columns.index(col)-1]
                    end = column_locations[input_columns.index(col)]
                    outlist[input_columns.index(col)] = line[start:end].strip()

                if not os.path.exists(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}"):
                    create_directory(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}",f"{Configuration['MAIN_DIRECTORY']}")
                Cube = fits.open(f"{Configuration['MAIN_DIRECTORY']}{basename}_cubelets/{basename}_{outlist[input_columns.index('id')]}_cube.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
                data = Cube[0].data
                hdr = Cube[0].header
                if hdr['NAXIS'] == 4:
                    data = data[0,:,:,:]
                    del hdr['NAXIS4']
                    hdr['NAXIS'] = 3
                # clean the header
                hdr = clean_header(Configuration,hdr,debug=debug)
                hdr['FATNOISE'] = np.mean([np.nanstd(data[0,:,:]),np.nanstd(data[-1,:,:])])
                if hdr['CDELT3'] < -1:
                    raise InputError(f"Your velocity axis is declining this won't work. exiting")
                fits.writeto(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/{basename}_{outlist[input_columns.index('id')]}_FAT.fits",data,hdr,overwrite=True)
                create_directory(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/Sofia_Output",f"{Configuration['MAIN_DIRECTORY']}")
                if not Configuration['SOFIA_DIR']:
                    Configuration['SOFIA_DIR']=f"{Configuration['MAIN_DIRECTORY']}{basename}_cubelets/"
                Configuration['FITTING_DIR']=f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/"
                Configuration['SOFIA_BASENAME'] = f"{basename}_{outlist[input_columns.index('id')]}"
                Configuration['BASE_NAME'] =  f"{basename}_{outlist[input_columns.index('id')]}_FAT"
                copy_homemade_sofia(Configuration,no_cat=True,debug=False)
                Catalogue['DISTANCE'].append(float(-1))
                Catalogue['NUMBER'].append(f"{outlist[input_columns.index('id')]}")
                Catalogue['DIRECTORYNAME'].append(f"{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}")
                Catalogue['CUBENAME'].append(f"{basename}_{outlist[input_columns.index('id')]}")
                with open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_cat.txt",'w') as cat:
                    for hline in headerlines:
                        cat.write(f"{hline}\n")
                    cat.write(line)
    Configuration['SOFIA_BASENAME'] = None
    Configuration['BASE_NAME'] = 'Unset'
    Configuration['FITTING_DIR'] = 'Unset'
    return Catalogue
sofia_input_catalogue.__doc__=f'''
NAME:
    sofia_input_catalogue

 PURPOSE:
     The input catalogue is a sofia 2 catalogue that should be arranged into directories and fitted.

 CATEGORY:
    read_functions

 INPUTS:
    filename = the name of the sofia catalogue

 OPTIONAL INPUTS:

 OUTPUTS:
    catalogue = the sofia catalogue configured to work as a FAT catalogue

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def sofia_template(debug = False):
    with import_res.open_text(templates, 'sofia_template.par') as tmp:
        template = tmp.readlines()
    result = {}
    counter = 0
    counter2 = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip())
        if key == '':
            result[f"EMPTY{counter}"] = line
            counter += 1
        elif key[0] == '#':
            result[f"HASH{counter2}"] = line
            counter2 += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result
sofia_template.__doc__ ='''
 NAME:
    sofia_template

 PURPOSE:
    Read the sofia2 file in Templates into a dictionary

 CATEGORY:
    read_functions

 INPUTS:

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    split, strip, open

 NOTE:
'''

def tirific_template(filename = '', debug = False):
    if filename == '':
        with import_res.open_text(templates, 'template.def') as tmp:
            template = tmp.readlines()
    elif filename == 'Installation_Check':
        from pyFAT_astro import Installation_Check as IC
        with import_res.open_text(IC, 'ModelInput.def') as tmp:
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

    debug =False

 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      split, strip, open

 NOTE:
'''
