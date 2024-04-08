# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in FAT to read input files

from pyFAT_astro.Support.fits_functions import check_mask,create_fat_cube
from pyFAT_astro.Support.fat_errors import BadCatalogueError, InputError, \
    BadSourceError
from pyFAT_astro.Support.log_functions import print_log,update_statistics
from pyFAT_astro.Support.modify_template import fix_sbr

import pyFAT_astro.Support.support_functions as sf

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
import copy
import numpy as np
try:
    from importlib.resources import open_text as pack_open_txt
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    from importlib_resources import open_text as pack_open_txt



from pyFAT_astro import Templates as templates
#Function to read a FAT input Catalogue


def catalogue(filename,split_char='|'):
    Catalogue = sf.Proper_Dictionary({})
    with open(filename,'r') as tmpfile:
        firstline = tmpfile.readline()
        all_columns_check = False
        required_columns= ['ID','DISTANCE','DIRECTORYNAME','CUBENAME']

        while not all_columns_check:
            input_columns = [x.strip().upper() for x in firstline.split(split_char)]
            Catalogue['ENTRIES'] = ['ENTRIES']
            Catalogue['ENTRIES'].extend(input_columns)
            required_columns= ['ID','DISTANCE','DIRECTORYNAME','CUBENAME']

            for key in required_columns:
                if key not in Catalogue['ENTRIES']:
                    if split_char == '|':
                        print(f'Key {key} not found')
                        split_char=' '
                        all_columns_check = False
                        break
                    else:
                        raise BadCatalogueError(f'We can not find the column for {key} in your input catalogue')
                else:
                    all_columns_check = True
                    continue



        for key in input_columns:
            Catalogue[key] = []


        for line in tmpfile.readlines():
            input = [x.strip() for x  in line.split(split_char)]
            if len(input) == len(input_columns):
                for i,key in enumerate(input_columns):
                    if key == 'DISTANCE':
                        Catalogue[key].append(float(input[i]))
                    else:
                        Catalogue[key].append(input[i])
            else:
                print(f'READ_CATALOGUE: Your line "{line}" in the input catalogue does not have correct number of columns, skipping it')

    #if 'NUMBER' in Catalogue['ENTRIES']:
    #    Catalogue['NUMBER'] = np.array(Catalogue['NUMBER'],dtype=int)

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


 OUTPUTS:
    Catalogue = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_edge_limits(xmin,xmax,ymin,ymax,zmin,zmax,Configuration ,\
        beam_edge = 0.5, vel_edge = 0.5 ):
    diff = np.array([xmin,abs(xmax - Configuration['NAXES'][0]),
                     ymin,abs(ymax - Configuration['NAXES'][1])],dtype = float)
    print_log(f'''CHECK_EDGE_LIMIT: We find these differences and this edge size {beam_edge*Configuration['BEAM_IN_PIXELS'][0]}
{'':8s} diff  = {diff}
''',Configuration, case=['debug_start'])
    if np.where(diff < beam_edge*Configuration['BEAM_IN_PIXELS'][0])[0].size:
        return True
    diff = np.array([zmin,abs(zmax-Configuration['NAXES'][2])],dtype=float)
    print_log(f'''CHECK_EDGE_LIMIT: And for velocity edge =  {vel_edge}
{'':8s} diff  = {diff}
''',Configuration, case=['debug_add'])
    if np.where(diff < vel_edge)[0].size:
        print_log(f"CHECK_EDGE_LIMIT: On the edge. \n",Configuration,\
            case=['verbose'])
        return True
    else:
        print_log(f"CHECK_EDGE_LIMIT: Off the edge. \n",Configuration,\
            case=['verbose'])
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
def check_source_brightness(Configuration,Fits_Files,moment= False):

    if not moment:
        Cube = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"\
            ,uint = False, do_not_scale_image_data=True,ignore_blank = True,\
             output_verify= 'ignore')
        Mask = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MASK']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, \
            output_verify= 'ignore')
        data = Cube[0].data[Mask[0].data > 0.5]
        snr = np.sort(data/Configuration['NOISE'])[::-1]
        checking= 'Cube'
        Cube.close()
        Mask.close()
        case = ['debug_add']
        fact = 1.
        fact_mean = 1.
    else:
        Moment = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT0']}")
    #model_mom0 = sf.remove_inhomogeneities(Configuration,model_mom0,inclination=float(current[1][0]), pa = float(current[2][0]), center = [current[3][0],current[4][0]])
        Channels = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['CHANNEL_MAP']}")
        noisemap = np.sqrt(Channels[0].data)*Configuration['NOISE']*\
            Configuration['CHANNEL_WIDTH']
        data = Moment[0].data
        snr= np.sort(Moment[0].data[Moment[0].data > 0]/\
            noisemap[Moment[0].data > 0])[::-1]
        checking= 'Moment 0 Map'
        case = ['main','screen']
        fact=2.
        fact_mean =  3./Configuration['SOURCE_MEAN_SNR']
    #The following should only be run if we ran sofia
    Too_Faint = False
    #Check that the source is bright enough
    Max_SNR = snr[int(len(snr)*Configuration['SOURCE_MAX_FRACTION']*fact)]

    if Max_SNR < Configuration['SOURCE_MAX_SNR']:
        print_log(f'''CHECK_SOURCE_BRIGHTNESS: The SNR of brightest {Configuration['SOURCE_MAX_FRACTION']*100.*fact}% of selected pixels does not always exceed {Configuration['SOURCE_MAX_SNR']} in the {checking}.
At the cutoff the SNR = {Max_SNR}.
''', Configuration,case = case)
        Too_Faint = True
            #raise BadSourceError(f"The SNR of brightest {Configuration['SOURCE_MAX_FRACTION']*100.}% of selected pixels does not always exceed {Configuration['SOURCE_MAX_SNR']}. Aborting.")
    else:
        print_log(f'''CHECK_SOURCE_BRIGHTNESS: The SNR of brightest {Configuration['SOURCE_MAX_FRACTION']*100.*fact}% of selected pixels always exceeds {Configuration['SOURCE_MAX_SNR']} in the {checking}.
At the cutoff the SNR = {Max_SNR}.
''', Configuration)
    Mean_SNR = np.nanmean(snr)
    if Mean_SNR < Configuration['SOURCE_MEAN_SNR']*fact_mean:
        print_log(f'''CHECK_SOURCE_BRIGHTNESS: The mean SNR of the pixels in the mask is {Mean_SNR} in the {checking}.
This is less than {Configuration['SOURCE_MEAN_SNR']*fact_mean} and therefore deemed too faint.
''', Configuration,case= case)
        Too_Faint = True
            #raise BadSourceError(f"The mean SNR of the pixels in the mask is {Mean_SNR}. This is too faint.")
    else:
        print_log(f'''CHECK_SOURCE_BRIGHTNESS: The mean SNR of the pixels in the mask is {Mean_SNR} in the {checking}.
''', Configuration)

    return Too_Faint,Max_SNR

check_source_brightness.__doc__ =f'''
 NAME:
    check_source_brightness

 PURPOSE:
    Check the cube or the moment 0 map whether the source is bright enogh the process

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Fat Configuration
    Fits_Files = The fits files used


 OPTIONAL INPUTS:
    moments  = False
        flag to check the moment 0 map instead of the cube
 OUTPUTS:
    Boolean which is true if the source is too faint

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    abs, np.where, np.array, print_log,

 NOTE:
'''



def extract_vrot(Configuration,map ,angle,center):
    print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
{'':8s} PA= {angle}
{'':8s} center= {center}
''',Configuration, case=['debug_start'])
   
    maj_profile,maj_axis,maj_resolution = sf.get_profile(Configuration,map,angle,center=center)
   
    # We should base extracting the RC on where the profile is negative and positive to avoid mistakes in the ceneter coming through
    neg_index = np.where(maj_profile < 0.)[0]
    pos_index = np.where(maj_profile > 0.)[0]

    print_log(f'''EXTRACT_VROT: The resolution on the extracted axis
{'':8s} resolution = {maj_resolution}
''',Configuration,case= ['debug_add'])
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
    print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
''',Configuration,case= ['debug_add'])
    ring_size_req = Configuration['BEAM_IN_PIXELS'][0]/maj_resolution
    print_log(f'''EXTRACT_VROT: We need a rings size of
{'':8s} ringsize= {ring_size_req}
{'':8s} because bmaj in pixels  ={Configuration['BEAM_IN_PIXELS'][0]}  and the resolution of the profile = {maj_resolution} pixels
''',Configuration,case= ['debug_add'])
    profile = np.array(avg_profile[0::int(ring_size_req)],dtype=float)

    try:
        profile[np.isnan(profile)] = profile[~np.isnan(profile)][-1]
    except IndexError:
        profile = []

    if len(profile) < 1:
        if len(avg_profile) > 0.:
            profile = np.array([np.max(avg_profile)],dtype=float)
        else:
            profile = np.array([Configuration['CHANNEL_WIDTH']*2.],dtype=float)
    print_log(f'''EXTRACT_VROT: Unlimited profile
{'':8s} unlimited  RC= {profile}
''',Configuration,case= ['debug_add'])
    profile[profile > 300.] = 300.
    profile[profile < Configuration['CHANNEL_WIDTH']*2.] = Configuration['CHANNEL_WIDTH']*2.


    #profile[0] = 0.
    print_log(f'''EXTRACT_VROT: Constructing the final RC
{'':8s} initial RC= {profile}
{'':8s} at step width= {ring_size_req}
''',Configuration,case= ['debug_add'])
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

    fit_type = 'Undefined'

 OUTPUTS:
    profile = The RC at given ring locations averaged over the approaching and receding side

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_DHI(Configuration, Model='Finalmodel'):
    #Get the sbrs
    radi,sbr,sbr_2,systemic = sf.load_tirific(Configuration,\
        f"{Configuration['FITTING_DIR']}{Model}/{Model}.def",\
        Variables = ['RADI','SBR','SBR_2','VSYS'],array=True)
    #convert to solar_mass/pc^2
    sbr_msolar = sf.columndensity(Configuration,sbr*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = sf.columndensity(Configuration,sbr_2*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    # interpolate these to ~0.1 beam steps
    new_radii = np.linspace(0,radi[-1],int(radi[-1]/(0.1*Configuration['BEAM'][0])+1))
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


    Model = 'Finalmodel'
    location of the def file to get DHI from. it should be in the fitting dir in the {{Model}}/{{Model}}.def

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_totflux(Configuration,map_name):
    image = fits.open(f"{Configuration['FITTING_DIR']}{map_name}")
    #We are taking these from the moment map so we have to divide out the km/s
    flux = float(np.nansum(image[0].data)/Configuration['BEAM_IN_PIXELS'][2]/Configuration['CHANNEL_WIDTH'])
    #Should this not have an additional channel width parameter
    error = np.sqrt((np.where(image[0].data> 0.)[0].size)/Configuration['BEAM_IN_PIXELS'][2])*Configuration['NOISE']
    image.close()
    return [flux,error]
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


 OUTPUTS:
    total flux in the map in Jy*km/s

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# Function to get the PA and inclination from the moment 0 for initial estimates
def guess_orientation(Configuration,Fits_Files, vsys = -1 ,center = None, \
            smooth = False):
    #open the moment 0

    print_log(f'''GUESS_ORIENTATION: starting extraction of initial parameters.
''',Configuration, case=['debug_start','verbose'])
    update_statistics(Configuration, message= "Starting the guess orientation run")

    Image = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT0']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = Image[0].data
    hdr = Image[0].header
    mom0 = copy.deepcopy(Image)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        mom0_wcs = WCS(hdr)
    Image.close()

    if not center:
        center = [hdr['NAXIS1']/2.-1,hdr['NAXIS2']/2.-1]
    Image = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['CHANNEL_MAP']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    noise_map = np.sqrt(Image[0].data)*Configuration['NOISE']*Configuration['CHANNEL_WIDTH']

    #these guesses get thrown off for very large galaxies so smooth when we have those:
    if smooth:
        sigma = sf.set_limits(Configuration['BEAM_IN_PIXELS'][0]*Configuration['MAX_SIZE_IN_BEAMS']/12.*2.,Configuration['BEAM_IN_PIXELS'][0],Configuration['BEAM_IN_PIXELS'][0]*5)
        print_log(f'''GUESS_ORIENTATION: We are smoothing the maps with {sigma} pixels
''',Configuration,case= ['debug_add'])
        tmp =  ndimage.gaussian_filter(map, sigma=(sigma, sigma), order=0)
        tmp[map <= 0.] = 0.
        map = copy.deepcopy(tmp)
        if Configuration['DEBUG']:
            fits.writeto(f"{Configuration['LOG_DIRECTORY']}smooth_mom_map.fits",map,hdr,overwrite = True)
        tmp = ndimage.gaussian_filter(noise_map, sigma=(sigma, sigma), order=0)
        tmp[noise_map <= 0.] = 0.
        noise_map = copy.deepcopy(tmp)
        if Configuration['DEBUG']:
            fits.writeto(f"{Configuration['LOG_DIRECTORY']}smooth_noise_map.fits",noise_map,hdr,overwrite = True)


    SNR = np.nanmean(map[noise_map > 0.]/noise_map[noise_map > 0.])
    Configuration['SNR'] = SNR
    noise_hdr = Image[0].header
    Image.close()
    noise_map [0. > map ] =0.
    median_noise_in_map = np.nanmedian(noise_map[noise_map > 0.])
    minimum_noise_in_map = np.nanmin(noise_map[noise_map > 0.])

    # if we have extended low level wings that are significant we do not really want to include then
    map[0.5*minimum_noise_in_map > noise_map] = 0.
    #Also remove negative values
    #map[0. > map ] =0


    mom0[0].data= map


    scale_factor = sf.set_limits(SNR/3.*minimum_noise_in_map/median_noise_in_map, 0.05, 1.)
    print_log(f'''GUESS_ORIENTATION: We find SNR = {SNR} and a scale factor {scale_factor} and the noise median {median_noise_in_map}
{'':8s} minimum {minimum_noise_in_map}
''',Configuration,case= ['debug_add'])
    '''
    if np.mean(Configuration['SIZE_IN_BEAMS']) >10:
        beam_check=[Configuration['BEAM_IN_PIXELS'][0],Configuration['BEAM_IN_PIXELS'][0]/2.]
    else:
    '''
    beam_check=[Configuration['BEAM_IN_PIXELS'][0]/2.]
    center_stable = False
    checked_center = False
    center_counter = 0.
    original_center = copy.deepcopy(center)
    print_log(f'''GUESS_ORIENTATION: Looking for the center, pa and inclination
''',Configuration, case= ['verbose'])

    update_statistics(Configuration, message= "Starting the initial search for the pa, inclination and center.")

    while not center_stable:
        inclination_av, pa_av, maj_extent_av =\
            sf.get_inclination_pa(Configuration, mom0, center, cutoff = scale_factor* median_noise_in_map)
        inclination_av = [inclination_av]
        int_weight = [2.]
        pa_av = [pa_av]
        maj_extent_av = [maj_extent_av]
        print_log(f'''GUESS_ORIENTATION: From the the initial guess with center {center}.
{'':8s} We get pa = {pa_av}, inclination = {inclination_av}, maj_extent_av {maj_extent_av}
''',Configuration,case= ['debug_add'])
        for mod in beam_check:

            for i in [[-1,-1],[-1,1],[1,-1],[1,1]]:
                center_tmp = [center[0]+mod*i[0],center[1]+mod*i[1]]

                print_log(f'''GUESS_ORIENTATION: Checking at location RA = {center_tmp[0]} pix, DEC = {center_tmp[1]} pix
''',Configuration,case= ['debug_add'])

                inclination_tmp, pa_tmp, maj_extent_tmp= \
                    sf.get_inclination_pa(Configuration, mom0, center_tmp,\
                    cutoff = scale_factor* median_noise_in_map,)
                inclination_av.append(inclination_tmp)
                pa_av.append(pa_tmp)
                int_weight.append(mod/beam_check[0]*0.25)
                maj_extent_av.append(maj_extent_tmp)

        if np.isnan(np.array(inclination_av,dtype=float)).all():
            print_log(f'''GUESS_ORIENTATION: We are unable to find an  inclination.
''',Configuration,case= ['main','screen'])
            raise BadSourceError(f'We are unable to find an initial inclination.')

        int_weight = np.array(int_weight)
        weight = np.array([1./x[1] for x in inclination_av],dtype= float)*int_weight

        inclination = np.array([np.nansum(np.array([x[0] for x in inclination_av],dtype=float)*weight)/np.nansum(weight),\
                                np.nansum(np.array([x[1] for x in inclination_av],dtype=float)*weight)/np.nansum(weight)],dtype=float)


        if np.isnan(np.array(inclination_av,dtype=float)).all():
            print_log(f'''GUESS_ORIENTATION: We are unable to find an PA.
''',Configuration,case=['main','screen'])
            raise BadSourceError(f'We are unable to find an initial PA.')
        weight = np.array([1./x[1] for x in pa_av],dtype= float)
        print_log(f'''GUESS_ORIENTATION: We find these pa and inclination
{'':8s} pa = {' '.join([f'{float(x[0]):.2f}' for x in pa_av])}
{'':8s} inclination = {' '.join([f'{float(x[0]):.2f}' for x in inclination_av])}
{'':8s} int_weights = {' '.join([f'{float(x):.2f}' for x in int_weight])}
{'':8s} weights = {' '.join([f'{float(x):.2f}' for x in weight])}
''',Configuration,case= ['debug_add'])
        pa = np.array([np.nansum(np.array([x[0] for x in pa_av],dtype=float)*weight)/np.nansum(weight),\
                                np.nansum(np.array([x[1] for x in pa_av],dtype=float)*weight)/np.nansum(weight)],dtype=float)

        maj_extent= np.nansum(maj_extent_av*weight)/np.nansum(weight)

        if center_counter == 0:
            original_pa = copy.deepcopy(pa)
            original_maj_extent = copy.deepcopy(maj_extent)
            original_inclination = copy.deepcopy(inclination)

        if center_counter > 0 and not any([np.isfinite(maj_extent),np.isfinite(pa[0]),np.isfinite(inclination[0])]):
            pa=original_pa
            inclination=original_inclination
            maj_extent= original_maj_extent
        #    center =original_center
        #    center_stable=True
        #    break
        # For very small galaxies we do not want to correct the extend
        if maj_extent/Configuration['BEAM'][0]/3600. > 3.:
            maj_extent = maj_extent+(Configuration['BEAM'][0]/3600.*0.2/scale_factor)

        print_log(f'''GUESS_ORIENTATION: From the maps we find
{'':8s} inclination = {inclination}
{'':8s} pa = {pa}
{'':8s} size in beams = {maj_extent/(Configuration['BEAM'][0]/3600.)}
''',Configuration,case= ['debug_add'])
                #map[3*minimum_noise_in_map > noise_map] = 0.
        # From these estimates we also get an initial SBR
        maj_profile,maj_axis,maj_resolution = sf.get_profile(Configuration,map, pa[0],center)
     
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
        print_log(f'''GUESS_ORIENTATION:'BMAJ in pixels, center of the profile, current center, difference between pos and neg
{'':8s}{Configuration['BEAM_IN_PIXELS'][0]} {center_of_profile} {center} {diff}
''',Configuration,case= ['debug_add'])

        # if the center of the profile is more than half a beam off from the Sofia center let's see which on provides a more symmetric profile
        #if False:
        if (abs(center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution) > Configuration['BEAM_IN_PIXELS'][0]*0.5 \
            or abs(center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution) > Configuration['BEAM_IN_PIXELS'][0]*0.5) and SNR > 3. and not checked_center:
            print_log(f'''GUESS_ORIENTATION: The SoFiA center and that of the SBR profile are separated by more than half a beam.
{'':8s}GUESS_ORIENTATION: Determining the more symmetric profile.
''',Configuration,case= ['debug_add'])
            # let's check against a central absorption

            cube = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}",\
                    uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
            buffer = int(round(np.mean(Configuration['BEAM_IN_PIXELS'][:2])/2.))
            central = cube[0].data[:,int(round(center[1]-buffer)):int(round(center[1]+buffer)),int(round(center[0]-buffer)):int(round(center[0]+buffer))]

            if np.count_nonzero(np.isnan(central))/central.size < 0.1:
                center,checked_center,center_stable = sf.get_new_center(Configuration,map,inclination=inclination[0],pa=pa[0],noise=median_noise_in_map)
                center_counter += 1
                print_log(f'''GUESS_ORIENTATION: We calculated a new center {center}.
''',Configuration,case= ['debug_add'])
            else:
                center_stable = True
                print_log(f'''GUESS_ORIENTATION: There appears to be a central absorption source. We are relying on sofia
''',Configuration,case= ['debug_add'])

        else:
            print_log(f'''GUESS_ORIENTATION: The previous center and that of the SBR profile are not separated by more than half a beam.
{'':8s}GUESS_ORIENTATION: Keeping the last center
''',Configuration,case= ['debug_add'])
            center_stable = True

    print_log(f'''GUESS_ORIENTATION: Looking for the Initial surface brightness profile.
''',Configuration, case= ['verbose'])
    update_statistics(Configuration, message= "Starting the initial search for the SBR and VROT.")
    Configuration['SIZE_IN_BEAMS'] = np.full(2,sf.set_limits(maj_extent/(Configuration['BEAM'][0]/3600.),1.0,Configuration['MAX_SIZE_IN_BEAMS']))

    ring_size_req = Configuration['BEAM_IN_PIXELS'][0]/maj_resolution*Configuration['RING_SIZE']

    SBR_initial = avg_profile[0::int(ring_size_req)]/(np.pi*Configuration['BEAM'][0]*Configuration['BEAM'][1]/(4.*np.log(2.))) # Jy*km/s
    #deproject the SBR
    SBR_initial = SBR_initial*np.cos(np.radians(inclination[0]))
    SBR_initial =np.hstack((SBR_initial[0],SBR_initial,SBR_initial[-1]))
    # unclear why we have this extra correction on the inner 3 points
    SBR_initial[0:3] = SBR_initial[0:3] * (1.2 -float(inclination[0])/90.)
    # make sure it is the right size
    number_of_rings = sf.calc_rings(Configuration)
    SBR_initial = list(SBR_initial)
    while len(SBR_initial) < number_of_rings:
         SBR_initial.append(SBR_initial[-1])
    while len(SBR_initial) > number_of_rings:
         SBR_initial.pop()
    SBR_initial = np.array(SBR_initial,dtype=float)
    #We need to know which is the approaching side and which is receding


    print_log(f'''GUESS_ORIENTATION: Looking for the Initial Rotation Curve.
''',Configuration,case= ['verbose'])
    Image = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT1']}",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = copy.deepcopy(Image[0].data)

    hdr = Image[0].header
    print_log(f'''GUESS_ORIENTATION: This is the amount of values we find initially {len(map[noise_map > 0.])}
''',Configuration,case= ['debug_add'])
    #Image.close()
    #map[3*minimum_noise_in_map > noise_map] = float('NaN')

    map[3.*minimum_noise_in_map > noise_map] = float('NaN')
    print_log(f'''GUESS_ORIENTATION: This is the amount of values we find after blanking low SNR {len(map[~np.isnan(map)])}
''',Configuration,case= ['debug_add'])
    if vsys == -1 or center_counter > 0.:
        #As python is utterly moronic the center goes in back wards to the map
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore",message="Mean of empty slice"\
                            ,category=RuntimeWarning)
            map_vsys = np.nanmean(map[int(round(center[1]-buffer)):int(round(center[1]+buffer)),int(round(center[0]-buffer)):int(round(center[0]+buffer))])
        if np.mean(Configuration['SIZE_IN_BEAMS']) < 10.:
            map_vsys = (vsys+map_vsys)/2.
    else:
        map_vsys = vsys


    if len(map[~np.isnan(map)]) < 10 or np.isnan(map_vsys):
        no_values  = True
        noise_level = 2.5
        sigma = [0.5,0.5]
        while no_values:
            print_log(f'''GUESS_ORIENTATION: We smooth the velocity field to use a lower Noise threshold
''',Configuration,case= ['debug_add'])
            tmp = copy.deepcopy(Image[0].data)
            tmp =  ndimage.gaussian_filter(Image[0].data, sigma=(sigma[1], sigma[0]), order=0)
            tmp[noise_level*minimum_noise_in_map > noise_map] =  float('NaN')
            #We need the cnter to be found else we get a crash
            if vsys == -1 or center_counter > 0.:
                #As python is utterly moronic the center goes in back wards to the map
                with warnings.catch_warnings():
                    warnings.filterwarnings("error",message="Mean of empty slice"\
                            ,category=RuntimeWarning)
                    try:
                        map_vsys = np.nanmean(tmp[int(round(center[1]-buffer)):int(round(center[1]+buffer)),int(round(center[0]-buffer)):int(round(center[0]+buffer))])
                    except RuntimeWarning:
                        map_vsys = vsys

                if np.mean(Configuration['SIZE_IN_BEAMS']) < 10.:
                    map_vsys = (vsys+map_vsys)/2.


            print_log(f'''GUESS_ORIENTATION: We used the threshold {noise_level} and sigma = {sigma}
{'':8s}We find the len {len(tmp[~np.isnan(tmp)])} and the central velocity {map_vsys}
''',Configuration,case= ['debug_add'])
            if len(tmp[~np.isnan(tmp)]) < 10 or np.isnan(map_vsys):
                sigma = [sigma[0]+0.5,sigma[0]+0.5]
                noise_level -= 0.5
            else:
                no_values = False
                map = tmp
            if noise_level < 1. and no_values:
                VROT_initial = [0]
                SBR_initial = [0]
                Image.close()
                return np.array(pa,dtype=float),np.array(inclination,dtype=float),SBR_initial,center[0],center[1],float('NaN'),VROT_initial

    Image.close()
    #we want to make sure the SBr is gaussian
    # We want to make this is sure this is a decent profile and thus we fit it with the gaussian
    #Configuration['SIZE_IN_BEAMS'] = len(SBR_initial)-1
    format = sf.set_format('SBR')
    Temp_template = {'SBR':' '.join([f'{float(x):{format}}' for x in SBR_initial])\
                        ,'SBR_2':' '.join([f'{float(x):{format}}' for x in SBR_initial]),\
                        'VSYS':' '.join([f'{float(x):{sf.set_format("VSYS")}}' for x in [vsys,vsys]])}
    SBR_initial = fix_sbr(Configuration,Temp_template,\
                        smooth = True,initial=True)
    SBR_initial = np.array([float(x) for x in Temp_template['SBR'].split()],dtype=float)

    noise_map = []

    # First we look for the kinematics center

    vel_pa = sf.get_vel_pa(Configuration,map,center=center)

    maj_profile,maj_axis,maj_resolution = sf.get_profile(Configuration,map,pa[0], center=center)
    zeros = np.where(maj_profile == 0.)[0]
    maj_profile[zeros] = float('NaN')
    if all(np.isnan(maj_profile)):
        print_log(f'''GUESS_ORIENTATION: The RC extracted from the VF is all NaN's, this means something has gone very wrong.
{'':8s} Raising an error.
''' , Configuration,case= ['main','screen'])
        VROT_initial = [0]
        SBR_initial = [0]
        Image.close()
        return np.array(pa,dtype=float),np.array(inclination,dtype=float),SBR_initial,center[0],center[1],float('NaN'),VROT_initial

    loc_max = np.mean(maj_axis[np.where(maj_profile == np.nanmax(maj_profile))[0]])
    loc_min = np.mean(maj_axis[np.where(maj_profile == np.nanmin(maj_profile))[0]])

    if loc_max > loc_min:
        pa[0] = pa[0]+180
        print_log(f'''GUESS_ORIENTATION: We have modified the pa by 180 deg as we found the maximum velocity west of the minimum.
''' , Configuration)

    print_log(f'''GUESS_ORIENTATION: this is the pa {pa} and the vel_pa {vel_pa}
''' , Configuration,case= ['debug_add'])

    if abs(pa[0]-vel_pa[0]) > 300.:
        if pa[0] > 180.:
            vel_pa[0] += 360.
        else:
            vel_pa[0] -= 360.
    if abs(pa[0]-vel_pa[0]) > 25. and ~np.isnan(vel_pa[0])  :
        # if the vel pa has a ridiculous error we use the normal pa
        morph_pa=copy.deepcopy(pa)
        if vel_pa[1] < 60.:
            pa=vel_pa
        #if the galaxy has a low inclination we do allow an error of the difference
        if inclination[0] < 35.:
            pa[1]=abs(morph_pa[0]-vel_pa[0])
    else:
        if ~np.isnan(vel_pa[0]):
            pa = [np.nansum([vel_pa[0]/vel_pa[1],\
                pa[0]/pa[1]])/np.nansum([1./vel_pa[1],1./pa[1]]),\
                2.*1./np.sqrt(np.nansum([1./vel_pa[1],1./pa[1]]))]
    print_log(f'''GUESS_ORIENTATION: this is the final pa {pa}
''' , Configuration,case= ['debug_add'])

    if smooth:
        buffer = int(Configuration['BEAM_IN_PIXELS'][0]*5.)
    else:
        buffer = int(round(np.mean(Configuration['BEAM_IN_PIXELS'][:2])/2.))
    print_log(f'''GUESS_ORIENTATION: We start with vsys = {vsys:.2f} km/s
{'':8s}GUESS_ORIENTATION:We subtract {map_vsys:.2f} km/s from the moment 1 map to get the VROT
''' , Configuration,case= ['debug_add'])

    map = map  - map_vsys
    VROT_initial = extract_vrot(Configuration,map ,pa[0],center)
    min_RC_length= len(VROT_initial)
    if pa[1] < 10.:
        for x in [pa[0]-pa[1],pa[0]-pa[1]/2.,pa[0]+pa[1]/2.,pa[0]+pa[1]]:
            tmp  = extract_vrot(Configuration,map ,x,center)
            if len(tmp) > 0.:
                RC_length = len(tmp)
                if RC_length < min_RC_length:
                    min_RC_length = RC_length
                print_log(f'''GUESS_ORIENTATION: We found the lengths for the angled rotation curves:
{'':8s}GUESS_ORIENTATION: RC length = {RC_length}, min_RC length = {min_RC_length}, Vrot ini = {len(VROT_initial)}, tmnp = {len(tmp)}
{'':8s}GUESS_ORIENTATION: Vrot {VROT_initial}
{'':8s}GUESS_ORIENTATION: tmp {tmp}
''',Configuration,case= ['debug_add'])
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


    print_log(f'''GUESS_ORIENTATION: We found the following initial rotation curve:
{'':8s}GUESS_ORIENTATION: RC = {VROT_initial}
''',Configuration,case= ['debug_add'])
    return np.array(pa,dtype=float),np.array(inclination,dtype=float),SBR_initial,center[0],center[1],map_vsys,VROT_initial
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
def load_basicinfo(Configuration,filename, Variables = None, unpack = True):
    if Variables is None:
        Variables = ['RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask','Tot FLux','D_HI','Distance','HI_Mass' ,'D_HI_kpc' ]
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
            RA, DEC = sf.convertRADEC(Configuration,tmp_str[0],tmp_str2[0],invert=True)
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



 OUTPUTS:
    outputarray = Array with values of all requested parameters

 OPTIONAL OUTPUTS:
       -

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def read_cube(Configuration,cube):
    cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}{cube}")
    Configuration['NOISE'] = cube_hdr['FATNOISE']
    Configuration['CHANNEL_WIDTH'] = cube_hdr['CDELT3']/1000.
    Configuration['PIXEL_SIZE'] = np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])])
    Configuration['NAXES'] = [cube_hdr['NAXIS1'],cube_hdr['NAXIS2'], cube_hdr['NAXIS3']]

    if np.sum(Configuration['VROT_INPUT_BOUNDARY']) == 0.:
        sf.set_boundaries(Configuration,'VROT',Configuration['CHANNEL_WIDTH'],600.,input=True)
    if np.sum(Configuration['SDIS_INPUT_BOUNDARY']) == 0.:
        sf.set_boundaries(Configuration,'SDIS',Configuration['CHANNEL_WIDTH']/3.,25.,input=True)

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

read_cube.__doc__ =f'''
 NAME:
    read_cube

 PURPOSE:
    Load values from the cube header into the Configuration

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    cube = cube to read values from

 OPTIONAL INPUTS:



 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:

'''
def check_parameters(Configuration,Variables,input_columns):
     # check that we found all parameters
    velocity ='v_app'
    for value in Variables:
        trig = False
        if value.lower() in input_columns:
            continue
        elif value.lower() == 'v_app':
            if 'v_rad' in input_columns:
                Variables[Variables.index('v_app')]='v_rad'
                velocity = 'v_rad'
                continue
            elif  'v_opt' in input_columns:
                Variables[Variables.index('v_app')]='v_opt'
                velocity = 'v_opt'
                continue
            else:
                trig = True    
        else:
            trig = True

        if trig:
            print_log(f'''SOFIA_CATALOGUE: We cannot find the required column for {value} in the sofia catalogue.
{"":8s}SOFIA_CATALOGUE: This can happen because a) you have tampered with the sofiainput.txt file in the Support directory,
{"":8s}SOFIA_CATALOGUE: b) you are using an updated version of SoFiA2 where the names have changed and FAT is not yet updated.'
{"":8s}SOFIA_CATALOGUE:    In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
{"":8s}SOFIA_CATALOGUE: c) You are using pre processed SoFiA output of your own and do not have all the output'
{"":8s}SOFIA_CATALOGUE:    Required output is {','.join(Variables)})
''',Configuration,case= ['main','screen'])
            raise BadCatalogueError("SOFIA_CATALOGUE: The required columns could not be found in the sofia catalogue.")
    return velocity    
check_parameters.__doc__ =f'''
 NAME:
    check_parameters(Variables,input_columns)

 PURPOSE:
    check wether all  variables are in the sofia catalogue
 
 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Variable = Reguired variable
    input_columns = found columns
    

 OPTIONAL INPUTS:



 OUTPUTS:
    velocity = the velocity found in the catalogue
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:

'''


# function to read the sofia catalogue
def sofia_catalogue(Configuration,Fits_Files, Variables = None ,no_edge_limit=False):
    strip_npix= False
    if Variables is None:
        Variables =['id','x','x_min','x_max','y','y_min','y_max','z','z_min',\
                    'z_max','ra','dec','v_app','f_sum','kin_pa','w50',\
                    'err_f_sum','err_x','err_y','err_z','rms','n_pix']
        strip_npix =True
    print_log(f'''SOFIA_CATALOGUE: Reading the source from the catalogue.
''',Configuration,case= ['debug_start'])
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
                    velocity = check_parameters(Configuration,Variables,input_columns) 
                    
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
        print_log(f'''SOFIA_CATALOGUE: Multiple sources were found we will try to select the correct one.
''',Configuration,case= ['debug_add'])
        if Configuration['SOFIA_RAN'] and not no_edge_limit:
            found = False
            beam_edge=2.
            if Configuration['VEL_SMOOTH_EXTENDED'] or Configuration['CHANNEL_DEPENDENCY'].lower() == 'hanning':
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
                                    Configuration,beam_edge = beam_edge, vel_edge= vel_edge)
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

                            if source_size/cube > 0.5:
                                print_log(f'''SOFIA_CATALOGUE: We discarded a very large source, so we will restore is and try for that.
!!!!!!!!!!!!!!!!!!!!!!!!! This means your original cube is in principle too small!!!!!!!!!!!!!!!!!!!!!!!!
''',Configuration)
                                many_sources[Variables.index('f_sum')][i]=outlist[Variables.index('f_sum')][i]
                        if np.nansum(many_sources[Variables.index('f_sum')]) == 0.:
                            print_log(f'''SOFIA_CATALOGUE:The found sources are too close to the edges of the cube. And not large enough to warrant trying them.
{'':8s} The edge limits were {beam_edge} beams spatially and {vel_edge} channels.
''',Configuration,case= ['main','screen'])
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
                    print_log(f'''SOFIA_CATALOGUE: We discarded a very bright source, let's check wether it satisfies our minimum boundaries.
''',Configuration,case= ['debug_add'])
                    index = np.where(np.nanmax(no_edge_fluxes) == no_edge_fluxes)[0][0]
                    edge = check_edge_limits(float(outlist[Variables.index('x_min')][index]),
                                    float(outlist[Variables.index('x_max')][index]),
                                    float(outlist[Variables.index('y_min')][index]),
                                    float(outlist[Variables.index('y_max')][index]),
                                    float(outlist[Variables.index('z_min')][index]),
                                    float(outlist[Variables.index('z_max')][index]),
                                    Configuration,vel_edge=min_vel_edge)
                    cube= float(Configuration['NAXES'][0])*float(Configuration['NAXES'][1])
                    source_size=(float(many_sources[Variables.index('y_max')][index])-float(many_sources[Variables.index('y_min')][index]))* \
                                (float(many_sources[Variables.index('x_max')][index])-float(many_sources[Variables.index('x_min')][index]))

                    if edge and source_size/cube < 0.5:
                        print_log(f'''SOFIA_CATALOGUE: The bright source is very close to limits and not more than half a channel in size.
''',Configuration)

                        if float(outlist[Variables.index('rms')][index])*1e6 < float(outlist[Variables.index('f_sum')][index]):
                            print_log(f'''SOFIA_CATALOGUE: There appears to be no noise in this cube. restoring the source.
''',Configuration)
                            many_sources  = copy.deepcopy(outlist)
                            fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
                    else:
                        if edge and source_size/cube >= 0.5:
                            print_log(f'''SOFIA_CATALOGUE: The bright source is very close to limits but spans more than half a channel so we are restoring it.
In principle your cube is too small.
''',Configuration)
                        else:
                            print_log(f'''SOFIA_CATALOGUE: The bright source is acceptable, restoring its flux.
''',Configuration)
                        many_sources  = copy.deepcopy(outlist)
                        fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
            print_log(f'''SOFIA_CATALOGUE: after checking edges we find these fluxes:
{'':8s}{many_sources[Variables.index('f_sum')]}
''',Configuration,case= ['debug_add'])
        else:
            #If we did not run SoFia or have no_edge_limit we simply assume we want to fit the brightest source in the cube
            many_sources  = copy.deepcopy(outlist)
            fluxes = np.array(outlist[Variables.index('f_sum')],dtype =float)
        outlist = []
        #We want the source with the most total flux.
        index = np.where(np.nanmax(fluxes) == fluxes)[0][0]
        print_log(f'''SOFIA_CATALOGUE: We select the {index} source of this list.
''',Configuration)
        outlist = [x[index] for x in many_sources]
        # check that our mask has the selected source
    else:
        outlist = [x[0] for x in outlist]
    if 'n_pix' in Variables:
        #If we have the amount of pixels in the mask we calculate a SNR in Jy/beam
        SNR=outlist[Variables.index('f_sum')]*Configuration['BEAM_IN_PIXELS'][2]/\
            (Configuration['CHANNEL_WIDTH']*float(outlist[Variables.index('n_pix')])*\
            Configuration['NOISE'])
    else:
        #if we don't know the amount of pixels we just leave the mask unmodified
        SNR = 5
    check_mask(Configuration,outlist[Variables.index('id')],Fits_Files,SNR=SNR)
    print_log(f'''SOFIA_CATALOGUE: we found these values:
{chr(10).join([f'{"":8s}{x} = {y}' for x,y in zip(Variables,outlist)])}
''',Configuration,case= ['debug_add'])
    if strip_npix:
        outlist.pop(Variables.index('n_pix'))
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



 OUTPUTS:
    outlist = list with the requested values for the brightest source in the catalogue provided it is not on the edge of the cube.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def sofia_input_catalogue(Configuration):
    Catalogue = sf.Proper_Dictionary({})

    Catalogue = {'ENTRIES':  ['ENTRIES','DISTANCE','ID','DIRECTORYNAME','CUBENAME'],\
                'DISTANCE': [],'DIRECTORYNAME': [],\
                'ID': [], 'CUBENAME': []}
    basename = Configuration['CATALOGUE'].split('_cat.txt')[0]
    Variables =['id','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra',\
                    'dec','v_app','f_sum','kin_pa','w50','err_f_sum','err_x','err_y','err_z']
    headerlines=[]
    if not Configuration['SOFIA_DIR']:
        Configuration['SOFIA_DIR']=f"{Configuration['MAIN_DIRECTORY']}{basename}_cubelets/"

    #Read the sofia catalogue
    with open(Configuration['CATALOGUE']) as sof_cat:
        for line in sof_cat.readlines():
            tmp =line.split()

            if line.strip() == '' or line.strip() == '#':
                pass
            elif tmp[0] == '#' and len(tmp) > 1:
                headerlines.append(line)
                if tmp[1].strip().lower() in ['name']:
                    input_columns = [x.strip() for x in tmp[1:]]
                    column_locations = []
                    for col in input_columns:
                        column_locations.append(line.find(col)+len(col))
                    # check that we found all parameters
                    velocity = check_parameters(Configuration,Variables,input_columns)                   
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
                    sf.create_directory(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}",f"{Configuration['MAIN_DIRECTORY']}")

                if 'create_fat_cube' in Configuration['FITTING_STAGES']:
                    hdr = create_fat_cube(Configuration,sofia_catalogue=True,name=basename,id = outlist[input_columns.index('id')])

                else:
                    Cube = fits.open(f"{Configuration['SOFIA_DIR']}{basename}_{outlist[input_columns.index('id')]}_cube.fits",uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
                    data = Cube[0].data
                    hdr = Cube[0].header
                    if hdr['NAXIS'] == 4:
                        data = data[0,:,:,:]
                        del hdr['NAXIS4']
                        hdr['NAXIS'] = 3
                    # clean the header
                    hdr = sf.clean_header(Configuration,hdr)
                    channel=int(0)
                    while np.isnan(data[channel,:,:]).all():
                        channel+=1
                    low_chan =np.nanstd(data[channel,:,:])

                    channel=int(-1)
                    while np.isnan(data[channel,:,:]).all():
                        channel-=1
                    high_chan = np.nanstd(data[channel,:,:])

                    hdr['FATNOISE'] = np.mean([low_chan,high_chan])

                    if hdr['CDELT3'] < -1:
                        raise InputError(f"Your velocity axis is declining this won't work. exiting")
                    fits.writeto(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/{basename}_{outlist[input_columns.index('id')]}_FAT.fits",\
                                 data,hdr,overwrite=True)
                    #Translate the big cube to parameters in this cube.
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    coordinate_frame = WCS(hdr)
                    x,y,z =coordinate_frame.wcs_world2pix(float(outlist[input_columns.index('ra')]),\
                                                    float(outlist[input_columns.index('dec')]), \
                                                    float(outlist[input_columns.index(velocity)]), 1.)
                shift = [int(float(outlist[input_columns.index('x')]) - x), \
                        int(float(outlist[input_columns.index('y')]) - y),\
                        int(float(outlist[input_columns.index('z')]) - z)]

                for i,coord in enumerate(['x','y','z']):
                    for chang in ['','_min','_max','_peak']:
                        val = f'{coord}{chang}'
                        if chang == '':
                            outlist[input_columns.index(val)] = f'{float(outlist[input_columns.index(val)])-shift[i]:.6f}'
                        else:
                            outlist[input_columns.index(val)] = f'{int(outlist[input_columns.index(val)])-shift[i]:d}'


                    
                sf.create_directory(f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/Sofia_Output",f"{Configuration['MAIN_DIRECTORY']}")

                Configuration['FITTING_DIR']=f"{Configuration['MAIN_DIRECTORY']}{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}/"
                Configuration['SOFIA_BASENAME'] = f"{basename}_{outlist[input_columns.index('id')]}"
                Configuration['BASE_NAME'] =  f"{basename}_{outlist[input_columns.index('id')]}_FAT"
                sf.copy_homemade_sofia(Configuration,no_cat=True)
                Catalogue['DISTANCE'].append(float(-1))
                Catalogue['ID'].append(f"{outlist[input_columns.index('id')]}")
                Catalogue['DIRECTORYNAME'].append(f"{basename}_FAT_cubelets/{basename}_{outlist[input_columns.index('id')]}")
                Catalogue['CUBENAME'].append(f"{basename}_{outlist[input_columns.index('id')]}")
                with open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_cat.txt",'w') as cat:
                    for hline in headerlines:
                        cat.write(f"{hline}")
                    comline=''
                    for col in input_columns:
                        if input_columns.index(col) == 0:
                            start = 0
                        else:
                            start = column_locations[input_columns.index(col)-1]
                        end = column_locations[input_columns.index(col)]
                        strlenfor = f'>{int(end-start)}s'
                        comline = f'{comline}{outlist[input_columns.index(col)]:{strlenfor}}'
                    cat.write(comline)

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

def sofia_template():
    with pack_open_txt(templates, 'sofia_template.par') as tmp:
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


 OUTPUTS:
    result = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    split, strip, open

 NOTE:
'''
