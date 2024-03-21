# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.
from astropy.io import fits
from astropy.wcs import WCS
from datetime import datetime
from shutil import copyfile
from pyFAT_astro.Support.fat_errors import FunctionCallError,BadCubeError\
    ,BadMaskError,InputError
from pyFAT_astro.Support.log_functions import print_log
import pyFAT_astro.Support.support_functions as sf


from make_moments.functions import moments

#from pyFAT_astro.Support.read_functions import obtain_border_pix
from scipy import ndimage
import numpy as np
import copy
import warnings
import os

#Check that the mask only contains the selected sources
def check_mask(Configuration,id,Fits_Files,SNR= 5.):
    print_log(f'''CHECK_MASK: Checking the mask to contain only the correct source.
and to have decent signal noise
''',Configuration, case =['debug_start','verbose'])
    if Configuration['DEBUG']:
        copyfile(f"{Configuration['FITTING_DIR']}/{Fits_Files['MASK']}",f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_Original_Mask.fits")
    mask = fits.open(f"{Configuration['FITTING_DIR']}/{Fits_Files['MASK']}",uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    mask = check_mask_sources(Configuration,id,mask)

    if SNR < 3.:
        mask = smooth_mask(Configuration,mask)
    fits.writeto(f"{Configuration['FITTING_DIR']}{Fits_Files['MASK']}",mask[0].data,\
            mask[0].header, overwrite = True)
    # to ensure compatible units and calculations with th models we make the maps ourselves
    del mask[0].header['C*3']
    chan_map = np.array(np.nansum(mask[0].data,axis=0),dtype =float)
    fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_chan.fits",\
        chan_map,header=mask[0].header,overwrite = True)
    mask.close()
    if Configuration['SOFIA_RAN']:
        print_log('CHECK_MASK: Creating Sofia Moments',Configuration,case=['verbose'])
        messages = moments(filename =  f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}",\
                        mask =  f"{Configuration['FITTING_DIR']}/{Fits_Files['MASK']}",\
                        overwrite = True, map_velocity_unit= 'km/s',\
                        debug = Configuration['DEBUG'], log=True,\
                        output_directory =   f"{Configuration['FITTING_DIR']}/Sofia_Output",\
                        output_name = f"{Configuration['BASE_NAME']}")
        print_log(messages,Configuration,case=['verbose'])
        
    
    #moments()
        #make_moments(Configuration, Fits_Files,fit_type='Generic_Initialize',vel_unit = 'm/s')

check_mask.__doc__ =f'''
 NAME:
    check_mask

 PURPOSE:
    Check that the mask only contains the source we want and if that source
    has low SNR that it is smoothed with a beam

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    id = the SoFiA ID of the source we are interested in
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    fit_type = 'Undefined'

 OUTPUTS:
    moment maps and the channel map

 OPTIONAL OUTPUTS:
    a cleaned mask

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def check_mask_sources(Configuration,id,mask):
    if float(id) not in mask[0].data:
        print_log(f'''CHECK_MASK: We cannot find the selected source in the mask. This will lead to errors. Aborting the fit.
    ''',Configuration,case=['main','screen'])
        BadMaskError(f" We can not find the sofia source id in the mask.")
    else:
        data = copy.deepcopy(mask[0].data)
        data[data != float(id)] = 0.
        diff = data-mask[0].data
        neg_index = np.where(diff < 0.)[0]
        if len(neg_index) != 0:
            print_log(f'''CHECK_MASK: The initial mask had more than a single source. redoing the mask.
    ''',Configuration, case = ['verbose'])
            mask[0].data = data


    mask[0].data[mask[0].data> 0.5] = 1.

    mask[0].header['BITPIX'] = -32
    return mask
check_mask_sources.__doc__ =f'''
 NAME:
    check_mask

 PURPOSE:
    Check that the mask only contains the source we want and not parts of other nearby sources.
    Subsequently make the channel_map and the moment maps.

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    id = the SoFiA ID of the source we are interested in
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    fit_type = 'Undefined'

 OUTPUTS:
    moment maps and the channel map

 OPTIONAL OUTPUTS:
    a cleaned mask

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
# Create a cube suitable for FAT
def create_fat_cube(Configuration, Fits_Files = None,sofia_catalogue=False,\
        id='No default',name='No default'):
    Configuration['PREPARATION_TIME'][0] = datetime.now()
    #First get the cubes
    if sofia_catalogue:
        Cube = fits.open(f"{Configuration['SOFIA_DIR']}{name}_{id}_cube.fits",\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    else:
        Cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['ORIGINAL_CUBE'],\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    data = Cube[0].data
    hdr = Cube[0].header
    if hdr['NAXIS'] == 4:
        data = data[0,:,:,:]
        del hdr['NAXIS4']
        hdr['NAXIS'] = 3
    # clean the header
    hdr = sf.clean_header(Configuration,hdr)
    data = prep_cube(Configuration,hdr,data)
    # and write our new cube

  
    if sofia_catalogue:
        print_log(f'''CREATE_FAT_CUBE: We are writing a FAT modified cube to be used for the fitting. This cube is called {name}_{id}/{name}_{id}_FAT.fits
''',Configuration)
        try:
            fits.writeto(f"{Configuration['MAIN_DIRECTORY']}{name}_FAT_cubelets/{name}_{id}/{name}_{id}_FAT.fits",data,hdr,overwrite=Configuration['SOFIA_OVERWRITE'])
        except OSError:
            raise InputError('pyFAT has already setup this directory, if you want to overwrite the bcubelets please set advanced.sofia_overwrite to True in the configuration')

    else:
        print_log(f'''CREATE_FAT_CUBE: We are writing a FAT modified cube to be used for the fitting. This cube is called {Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE']}
''',Configuration)
        fits.writeto(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'],data,hdr)
    # Release the arrays

    Cube.close()

    data = []
    if sofia_catalogue:
        return hdr
    hdr = []
    Configuration['PREPARATION_TIME'][1] = datetime.now()

create_fat_cube.__doc__ =f'''
 NAME:
    create_fat_cube

 PURPOSE:
    As we do not want to work from the original cube
    this function copies that cube into a new cube
    that we can then happily modify and work with.
    To keep disk space overhead this cube might later on also
    be cut down to size.

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:

 OUTPUTS:
    The FAT input cube

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def cut_cubes(Configuration, Fits_Files, galaxy_box):
    print_log(f'''CUT_CUBES: Starting to cut the cube_size.
''', Configuration, case= ['debug_start'])
    cube_edge= [6.,5.*round(Configuration['BEAM_IN_PIXELS'][0]),5.*round(Configuration['BEAM_IN_PIXELS'][0])]
    cube_size= []
    new_cube = []
    for i in [2,1,0]:
        cube_size.append(Configuration['NAXES'][i])
        new_cube.append([0,Configuration['NAXES'][i]])
    new_cube = np.array(new_cube,dtype=int)
    cut  = False
    for i,limit in enumerate(galaxy_box):
        if limit[0] > cube_edge[i]:
            cut = True
            new_cube[i,0] = limit[0]-int(cube_edge[i])
        if limit[1] < cube_size[i] - cube_edge[i]:
            cut = True
            new_cube[i,1] = limit[1]+int(cube_edge[i])
    if 'restart_fitting' in [x.lower() for x in Configuration['FITTING_STAGES']]:
        if 'create_fat_cube' in [x.lower() for x in Configuration['FITTING_STAGES']]:
            files_to_cut = [Fits_Files['FITTING_CUBE']]
        else:
            cut = False
    else:
        files_to_cut = [Fits_Files['FITTING_CUBE'],Fits_Files['MASK'],\
                        Fits_Files['MOMENT0'],\
                        Fits_Files['MOMENT1'],\
                        Fits_Files['MOMENT2'],\
                        Fits_Files['CHANNEL_MAP'],\
                        ]
    if cut:
        print_log(f'''CUT_CUBES: Your input cube is significantly larger than the detected source.
{"":8s}CUT_CUBES: we will cut to x-axis = [{new_cube[2,0]},{new_cube[2,1]}] y-axis = [{new_cube[1,0]},{new_cube[1,1]}]
{"":8s}CUT_CUBES: z-axis = [{new_cube[2,0]},{new_cube[2,1]}].
{"":8s}CUT_CUBES: We will cut the following files:
{"":8s}CUT_CUBES: {', '.join(files_to_cut)}
''', Configuration)

        for file in files_to_cut:
            print_log(f'''CUT_CUBES: We are cutting the file {file}
''', Configuration)
            if os.path.exists(f"{Configuration['FITTING_DIR']}{file}"):
                cutout_cube(Configuration,file,new_cube)
            else:
                if file == Fits_Files['CHANNEL_MAP']:
                    pass
                else:
                    raise FunctionCallError(f'We are trying to cut {file} but it does not exist')


    #We want to check if the cube has a decent number of pixels per beam.
    if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['OPTIMIZED_CUBE']}"):
        optimized_cube(Configuration, Fits_Files)
    return new_cube
cut_cubes.__doc__ =f'''
 NAME:
    cut_cubes

 PURPOSE:
    Cut all FAT related output back in size to a system that fits snugly around the sofia detection

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    galaxy_box = array that contains the new size as
                 [[z_min,z_max],[y_min,y_max], [x_min,x_max]]
                 adhering to fits' idiotic way of reading fits files

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    cutout_cube, optimized_cube

 NOTE: The cut cubes will overwrite the _FAT products
       down to size and overwrite the existing files,
       User supplied SoFIA products are hence copied as fat products.
       To keep them safe.
'''

def cutout_cube(Configuration,filename,sub_cube):
    Cube = fits.open(Configuration['FITTING_DIR']+filename,uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    hdr = Cube[0].header

    if hdr['NAXIS'] == 3:
        data = Cube[0].data[sub_cube[0,0]:sub_cube[0,1],sub_cube[1,0]:sub_cube[1,1],sub_cube[2,0]:sub_cube[2,1]]
        hdr['NAXIS1'] = sub_cube[2,1]-sub_cube[2,0]
        hdr['NAXIS2'] = sub_cube[1,1]-sub_cube[1,0]
        hdr['NAXIS3'] = sub_cube[0,1]-sub_cube[0,0]
        hdr['CRPIX1'] = hdr['CRPIX1'] -sub_cube[2,0]
        hdr['CRPIX2'] = hdr['CRPIX2'] -sub_cube[1,0]
        hdr['CRPIX3'] = hdr['CRPIX3'] -sub_cube[0,0]
        #Only update when cutting the cube
        Configuration['NAXES'] =[ hdr['NAXIS1'],hdr['NAXIS2'],hdr['NAXIS3']]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coordinate_frame = WCS(hdr)
            xlow,ylow,zlow = coordinate_frame.wcs_pix2world(1,1,1., 1.)
            xhigh,yhigh,zhigh = coordinate_frame.wcs_pix2world(*Configuration['NAXES'], 1.)
            xlim = np.sort([xlow,xhigh])
            ylim = np.sort([ylow,yhigh])
            zlim =np.sort([zlow,zhigh])/1000.
            sf.set_boundaries(Configuration,'VSYS',*zlim,input=True)
            sf.set_boundaries(Configuration,'XPOS',*xlim,input=True)
            sf.set_boundaries(Configuration,'YPOS',*ylim,input=True)

    elif hdr['NAXIS'] == 2:
        data = Cube[0].data[sub_cube[1,0]:sub_cube[1,1],sub_cube[2,0]:sub_cube[2,1]]
        hdr['NAXIS1'] = sub_cube[2,1]-sub_cube[2,0]
        hdr['NAXIS2'] = sub_cube[1,1]-sub_cube[1,0]
        hdr['CRPIX1'] = hdr['CRPIX1'] -sub_cube[2,0]
        hdr['CRPIX2'] = hdr['CRPIX2'] -sub_cube[1,0]

    Cube.close()
    fits.writeto(Configuration['FITTING_DIR']+filename,data,hdr,overwrite = True)
cutout_cube.__doc__ =f'''
 NAME:
    cutout_cube

 PURPOSE:
    Cut filename back to the size of subcube, update the header and write back to disk.

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    filename = name of the cube to be cut
    sub_cube = array that contains the new size as
                [[z_min,z_max],[y_min,y_max], [x_min,x_max]]
                adhering to fits' idiotic way of reading fits files.

 OPTIONAL INPUTS:

 OUTPUTS:
    the cut cube is written to disk.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def prep_cube(Configuration,hdr,data):
    print_log( f'''PREPROCESSING: starting the preprocessing of the cube
''',Configuration, case = ['debug_start'])

    if hdr['CDELT3'] < -1:
        print_log( f'''PREPROCESSING: Your velocity axis is declining with increasing channels
{"":8s}PREPROCESSING: We reversed the velocity axis.
''', Configuration)
        hdr['CDELT3'] = abs(hdr['CDELT3'])
        hdr['CRPIX3'] = hdr['NAXIS3']-hdr['CRPIX3']+1
        data = data[::-1,:,:]
    #Check for zeros
    if np.where(data == 0.):
        print_log(f'''PREPROCESSING: Your cube contains values exactly 0. If this is padding these should be blanks.
{"":8s}PREPROCESSING: We have changed them to blanks.
''', Configuration)
        data[np.where(data == 0.)] = float('NaN')

    # check for blank channels and noise statistics
    cube_ok = False
    prev_first_comparison = 0.
    prev_last_comparison = 0.
    times_cut_first_channel = 0
    times_cut_last_channel = 0
    corner_box_size = int(np.floor(2*hdr['BMAJ']/abs(hdr['CDELT1'])))
    if corner_box_size > np.mean([hdr['NAXIS1'],hdr['NAXIS2']])/10.:
        corner_box_size=int(np.floor(np.mean([hdr['NAXIS1'],hdr['NAXIS2']])/10.))
    while not cube_ok:
        while np.isnan(data[0,:,:]).all():
            data=data[1:,:,:]
            hdr['NAXIS3'] = hdr['NAXIS3']-1
            hdr['CRPIX3'] = hdr['CRPIX3'] - 1
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has too many blanked channels.
''', Configuration,case = ['main', 'screen'])
                raise BadCubeError("The cube has too many blanked channels")
            log_statement = f'''PREPROCESSING: We are cutting the cube as the first channel is completely blank.
'''
            print_log(f'''PREPROCESSING: We are cutting the cube as the first channel is completely blank.
''', Configuration)
        while np.isnan(data[-1,:,:]).all():
            data=data[:-1,:,:]
            hdr['NAXIS3'] = hdr['NAXIS3']-1
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has too many blanked channels.
''', Configuration,case = ['main', 'screen'])
                raise BadCubeError("The cube has too many blanked channels")
            print_log(f'''PREPROCESSING: We are cutting the cube as the last channel is completely blank.
''', Configuration)
        #Then check the noise statistics
        noise_first_channel = np.nanstd(data[0,:,:])
        noise_last_channel = np.nanstd(data[-1,:,:])
        noise_bottom_right =  np.nanstd(data[:,:corner_box_size,-corner_box_size:])
        noise_top_right =  np.nanstd(data[:,-corner_box_size:,-corner_box_size:])
        noise_bottom_left =  np.nanstd(data[:,:corner_box_size,:corner_box_size])
        noise_top_left =  np.nanstd(data[:,-corner_box_size:,:corner_box_size])
        noise_corner = np.nanmean([noise_top_right,noise_bottom_right,noise_bottom_left,noise_top_left])
        channel_noise = np.nanmean([noise_first_channel,noise_last_channel])
        if ~np.isfinite(noise_corner):
            noise_corner = copy.deepcopy(channel_noise)
        difference = abs((noise_first_channel-noise_last_channel)/noise_first_channel)
        if noise_corner/channel_noise >1.5 and difference < 0.2:
            noise_corner = copy.deepcopy(channel_noise)
        difference2 = abs(channel_noise-noise_corner)/channel_noise
        if difference < 0.2 and np.isfinite(difference) and difference2 < 0.25:
            cube_ok = True
        else:
            print_log(f'''PREPROCESSING: We are cutting the cube as clearly the noise statistics are off.
{"":8s}PREPROCESSING: Noise in the first channel is {noise_first_channel}
{"":8s}PREPROCESSING: Noise in the last channel is {noise_last_channel}
{"":8s}PREPROCESSING: Noise in the corners is {noise_corner}
''',Configuration,case= ['verbose'])
            first_comparison = abs(noise_first_channel-noise_corner)/noise_corner
            last_comparison = abs(noise_last_channel-noise_corner)/noise_corner
            if prev_first_comparison == 0:
                prev_first_comparison = first_comparison
            if prev_last_comparison == 0.:
                prev_last_comparison = last_comparison
            if (first_comparison > last_comparison and \
               abs((first_comparison-prev_first_comparison)/first_comparison) >= abs((last_comparison-prev_last_comparison)/last_comparison) and \
               times_cut_first_channel < 0.1*hdr['NAXIS3']) or \
               ~np.isfinite(noise_first_channel):
                prev_first_comparison = first_comparison
                times_cut_first_channel += 1
                data=data[1:,:,:]
                hdr['CRPIX3'] = hdr['CRPIX3'] - 1
                hdr['NAXIS3'] = hdr['NAXIS3'] - 1
            else:
                if times_cut_last_channel < 0.1*hdr['CRPIX3'] or \
                    ~np.isfinite(noise_last_channel):
                    prev_last_comparison = last_comparison
                    times_cut_last_channel += 1
                    data=data[:-1,:,:]
                    hdr['NAXIS3'] = hdr['NAXIS3'] - 1
                else:
                    cube_ok = True

            if times_cut_first_channel >= 8 and times_cut_last_channel >= 8:
                cube_ok = True
                print_log( f'''PREPROCESSING: This cube has non-uniform noise statistics.'
''',Configuration)
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has noise statistics that cannot be dealt with.
''',Configuration,case = ['main','screen'])
                raise BadCubeError('The Cube has noise statistics that cannot be dealt with')
    if ~np.isfinite(noise_corner):
        print_log(f'''PREPROCESSING: This cube has noise statistics that cannot be dealt with.
''',Configuration,case = ['main','screen'])
        raise BadCubeError('The Cube has noise statistics that cannot be dealt with')
    #hdr['FATNOISE'] = noise_corner
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low_noise_indices = np.where(data < -10*noise_corner)
    if len(low_noise_indices) > 0.0001*hdr['NAXIS1']*hdr['NAXIS2']*hdr['NAXIS3']:
        data[low_noise_indices] = float('NaN')
        print_log(f'''PREPROCESSING: Your cube had a significant amount of values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.
{"":8s}PREPROCESSING: We blanked these values.
''',Configuration)

    #Let's make sure that a central absorption source is treated the same regardless od the size of the cube
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        buffer=5.*hdr['BMAJ']/np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])
        central = data[:,int(hdr['NAXIS2']/2.-buffer):int(hdr['NAXIS2']/2.+buffer),\
                        int(hdr['NAXIS1']/2.-buffer):int(hdr['NAXIS1']/2.+buffer)]

        central_noise_indices = np.where(central < -10*noise_corner)
    if len(central_noise_indices) > 0.0001*4*buffer*hdr['NAXIS3']:
        central[central_noise_indices] = float('NaN')
        data[:,int(hdr['NAXIS2']/2.-buffer):int(hdr['NAXIS2']/2.+buffer),\
                        int(hdr['NAXIS1']/2.-buffer):int(hdr['NAXIS1']/2.+buffer)] =\
                        central
        print_log(f'''PREPROCESSING: Your cube had a significant amount of values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.
{"":8s}PREPROCESSING: We blanked these values.
''',Configuration)


    # Check whether any channels got fully blanked
    blanked_channels = []
    for z in range(hdr['NAXIS3']):
        if np.isnan(data[z,:,:]).all():
            blanked_channels.append(str(z))
    if len(blanked_channels) == 0.:
        blanked_channels = ['-1']
    hdr['BL_CHAN'] = ','.join(blanked_channels)
    #Previously SoFiA was not dealing with blanks properly and the statement went here
    #finally we want to get the noise level from the negative in the final cube

    new_noise = np.nanstd(np.hstack([data[data<0],-1.*data[data<0]]))
    #If we have noiseles cubes this will be far to low
    if abs(new_noise/np.mean([channel_noise,noise_corner])) > 4. \
        or abs(new_noise/np.mean([channel_noise,noise_corner])) < 0.25:
            print_log(f'''PREPROCESSING: There is something odd in the noise statistics of your cube. We are not using the nehgative values
''',Configuration)
            #new_noise = np.mean([channel_noise,noise_corner])
    hdr['FATNOISE'] = new_noise

    return data
prep_cube.__doc__ =f'''
 NAME:
    prep_cube

 PURPOSE:
    Clean up the cube  and make sure it has all the right
    characteristics and blanks

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    hdr = header to be cleaned
    data = data array to be cleaned

 OPTIONAL INPUTS:

 OUTPUTS:
    the updated header

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Create an optimized cube if required
def optimized_cube(Configuration,Fits_Files):
    pix_per_beam = Configuration['BEAM_IN_PIXELS'][1]
    if pix_per_beam > Configuration['OPT_PIXEL_BEAM']:
        cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'])
        data = cube[0].data
        hdr = cube[0].header
        if f"{abs(hdr['CDELT1']):.16f}" != f"{abs(hdr['CDELT2']):.16f}":
            print_log(f'''OPTIMIZED_CUBE: Your input cube does not have square pixels.
{"":8s}OPTIMIZED_CUBE: FAT cannot optimize your cube.
''', Configuration)
        cube.close()
        required_cdelt = hdr['BMIN']/int(Configuration['OPT_PIXEL_BEAM'])
        ratio = required_cdelt/abs(hdr['CDELT2'])
        opt_data,opt_hdr = regrid_cube(data, hdr, ratio)

        fits.writeto(Configuration['FITTING_DIR']+Fits_Files['OPTIMIZED_CUBE'], opt_data,opt_hdr)

        Configuration['OPTIMIZED'] = True
        print_log(f'''OPTIMIZED_CUBE: Your input cube has {pix_per_beam} pixels along the minor FWHM.
{"":8s}OPTIMIZED_CUBE: You requested { Configuration['OPT_PIXEL_BEAM']} therefore we regridded the cube into a new cube.
{"":8s}OPTIMIZED_CUBE: We are using the pixel size of {required_cdelt}.
''', Configuration)

    else:
        print_log(f'''OPTIMIZED_CUBE: Your input cube has {pix_per_beam} pixels along the minor FWHM.
{"":8s}OPTIMIZED_CUBE: You requested {Configuration['OPT_PIXEL_BEAM']} but we cannot improve the resolution.
{"":8s}OPTIMIZED_CUBE: We are using the pixel size of the original cube.
''', Configuration)
optimized_cube.__doc__ =f'''
 NAME:
    optimized_cube

 PURPOSE:
    Check the amount of pixels that are in the minor axis FWHM,
    if more than OPT_PIXELBEAM regrid the cube to less pixels
    and write FAT_opt cube and set Configuration['OPTIMIZED'] = True

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames

 OPTIONAL INPUTS:

 OUTPUTS:
    An optimized cube is created and written to disk.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def regrid_cube(data,hdr,ratio):
    regrid_hdr = copy.deepcopy(hdr)
    # First get the shape of the data
    shape = np.array(data.shape, dtype=float)
    #The new shape has to be integers
    new_shape = [int(x/ratio) for x in shape]
    # Which means our real ratio is
    real_ratio = shape/new_shape
    #* np.ceil(shape / ratio).astype(int)
    new_shape[0] = shape[0]
    # Create the zero-padded array and assign it with the old density
    regrid_data = regridder(data,new_shape)
    regrid_hdr['CRPIX1'] =   (regrid_hdr['CRPIX1']-1)/real_ratio[2]+1
    regrid_hdr['CRPIX2'] =   (regrid_hdr['CRPIX2']-1)/real_ratio[1]+1
    wcs_found = False
    try:
        regrid_hdr['CDELT1'] =   regrid_hdr['CDELT1']*real_ratio[2]
        regrid_hdr['CDELT2'] =   regrid_hdr['CDELT2']*real_ratio[1]
        wcs_found = True
    except KeyError:
        print("No CDELT found")
    try:
        regrid_hdr['CD1_1'] =   regrid_hdr['CD1_1']*real_ratio[2]
        regrid_hdr['CD2_2'] =   regrid_hdr['CD2_2']*real_ratio[1]
        wcs_found = True
    except KeyError:
        if not wcs_found:
            print("No CD corr matrix found")
    try:
        regrid_hdr['CD1_2'] =   regrid_hdr['CD1_2']*real_ratio[2]
        regrid_hdr['CD2_1'] =   regrid_hdr['CD2_1']*real_ratio[1]
        wcs_found = True
    except KeyError:
        if not wcs_found:
            print("No CD cross-corr matrix found")

    regrid_hdr['NAXIS1'] =   new_shape[2]
    regrid_hdr['NAXIS2'] =   new_shape[1]
    return regrid_data,regrid_hdr
regrid_cube.__doc__ =f'''
 NAME:
    regrid_cube

 PURPOSE:
    Regrid a cube to an arbitrary smaller cube in the spatial plane with larger pixels

 CATEGORY:
    fits_functions

 INPUTS:
    data = the array containing the input cube.
    hdr = the header data
    ratio = the requested ratio of old/new

 OPTIONAL INPUTS:

 OUTPUTS:
    regriddata = regridded data array
    regrid_hdr = header corresponding to the regridded data

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    copy.deepcopy, regridder, np.array

 NOTE: This might not work well with extended projections.
'''

def regridder(oldarray, newshape):
    oldshape = np.array(oldarray.shape)
    newshape = np.array(newshape, dtype=float)
    ratios = oldshape/newshape
        # calculate new dims
    nslices = [ slice(0,j) for j in list(newshape) ]
    #make a list with new coord
    new_coordinates = np.mgrid[nslices]
    #scale the new coordinates
    for i in range(len(ratios)):
        new_coordinates[i] *= ratios[i]
    #create our regridded array
    newarray = ndimage.map_coordinates(oldarray, new_coordinates,order=1)

    return newarray
regridder.__doc__ =f'''
 NAME:
regridder
 PURPOSE:
Regrid an array into a new shape through the ndimage module
 CATEGORY:
    fits_functions

 INPUTS:
    oldarray = the larger array
    newshape = the new shape that is requested

 OPTIONAL INPUTS:

 OUTPUTS:
    newarray = regridded array

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    scipy.ndimage.map_coordinates, np.array, np.mgrid

 NOTE:
'''

def smooth_mask(Configuration,mask):
    print_log(f'''SMOOTH_MASK: Starting smooth mask
''', Configuration)

    # First we smooth the mask
    mask_data = np.array(mask[0].data,dtype=float)
    mask_data = ndimage.gaussian_filter(mask_data,sigma=(0.,\
        Configuration['BEAM_IN_PIXELS'][1] / np.sqrt(8 * np.log(2)),\
        Configuration['BEAM_IN_PIXELS'][0] / np.sqrt(8 * np.log(2))),order=0)
    mask_data[mask_data > 0.05] = 1.
    mask_data[mask_data < 0.05] = 0.
    mask[0].data=mask_data
    return mask
    #reaplly it to get the moments


    # If we have a lot of pixels and stuff we need smooth more
    #So we take the inverse of the initial factor with a standard of 0.75 minimum of 0.5 and a maximum of 1.25

smooth_mask.__doc__ =f'''
 NAME:
    smooth_mask
 PURPOSE:
    smooth the mask with a beam
 CATEGORY:
    fits_functions

 INPUTS:

 OPTIONAL INPUTS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:

 NOTE:
'''
