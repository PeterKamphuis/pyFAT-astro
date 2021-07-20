# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.
from astropy.io import fits
from astropy.wcs import WCS
from pyFAT_astro.Support.support_functions import linenumber,print_log,set_limits,clean_header
from pyFAT_astro.Support.read_functions import obtain_border_pix
from scipy import ndimage
import numpy as np
import copy
import warnings
import os

class BadHeaderError(Exception):
    pass
class BadCubeError(Exception):
    pass
class BadMaskError(Exception):
    pass
#Check that the mask only contains the selected sources
def check_mask(Configuration,id,Fits_Files,debug=False):
    if debug:
        print_log(f'''CHECK_MASK: Checking the mask to contain only the correct source.
''',Configuration['OUTPUTLOG'],debug=True)
    mask = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}",uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')

    if float(id) not in mask[0].data:
        print_log(f'''CHECK_MASK: We cannot find the selected source in the mask. This will lead to errors. Aborting the fit.
''',Configuration['OUTPUTLOG'],screen = True,debug=debug)
        BadMaskError(f" We can not find the sofia source id in the mask.")
    else:
        data = copy.deepcopy(mask[0].data)
        data[data != float(id)] = 0.
        diff = data-mask[0].data
        neg_index = np.where(diff < 0.)[0]
        if neg_index.shape:
            if debug:
                print_log(f'''CHECK_MASK: The initial mask had more than a single source. redoing the mask.
''',Configuration['OUTPUTLOG'],screen = True)
            mask[0].data = data
            fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}",mask[0].data,mask[0].header, overwrite = True)
# to ensure compatible units and calculations with th models we make the maps ourselves
    del mask[0].header['C*3']
    mask[0].data[mask[0].data> 0.5] = 1.
    chan_map = np.array(np.nansum(mask[0].data,axis=0),dtype =float)
    mask[0].header['BITPIX'] = -32
    fits.writeto(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Configuration['BASE_NAME']}_chan.fits",chan_map,header=mask[0].header,overwrite = True)
    mask.close()
    if Configuration['SOFIA_RAN']:
        make_moments(Configuration, Fits_Files,fit_type='Generic_Initialize',vel_unit = 'm/s',debug=debug)
check_mask.__doc__ =f'''
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
    debug = False
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
def create_fat_cube(Configuration, Fits_Files, debug = False):
    #First get the cubes
    Cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['ORIGINAL_CUBE'],uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    data = Cube[0].data
    hdr = Cube[0].header
    if hdr['NAXIS'] == 4:
        data = data[0,:,:,:]
        del hdr['NAXIS4']
        hdr['NAXIS'] = 3
    # clean the header
    hdr = clean_header(Configuration,hdr,debug=debug)
    data = prep_cube(Configuration,hdr,data,debug=debug)
    # and write our new cube
    log_statement = f'''CREATE_FAT_CUBE: We are writing a FAT modfied cube to be used for the fitting. This cube is called {Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE']}
'''
    print_log(log_statement,Configuration["OUTPUTLOG"])
    fits.writeto(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'],data,hdr)
    # Release the arrays

    Cube.close()
    data = []
    hdr = []
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
    debug = False

 OUTPUTS:
    The FAT input cube

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def cut_cubes(Configuration, Fits_Files, galaxy_box, debug = False):

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


    if cut:
        files_to_cut = [Fits_Files['FITTING_CUBE'],'Sofia_Output/'+Fits_Files['MASK'],\
                        'Sofia_Output/'+Fits_Files['MOMENT0'],\
                        'Sofia_Output/'+Fits_Files['MOMENT1'],\
                        'Sofia_Output/'+Fits_Files['MOMENT2'],\
                        'Sofia_Output/'+Fits_Files['CHANNEL_MAP'],\
                        ]


        print_log(f'''CUT_CUBES: Your input cube is significantly larger than the detected source.
{"":8s}CUT_CUBES: we will cut to x-axis = [{new_cube[2,0]},{new_cube[2,1]}] y-axis = [{new_cube[1,0]},{new_cube[1,1]}]
{"":8s}CUT_CUBES: z-axis = [{new_cube[2,0]},{new_cube[2,1]}].
{"":8s}CUT_CUBES: We will cut the following files:
{"":8s}CUT_CUBES: {', '.join(files_to_cut)}
''', Configuration['OUTPUTLOG'])

        for file in files_to_cut:
            cutout_cube(Configuration,file,new_cube,debug=debug)

    #We want to check if the cube has a decent number of pixels per beam.
    if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['OPTIMIZED_CUBE']}"):
        optimized_cube(Configuration, Fits_Files,debug=debug)
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
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    cutout_cube, optimized_cube

 NOTE: The cut cubes will overwrite the _FAT products
       down to size and overwrite the existing files,
       User supplied SoFIA products are hence copied as fat products.
       To keep them safe.
'''

def cutout_cube(Configuration,filename,sub_cube, debug = False):
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
            Configuration['NAXES_LIMITS'] = [np.sort([xlow,xhigh]),np.sort([ylow,yhigh]),np.sort([zlow,zhigh])/1000.]
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
    debug = False

 OUTPUTS:
    the cut cube is written to disk.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# Extract a PV-Diagrams
def extract_pv(Configuration,cube_in,angle,center=[-1,-1,-1],finalsize=[-1,-1],convert=-1, debug = False):
    if debug:
        print_log(f'''EXTRACT_PV: We are the extraction of a PV-Diagram
{'':8s} PA = {angle}
{'':8s} center = {center}
{'':8s} finalsize = {finalsize}
{'':8s} convert = {convert}
''', Configuration['OUTPUTLOG'], debug =True)


    cube = copy.deepcopy(cube_in)
    hdr = copy.deepcopy(cube[0].header)
    TwoD_hdr= copy.deepcopy(cube[0].header)
    data = copy.deepcopy(cube[0].data)
    #Because astro py is even dumber than Python
    try:
        if hdr['CUNIT3'].lower() == 'km/s':
            hdr['CUNIT3'] = 'm/s'
            hdr['CDELT3'] = hdr['CDELT3']*1000.
            hdr['CRVAL3'] = hdr['CRVAL3']*1000.
        elif hdr['CUNIT3'].lower() == 'm/s':
            hdr['CUNIT3'] = 'm/s'
    except KeyError:
        hdr['CUNIT3'] = 'm/s'
    if center[0] == -1:
        center = [hdr['CRVAL1'],hdr['CRVAL2'],hdr['CRVAL3']]
        xcenter,ycenter,zcenter = hdr['CRPIX1'],hdr['CRPIX2'],hdr['CRPIX3']
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            coordinate_frame = WCS(hdr)
        xcenter,ycenter,zcenter = coordinate_frame.wcs_world2pix(center[0], center[1], center[2], 1.)

    nz, ny, nx = data.shape
    if finalsize[0] != -1:
        if finalsize[1] >nz:
            finalsize[1] = nz
        if finalsize[0] > nx:
            finalsize[0] = nx
    # if the center is not set assume the crval values
    if debug:
        print_log(f'''EXTRACT_PV: The shape of the output
{'':8s} nz = {nz}
{'':8s} ny = {ny}
{'':8s} nx = {nx}
''', Configuration['OUTPUTLOG'])
    x1,x2,y1,y2 = obtain_border_pix(Configuration,angle,[xcenter,ycenter],debug=debug)
    linex,liney,linez = np.linspace(x1,x2,nx), np.linspace(y1,y2,nx), np.linspace(0,nz-1,nz)
    #This only works when ny == nx hence nx is used in liney
    new_coordinates = np.array([(z,y,x)
                        for z in linez
                        for y,x in zip(liney,linex)
                        ],dtype=float).transpose().reshape((-1,nz,nx))
    #spatial_resolution = abs((abs(x2-x1)/nx)*np.sin(np.radians(angle)))+abs(abs(y2-y1)/ny*np.cos(np.radians(angle)))
    PV = ndimage.map_coordinates(data, new_coordinates,order=1)
    if hdr['CDELT1'] < 0:
        PV = PV[:,::-1]

    if finalsize[0] == -1:
        # then lets update the header
        # As python is stupid making a simple copy will mean that these changes are still applied to hudulist
        TwoD_hdr['NAXIS2'] = nz
        TwoD_hdr['NAXIS1'] = nx

        TwoD_hdr['CRPIX2'] = hdr['CRPIX3']
        if convert !=-1:
            TwoD_hdr['CRVAL2'] = hdr['CRVAL3']/convert
        else:
            TwoD_hdr['CRVAL2'] = hdr['CRVAL3']
        TwoD_hdr['CRPIX1'] = xcenter+1
    else:
        zstart = set_limits(int(zcenter-finalsize[1]/2.),0,int(nz))
        zend = set_limits(int(zcenter+finalsize[1]/2.),0,int(nz))
        xstart = set_limits(int(xcenter-finalsize[0]/2.),0,int(nx))
        xend = set_limits(int(xcenter+finalsize[0]/2.),0,int(nx))
        PV =  PV[zstart:zend, xstart:xend]
        TwoD_hdr['NAXIS2'] = int(finalsize[1])
        TwoD_hdr['NAXIS1'] = int(finalsize[0])
        TwoD_hdr['CRPIX2'] = hdr['CRPIX3']-int(nz/2.-finalsize[1]/2.)
        if convert !=-1:
            TwoD_hdr['CRVAL2'] = hdr['CRVAL3']/convert
        else:
            TwoD_hdr['CRVAL2'] = hdr['CRVAL3']

        TwoD_hdr['CRPIX1'] = int(finalsize[0]/2.)+1
    if convert !=-1:
        TwoD_hdr['CDELT2'] = hdr['CDELT3']/convert
    else:
        TwoD_hdr['CDELT2'] = hdr['CDELT3']
    TwoD_hdr['CTYPE2'] = hdr['CTYPE3']
    try:
        if hdr['CUNIT3'].lower() == 'm/s' and convert == -1:
            TwoD_hdr['CDELT2'] = hdr['CDELT3']/1000.
            TwoD_hdr['CRVAL2'] = hdr['CRVAL3']/1000.
            TwoD_hdr['CUNIT2'] = 'km/s'
            del (TwoD_hdr['CUNIT3'])
        elif  convert != -1:
            del (TwoD_hdr['CUNIT3'])
            del (TwoD_hdr['CUNIT2'])
        else:
            TwoD_hdr['CUNIT2'] = hdr['CUNIT3']
            del (TwoD_hdr['CUNIT3'])
    except:
        print_log(f'''EXTRACT_PV: We could not find units in the header for the 3rd axis.
''', Configuration['OUTPUTLOG'])
    del (TwoD_hdr['CRPIX3'])
    del (TwoD_hdr['CRVAL3'])
    del (TwoD_hdr['CDELT3'])
    del (TwoD_hdr['CTYPE3'])

    del (TwoD_hdr['NAXIS3'])
    TwoD_hdr['CRVAL1'] = 0.
    #Because we used nx in the linspace for liney we also use it here
    TwoD_hdr['CDELT1'] = np.sqrt(((x2-x1)*abs(hdr['CDELT1'])/nx)**2+((y2-y1)*abs(hdr['CDELT2'])/nx)**2)*3600.

    TwoD_hdr['CTYPE1'] = 'OFFSET'
    TwoD_hdr['CUNIT1'] = 'ARCSEC'
    TwoD_hdr['HISTORY'] = f'EXTRACT_PV: PV diagram extracted with angle {angle} and center {center}'
    # Then we change the cube and rteturn the PV construct
    cube[0].header = TwoD_hdr
    cube[0].data = PV

    return cube
extract_pv.__doc__ = '''
 NAME:
    extract_pv

 PURPOSE:
    extract a PV diagram from a cube object. Angle is the PA and center the central location. The profile is extracted over the full length of the cube and afterwards cut back to the finalsize.

 CATEGORY:
     fits_functions

 INPUTS:
    Configuration = Standard FAT Configuration
    cube_in = is a fits cube object
    angle = Pa of the slice

 OPTIONAL INPUTS:
    center = [-1,-1,-1]
    the central location of the slice in WCS map_coordinates [RA,DEC,VSYS], default is the CRVAL values in the header

    finalsize = [-1,-1,-1]
    final size of the PV-diagram in pixels, default is no cutting

    convert=-1
    conversion factor for velocity axis, default no conversion

    debug = False

 KEYWORD PARAMETERS:

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def prep_cube(Configuration,hdr,data, debug = False):
    if debug:
        print_log( f'''PREPROCESSING: starting the preprocessing of the cube
''',Configuration['OUTPUTLOG'],debug = True)

    if hdr['CDELT3'] < -1:
        print_log( f'''PREPROCESSING: Your velocity axis is declining with increasing channels
{"":8s}PREPROCESSING: We reversed the velocity axis.
''', Configuration['OUTPUTLOG'])
        hdr['CDELT3'] = abs(hdr['CDELT3'])
        hdr['CRPIX3'] = hdr['NAXIS3']-hdr['CRPIX3']+1
        data = data[::-1,:,:]
    #Check for zeros
    if np.where(data == 0.):
        print_log(f'''PREPROCESSING: Your cube contains values exactly 0. If this is padding these should be blanks.
{"":8s}PREPROCESSING: We have changed them to blanks.
''', Configuration['OUTPUTLOG'])
        data[np.where(data == 0.)] = float('NaN')

    # check for blank channels and noise statistics
    cube_ok = False
    prev_first_comparison = 0.
    prev_last_comparison = 0.
    times_cut_first_channel = 0
    times_cut_last_channel = 0
    while not cube_ok:
        while np.isnan(data[0,:,:]).all():
            data=data[1:,:,:]
            hdr['NAXIS3'] = hdr['NAXIS3']-1
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has too many blanked channels.
''', Configuration['OUTPUTLOG'],screen=True)
                raise BadCubeError("The cube has too many blanked channels")
            log_statement = f'''PREPROCESSING: We are cutting the cube as the first channel is completely blank.
'''
            print_log(f'''PREPROCESSING: We are cutting the cube as the first channel is completely blank.
''', Configuration['OUTPUTLOG'])
        while np.isnan(data[-1,:,:]).all():
            data=data[:-1,:,:]
            hdr['NAXIS3'] = hdr['NAXIS3']-1
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has too many blanked channels.
''', Configuration['OUTPUTLOG'],screen = True)
                raise BadCubeError("The cube has too many blanked channels")
            print_log(f'''PREPROCESSING: We are cutting the cube as the last channel is completely blank.
''', Configuration['OUTPUTLOG'])
        #Then check the noise statistics
        noise_first_channel = np.nanstd(data[0,:,:])
        noise_last_channel = np.nanstd(data[-1,:,:])
        noise_bottom_right =  np.nanstd(data[:,:6,-6:])
        noise_top_right =  np.nanstd(data[:,-6:,-6:])
        noise_bottom_left =  np.nanstd(data[:,:6,:6])
        noise_top_left =  np.nanstd(data[:,-6:,:6])
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
''',Configuration['OUTPUTLOG'])
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
''',Configuration['OUTPUTLOG'])
            if hdr['NAXIS3'] < 5:
                print_log(f'''PREPROCESSING: This cube has noise statistics that cannot be dealt with.
''',Configuration['OUTPUTLOG'],screen=True)
                raise BadCubeError('The Cube has noise statistics that cannot be dealt with')
    if ~np.isfinite(noise_corner):
        print_log(f'''PREPROCESSING: This cube has noise statistics that cannot be dealt with.
''',Configuration['OUTPUTLOG'],screen=True)
        raise BadCubeError('The Cube has noise statistics that cannot be dealt with')
    hdr['FATNOISE'] = noise_corner
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low_noise_indices = np.where(data < -10*noise_corner)
    if len(low_noise_indices) > 0:
        data[low_noise_indices] = float('NaN')
        print_log(f'''PREPROCESSING: Your cube had values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.
{"":8s}PREPROCESSING: We blanked these values.
''',Configuration['OUTPUTLOG'])

    # Check whether any channels got fully blanked
    blanked_channels = []
    for z in range(hdr['NAXIS3']):
        if np.isnan(data[z,:,:]).all():
            blanked_channels.append(str(z))
    if len(blanked_channels) == 0.:
        blanked_channels = ['-1']
    hdr['BL_CHAN'] = ','.join(blanked_channels)
    #Previously SoFiA was not dealing with blanks properly and the statement went here

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
    debug = False

 OUTPUTS:
    the updated header

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Create an optimized cube if required
def optimized_cube(Configuration,Fits_Files, debug = False):
    pix_per_beam = Configuration['BEAM_IN_PIXELS'][1]
    if pix_per_beam > Configuration['OPT_PIXEL_BEAM']:
        cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'])
        data = cube[0].data
        hdr = cube[0].header
        if f"{abs(hdr['CDELT1']):.16f}" != f"{abs(hdr['CDELT2']):.16f}":
            print_log(f'''OPTIMIZED_CUBE: Your input cube does not have square pixels.
{"":8s}OPTIMIZED_CUBE: FAT cannot optimize your cube.
''', Configuration['OUTPUTLOG'])
        cube.close()
        required_cdelt = hdr['BMIN']/int(Configuration['OPT_PIXEL_BEAM'])
        ratio = required_cdelt/abs(hdr['CDELT2'])
        opt_data,opt_hdr = regrid_cube(data, hdr, ratio)

        fits.writeto(Configuration['FITTING_DIR']+Fits_Files['OPTIMIZED_CUBE'], opt_data,opt_hdr)

        Configuration['OPTIMIZED'] = True
        print_log(f'''OPTIMIZED_CUBE: Your input cube has {pix_per_beam} pixels along the minor FWHM.
{"":8s}OPTIMIZED_CUBE: You requested { Configuration['OPT_PIXEL_BEAM']} therefore we regridded the cube into a new cube.
{"":8s}OPTIMIZED_CUBE: We are using the pixel size of {required_cdelt}.
''', Configuration['OUTPUTLOG'])

    else:
        print_log(f'''OPTIMIZED_CUBE: Your input cube has {pix_per_beam} pixels along the minor FWHM.
{"":8s}OPTIMIZED_CUBE: You requested {Configuration['OPT_PIXEL_BEAM']} but we cannot improve the resolution.
{"":8s}OPTIMIZED_CUBE: We are using the pixel size of the original cube.
''', Configuration['OUTPUTLOG'])
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
    debug = False

 OUTPUTS:
    An optimized cube is created and written to disk.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def make_moments(Configuration,Fits_Files,fit_type = 'Undefined',
                 moments = [0,1,2],overwrite = False, level=None,
                  vel_unit= None, debug = False):
    if fit_type == 'Generic_Initialize':
        filename = f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}"
        basename = f"{Configuration['BASE_NAME']}"
        directory = f"{Configuration['FITTING_DIR']}/Sofia_Output/"
        mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}"
    elif fit_type == 'Generic_Final':
        filename = f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits"
        basename = 'Finalmodel'
        directory = f"{Configuration['FITTING_DIR']}/Finalmodel/"
        mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}"
    else:
        filename = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.fits"
        basename = fit_type
        directory = f"{Configuration['FITTING_DIR']}{fit_type}"
        if not level:
            mask_cube = f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MASK']}"

    cube = fits.open(filename)
    if vel_unit:
        cube[0].header['CUNIT3'] = vel_unit
    if mask_cube:
        mask = fits.open(mask_cube)
        with np.errstate(invalid='ignore', divide='ignore'):
            cube[0].data[mask[0].data < 0.5] = float('NaN')
        mask.close()
    else:
        if not level:
            level = 3.*np.mean([np.nanstd(cube[0].data[0:2,:,:]),np.nanstd(cube[0].data[-3:-1,:,:])])
        with np.errstate(invalid='ignore', divide='ignore'):
            cube[0].data[cube[0].data < level] = float('NaN')
    try:
        if cube[0].header['CUNIT3'].lower().strip() == 'm/s':
            print_log(f"We convert your m/s to km/s", Configuration['OUTPUTLOG'])
            cube[0].header['CUNIT3'] = 'km/s'
            cube[0].header['CDELT3'] = cube[0].header['CDELT3']/1000.
            cube[0].header['CRVAL3'] = cube[0].header['CRVAL3']/1000.
        elif cube[0].header['CUNIT3'].lower().strip() == 'km/s':
            pass
        else:
            print_log(f"Your Velocity unit {cube[0].header['CUNIT3']} is weird. Your units could be off", Configuration['OUTPUTLOG'])
    except KeyError:
        print_log(f"Your CUNIT3 is missing, that is bad practice. We'll add a blank one but we're not guessing the value", Configuration['OUTPUTLOG'])
        cube[0].header['CUNIT3'] = 'Unknown'
    #Make a 2D header to use
    hdr2D = copy.deepcopy(cube[0].header)
    hdr2D.remove('NAXIS3')
    hdr2D['NAXIS'] = 2
    # removing the third axis means we cannot correct for varying platescale, Sofia does so this is and issue so let's not do this
    hdr2D.remove('CDELT3')
    hdr2D.remove('CTYPE3')
    hdr2D.remove('CUNIT3')
    hdr2D.remove('CRPIX3')
    hdr2D.remove('CRVAL3')

    # we need a moment 0 for the moment 2 as well
    if 0 in moments:
        hdr2D['BUNIT'] = f"{cube[0].header['BUNIT']}*{cube[0].header['CUNIT3']}"
        moment0 = np.nansum(cube[0].data, axis=0) * cube[0].header['CDELT3']
        moment0[np.invert(np.isfinite(moment0))] = float('NaN')
        hdr2D['DATAMAX'] = np.nanmax(moment0)
        hdr2D['DATAMIN'] = np.nanmin(moment0)
        fits.writeto(f"{directory}/{basename}_mom0.fits",moment0,hdr2D,overwrite = overwrite)
    if 1 in moments or 2 in moments:
        zaxis = cube[0].header['CRVAL3'] + (np.arange(cube[0].header['NAXIS3'])+1 \
              - cube[0].header['CRPIX3']) * cube[0].header['CDELT3']
        c=np.transpose(np.resize(zaxis,[cube[0].header['NAXIS1'],cube[0].header['NAXIS2'],len(zaxis)]),(2,1,0))
        hdr2D['BUNIT'] = f"{cube[0].header['CUNIT3']}"
        # Remember Python is stupid so z,y,x
        with np.errstate(invalid='ignore', divide='ignore'):
            moment1 = np.nansum(cube[0].data*c, axis=0)/ np.nansum(cube[0].data, axis=0)
        moment1[np.invert(np.isfinite(moment1))] = float('NaN')
        hdr2D['DATAMAX'] = np.nanmax(moment1)
        hdr2D['DATAMIN'] = np.nanmin(moment1)
        if 1 in moments:
            fits.writeto(f"{directory}/{basename}_mom1.fits",moment1,hdr2D,overwrite = overwrite)
        if 2 in moments:
            d = c - np.resize(moment1,[len(zaxis),cube[0].header['NAXIS2'],cube[0].header['NAXIS1']])
            with np.errstate(invalid='ignore', divide='ignore'):
                moment2 = np.sqrt(np.nansum(cube[0].data*d**2, axis=0)/ np.nansum(cube[0].data, axis=0))
            moment2[np.invert(np.isfinite(moment1))] = float('NaN')
            hdr2D['DATAMAX'] = np.nanmax(moment2)
            hdr2D['DATAMIN'] = np.nanmin(moment2)
            fits.writeto(f"{directory}/{basename}_mom2.fits",moment2,hdr2D,overwrite = overwrite)
    cube.close()

make_moments.__doc__ =f'''
 NAME:
    make_moments

 PURPOSE:
    Make the moment maps

 CATEGORY:
    fits_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames


 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of ftting

    moments = [0,1,2]
    moment maps to create

    overwrite = False
    overwrite existing maps

    level=None
    cutoff level to use, if set the mask will not be used

    vel_unit= none
    velocity unit of the input cube

 OUTPUTS:

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
    except:
        print("No CDELT found")
    try:
        regrid_hdr['CD1_1'] =   regrid_hdr['CD1_1']*real_ratio[2]
        regrid_hdr['CD2_2'] =   regrid_hdr['CD2_2']*real_ratio[1]
        wcs_found = True
    except:
        if not wcs_found:
            print("No CD corr matrix found")
    try:
        regrid_hdr['CD1_2'] =   regrid_hdr['CD1_2']*real_ratio[2]
        regrid_hdr['CD2_1'] =   regrid_hdr['CD2_1']*real_ratio[1]
        wcs_found = True
    except:
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
