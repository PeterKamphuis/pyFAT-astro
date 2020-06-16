#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from astropy.io import fits
from support_functions import linenumber,print_log
import numpy as np

class BadHeaderError(Exception):
    pass
class BadCubeError(Exception):
    pass

import warnings


# clean the header
def clean_header(hdr,log):
    if not 'EPOCH' in hdr:
        if 'EQUINOX' in hdr:
            log_statement = f'''CLEAN_HEADER: Your cube has no EPOCH keyword but we found EQUINOX.
{"":8s}We have set EPOCH to {hdr['EQUINOX']}
'''
            print_log(log_statement,log)
            hdr['EPOCH'] = hdr['EQUINOX']
            del hdr['EQUINOX']
        else:
            log_statement = f'''CLEAN_HEADER: Your cube has no EPOCH keyword
{"":8s}CLEAN_HEADER: We assumed J2000
'''
            print_log(log_statement,log)
            hdr['EPOCH'] = 2000.

    if hdr['CUNIT3'].upper() == 'HZ' or hdr['CTYPE3'].upper() == 'FREQ':
        print_log('CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.',log)
        raise BadHeaderError('The Cube has frequency as a velocity axis this is not supported')

    if not 'CUNIT3' in hdr:
        if hdr['CDELT3'] > 100:
            hdr['CUNIT3'] = 'm/s'
        else:
            hdr['CUNIT3'] = 'km/s'
        log_statement = f'''CLEAN_HEADER: Your header did not have a unit for the third axis, that is bad policy.
{"":8s} We have set it to {hdr['CUNIT3']}. Please ensure that is correct.'
'''
        print_log(log_statement,log)
    vel_types = ['VELO-HEL','VELO-LSR','FELO-HEL','FELO-LSR','VELO', 'VELOCITY']
    if hdr['CTYPE3'].upper() not in vel_types:
        if hdr['CTYPE3'].split('-')[0].upper() in ['RA','DEC']:
            log_statement = f'''CLEAN_HEADER: Your zaxis is a spatial axis not a velocity axis.
{"":8s}CLEAN_HEADER: Please arrange your cube logically
'''
            print_log(log_statement,log)
            raise BadHeaderError("The Cube's third axis is not a velocity axis")
        hdr['CTYPE3'] = 'VELO'
        log_statement = f'''CLEAN_HEADER: Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.
'''
        print_log(log_statement,log)

    if hdr['CUNIT3'].lower() == 'km/s':
        log_statement = f'''CLEAN_HEADER: The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'
'''
        print_log(log_statement,log)
        hdr['CUNIT3'] = 'm/s'
        hdr['CDELT3'] = hdr['CDELT3']*1000.
        hdr['CRVAL3'] = hdr['CRVAL3']*1000.
    # Check for the beam
    if not 'BMAJ' in hdr:
        if 'BMMAJ' in hdr:
            hdr['BMAJ']= hdr['BMMAJ']/3600.
        else:
            found = False
            for line in hdr['HISTORY']:
                tmp = [x.strip().upper() for x in line.split()]
                if 'BMAJ=' in tmp:
                    hdr['BMAJ'] = tmp[tmp.index('BMAJ=') + 1]
                    found = True
                if 'BMIN=' in tmp:
                    hdr['BMIN'] = tmp[tmp.index('BMIN=') + 1]
                if 'BPA=' in tmp:
                    hdr['BPA'] = tmp[tmp.index('BPA=') + 1]
                if found:
                    break
            if not found:
                log_statement = f'''CLEAN_HEADER: WE CANNOT FIND THE MAJOR AXIS FWHM IN THE HEADER
'''
                print_log(log_statement,log)
                raise BadHeaderError("The Cube has no major axis FWHM in the header.")
    if not 'CTYPE1' in hdr or not 'CTYPE2' in hdr:
        log_statement = f'''CLEAN_HEADER: Your spatial axes have no ctype. this can lead to errors.
'''
        print_log(log_statement,log)
        raise BadHeaderError("The Cube header has no ctypes.")

    if hdr['CTYPE1'].split('-')[0].upper() in ['DEC']:
        log_statement = f'''CLEAN_HEADER: !!!!!!!!!!Your declination is in the first axis, that's ignorant !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
'''
        print_log(log_statement,log)
        raise BadHeaderError("Your spatial axes are reversed")
    if hdr['CTYPE2'].split('-')[0].upper() in ['RA']:
        log_statement = f'''CLEAN_HEADER: !!!!!!!!!!Your right ascension is on the second axis, that's ignorant !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
'''
        print_log(log_statement,log)
        raise BadHeaderError("Your spatial axes are reversed")

    if not 'BMIN' in hdr:
        if 'BMMIN' in hdr:
            hdr['BMIN']= hdr['BMMIN']/3600.
        else:
            log_statement = f'''CLEAN_HEADER: We cannot find the minor axis FWHM. Assuming a circular beam.
'''
            print_log(log_statement,log)
            hdr['BMIN'] = hdr['BMAJ']

    if len(hdr['HISTORY']) > 10:
        del hdr['HISTORY']
        log_statement = f'''CLEAN_HEADER: Your cube has a significant history attached we are removing it for easier interpretation.
'''
        print_log(log_statement,log)

    if abs(hdr['BMAJ']/hdr['CDELT1']) < 2:
        log_statement = f'''CLEAN_HEADER: !!!!!!!!!!Your cube has less than two pixels per beam major axis.!!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!           This will lead to bad results.              !!!!!!!!!!!!!!!!'
'''
        print_log(log_statement,log)

    if abs(hdr['BMAJ']/hdr['CDELT1']) > hdr['NAXIS1']:
        log_statement = f'''CLEAN_HEADER: !!!!!!!!!!Your cube is smaller than the beam major axis. !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
'''
        print_log(log_statement,log)
        raise BadHeaderError("Your cube is too small for your beam")

clean_header.__doc__ = '''
;+
; NAME:
;       CLEAN_HEADER
;
; PURPOSE:
;       Clean up the cube header and make sure it has all the right
;       variables that we require in the process of fitting
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       CLEAN_HEADER,header,log=log
;
;
; INPUTS:
;       hdr = the header of the cube
;      log = the logging file
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;       clean_header = a header that is ok.
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;
;
; EXAMPLE:

;
; NOTE:
;
;-
'''
# Create a cube suitable for FAT
def create_fat_cube(Configuration, Fits_Files):
    #First get the cubes
    Cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['ORIGINAL_CUBE'],uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    data = Cube[0].data
    hdr = Cube[0].header
    if hdr['NAXIS'] == 4:
        data = data[0,:,:,:]
        del hdr['NAXIS4']
        hdr['NAXIS'] = 3
    # clean the header
    clean_header(hdr,Configuration["OUTPUTLOG"])
    data = prep_cube(hdr,data,Configuration["OUTPUTLOG"])

    # and write our new cube
    log_statement = f'''CREATE_FAT_CUBE: We are writing a FAT modfied cube to be used for the fitting. This cube is called {Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE']}
'''
    print_log(log_statement,Configuration["OUTPUTLOG"])
    fits.writeto(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'],data,hdr)
    # Release the arrays
    Cube.close()
    data = []
    hdr = []


create_fat_cube.__doc__ = '''

;+
; NAME:
;       create_fat_cube
;
; PURPOSE:
;       As we do not want to work from the original cube this function copies that cube into a new cube that we can then happily modify and work with
;       To keep disk space overhead this cube might later on also be cut down to size
;
; CATEGORY:
;       Fits
;
;
; INPUTS:
;       dir =  directory where the fitting is done
;  cubename = name of the original cube
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          beam = the beam recorded in the header. If no beam is
;          present it will return [NaN,Nan] and FAT will abort
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       astropy.io.fits
;
; EXAMPLE:
;
;
'''

#preprocess the cube

def prep_cube(hdr,data,log):

    if hdr['CDELT3'] < -1:
        log_statement = f'''PREPROCESSING: Your velocity axis is declining with increasing channels
{"":8s}PREPROCESSING: We reversed the velocity axis.
'''
        print_log(log_statement, log)
        hdr['CDELT3'] = abs(hdr['CDELT3'])
        hdr['CRPIX3'] = hdr['NAXIS3']-hdr['CRPIX3']+1
        data = data[::-1,:,:]
    #Check for zeros
    if np.where(data == 0.):
        log_statement = f'''PREPROCESSING: Your cube contains values exactly 0. If this is padding these should be blanks.
{"":8s}PREPROCESSING: We have changed them to blanks.
'''
        print_log(log_statement, log)
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
                log_statement = f'''PREPROCESSING: This cube has too many blanked channels.
'''
                print_log(log_statement, log)
                raise BadCubeError("The cube has too many blanked channels")
            log_statement = f'''PREPROCESSING: We are cutting the cube as the first channel is completely blank.
'''
            print_log(log_statement, log)
        while np.isnan(data[-1,:,:]).all():
            data=data[:-1,:,:]
            hdr['NAXIS3'] = hdr['NAXIS3']-1
            if hdr['NAXIS3'] < 5:
                log_statement = f'''PREPROCESSING: This cube has too many blanked channels.
'''
                print_log(log_statement, log)
                raise BadCubeError("The cube has too many blanked channels")
            log_statement = f'''PREPROCESSING: We are cutting the cube as the last channel is completely blank.
'''
            print_log(log_statement, log)
        #Then check the noise statistics
        noise_first_channel = np.nanstd(data[0,:,:])
        noise_last_channel = np.nanstd(data[-1,:,:])
        noise_bottom_right =  np.nanstd(data[:,:6,-6:])
        noise_top_right =  np.nanstd(data[:,-6:,-6:])
        noise_bottom_left =  np.nanstd(data[:,:6,:6])
        noise_top_left =  np.nanstd(data[:,-6:,:6])
        noise_corner = np.mean([noise_top_right,noise_bottom_right,noise_bottom_left,noise_top_left])
        channel_noise = np.mean([noise_first_channel,noise_last_channel])
        difference = abs((noise_first_channel-noise_last_channel)/noise_first_channel)
        if noise_corner/channel_noise >1.5 and difference < 0.2:
            noise_corner = copy.deepcopy(channel_noise)
        difference2 = abs(channel_noise-noise_corner)/channel_noise
        if difference < 0.2 and np.isfinite(difference) and difference2 < 0.25:
            cube_ok = True
        else:
            log_statement = f'''PREPROCESSING: We are cutting the cube as clearly the noise statistics are off.
{"":8s}PREPROCESSING: Noise in the first channel is {noise_first_channel}
{"":8s}PREPROCESSING: Noise in the last channel is {noise_last_channel}
{"":8s}PREPROCESSING: Noise in the corners is {noise_corner}
'''
            print_log(log_statement,log)
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
                    ~np.finite(noise_last_channel):
                    prev_last_comparison = last_comparison
                    times_cut_last_channel += 1
                    data=data[:-1,:,:]
                    hdr['NAXIS3'] = hdr['NAXIS3'] - 1
                else:
                    cube_ok = True

            if times_cut_first_channel >= 8 and times_cut_last_channel >= 8:
                cube_ok = True
                log_statement = f'''  PREPROCESSING: This cube has non-uniform noise statistics.'
'''
                print_log(log_statement,log)
            if hdr['NAXIS3'] < 5:
                log_statement = f'''PREPROCESSING: This cube has noise statistics that cannot be dealt with.
'''
                print_log(log_statement,log)
                raise BadCubeError('The Cube has noise statistics that cannot be dealt with')

    hdr['FATNOISE'] = noise_corner
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        low_noise_indices = np.where(data < -10*noise_corner)
    if len(low_noise_indices) > 0:
        data[low_noise_indices] = float('NaN')
        log_statement=f'''PREPROCESSING: Your cube had values below -10*sigma. If you do not have a central absorption source there is something seriously wrong with the cube.
{"":8s}PREPROCESSING: We blanked these values.
'''
        print_log(log_statement,log)

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

prep_cube.__doc__ = '''
;+
; NAME:
;       prep_cube
;
; PURPOSE:
;       Clean up the cube  and make sure it has all the right
;       characteristics and blanks
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;         cube = the array containing the input cube.
;          log = the logging file.
;
; OPTIONAL INPUTS:
;
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
;       np.where(), np.finite(), np.isnan()
;
; EXAMPLE:
'''
