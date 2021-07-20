# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from inspect import getframeinfo,stack

from scipy.optimize import curve_fit, OptimizeWarning
from scipy import ndimage
from astropy.wcs import WCS
from astropy.io import fits
from dataclasses import  asdict
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res
import matplotlib.pyplot as plt
import os
import sys
import inspect
import signal
import time
import traceback
import numpy as np
import copy
import warnings
import re
import subprocess
class SupportRunError(Exception):
    pass
class SmallSourceError(Exception):
    pass
class FileNotFoundError(Exception):
    pass
class TirificKillError(Exception):
    pass
class InputError(Exception):
    pass
class ProgramError(Exception):
    pass

# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def __setitem__(self, key, value):
        if key not in self:
            # If it is a new item we only allow it if it is not Configuration or Original_Cube or if we are in setup_configuration
            function,variable,empty = traceback.format_stack()[-2].split('\n')
            function = function.split()[-1].strip()
            variable = variable.split('[')[0].strip()
            if variable == 'Original_Configuration' or variable == 'Configuration':
                if function != 'setup_configuration':
                    raise ProgramError("FAT does not allow additional values to the Configuration outside the setup_configuration in support_functions.")
        OrderedDict.__setitem__(self,key, value)
    #    print("what habbens now")
    def insert(self, existing_key, new_key, key_value):
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

def calc_rings(Configuration,size_in_beams = 0., ring_size  = 0.,debug=False):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    if size_in_beams == 0.:
        size_in_beams = Configuration['SIZE_IN_BEAMS']
    if debug:
        print_log(f'''CALC_RINGS: Calculating the number of rings in the model.
{'':8s} size in beams = {size_in_beams}
{'':8s} ring_size = {ring_size}
{'':8s} the maximum amount of rings = {Configuration['MAX_SIZE_IN_BEAMS']}
''',Configuration['OUTPUTLOG'],debug = True)
    est_rings = round((size_in_beams)/(ring_size)+2.)
    #if est_rings > 20 and Configuration['MAX_SIZE_IN_BEAMS'] > 25:
    if est_rings > 20:
        Configuration['OUTER_RINGS_DOUBLED'] = True
        # 0.5 should always be rounding up else we get problemsin the reverse
        no_rings = 2+10+np.floor((est_rings-12)/2.+0.5)
    else:
        Configuration['OUTER_RINGS_DOUBLED'] = False
        no_rings = est_rings
    if debug:
        print_log(f'''CALC_RINGS: The model will have {no_rings} Rings.
{'':8s} We estimated {est_rings} and double outer rings is {Configuration['OUTER_RINGS_DOUBLED']}
''',Configuration['OUTPUTLOG'])
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
    debug = False

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


def check_sofia(Configuration,Fits_Files,debug=False):
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
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
# clean the header
def clean_header(Configuration,hdr_in,two_dim=False,mask_file=False, debug = False):
    hdr = copy.deepcopy(hdr_in)
    if debug:
        print_log(f'''CLEAN_HEADER: Starting to clean the header.
''',Configuration['OUTPUTLOG'],debug=True)
    keywords = ['CDELT','CUNIT','CRPIX','CRVAL','CTYPE']
    for key in keywords:
        try:
            del hdr[f'{key}4']
        except:
            pass
    if not two_dim:
        hdr['NAXIS'] = 3
        if not 'CUNIT3' in hdr:
            if hdr['CDELT3'] > 500:
                hdr['CUNIT3'] = 'm/s'
            else:
                hdr['CUNIT3'] = 'km/s'
            print_log(f'''CLEAN_HEADER: Your header did not have a unit for the third axis, that is bad policy.
    {"":8s} We have set it to {hdr['CUNIT3']}. Please ensure that is correct.'
    ''',Configuration['OUTPUTLOG'])
        if hdr['CUNIT3'].upper() == 'HZ' or hdr['CTYPE3'].upper() == 'FREQ':
            print_log('CLEAN_HEADER: FREQUENCY IS NOT A SUPPORTED VELOCITY AXIS.', Configuration['OUTPUTLOG'],screen=True)
            raise BadHeaderError('The Cube has frequency as a velocity axis this is not supported')

        vel_types = ['VELO-HEL','VELO-LSR','VELO', 'VELOCITY']
        if hdr['CTYPE3'].upper() not in vel_types:
            if hdr['CTYPE3'].split('-')[0].upper() in ['RA','DEC']:
                print_log(f'''CLEAN_HEADER: Your zaxis is a spatial axis not a velocity axis.
    {"":8s}CLEAN_HEADER: Please arrange your cube logically
    ''',Configuration['OUTPUTLOG'],screen=True)
                raise BadHeaderError("The Cube's third axis is not a velocity axis")
            hdr['CTYPE3'] = 'VELO'
            print_log(f'''CLEAN_HEADER: Your velocity projection is not standard. The keyword is changed to VELO (relativistic definition). This might be dangerous.
    ''',Configuration['OUTPUTLOG'])

        if hdr['CUNIT3'].lower() == 'km/s':
            print_log( f'''CLEAN_HEADER: The channels in your input cube are in km/s. This sometimes leads to problems with wcs lib, hence we change it to m/s.'
    ''',Configuration['OUTPUTLOG'])
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
''',Configuration['OUTPUTLOG'])
            hdr['EPOCH'] = hdr['EQUINOX']
            del hdr['EQUINOX']
        else:
            print_log(f'''CLEAN_HEADER: Your cube has no EPOCH keyword
{"":8s}CLEAN_HEADER: We assumed J2000
''',Configuration['OUTPUTLOG'])
            hdr['EPOCH'] = 2000.



    if f'PC01_01' in hdr:
        del hdr['PC0*_0*']


    # Check for the beam
    if not 'BMAJ' in hdr and not mask_file:
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
                print_log(f'''CLEAN_HEADER: WE CANNOT FIND THE MAJOR AXIS FWHM IN THE HEADER
''',Configuration['OUTPUTLOG'],screen=True)
                raise BadHeaderError("The Cube has no major axis FWHM in the header.")
    if not 'CTYPE1' in hdr or not 'CTYPE2' in hdr:
        print_log(f'''CLEAN_HEADER: Your spatial axes have no ctype. this can lead to errors.
''',Configuration['OUTPUTLOG'],screen=True)
        raise BadHeaderError("The Cube header has no ctypes.")

    if hdr['CTYPE1'].split('-')[0].upper() in ['DEC']:
        print_log(f'''CLEAN_HEADER: !!!!!!!!!!Your declination is in the first axis. !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
''',Configuration['OUTPUTLOG'],screen = True)
        raise BadHeaderError("Your spatial axes are reversed")
    if hdr['CTYPE2'].split('-')[0].upper() in ['RA']:
        print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your right ascension is on the second axis. !!!!!!!!!!!!!!!!!
{"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
''',Configuration['OUTPUTLOG'],screen = True)
        raise BadHeaderError("Your spatial axes are reversed")
    if hdr['CRVAL1'] < 0.:
        print_log(f'''CLEAN_HEADER: your RA crval is negative, this can lead to errors. Adding 360. deg
''',Configuration['OUTPUTLOG'])
        hdr['CRVAL1'] = hdr['CRVAL1']+360.
    if not 'BMIN' in hdr and not mask_file:
        if 'BMMIN' in hdr:
            hdr['BMIN']= hdr['BMMIN']/3600.
        else:
            print_log(f'''CLEAN_HEADER: We cannot find the minor axis FWHM. Assuming a circular beam.
''',Configuration['OUTPUTLOG'])
            hdr['BMIN'] = hdr['BMAJ']
    if not 'BPA' in hdr and not mask_file:
            print_log(f'''CLEAN_HEADER: We cannot find the Beam PA assuming it to be 0
''',Configuration['OUTPUTLOG'])
            hdr['BPA'] = 0.

    try:
        if len(hdr['HISTORY']) > 10:
            del hdr['HISTORY']
            print_log( f'''CLEAN_HEADER: Your cube has a significant history attached we are removing it for easier interpretation.
''',Configuration['OUTPUTLOG'])
    except KeyError:
        pass
    try:
        if abs(hdr['BMAJ']/hdr['CDELT1']) < 2:
            print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your cube has less than two pixels per beam major axis.!!!!!!!!!!!!!!!!!
    {"":8s}CLEAN_HEADER: !!!!!!!!!!           This will lead to bad results.              !!!!!!!!!!!!!!!!'
    ''',Configuration['OUTPUTLOG'])

        if abs(hdr['BMAJ']/hdr['CDELT1']) > hdr['NAXIS1']:
            print_log( f'''CLEAN_HEADER: !!!!!!!!!!Your cube is smaller than the beam major axis. !!!!!!!!!!!!!!!!!
    {"":8s}CLEAN_HEADER: !!!!!!!!!!         This will not work.          !!!!!!!!!!!!!!!!
    ''',Configuration['OUTPUTLOG'])
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
    debug = False

 OUTPUTS:
    the updated header

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def columndensity(Configuration,levels,systemic = 100.,beam=[-1.,-1.],channel_width=-1.,column= False,arcsquare=False,solar_mass_input =False,solar_mass_output=False, debug = False):

    if beam[0] == -1:
        if not arcsquare:
            beam = Configuration['BEAM'][:2]
    if channel_width == -1:
        if arcsquare:
            channel_width = 1.
        else:
            channel_width = Configuration['CHANNEL_WIDTH']

    if debug:
        print_log(f'''COLUMNDENSITY: Starting conversion from the following input.
{'':8s}Levels = {levels}
{'':8s}Beam = {beam}
{'':8s}channel_width = {channel_width}
''',Configuration['OUTPUTLOG'],debug =True)
    beam=np.array(beam,dtype=float)
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg
    if debug:
        print_log(f'''COLUMNDENSITY: We have the following input for calculating the columns.
{'':8s}COLUMNDENSITY: level = {levels}, channel_width = {channel_width}, beam = {beam}, systemic = {systemic})
''',Configuration['OUTPUTLOG'])
    if systemic > 10000:
        systemic = systemic/1000.
    f = f0 * (1 - (systemic / c)) #Systemic frequency
    if arcsquare:
        HIconv = 605.7383 * 1.823E18 * (2. *np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass_input:
                levels=levels*solarmass/(mHI*pc**2)
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
    if ~column and solar_mass_input:
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
    debug = False

    systemic = 100.
    the systemic velocity of the source

    beam  = [-1.,-1.]
    the FWHM of the beam in arcsec, if unset taken from Configuration

    channelwidth = -1. width of a channel in km/s
    channelwidth of the observation if unset taken from Configuration

    column = false
    if True input is columndensities else in mJy/beam

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

 NOTE:
'''

        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(Configuration,RAin,DECin,invert=False, colon=False,debug = False):
    if debug:
        print_log(f'''CONVERTRADEC: Starting conversion from the following input.
    {'':8s}RA = {RAin}
    {'':8s}DEC = {DECin}
''',Configuration['OUTPUTLOG'],debug =True)
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
    debug = False

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
def convertskyangle(Configuration, angle, distance=1., unit='arcsec', distance_unit='Mpc', physical=False,debug = False):
    if debug:
        print_log(f'''CONVERTSKYANGLE: Starting conversion from the following input.
    {'':8s}Angle = {angle}
    {'':8s}Distance = {distance}
''',Configuration['OUTPUTLOG'],debug =True)
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
        print_log('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.\n',Configuration['OUTPUTLOG'],screen=True)
        print_log('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n',Configuration['OUTPUTLOG'],screen=True)
        raise SupportRunError('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.')
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print_log('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n',Configuration['OUTPUTLOG'],screen=True)
            print_log('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n',Configuration['OUTPUTLOG'],screen=True)
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')


        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle / (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle * (10 ** 3)
        else:
            print_log('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n',Configuration['OUTPUTLOG'],screen=True)
            print_log('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n',Configuration['OUTPUTLOG'],screen=True)
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
    debug = False

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


def copy_homemade_sofia(Configuration,no_cat=False,debug=False):
    create_directory('Sofia_Output',Configuration['FITTING_DIR'])
    files =['_mask.fits','_mom0.fits','_mom1.fits','_chan.fits','_mom2.fits']
    if not no_cat:
        files.append('_cat.txt')
    for file in files:
        os.system(f'''cp {Configuration['SOFIA_DIR']}{Configuration['SOFIA_BASENAME']}{file} {Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']}{file}''')
        if not os.path.exists(f"{Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']+file}"):
            if file in ['_mask.fits','_cat.txt']:
                print_log(f'''COPY_HOMEMADE_SOFIA: Something went wrong copying the file  {Configuration['SOFIA_DIR']}{Configuration['SOFIA_BASENAME']+file}
{'':8s} to the file {Configuration['FITTING_DIR']}Sofia_Output/{Configuration['BASE_NAME']+file}.
{'':8s}We are aborting this fit as we cannot make it without Sofia input.
''',Configuration['OUTPUTLOG'],screen=True,debug=debug )
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
            two_dim = True
            hdr = clean_header(Configuration,hdr_in,two_dim=two_dim,mask_file=mask_file,debug=debug)
            update = False
            for key in hdr:
                try:
                    if hdr[key] != hdr_in[key]:
                        update =True
                except:
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
    debug = False
    check = False
    if set to true FAT merely checks for the existings of the sofia files where they expect them
 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''


def create_directory(directory,base_directory,debug=False):
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
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:
    The requested directory is created but only if it does not yet exist

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def deproject(Configuration,map,angle, center = 0., invert = False,debug=False):
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
    debug = False

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

def finish_current_run(Configuration,current_run,debug=False):
    print_log(f"FINISH_CURRENT_RUN: Is Tirific Running? {Configuration['TIRIFIC_RUNNING']}. \n",Configuration['OUTPUTLOG'],debug=debug)
    if Configuration['TIRIFIC_RUNNING']:
        try:
            current_run.stdout.close()
            current_run.stderr.close()
        except:
            pass
        try:
            os.kill(Configuration['TIRIFIC_PID'], signal.SIGKILL)
            print_log(f"FINISH_CURRENT_RUN: We killed PID = {Configuration['TIRIFIC_PID']}. \n",Configuration['OUTPUTLOG'])
        except:
            try:
                current_run.kill()
                print_log(f"FINISH_CURRENT_RUN: We killed the current run although we failed on the PID = {Configuration['TIRIFIC_PID']}. \n",\
                          Configuration['OUTPUTLOG'])
            except AttributeError:
                print_log(f"FINISH_CURRENT_RUN: We failed to kill the current run with PID {Configuration['TIRIFIC_PID']} even though we have tirific running",Configuration['OUTPUTLOG'],screen=True)
                raise TirificKillError('FINISH_CURRENT_RUN: Despite having an initialized tirific we could not kill it. This should not happen.')
        Configuration['TIRIFIC_RUNNING'] = False
        Configuration['TIRIFIC_PID'] = 'Not Initialized'
    else:
        print_log(f"FINISH_CURRENT_RUN: No run is initialized. \n",Configuration['OUTPUTLOG'])

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
    debug = False

 OUTPUTS:
    kills tirific if initialized or raises an error when it fails to do so while tirific is running

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def fit_sine(Configuration,x,y,debug = False):
    if debug:
        print_log(f'''FIT_SINE: Starting to fit a sin.
{'':8s}x = {x}
{'':8s}y = {y}
{'':8s} x size = {x.size} y size = {y.size}
''', Configuration['OUTPUTLOG'],debug =True)
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

    #print(est_width)
    est_amp = np.mean(y)
    peak_location = np.where(y == np.nanmax(y))[0]
    if peak_location.size > 1:
        peak_location = int(peak_location[0])

    est_center = float(x[peak_location])

    est_width=est_width*(2.*np.pi/180.)
    #print(est_peak,est_center,est_width,est_amp)
    with warnings.catch_warnings():
        warnings.simplefilter("error", OptimizeWarning)
        try:
            sin_parameters, sin_covariance = curve_fit(sine, x[~np.isnan(y)], y[~np.isnan(y)],p0=[est_peak,est_center,est_width,est_amp])
        except RuntimeError:
            sin_parameters = [float('NaN') for z in range(4)]
        except OptimizeWarning:
            if debug:
                print_log(f'''FIT_SINE: the covariance could not be estimated. Using initial estimates
{'':8s}peak est = {est_peak}
{'':8s}center est = {est_center}
{'':8s}width est = {est_width}
{'':8s}amp est = {est_amp}
''', Configuration['OUTPUTLOG'],debug =True)
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

   debug = False

OUTPUTS:
   ratios
   the fitted sin profile or when the  width is too small or to wide the original profiles

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:
'''


def fit_gaussian(Configuration,x,y, covariance = False,errors = None, debug = False):
    if debug:
        print_log(f'''FIT_GAUSSIAN: Starting to fit a Gaussian.
{'':8s}x = {x}
{'':8s}y = {y}
''', Configuration['OUTPUTLOG'],debug =True)
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
    gauss_parameters, gauss_covariance = curve_fit(gaussian_function, x, y,p0=[est_peak,est_center,est_sigma],sigma= errors,absolute_sigma= absolute_sigma)
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

    debug = False

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

def get_from_template(Configuration,Tirific_Template,Variables, debug = False):
    out = []
    if debug:
        print_log(f'''GET_FROM_TEMPLATE: Trying to get the following profiles {Variables}
''',Configuration['OUTPUTLOG'] ,debug= True)
    for key in Variables:
        out.append([float(x) for x  in Tirific_Template[key].split()])
    #Because lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    if debug:
        print_log(f'''GET_FROM_TEMPLATE: We extracted the following profiles from the Template.
{'':8s}GET_FROM_TEMPLATE: {out}
''',Configuration['OUTPUTLOG'])
    #Beware that lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    # However if you make a np. array from it make sure that you specify float  or have lists of the same length else you get an array of lists which behave just as dumb
    return out

get_from_template.__doc__ =f'''
 NAME:
    get_from_template
 PURPOSE:
    Return a specified list of prameters from the template. This puts the template values in a list !!!!!!!!

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template =  Standard tirific template
    Variables = parameters to be extracted

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    list with the values of the requested variables

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    This is very similar to load_template in read_functions but this returns an
    unchecked list whereas load_template return a np array with length NUR,
    the latter can thus can include zeros
'''

def get_inclination_pa(Configuration, Image, center, cutoff = 0., debug = False):
    map = copy.deepcopy(Image[0].data)
    for i in [0,1]:
        if debug:
            print_log(f'''GET_INCLINATION_PA: Doing iteration {i} in the estimates
''',Configuration['OUTPUTLOG'])
        # now we need to get profiles under many angles let's say 100
        #extract the profiles under a set of angles
        angles = np.linspace(0, 360, 180)
        ratios, maj_extent = obtain_ratios(Configuration,map, center, angles,noise = cutoff,debug=debug)
        if np.any(np.isnan(ratios)):
            return [float('NaN'),float('NaN')],  [float('NaN'),float('NaN')],float('NaN')
        sin_ratios,sin_parameters = fit_sine(Configuration,angles,ratios,debug=debug)
        if np.any(np.isnan(sin_parameters)):
            return [float('NaN'),float('NaN')],  [float('NaN'),float('NaN')],float('NaN')

        '''
        import matplotlib
        matplotlib.use('TkAgg')
        fig = plt.figure()
        plt.plot(angles,ratios)
        plt.plot(angles,sin_ratios)
        plt.show()
        '''
        ratios=sin_ratios
        if debug:
            if i == 0:
                print_log(f'''GET_INCLINATION_PA: We initially find radius of {maj_extent*3600./(Configuration['BEAM'][0])} beams.
''',Configuration['OUTPUTLOG'], debug = True)
                print_log(f'''GET_INCLINATION_PA: We initially find the ratios:
{'':8s} ratios = {ratios}
''',Configuration['OUTPUTLOG'])
            else:
                print_log(f'''GET_INCLINATION_PA: From the cleaned map we find radius of {maj_extent*3600./Configuration['BEAM'][0]} beams.
''',Configuration['OUTPUTLOG'])
                print_log(f'''GET_INCLINATION_PA: We  find these ratios from the cleaned map:
{'':8s} ratios = {ratios}
''',Configuration['OUTPUTLOG'])
        max_index = np.where(ratios == np.nanmax(ratios))[0]
        if max_index.size > 1:
            max_index =int(max_index[0])
        min_index = np.where(ratios == np.nanmin(ratios))[0]
        if min_index.size > 1:
            min_index =int(min_index[0])
        max_index = set_limits(max_index,2,177)
        min_index = set_limits(min_index,2,177)
        if debug:
            if i == 0:
                print_log(f'''GET_INCLINATION_PA: We initially find these indeces min {min_index } {angles[min_index]} max {max_index} {angles[max_index]}.
''',Configuration['OUTPUTLOG'])
            else:
                print_log(f'''GET_INCLINATION_PA: From the cleaned map we find these indeces min {min_index }  {angles[min_index]} max {max_index} {angles[max_index]}.
''',Configuration['OUTPUTLOG'])
        #get a 10% bracket

        tenp_max_index = np.where(ratios > np.nanmax(ratios)*0.9)[0]
        tenp_min_index = np.where(ratios < np.nanmin(ratios)*1.1)[0]

        if tenp_max_index.size <= 1 and 2 <= max_index <=177 :
            tenp_max_index= [max_index-2,max_index+2]

        if tenp_min_index.size <= 1:
            tenp_min_index= [min_index-2,min_index+2]
        if angles[min_index]-90 > 0.:
            if angles[max_index] > 165:
                pa = float(np.nanmean(np.array([angles[min_index]+90,angles[max_index]],dtype=float)))
            else:
                pa = float(np.nanmean(np.array([angles[min_index]-90,angles[max_index]],dtype=float)))
        else:
            if angles[max_index] < 15:
                pa = float(np.nanmean(np.array([angles[min_index]-90,angles[max_index]],dtype=float)))
            else:
                pa = float(np.nanmean(np.array([angles[min_index]+90,angles[max_index]],dtype=float)))
        if 180. < pa:
            pa = pa -180
        if debug:
            if i == 0:
                print_log(f'''GET_INCLINATION_PA: The initial final pa = {pa}.
''',Configuration['OUTPUTLOG'])
            else:
                print_log(f'''GET_INCLINATION_PA: The pa from the cleaned map  = {pa}.
''',Configuration['OUTPUTLOG'])
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
        if debug:
            if i == 0:
                print_log(f'''GET_INCLINATION_PA: The initial inclination = {inclination}.
''',Configuration['OUTPUTLOG'])
            else:
                print_log(f'''GET_INCLINATION_PA: From the cleaned map we find inclination = {inclination}.
''',Configuration['OUTPUTLOG'])
        # this leads to trouble for small sources due to uncertain PA and inclination estimates
        if inclination < 70. and maj_extent*3600./Configuration['BEAM'][0] > 4:
            Image_clean = remove_inhomogeneities(Configuration,Image,inclination=inclination, pa = pa,iteration=i, center = center,WCS_center = False, debug=debug)
            map = Image_clean[0].data
            Image_clean.close()
        else:
            break
    inclination_error = inclination_error/np.sin(np.radians(inclination))
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
    debug = False

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





def check_tiltogram(Configuration, tiltogram,inner_min=3,debug=False):
    for i in [0,1]:
        theta_inner = np.array([np.mean(tiltogram[i,0:x+1,0:x+1]) for x in range(tiltogram.shape[1])],dtype=float)
        theta_mutual = np.array([np.mean(tiltogram[i,x:,0:x+1]) for x in range(tiltogram.shape[1])],dtype=float)
        #theta_outer = np.array([np.mean(tiltogram[i,x+1:,x+1:]) for x in range(tiltogram.shape[1])],dtype=float)
        #And then we want to apply the rules
        # (i) the difference between theta_inner and theta_mutual is larger than the differences observed at other radii
        # (ii) theta_inner < 5 deg
        # (iii) thetamut > 15 deg
        diff = np.array(abs(theta_inner-theta_mutual),dtype = float)
        if debug:
            print_log(f'''CHECK_TILTOGRAM: Checking the tiltogram in side {i}.
{'':8s} Diff = {diff}
''',Configuration['OUTPUTLOG'])
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
   debug = False
   inner_min = minimum set of inner rings that should be fixed

OUTPUTS:
   Updates Configuration['INNER_FIX']

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
   Unspecified

NOTE:

'''
def find_program(name,search):
    found = False
    while not found:
        try:
            run = subprocess.Popen([name], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
            run.stdout.close()
            run.stderr.close()
            os.kill(run.pid, signal.SIGKILL)
            found = True
        except:
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

# Function to get the amount of inner rings to fix
def get_inner_fix(Configuration,Tirific_Template, debug =False):
    if debug:
        print_log(f'''GET_INNER_FIX: Attempting to get the inner rings to be fixed.
''',Configuration['OUTPUTLOG'], debug = True)
    sbr_av = np.array([(float(x)+float(y))/2. for x,y  in zip(Tirific_Template['SBR'].split(),Tirific_Template['SBR_2'].split())],dtype = float)
    column_levels = columndensity(Configuration,sbr_av, arcsquare = True, debug = debug)
    column_levels[0]= 1e21
    tmp = np.where(column_levels > 1e20)[0]
    if Configuration['OUTER_RINGS_DOUBLED']:
        inner_min =set_limits(Configuration['NO_RINGS']/3.,5.,11)
    else:
        inner_min =set_limits(Configuration['NO_RINGS']/5.,4.,Configuration['NO_RINGS']*0.7)
    inner_min = int(set_limits(np.floor(tmp[-1]/1.5-1), inner_min, Configuration['NO_RINGS']*0.9))


    tiltogram = make_tiltogram(Configuration,Tirific_Template,debug=debug)
    check_tiltogram(Configuration,tiltogram,inner_min=inner_min,debug=debug)
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
    debug = False

 OUTPUTS:
    Update INNER_FIX in the configuration

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    !!!!!!!!!!!!!This appears to currently not be working well.
'''

def get_kinematical_center(Configuration,map,angle,center= [0.,0],debug=False ):
    if np.sum(center) == 0.:
        center = [len(map[0,:])/2.,len(map[:,0])/2.]
    if angle > 180.:
        angle= angle-180.

    if angle < 90:
        angle = angle+90.
    else:
        angle= angle-90.
    buffer = int(round(np.mean(Configuration['BEAM_IN_PIXELS'][:2])/3.))
    min_axis_prof, min_axis, min_res = get_profile(Configuration,map,angle,center= center,debug=debug)
    found_vsys = [np.nanmean(min_axis_prof),0,0]
    found_diff=  abs(np.nanmin(min_axis_prof)-np.nanmax(min_axis_prof))+ np.nansum([abs(x-found_vsys[0]) for x in min_axis_prof])
    for x in range(-buffer,buffer):
        for y in range(-buffer,buffer):
                var_center = [int(round(center[0]+x)),int(round(center[1]+y))]
                min_axis_prof, min_axis, min_res = get_profile(Configuration,map,angle,center= var_center,debug=debug)
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
    debug = False
    center= [0.,0.]
    default is the cenetr of the map given in pixels

 OUTPUTS:
    profile, the corresponding axis in pixel size, resolution of a step on the axis

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def get_profile(Configuration,map,angle,center= [0.,0.],debug=False):
    if np.sum(center) == 0.:
        center = [len(map[0,:])/2.,len(map[:,0])/2.]
    x1,x2,y1,y2 = obtain_border_pix(Configuration,angle,center,debug=debug)
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
    debug = False
    center= [0.,0.]
    default is the cenetr of the map given in pixels

 OUTPUTS:
    profile, the corresponding axis in pixel size, resolution of a step on the axis

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_ring_weights(Configuration,Tirific_Template,debug = False):
    if debug:
        print_log(f'''GET_RING_WEIGTHS: Getting the importance of the rings in terms of SBR.
''',Configuration['OUTPUTLOG'], debug = True)
    sbr = np.array(get_from_template(Configuration,Tirific_Template, ["SBR",f"SBR_2"]),dtype=float)
    systemic = np.array(get_from_template(Configuration,Tirific_Template, ["VSYS"]),dtype=float)
    systemic = systemic[0,0]
    radii,cut_off_limits = sbr_limits(Configuration, systemic= systemic , debug = debug)
    weights= [[],[]]
    for i in [0,1]:
        weights[i] = [set_limits(x/y,0.1,10.) for x,y in zip(sbr[i],cut_off_limits)]
        weights[i] = weights[i]/np.nanmax(weights[i])
        weights[i][0:2] = np.nanmin(weights[i])
        weights[i] = [set_limits(x,0.1,1.) for x in weights[i]]
    if debug:
        print_log(f'''GET_RING_WEIGTHS: Obtained the following weights.
{'':8s}{weights}
''',Configuration['OUTPUTLOG'])
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
    debug = False

 OUTPUTS:
    numpy array with the weight normalized to the the maximum.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    Weight 1 is most important, weight 0. least important.
    Errors should be divided by these weights to reflect the importance
'''

def get_usage_statistics(Configuration,process_id, debug = False):
    result = subprocess.check_output(['top',f'-p {process_id}','-d 1','-n 1'])
    #result = subprocess.check_output(['ps','u'])
    lines = result.decode('utf8').split('\n')
    column_names = [x.upper() for x in lines[6].strip().split()]
    if debug:
        print_log(f'''{'':8s}GET_usage_statistics: We extracted the following column names {column_names}
''',Configuration['OUTPUTLOG'],debug=True)
    CPU = float(0.)
    mem=float(0.)
    column_var = [x for x in lines[7].strip().split()]
    if debug:
        print_log(f'''{'':8s}GET_usage_statistics: We extracted the following variables {column_var}
''',Configuration['OUTPUTLOG'])
    try:
        if int(column_var[column_names.index('PID')]) == int(process_id):
            CPU = float(column_var[column_names.index('%CPU')])
            mem = float(column_var[column_names.index('RES')])/1024**2
    except:
        #if the PID is not numeric it got merged with the crap in shiftcentercounter
        try:
            if column_var[column_names.index('COMMAND')-1] == 'tirific':
                CPU = float(column_var[column_names.index('%CPU')-1])
                mem = float(column_var[column_names.index('RES')-1])/1024**2
        except:
            pass
    return CPU,mem

get_usage_statistics.__doc__ =f'''
 NAME:
    get_usage_statistics
 PURPOSE:
    use top to get the current CPU and memory usage of tirific

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    process_id = process id of the tirific currently running.

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    CPU = The current CPU usage
    mem = current memory usage in Mb

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_vel_pa(Configuration,velocity_field,center= [0.,0.], debug =False):
    if debug:
        print_log(f'''GET_VEL_PA: This is the center we use {center}
''',Configuration['OUTPUTLOG'],debug = True)
    # because python is stupid the ceneter should be reversed to match the map
    center.reverse()
    if np.sum(center) == 0.:
        center = [x/2. for x in velocity_field.shape]

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

    if debug:
        print_log(f'''GET_VEL_PA: This is the location of the maximum {max_pos} and minimum {min_pos}
''',Configuration['OUTPUTLOG'])
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

    if np.degrees(abs(pa-pa_from_max)) > 170. or   np.degrees(abs(pa-pa_from_min)) > 170:
        pa = pa
    else:
        pa = np.nanmean([pa,pa_from_max,pa_from_min])
    center.reverse()
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
    debug = False
    center = center of the galaxy

 OUTPUTS:
    mean and standard deviation of the pa from minimum value to center, maximum to center and minimum to maximum

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

# A simple function to return the line numbers in the stack from where the functions are called
def linenumber(debug=False):
    line = []
    for key in stack():
        if key[1] == 'main.py':
            break
        if key[3] != 'linenumber' and key[3] != 'print_log':
            file = key[1].split('/')
            line.append(f"In the function {key[3]} at line {key[2]} in file {file[-1]}")
    if len(line) > 0:
        if debug:
            line = ', '.join(line)+f'\n{"":8s}'
        else:
            line = f'{"":8s}'
    else:
        for key in stack():
            if key[1] == 'main.py':
                line = f"{'('+str(key[2])+')':8s}"
                break
    return line

linenumber.__doc__ =f'''
 NAME:
    linenumber

 PURPOSE:
    get the line number of the print statement in the main. Not sure how well this is currently working

 CATEGORY:
    support_functions

 INPUTS:

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    the line number of the print statement

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
    If debug = True the full stack of the line print will be given, in principle
    the first debug message in every function should set this to true and later messages not.
    !!!!Not sure whether currently the linenumber is produced due to the restructuring.
'''
def make_tiltogram(Configuration,Tirific_Template,debug =False):
    if debug:
        print_log(f'''MAKE_TILTOGRAM: Starting tiltogram.
''',Configuration['OUTPUTLOG'])
    pa_incl = np.array(get_from_template(Configuration,Tirific_Template,Variables=['PA','PA_2','INCL','INCL_2']),dtype=float)
    sbr = np.array(get_from_template(Configuration,Tirific_Template, ["SBR",f"SBR_2"]),dtype=float)
    systemic = np.array(get_from_template(Configuration,Tirific_Template, ["VSYS"]),dtype=float)
    systemic = systemic[0,0]
    radii,cut_off_limits = sbr_limits(Configuration, systemic= systemic , debug = debug)
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
            if debug:
                print_log(f'''MAKE_TILTOGRAM: For the cartesian coordinates we find in side {i}
{'':8s} x = {x[i]}
{'':8s} y = {y[i]}
{'':8s} z = {z[i]}
''',Configuration['OUTPUTLOG'])

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
    debug = False

 OUTPUTS:
    a multidimensionals array containig the tiltograms for both sides

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:

'''

def obtain_border_pix(Configuration,angle,center, debug = False):
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
    debug = False

 OUTPUTS:
    x,y
    pixel locations of how the line runs through the map

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    np.tan,np.radians,.reverse()

 NOTE:
'''

def obtain_ratios(Configuration, map, center, angles, noise = 0. ,debug = False):
    ratios = []
    max_extent = 0.
    for angle in angles:
        #major axis
        maj_profile,maj_axis,maj_resolution = get_profile(Configuration,map,angle,center=center,debug=debug)
        tmp = np.where(maj_profile > noise)[0]
        #gauss =fit_gaussian(maj_axis[tmp], maj_profile[tmp])
        #maj_gaus = gaussian_function(maj_profile, *gauss)
        #tmp = np.where(maj_gaus > 0.2*np.nanmax(maj_gaus))[0]
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
            min_profile,min_axis,min_resolution = get_profile(Configuration,map,angle+90,center=center,debug=debug)
        else:
            min_profile,min_axis,min_resolution = get_profile(Configuration,map,angle-90,center=center,debug=debug)
        tmp = np.where(min_profile > noise)[0]
        #gauss =fit_gaussian(min_axis[tmp], maj_profile[tmp])
        #min_gaus = gaussian_function(min_profile,*gauss)
        #tmp = np.where(min_gaus > 0.2*np.nanmax(min_gaus))[0]
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

       debug = False

 OUTPUTS:
         ratios =  ratios corresponding to the input angles
         max_extent = the maximum extend of any profile analysed (in degree)

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
      np.mean, ndimage.map_coordinates, np.linspace, np.vstack, np.cos, np.radians, np.sin

 NOTE:
'''

def print_log(log_statement,log, screen = False,debug = False):
    log_statement = f"{linenumber(debug=debug)}{log_statement}"
    if screen or not log:
        print(log_statement)
    if log:
        with open(log,'a') as log_file:
            log_file.write(log_statement)

print_log.__doc__ =f'''
 NAME:
    print_log
 PURPOSE:
    Print statements to log if existent and screen if Requested
 CATEGORY:
    support_functions

 INPUTS:
    log_statement = statement to be printed
    log = log to print to, can be None

 OPTIONAL INPUTS:
    debug = False

    screen = False
    also print the statement to the screen

 OUTPUTS:
    line in the log or on the screen

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    linenumber, .write

 NOTE:
    If the log is None messages are printed to the screen.
    This is useful for testing functions.
'''

def remove_inhomogeneities(Configuration,fits_map_in,inclination=30., pa = 90. , center = [0.,0.],WCS_center = True, iteration= 0 , debug=False):
    fits_map = copy.deepcopy(fits_map_in)
    if debug:
        print_log(f'''REMOVE_INHOMOGENEITIES: These are the values we get as input
{'':8s}Inclination = {inclination}
{'':8s}pa = {pa}
{'':8s}center = {center}, WCS = {WCS_center}
{'':8s}map shape = {np.shape(fits_map[0].data)}
''',Configuration['OUTPUTLOG'],debug = True)
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
    rot_map = rotateImage(Configuration,map,pa-90,[x,y],debug=debug)
    # deproject
    dep_map = deproject(Configuration,copy.deepcopy(rot_map),inclination,center = y, debug=debug)

    if debug:
        fits.writeto(f"{Configuration['FITTING_DIR']}rot_map_{int(iteration)}.fits",rot_map,fits_map[0].header,overwrite = True)
        fits.writeto(f"{Configuration['FITTING_DIR']}dep_map_{int(iteration)}.fits",dep_map,fits_map[0].header,overwrite = True)

    angles = np.linspace(5.,360.,71)
    minimum_map = copy.deepcopy(dep_map)
    for angle in angles:
        rot_dep_map =  rotateImage(Configuration,copy.deepcopy(dep_map),angle,[x,y],debug=debug)

        #tmp = np.where(rot_dep_map < minimum_map)[0]
        minimum_map[rot_dep_map < minimum_map] =rot_dep_map[rot_dep_map < minimum_map]
    clean_map = rotateImage(Configuration,deproject(Configuration,copy.deepcopy(minimum_map),inclination,center = y,invert= True,debug=debug),-1*(pa-90),[x,y],debug=debug)

    if debug:
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
    debug = False

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

def rename_fit_products(Configuration,stage = 'initial', fit_type='Undefined',debug = False):
    extensions = ['def','log','ps','fits']
    for filetype in extensions:
        if filetype == 'log':
            if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                os.system(f"cp {Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype} {Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev.{filetype} ")
        else:
            if filetype == 'def':
                if stage in ['run_ec','run_os','run_cc']:
                    Loopnr = f"Iteration_{Configuration['ITERATIONS']-1}"
                elif stage in ['after_cc','after_ec','after_os'] :
                    Loopnr = f"Iteration_{Configuration['ITERATIONS']}"
                elif fit_type == Configuration['USED_FITTING'] and stage in ['final_os']:
                    Loopnr = f"Smoothed_1"
                else:
                    Loopnr = 'Output_before_'+stage
                if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                    os.system(f"mv {Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype} {Configuration['FITTING_DIR']}{fit_type}/{fit_type}_{Loopnr}.{filetype}")

            elif os.path.exists(f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype}"):
                os.system(f"mv {Configuration['FITTING_DIR']}{fit_type}/{fit_type}.{filetype} {Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Prev.{filetype}")

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
    debug = False

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
def rotateImage(Configuration,image, angle, pivot, debug = False):
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
    debug = False

 OUTPUTS:
    The rotated image is returned

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Function to read FAT configuration file into a dictionary
def setup_configuration(cfg):
    if cfg.installation_check:
        cfg.fitting.fixed_parameters=['INCL','PA','SDIS']
        cfg.advanced.max_iterations= 1
        cfg.advanced.loops= 1
        cfg.input.main_directory= f'{cfg.input.main_directory}/FAT_Installation_Check/'
        cfg.output.output_quantity = 0
        cfg.fitting.fitting_stages = ['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC']
        cfg.cube_name = 'NGC_2903.fits'
        test_files = ['NGC_2903.fits','ModelInput.def']
        if not os.path.isdir(cfg.input.main_directory):
            os.mkdir(cfg.input.main_directory)
        else:
            for file in test_files:
                try:
                    os.remove(test_dir+file)
                except:
                    pass
        my_resources = import_res.files('pyFAT_astro.Installation_Check')
        for file in test_files:
            data = (my_resources / file).read_bytes()
            with open(cfg.input.main_directory+file,'w+b') as tmp:
                tmp.write(data)
    if cfg.cube_name:
        cfg.input.catalogue = None


    Configuration = Proper_Dictionary({})

    for key in cfg._content:
        input_key = getattr(cfg,key)
        if str(type(input_key)) == 'omegaconf.dictconfig.DictConfig':
            for sub_key in input_key._content:
                if str(key) == 'output' and str(sub_key) == 'catalogue':
                    Configuration['OUTPUT_CATALOGUE'] =  getattr(input_key,sub_key)
                elif str(key) == 'fitting' and str(sub_key) == 'fixed_parameters':
                    value = getattr(input_key,sub_key)
                    for req_key in ['Z0','XPOS','YPOS','VSYS']:
                        if req_key not in value:
                            value.append(req_key)
                    Configuration[str(sub_key).upper()] =  [value,value]
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
                'SOFIA_RAN', #Check if we have ran Sofia
                'NO_RADEC'
                ]
#
    for key in boolean_keys:
        Configuration[key] = False

    other_keys={'ID_NR': 'Unset', # ID of the galaxy in the catalogue , set from the catalogue at start of loop
               'SUB_DIR': 'Unset', # Name of the directory in which galaxy resides, set from the catalogue at start of loop
               'FITTING_DIR': 'Unset', # Full path of the directory in which the fitting takes place, set at start of loop
               'BASE_NAME': 'Unset', #Basename for FAT products, typically {input_cube}_FAT, set at start of loop
               'LOG_DIR': 'Unset', #Directory to put log files from run, set at start of loop

               'PREP_END_TIME': 'Not completed',
               'START_TIME':'Not completed',
               'END_TIME':'Not completed',
               'OUTPUTLOG':'Not set yet',
               'RUN_COUNTER': 0,
               'ITERATIONS': 0,
               'CURRENT_STAGE': 'initial', #Current stage of the fitting process, set at switiching stages
               'USED_FITTING': None,
               'TIRIFIC_PID': 'Not Initialized', #Process ID of tirific that is running
               'FINAL_COMMENT': "This fitting stopped with an unregistered exit.",

               'MAX_SIZE_IN_BEAMS': 30, # The galaxy is not allowed to extend beyond this number of beams in radius, set in check_source
               'MIN_SIZE_IN_BEAMS': 0., # Minimum allowed radius in number of beams of the galaxy, set in check_source
               'SIZE_IN_BEAMS': 0, # The radius of the galaxy in number of beams, adapted after running Sofia
               'NO_RINGS': 0., # The number of rings in the fit,
               'LAST_RELIABLE_RINGS': [0.,0.], # Location of the rings where the SBR drops below the cutoff limits, adapted after every run. Should only be modified in check_size
               'LIMIT_MODIFIER': [1.], #Modifier for the cutoff limits based on the inclination , adapted after every run.
               'OLD_RINGS': [], # List to keep track of the ring sizes that have been fitted.

               'NO_POINTSOURCES': 0. , # Number of point sources, set in run_tirific

               'INNER_FIX': [4.,4.], #Number of rings that are fixed in the inner part for the INCL and PA, , adapted after every run in get_inner_fix in support_functions and for both sides
               'WARP_SLOPE': [0.,0.], #Ring numbers from which outwards the warping should be fitted as a slope, set in get_warp_slope in modify_template
               'OUTER_SLOPE_START': 1, # Ring number from where the RC is fitted as a slope
               'RC_UNRELIABLE': 1, # Ring number from where the RC values are set flat. Should only be set in check_size

               'NOISE': 0. , #Noise of the input cube in Jy/beam, set in read_cube
               'BEAM_IN_PIXELS': [0.,0.,0.], #FWHM BMAJ, BMIN in pixels and total number of pixels in beam area, set in main
               'BEAM': [0.,0.,0.], #  FWHM BMAJ, BMIN in arcsec and BPA, set in main
               'BEAM_AREA': 0., #BEAM_AREA in arcsec set in main
               'NAXES': [0.,0.,0.], #  Size of the cube in pixels x,y,z arranged like sane people not python, set in main
               'NAXES_LIMITS': [[0.,0.],[0.,0.],[0.,0.]], #  Size of the cube in degree and km/s,  x,y,z arranged like sane people not python, set in main updated in cut_cubes
               'MAX_ERROR': {}, #The maximum allowed errors for the parameters, set in main derived from cube
               'MIN_ERROR': {}, #The minumum allowed errors for the parameters, initially set in check_source but can be modified through out INCL,PA,SDSIS,Z0 errors change when the parameters is fixed or release
               'CHANNEL_WIDTH': 0., #Width of the channel in the cube in km/s, set in main derived from cube
               'PIXEL_SIZE': 0., #'Size of the pixels in degree'
               }
    #####!!!!!!!!!!!!!!!!!!!!!!!!!!!! The use of min_error is not full implemented yet!!!!!!!!!!!!!!!!!!!!!!!!!
    for key in other_keys:
        Configuration[key] = other_keys[key]

# The parameters that need boundary limits are set here
    boundary_limit_keys = ['PA','INCL', 'SDIS', 'Z0','VSYS','XPOS','YPOS']
    for key in boundary_limit_keys:
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
        while not os.path.exists(Configuration['CATALOGUE']):
            Configuration['CATALOGUE'] = input(f'''
Your input catalogue ({Configuration['CATALOGUE']}) does not exist.
Please provide the correct file name.
: ''')
    #Make sure there is only one Fit_ stage

    # Make sure all selected stages exist
    possible_stages = ['Fit_Tirific_OSC','Create_FAT_Cube','Run_Sofia','Existing_Sofia','Sofia_Catalogue']
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

    Configuration['FITTING_STAGES'] =approved_stages
    if 'sofia_catalogue' in Configuration['FITTING_STAGES'] and 'existing_sofia' in Configuration['FITTING_STAGES']:
        Configuration['FITTING_STAGES'].remove('existing_sofia')
    if ('sofia_catalogue' in Configuration['FITTING_STAGES'] or 'existing_sofia' in Configuration['FITTING_STAGES']) and 'run_sofia' in Configuration['FITTING_STAGES']:
        Configuration['FITTING_STAGES'].remove('run_sofia')
    if 'sofia_catalogue' in Configuration['FITTING_STAGES'] and 'create_fat_cube' in Configuration['FITTING_STAGES']:
        Configuration['FITTING_STAGES'].remove('create_fat_cube')

    if 'run_sofia' in Configuration['FITTING_STAGES'] and 'existing_sofia' in Configuration['FITTING_STAGES']:
        raise InputError(f"You are both providing existing sofia input and ask for sofia to be ran. This won't work exiting.")



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


    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['OUTPUT_QUANTITY'] == 5:
        Configuration['OUTPUT_QUANTITY'] = 4

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
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:
    Configuration = dictionary with the config file input

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: Only this function can add entries to Original_Configuration or Configuration
'''



def sbr_limits(Configuration, systemic= 100. , debug = False):
    radii = set_rings(Configuration,debug=debug)
    if debug:
        print_log(f'''SBR_LIMITS: Got {len(radii)} radii
''',Configuration['OUTPUTLOG'], debug=True)
    level = Configuration['NOISE']*1000
    noise_in_column = columndensity(Configuration,level,systemic = systemic)
    J2007col=9.61097e+19
    ratio=(noise_in_column/J2007col)**0.5
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
    if len(Configuration['LIMIT_MODIFIER']) == 1:

        sbr_ring_limits= sbr_ring_limits*Configuration['LIMIT_MODIFIER']
    else:
        mod_list = list(Configuration['LIMIT_MODIFIER'])
        while len(mod_list) < len(sbr_ring_limits):
            mod_list.append(Configuration['LIMIT_MODIFIER'][-1])
        Configuration['LIMIT_MODIFIER'] = np.array(mod_list,dtype=float)
        sbr_ring_limits=[x*y for x,y in zip(sbr_ring_limits,Configuration['LIMIT_MODIFIER'])]
    if debug:
        print_log(f'''SBR_LIMITS: Retrieved these radii and limits:
{'':8s}{radii}
{'':8s}{sbr_ring_limits}
''',Configuration['OUTPUTLOG'])
    return radii,sbr_ring_limits

sbr_limits.__doc__ =f'''
 NAME:
    sbr_limits

 PURPOSE:
    Calculate the sbr limits below which rings can not be trusted.

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

    systemic= 100.
    systemic velocity

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
def set_limit_modifier(Configuration,Inclination, debug= False):
    if debug:
        print_log(f'''SET_LIMIT_MODIFIER: Checking the limit modifier.
''', Configuration['OUTPUTLOG'], debug=True)
    if not Inclination.shape:
        Inclination = [Inclination]
    modifier_list = []
    # Correction because the noise applies to non-face on rings while the SBR is face on,correction is normalized to the average inclination
    # The square root is because of the square root between the noise in cube and the noise in J2007
    # The lower limit corresponds to a inclination 80 above which the cos becomes too steep
    for inc in Inclination:
        if inc > 90.:
            inc  = 180.-inc
        if inc < 0.:
            inc= abs(inc)
        modifier_list.append(set_limits(np.sqrt(np.cos(np.radians(inc))/np.cos(np.radians(60.))),0.75,2.))
    if Configuration['OUTER_RINGS_DOUBLED']:
        if len(modifier_list) > 10:
            modifier_list[10:]= np.sqrt(modifier_list[10:])
    Configuration['LIMIT_MODIFIER'] = np.array(modifier_list,dtype=float)
    print_log(f'''SET_LIMIT_MODIFIER: We updated the LIMIT_MODIFIER to {Configuration['LIMIT_MODIFIER']}.
''', Configuration['OUTPUTLOG'])

set_limit_modifier.__doc__ =f'''
 NAME:
    set_limit_modifier

 PURPOSE:
    Write to the Configuration a value or list of values to modify inclination dependent values. Foremost the sbr_limits

 CATEGORY:
    support_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Inclination = Inclination of the model, can be singular or for all rings.

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    'LIMIT_MODIFIER' is updated in the Configuration.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_limits(value,minv,maxv,debug = False):
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
    debug = False

 OUTPUTS:
    the limited Value

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_ring_size(Configuration, debug = False, size_in_beams = 0., check_set_rings = False):
    if size_in_beams == 0.:
        size_in_beams =  Configuration['SIZE_IN_BEAMS']
    ring_size = Configuration['RING_SIZE']
    no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
    if debug:
        print_log(f'''SET_RING_SIZE: Starting with the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
''', Configuration['OUTPUTLOG'],debug=True)

    while ring_size > 0.5 and  no_rings < 8.:
        previous_ringsize = ring_size
        ring_size = set_limits(ring_size/1.5,0.5,float('NaN'),debug=debug)
        no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
        print_log(f'''SET_RING_SIZE: Because we had less than four rings we have reduced the ring size from {previous_ringsize} to {ring_size}
''',Configuration['OUTPUTLOG'])

    while no_rings < Configuration['MINIMUM_RINGS'] and not size_in_beams >=  Configuration['MAX_SIZE_IN_BEAMS']:
        size_in_beams = set_limits(size_in_beams+1.*ring_size,1, Configuration['MAX_SIZE_IN_BEAMS'])
        no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
        print_log(f'''SET_RING_SIZE: The initial estimate is too small to fit adding a ring to it.
''',Configuration['OUTPUTLOG'])

    if check_set_rings:
        if debug:
            print_log(f'''SET_RING_SIZE: Setting the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
{'':8s}NO_RINGS = {no_rings}
''', Configuration['OUTPUTLOG'])
        return size_in_beams,ring_size,int(no_rings)
    else:
        if debug:
            print_log(f'''SET_RING_SIZE: Setting the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
{'':8s}NO_RINGS = {no_rings}
''', Configuration['OUTPUTLOG'])
        Configuration['NO_RINGS'] = int(no_rings)
        Configuration['SIZE_IN_BEAMS'] = size_in_beams
        Configuration['RING_SIZE'] = ring_size
        if Configuration['NO_RINGS'] < Configuration['MINIMUM_RINGS']:
            print_log(f'''SET_RING_SIZE: With a ring size of {Configuration['RING_SIZE']} we still only find {Configuration['NO_RINGS']}.
    {"":8s}SET_RING_SIZE: This is not enough for a fit.
    ''',Configuration['OUTPUTLOG'],screen=True)
            raise SmallSourceError('This source is too small to reliably fit.')

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
    debug = False

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

def set_rings(Configuration,ring_size = 0. , debug = False):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    if debug:
        print_log(f'''SET_RINGS: Starting with the following parameters.
{'':8s}RING_SIZE = {ring_size}
''', Configuration['OUTPUTLOG'],debug=True)
    no_rings = calc_rings(Configuration,debug=debug)
    if debug:
        print_log(f'''SET_RINGS: We find {no_rings} rings.
''', Configuration['OUTPUTLOG'])
    if Configuration['OUTER_RINGS_DOUBLED']:
        print_log(f'''SET_RINGS: This is a large galaxy. Therefore we use twice the ring size in the outer parts.
''',Configuration['OUTPUTLOG'],screen =True)
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
    if debug:
        if Configuration['OUTER_RINGS_DOUBLED']:
            req_outer_ring = 2.*Configuration['BEAM'][0]*ring_size
        else:
            req_outer_ring = Configuration['BEAM'][0]*ring_size
        print_log(f'''SET_RINGS: Got the following radii.
{'':8s}{radii}
{'':8s}We should have {Configuration['NO_RINGS']} or we incorrectly updated and should have {no_rings}
{'':8s}We have {len(radii)} rings.
{'':8s}The last ring should be around {Configuration['BEAM'][0]*Configuration['SIZE_IN_BEAMS']}
{'':8s}The rings should be size {Configuration['BEAM'][0]*ring_size} and outer rings {req_outer_ring}
{'':8s}They are {radii[3]-radii[2]} and {radii[-1]-radii[-2]}
''', Configuration['OUTPUTLOG'])
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
    debug = False

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

def sofia_output_exists(Configuration,Fits_Files, debug = False):
    if debug:
        print_log(f'''SOFIA_OUTPUT_EXISTS: Starting check
''', Configuration['OUTPUTLOG'],debug = True)
    req_files= ['MOMENT1','MOMENT0','MASK']
    for file in req_files:
        if os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]):
            continue
        else:
            log_statement = f"SOFIA_OUTPUT_EXISTS: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]} is not found."
            print_log(log_statement, Configuration['OUTPUTLOG'],screen =True)
            raise FileNotFoundError(log_statement)

    if not os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'):
        log_statement = f"SOFIA_OUTPUT_EXISTS: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'} is not found."
        print_log(log_statement, Configuration['OUTPUTLOG'],screen =True)
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
    debug = False

 OUTPUTS:
    Raises a FileNotFoundError if the files are not found.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
