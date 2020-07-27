#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from inspect import getframeinfo,stack
from scipy.optimize import curve_fit

import os
import traceback
import numpy as np
import copy
import warnings
import re

class SupportRunError(Exception):
    pass

# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
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
                    "----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")






#batch convert types
def convert_type(array, type = 'float'):
    if type =='int':
        return  [int(x) for x in array]
    elif type =='str':
        return  [str(x) for x in array]
    else:
        return  [float(x) for x in array]

convert_type.__doc__ = '''
;+
; NAME:
;       def convert_type(array, type = 'float'):
;
; PURPOSE:
;       Convert a list of variables from one type to another and be able to have them unpack into single varaiable again.
;
; CATEGORY:
;       Support
;
; INPUTS:
;     Array to convert
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''

#Function to convert column densities
# levels should be n mJy/beam when flux is given
def columndensity(levels,systemic = 100.,beam=[1.,1.],channel_width=1.,column= False,arcsquare=False,solar_mass_input =False,solar_mass_output=False):
    beam=np.array(beam)
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg

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
columndensity.__doc__ = '''
;+
; NAME:
;       columndensity(levels,systemic = 100.,beam=[1.,1.],channel_width=1.,column= False,arcsquare=False,solar_mass_input =False,solar_mass_output=False)
;
; PURPOSE:
;       Convert the various surface brightnesses to other values
;
; CATEGORY:
;       Support
;
; INPUTS:
;       levels = the values to convert
;       systemic = the systemic velocity of the source
;        beam  =the beam in arcse
;       channelwidth = width of a channel in km/s
;     column = if true input is columndensities
;
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(RAin,DECin,invert=False, colon=False):
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


# function for converting kpc to arcsec and vice versa

def convertskyangle(angle, distance=1., unit='arcsec', distance_unit='Mpc', physical=False):
    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.\n')
        print('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n')
        raise SupportRunError('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.')
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n')
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
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n')
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')

        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc


def gaussian_function(axis,peak,center,sigma):
    return peak*np.exp(-(axis-center)**2/(2*sigma**2))

def fit_gaussian(x,y, covariance = False):
    # First get some initial estimates
    est_peak = np.nanmax(y)
    est_center = float(x[np.where(y == est_peak)])
    est_sigma = np.nansum(y*(x-est_center)**2)/np.nansum(y)
    gauss_parameters, gauss_covariance = curve_fit(gaussian_function, x, y,p0=[est_peak,est_center,est_sigma])
    if covariance:
        return gauss_parameters, gauss_covariance
    else:
        return gauss_parameters

# A simple function to return the line numbers in the stack from where the functions are called
def linenumber(debug=False):
    line = []
    for key in stack():
        if key[1] == 'FAT.py':
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
            if key[1] == 'FAT.py':
                line = f"{'('+str(key[2])+')':8s}"
                break
    return line


def print_log(log_statement,log, screen = False,debug = False):
    log_statement = f"{linenumber(debug=debug)}{log_statement}"
    if screen or not log:
        print(log_statement)
    if log:
        with open(log,'a') as log_file:
            log_file.write(log_statement)

print_log.__doc__ = '''
;+
; NAME:
;       print_log
;
; PURPOSE:
;       Print statements to log if existent and screen if Requested
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       cleanup(Configuration)
;
; INPUTS:
;     Configuration = Structure that has the current fitting directory and expected steps
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
; MODULES CALLED:
;       os
;
; EXAMPLE:
;

'''

def rename_fit_products(Configuration,stage = 'initial', fit_stage='Undefined_Stage'):
    extensions = ['def','log','ps','fits']
    for filetype in extensions:
        if filetype == 'log':
            if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype}"):
                os.system(f"cp {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype} ")
        else:
            if filetype == 'def':
                if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype}"):
                    os.system(f"mv {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev2.{filetype} ")
            if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype}"):
                os.system(f"mv {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype}")


rename_fit_products.__doc__ = '''
;+
; NAME:
;       rename_fit_products(Configuration,stage)
;
; PURPOSE:
;       rename the tirific product from the previous stage.
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''


def sbr_limits(Configuration,hdr, systemic= 100.):
    radii = set_rings(Configuration, hdr)
    level = hdr['FATNOISE']*1000
    bm = [hdr['BMAJ']*3600.,hdr['BMIN']*3600.]
    noise_in_column = columndensity(level,beam = bm,systemic = systemic,channel_width=hdr['CDELT3']/1000.)
    J2007col=9.61097e+19
    ratio=(noise_in_column/J2007col)**0.5
    beamsolid=(np.pi*bm[0]*bm[1])/(4.*np.log(2.))
    ringarea= [0 if radii[0] == 0 else np.pi*((radii[0]+radii[1])/2.)**2]
    ringarea = np.hstack((ringarea,
                         [np.pi*(((y+z)/2.)**2-((y+x)/2.)**2) for x,y,z in zip(radii,radii[1:],radii[2:])],
                         [np.pi*((radii[-1]+0.5*(radii[-1]-radii[-2]))**2-((radii[-1]+radii[-2])/2.)**2)]
                         ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sbr_ring_limits=9e-4*(ringarea/beamsolid)**(-0.82)*ratio
    if ringarea[0] == 0.:
         sbr_ring_limits[0]=0.
    if len(Configuration['LIMIT_MODIFIER']) == 1:
        sbr_ring_limits= sbr_ring_limits*Configuration['LIMIT_MODIFIER']
    else:
        sbr_ring_limits=[x*y for x,y in zip(sbr_ring_limits,Configuration['LIMIT_MODIFIER'])]

    return radii,sbr_ring_limits


sbr_limits.__doc__ = '''
;+
; NAME:
;       sbr_limits(Configuration,hdr)
;
; PURPOSE:
;       Create the radii at which to evaluate  the model and a corresponding array that has the sbr reliability sbr_limits
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''

def set_limits(value,minv,maxv):
    if value < minv:
        return minv
    elif value > maxv:
        return maxv
    else:
        return value

set_limits.__doc__ = '''
;+
; NAME:
;       set_limits(value,min,max)
;
; PURPOSE:
;       Make sure Value is between min and max else set to min when smaller or max when larger.
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     value
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
#simple function keep track of how to modify the edge limits
def set_limit_modifier(Configuration,Inclination):

    if not Inclination.shape:
        Inclination = [Inclination]
    modifier_list = []
    for inc in Inclination:
        if 40 < inc < 50:
            modifier_list.append(set_limits(1.+(50-(inc)*0.05),1,2.5))
        elif inc < 40:
            modifier_list.append(set_limits(np.sin(np.radians(75.))/np.sin(np.radians(inc))),1.,2.5)
        else:
            modifier_list.append(1.)
    if Configuration['OUTER_RINGS_DOUBLED']:
        if len(modifier_list) > 10:
            modifier_list[10:]= np.sqrt(modifier_list[10:])

    Configuration['LIMIT_MODIFIER'] = np.array(modifier_list,dtype=float)

def set_rings(Configuration,hdr):
    bmaj = float(hdr['BMAJ'])*3600.
    if Configuration['OUTER_RINGS_DOUBLED']:
        radii = [0.,1./5.*bmaj]
        radii = np.hstack((radii,(np.linspace(bmaj*Configuration['RING_SIZE'],bmaj*10.*Configuration['RING_SIZE'], \
                                        10)+1./5*bmaj)))
        radii = np.hstack((radii,(np.linspace(bmaj*11.*Configuration['RING_SIZE'],bmaj \
                                      *Configuration['NO_RINGS']*Configuration['RING_SIZE'],int((Configuration['NO_RINGS']-9)/2.))+1./5*bmaj)))
    else:
        radii = [0.,1./5.*bmaj]
        radii = np.hstack((radii,(np.linspace(bmaj*Configuration['RING_SIZE'],bmaj*Configuration['NO_RINGS']*Configuration['RING_SIZE'], \
                                        Configuration['NO_RINGS'])+1./5.*bmaj)))
    return np.array(radii,dtype = float)


set_rings.__doc__ = '''
;+
; NAME:
;       set_rings(Configuration,hdr)
;
; PURPOSE:
;       make an array that chas all the ring radii in ARCSEC
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
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
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
