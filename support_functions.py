#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
import numpy as np # Used in convertskyangle and columndensity and
import copy # Used in columndensities, convertRADEC
import re #used in convertRADEC
from scipy import ndimage
from scipy.optimize import curve_fit # Used in Fit_Gaussian
import scipy.constants # used in sofia_pix2world
from astropy import wcs #used throughout
from astropy.io import fits
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.axes as maxes
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
        sys.exit()
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
            sys.exit()

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
            sys.exit()
        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc
#Function to extract PV-Diagrams from a datacube.
def extract_pv(cube_in,angle,center=[-1,-1,-1],finalsize=[-1,-1],convert=-1):
    cube = copy.deepcopy(cube_in)
    hdr = copy.deepcopy(cube[0].header)
    data = copy.deepcopy(cube[0].data)
    coordinate_frame = wcs.WCS(hdr)


    nz, ny, nx = data.shape

    if finalsize[0] != -1:
        if finalsize[0] >nz:
            finalsize[0] = nz
        if finalsize[1] >nx:
            finalsize[1] = nx
    # if the center is not set assume the crval values
    if center[0] == -1:
        center = [hdr['CRVAL1'],hdr['CRVAL2'],hdr['CRVAL3']]

    xcenter,ycenter,zcenter = coordinate_frame.wcs_world2pix(center[0], center[1], center[2], 0.)

    nextcube = rotateCube(data, angle - 90 + 180, [xcenter,ycenter])

    # Let's slice out the PV-Diagram
    #sliced = np.zeros((nz, nx), dtype=float)
    if finalsize[0] == -1:
        print("What is happening here")
        sliced = np.zeros((nz, nx), dtype=float)
        sliced = nextcube[:, int(ycenter), :]
        # then lets update the header
        # As python is stupid making a simple copy will mean that these changes are still applied to hudulist
        hdr['NAXIS2'] = nz
        hdr['NAXIS1'] = nx
        hdr['CRPIX2'] = hdr['CRPIX3']
        if convert !=-1:
            hdr['CRVAL2'] = hdr['CRVAL3']/convert
        else:
            hdr['CRVAL2'] = hdr['CRVAL3']
        hdr['CRPIX1'] = xcenter+1
    else:
        print("This is the attempted  final size of xv {}".format(finalsize))
        print("This is the achieved final size of xv {} {}".format(int(finalsize[0]),int(finalsize[1])))

        sliced = np.zeros((int(finalsize[0]),int(finalsize[1])),dtype=float)
        zstart = int(nz/2.-finalsize[0]/2.)
        zend = int(nz/2.+finalsize[0]/2.)
        xstart = int(xcenter-finalsize[1]/2.)
        xend = int(xcenter+finalsize[1]/2.)
        print("WTF {} {} {}".format(xstart,xend,nextcube.shape))
        sliced =  nextcube[zstart:zend, int(ycenter), xstart:xend]
        hdr['NAXIS2'] = int(finalsize[0])
        hdr['NAXIS1'] = int(finalsize[1])
        print(hdr['NAXIS1'])
        hdr['CRPIX2'] = hdr['CRPIX3']-int(nz/2.-finalsize[0]/2.)
        if convert !=-1:
            hdr['CRVAL2'] = hdr['CRVAL3']/convert
        else:
            hdr['CRVAL2'] = hdr['CRVAL3']

        hdr['CRPIX1'] = int(finalsize[1]/2.)+1
    if convert !=-1:
        hdr['CDELT2'] = hdr['CDELT3']/convert
    else:
        hdr['CDELT2'] = hdr['CDELT3']
    hdr['CTYPE2'] = hdr['CTYPE3']
    try:
        if hdr['CUNIT3'].lower() == 'm/s' and convert == -1:
            hdr['CDELT2'] = hdr['CDELT3']/1000.
            hdr['CRVAL2'] = hdr['CRVAL3']/1000.
            hdr['CUNIT2'] = 'km/s'
            del (hdr['CUNIT3'])
        elif  convert != -1:
            del (hdr['CUNIT3'])
            del (hdr['CUNIT2'])
        else:
            hdr['CUNIT2'] = hdr['CUNIT3']
            del (hdr['CUNIT3'])
    except:
        print("No units")
    del (hdr['CRPIX3'])
    del (hdr['CRVAL3'])
    del (hdr['CDELT3'])
    del (hdr['CTYPE3'])

    del (hdr['NAXIS3'])
    hdr['CRVAL1'] = 0.
    hdr['CDELT1'] = abs(hdr['CDELT1'] * 3600.)
    hdr['CTYPE1'] = 'OFFSET'
    hdr['CUNIT1'] = 'ARCSEC'
    # Then we change the cube and rteturn the PV construct
    print(hdr['NAXIS1'],sliced.shape)
    cube[0].header = hdr
    cube[0].data = sliced
    print(cube[0].header['NAXIS1'],cube[0].data.shape)
    return cube
#Function to fit a Gaussian to an array of invalues
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

#function to fit a plane to a numpy array image. Found in
def fit_plane(image):
 # Following https://stackoverflow.com/questions/35005386/fitting-a-plane-to-a-2d-array
    ysize,xsize = image.shape
    X1,X2 =np.mgrid[:ysize,:xsize]
    #Linear regression fitting
    X = np.hstack((np.reshape(X1, (ysize * xsize, 1)), np.reshape(X2, (ysize * xsize, 1))))
    X = np.hstack((np.ones((ysize * xsize, 1)), X))
    YY = np.reshape(image, (ysize * xsize, 1))
    theta = np.dot(np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()),YY)

    plane = np.reshape(np.dot(X, theta), (ysize, xsize));
    return plane

def nanfit_plane(image):
 # Following https://stackoverflow.com/questions/35005386/fitting-a-plane-to-a-2d-array
    ysize,xsize = image.shape
    X1,X2 =np.mgrid[:ysize,:xsize]
    #Linear regression fitting
    X = np.hstack((np.reshape(X1, (ysize * xsize, 1)), np.reshape(X2, (ysize * xsize, 1))))
    X = np.hstack((np.ones((ysize * xsize, 1)), X))
    YY = np.reshape(image, (ysize * xsize, 1))
    YY_m =np.ma.array(YY, mask= np.isnan(YY))
    XX = np.dot(np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose())
    theta = np.ma.dot(XX,YY_m)
    theta = np.array(theta,dtype=np.double)
    plane = np.reshape(np.dot(X, theta), (ysize, xsize));
    return plane


# Function to input a boolean answer
def get_bool(print_str="Please type True or False",default=True):
    invalid_input = True
    while invalid_input:
        inp = input(print_str)
        if inp == "":
            if default:
                return True
            else:
                return False
        elif inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes":
            return True
        elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no":
            return False
        else:
            print("Error: the answer must be true/false or yes/no.")

#Function to read in a simple catalogue
def get_catalog(filename, unpack =True):
    tmp = open(filename, 'r')
    unarranged = tmp.readlines()
    # Separate the keyword names
    for line in unarranged:
        var_concerned = str(line.split('=')[0].strip().lower())
        if var_concerned == 'catalogue':
            input=str(line.split('=')[1].strip().split('/')[-1])
        if var_concerned == 'outputcatalogue':
            result=str(line.split('=')[1].strip().split('/')[-1])

    if unpack:
        return input,result
    else:
        return [input,result]

#Function for loading the variables of a tirific def file into a set of variables to be used
def load_tirific(filename,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True ):
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    try:
        tmp = open(filename, 'r')
    except FileNotFoundError:
        print("LOAD_TIRIFIC: We could not find the file returning an empty array")
        outputarray = np.zeros((3, len(Variables)), dtype=float)
        outputarray[:,:] = float(0.)
        if unpack:
            return (*outputarray.T,)
        else:
            return outputarray
    numrings = [int(e.split('=')[1].strip()) for e in tmp.readlines() if e.split('=')[0].strip().upper() == 'NUR']
    tmp.seek(0)
    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
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
                varpos = np.where(var_concerned[2:] == Variables)[0]
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
        #load barollo output in the same way as tirfic
def load_barollo(filename,Variables = ['DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2'],
                 unpack = True ):

    tir_bar_translate = {
    'RADI': 'RAD',
    'DISTANCE' : 'DISTANCE',
    'NUR' : 'LENGTH',
    'VROT': 'VROT',
    'VROT_2': 'VROT',
    'VROT_ERR': 'NULL',
    'VROT_2_ERR': 'NULL',
    'Z0': 'Z0',
    'Z0_2': 'Z0',
    'SBR': 'SURFDENS',
    'SBR_ERR': 'ERR_SD',
    'SBR_2': 'SURFDENS',
    'SBR_2_ERR': 'ERR_SD',
    'SDIS': 'DISP',
    'SDIS_2': 'DISP',
    'SDIS_ERR': 'NULL',
    'SDIS_2_ERR': 'NULL',
    'XPOS': 'XPOS',
    'YPOS': 'YPOS',
    'VSYS': 'VSYS',
    'XPOS_2': 'XPOS',
    'YPOS_2': 'YPOS',
    'VSYS_2': 'VSYS',
    'PA': 'P.A.',
    'PA_2': 'P.A.',
    'INCL': 'INC',
    'INCL_2':'INC',
    'PA_ERR': 'NULL',
    'PA_2_ERR': 'NULL',
    'INCL_ERR': 'NULL',
    'INCL_2_ERR':'NULL',
    }
    Variables = np.array([e.upper() for e in Variables],dtype=str)

    try:
        tmp = open(filename, 'r')
    except FileNotFoundError:
        print("We cannot find that barollo file")
        outputarray = np.zeros((3, len(Variables)), dtype=float)
        outputarray[:,:] = float(0.)
        if unpack:
            return (*outputarray.T,)
        else:
            return outputarray
    unarranged = tmp.readlines()
    tmp.close()
    stripped_variables =  [f.split('_')[0] for f in Variables]
    if 'SBR' in stripped_variables:
        dir = '/'.join(filename.split('/')[:-1])
        print(dir)
        tmp = open(dir+'/densprof.txt', 'r')
        sbrfile = tmp.readlines()
        tmp.close()
        Var_inSBR = [f.strip() for f in sbrfile[13].split() if f != '#']

    numrings = len(unarranged)-1
    outputarray=np.zeros((numrings,len(Variables)),dtype=float)
    Var_inFile = [f.split('(')[0] for f in unarranged[0].split() if f != '']
    # Separate the keyword names
    for i,var in enumerate(Variables):
        bar_var = tir_bar_translate[var]
        if bar_var == 'ERR_SD' or bar_var == 'SURFDENS':
            for x in range(0,len(sbrfile)-16):
                tmp = [f for f in sbrfile[x+16].split() if f != '']
                outputarray[x,i] = float(tmp[Var_inSBR.index(bar_var)])
        elif bar_var == 'NULL':
            outputarray[:,i] = 0.
        else:
            for x in range(0,len(unarranged)-1):
                tmp = [f for f in unarranged[x+1].split() if f != '']
                outputarray[x,i] = float(tmp[Var_inFile.index(bar_var)])
    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray
# load the basic info file to get the Sofia FAT Initial_Estimates
def load_basicinfo(filename, Variables = ['RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask','Tot FLux','D_HI','Distance','HI_Mass' ,'D_HI' ], unpack = True):
    outputarray = np.zeros((2,len(Variables)), dtype=float)
    try:
        tmp = open(filename, 'r')
    except FileNotFoundError:
        print("That Basic info file does not exist. Returning empty array.")
        if unpack:
            return (*outputarray.T,)
        else:
            return outputarray
    fileline = tmp.readlines()
    Var_inFile = [f.strip() for f in fileline[3].split('  ') if f != '']
    invalues = [f.strip() for f in fileline[6].split('  ') if f != '']
    for i,var in enumerate(Variables):
        if var == 'RA' or var == 'DEC':
            tmp_str = invalues[Var_inFile.index('RA')].split('+/-')
            tmp_str2 = invalues[Var_inFile.index('DEC')].split('+/-')
            RA, DEC = convertRADEC(tmp_str[0],tmp_str2[0],invert=True)
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


# Function to plot fits files, Might not be working
def plot_fits(filetoplot, figure, contours=None, figure_coordinates=[0.1, 0.1, 0.8, 0.8], cmap='jet', log=False,
                  aspect=None, cbar = None):
    hdu = fits.open(filetoplot)[0]
    plt.rc('xtick', color='k', labelsize='medium', direction='in')
    plt.rc('ytick', color='k', labelsize='medium', direction='in')
    space = figure_coordinates[2] * 7.5
    wcsfound = False
    try:

        wcsin = wcs.WCS(hdu.header)
        fitsplt = figure.add_axes(figure_coordinates, projection=wcsin)
        wcsfound = True
        print("This happening now?")
        exit()
    except:
        #fitsplt = figure.add_axes(figure_coordinates)
        fitsplt= figure
        xaxis = [hdu.header['CRVAL1'] + (i - hdu.header['CRPIX1'] + 1) * (hdu.header['CDELT1']) for i in
                 range(hdu.header['NAXIS1'])]
        yaxis = [hdu.header['CRVAL2'] + (i - hdu.header['CRPIX2'] + 1) * (hdu.header['CDELT2']) for i in
                 range(hdu.header['NAXIS2'])]

    # we have to check for nan blanks and set them to 0
    tmp = np.isnan(hdu.data)
    nonzero = np.array(np.where(tmp == True))

    if nonzero.shape[1] < 1:
        nonzero = np.array(hdu.data[np.where(hdu.data != 0)])
    else:
        nonzero = np.array(hdu.data[np.where(tmp == False)])
    mini = np.min(nonzero)
    counter = 0
    while np.array(np.where(nonzero < mini)[0]).shape[0] < nonzero.shape[0] / 10.:
        mini = mini + abs(mini / 10.)
        counter += 1
        if counter > 100:
            break
    counter = 0
    maxi = np.max(nonzero)
    while np.array(np.where(nonzero > maxi)[0]).shape[0] < nonzero.shape[0] / 200.:
        maxi = maxi - abs(maxi / 50.)
        counter += 1
        if counter > 100:
            break

    if log:
        if mini < 0:
            print('Why are we not adjusting?')
            rows, cols = np.where(hdu.data > 1e-3)
            newarray = hdu.data[rows, cols]
            mini = np.min(newarray)
            counter = 0
            while np.array(np.where(newarray < mini)[0]).shape[0] < newarray.shape[0] / 10.:
                mini = mini + abs(mini / 10.)
                counter += 1
                if counter > 100:
                    break
        else:
            newarray = hdu.data
        print(mini)

        i = fitsplt.imshow(hdu.data, origin='lower', cmap=cmap, norm=colors.LogNorm(vmin=mini, vmax=maxi),
                           aspect=aspect)

    else:
        i = fitsplt.imshow(hdu.data, origin='lower', cmap=cmap, vmin=mini, vmax=maxi, aspect=aspect)
    if cbar.lower() == 'horizontal':
        divider = make_axes_locatable(fitsplt)
        cax = divider.append_axes("top", size="5%", pad=0.00, axes_class=maxes.Axes)
        cbar = plt.colorbar(i, cax=cax, orientation='horizontal')
        cax.xaxis.set_ticks_position('top')
        cbar.set_ticks([mini, maxi])
    #try:
    #    cbar = figure.colorbar(i)
    #except:
    #    pass

    if wcsfound:
        ra = fitsplt.coords[0]
        ra.set_major_formatter('hh:mm:ss')
        ra.set_ticks(number=round(space))
    else:
        plt.gca().set_xticks(range(len(xaxis))[0:-1:int(len(xaxis) / 5)])
        plt.gca().set_yticks(range(len(yaxis))[0:-1:int(len(yaxis) / 5)])
        plt.gca().set_xticklabels(['{:10.0f}'.format(i) for i in xaxis[0:-1:int(len(xaxis) / 5)]])
        plt.gca().set_yticklabels(['{:10.1f}'.format(i) for i in yaxis[0:-1:int(len(yaxis) / 5)]])
    fitsplt.set_xlabel(hdu.header['CTYPE1'])
    fitsplt.set_ylabel(hdu.header['CTYPE2'])
    if log:
        levels = [1E20, 2E20, 4E20, 8E20, 16E20, 32E20]
        try:
            beam = [hdu.header['BMAJ'] * 3600., hdu.header['BMIN'] * 3600.]
        except:
            beam = [1., 1.]
        levels = columndensity(levels, vwidth=1, beam=beam, ncolumn=True) / 1000.
    else:
        if mini > 0:
            maxvrot = (maxi - mini) / 2.
            velostep = int(((int(maxi - mini) - int(maxi - mini) / 10.) / 10.) * 1.5)
            levels = [(i) * velostep + mini for i in range(15)]

        else:

            sig1 = np.std(hdu.data[0:10, -10:-1])
            sig2 = np.std(hdu.data[0:10, 0:10])
            sig3 = np.std(hdu.data[-10:-1, -10:-1])
            sig4 = np.std(hdu.data[-10:-1, 0:10])
            sig = np.mean([sig1, sig2, sig3, sig4])
            if np.isnan(sig):
                sig = abs(mini) / 3.
            while sig > maxi:
                sig = sig / 2.
            levels = np.array([-2, -1, 1, 2, 4, 8, 16, 32, 64, 128]) * 1.5 * sig

    contourdata = fitsplt.contour(hdu.data, levels, colors='k', origin='lower', linewidths=0.75)
    if contours:

        hdumod = fits.open(contours)[0]
        if mini < 0:
            contour = fitsplt.contour(hdumod.data, levels, colors='r', origin='lower')
        else:
            try:
                contour = fitsplt.contour(hdumod.data, levels, colors='w', origin='lower')
            except UserWarning:
                contour = fitsplt.contour(hdumod.data, levels / 2., colors='w', origin='lower')

#Function to read simple input files that  use = as a separator between the required input and the values
def read_input_file(filename):
    tmpfile = open(filename, 'r')
    File = Proper_Dictionary({})
    unarranged = tmpfile.readlines()
    # Separate the keyword names
    for tmp in unarranged:
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
        File[tmp.split('=', 1)[0].strip().upper()] = tmp.rstrip()
    return File

# A function to regrid a fits image to different pixels size

def regrid_image(img,ratio):
    # First get the shape of the data
    shape = np.array(img[0].data.shape, dtype=float)
    # As this only work when shape is an exact multiple of ratio we will make a zero-paaded array that satisfies that condition.
    new_shape = ratio * np.ceil(shape / ratio).astype(int)

    # Create the zero-padded array and assign it with the old density
    new_img = np.zeros(new_shape)
    shape = np.array(img[0].data.shape, dtype=int)
    new_img[:shape[0], :shape[1]] = img[0].data

    # Now reshape and regrid to new image
    tmp = new_img.reshape((new_shape[0] // ratio, ratio,
                                   new_shape[1] // ratio, ratio))
    coarse_img = np.sum(tmp, axis=(1, 3))/(ratio**2)
    img[0].data = coarse_img
    img[0].header['CRPIX1'] =   np.ceil(img[0].header['CRPIX1']/ratio)
    img[0].header['CRPIX2'] =   np.ceil(img[0].header['CRPIX2']/ratio)
    try:
        img[0].header['CDELT1'] =   img[0].header['CDELT1']*ratio
        img[0].header['CDELT2'] =   img[0].header['CDELT2']*ratio
    except:
        print("No CDELT found")
    try:
        img[0].header['CD1_1'] =   img[0].header['CD1_1']*ratio
        img[0].header['CD2_2'] =   img[0].header['CD2_2']*ratio
    except:
        print("No CD corr matrix found")
    try:
        img[0].header['CD1_2'] =   img[0].header['CD1_2']*ratio
        img[0].header['CD2_1'] =   img[0].header['CD2_1']*ratio
    except:
        print("No CD cross-corr matrix found")

    img[0].header['NAXIS1'] =   new_shape[1]
    img[0].header['NAXIS2'] =   new_shape[0]
    return img
#function to rotate a cube without losing info
def rotateCube(Cube, angle, pivot):
    padX = [int(Cube.shape[2] - pivot[0]), int(pivot[0])]
    padY = [int(Cube.shape[1] - pivot[1]), int(pivot[1])]
    imgP = np.pad(Cube, [[0, 0], padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, angle, axes=(2, 1), reshape=False)
    return imgR[:, padY[0]: -padY[1], padX[0]: -padX[1]]

#function to align and scale an image based on the histogram of values
#where the image is treated such that the peak occurs at the peak value and the bottom corresponds to -3 sigma
#if bootom is unset the image is only shifted to match the peak value
def scale_and_align(img,peak=0.3,bottom='Empty'):
    # first we want to determine a range that includes the bulk of values but not extreme outliers
    std_dev=np.nanstd(img)
    bin_counts, bin_edges = np.histogram(img, bins = 10000, range=[-5*std_dev,10*std_dev])
    bin_centers = bin_edges[:-1]+abs(bin_edges[5] - bin_edges[6])/2.
    #locate our maximum
    max_loc =np.where( bin_counts == np.nanmax(bin_counts))[0][0]
    if max_loc > len(bin_counts)-100:
        positive_max_buffer = len(bin_counts)-max_loc-1
    else:
        positive_max_buffer = 100
    # we fit a GAUSSIAN to the negative sides
    gauss_parameters = fit_gaussian(bin_centers[0:max_loc+positive_max_buffer], bin_counts[0:max_loc+positive_max_buffer])
    #shift the image
    new_img = img -(gauss_parameters[1]-peak)
    #if we have a bottom we want to scale as swell
    if bottom != 'Empty':
        new_img = (new_img - peak)*abs((peak-bottom)/(gauss_parameters[1]-3*gauss_parameters[2]))+peak
    return new_img



#Function to convert pixels to coordinates in a correct manner (Ripped from SoFiA's wcs_coordinates)
#https://github.com/SoFiA-Admin/SoFiA/blob/master/sofia/wcs_coordinates.py
def sofia_pix2world(pixels,header,unpack =True):
    # Fix headers where "per second" is written "/S" instead of "/s"
    # (assuming they mean "per second" and not "per Siemens").
    if "cunit3" in header and "/S" in header["cunit3"]:
        err.warning("Converting '/S' to '/s' in CUNIT3.")
        header["cunit3"] = header["cunit3"].replace("/S", "/s")

    # Constrain the RA axis reference value CRVAL_ to be between 0 and 360 deg
    rafound = 0
    for kk in range(header["naxis"]):
        if header["ctype1"][:2] == "RA":
            rafound = 1
            break
    if rafound:
        if header["crval%i" % (kk + 1)] < 0:
            err.warning("Adding 360 deg to RA reference value.")
            header["crval%i" % (kk + 1)] += 360
        elif header["crval%i" % (kk + 1)] > 360:
            err.warning("Subtracting 360 deg from RA reference value.")
            header["crval%i" % (kk + 1)] -= 360

    if header['naxis'] == 2:
        wcsin = wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL])
        xy = np.array(pixels, dtype=float)
        objects = wcsin.wcs_pix2world([xy], 0)
    else:
        wcsin = wcs.WCS(header, naxis=[wcs.WCSSUB_CELESTIAL, wcs.WCSSUB_SPECTRAL])
        xyz = np.array(pixels, dtype=float)
        if "cellscal" in header and header["cellscal"] == "1/F":
            print(
                "CELLSCAL keyword with value of 1/F found.\n"
                "Will account for varying pixel scale in WCS coordinate calculation.")
            x0, y0 = header["crpix1"] - 1, header["crpix2"] - 1
            # Will calculate the pixscale factor of each channel as:
            # pixscale = ref_frequency / frequency
            if header["ctype3"] == "VELO-HEL":
                pixscale = (1 - header["crval3"] / scipy.constants.c) / (1 - (
                            ((xyz[2] + 1) - header["crpix3"]) * header["cdelt3"] + header[
                        "crval3"]) / scipy.constants.c)
            else:
                print(
                    "Cannot convert 3rd axis coordinates to frequency. Ignoring the effect of CELLSCAL = 1/F.")
                pixscale = 1.0
            xyz[0] = (xyz[0] - x0) * pixscale + x0
            xyz[1] = (xyz[1] - y0) * pixscale + y0
            print(xyz)
            objects =  wcsin.wcs_pix2world([xyz], 0)
    objects=np.array(objects,dtype=float)
    if unpack:
        return (*objects.T,)
    else:
        return objects
