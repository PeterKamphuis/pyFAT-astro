#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
# This is the python version of FAT
import sys
import os
import copy
import numpy as np
from optparse import OptionParser
import traceback
import warnings
from datetime import datetime
from astropy.io import fits
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Ellipse
    import matplotlib.axes as maxes
from astropy.io import fits
from astropy.wcs import WCS
import pyFAT_astro

import pyFAT_astro.Support.modify_template as mt
import pyFAT_astro.Support.read_functions as rf
import pyFAT_astro.Support.run_functions as runf
import pyFAT_astro.Support.write_functions as wf
import pyFAT_astro.Support.support_functions as sf
import pyFAT_astro.Support.fits_functions as ff
#import pyFAT_astro.Support.development_functions as df

homedir = '/home/peter/'
#homedir = '/Users/peter/'

Configuration = {'INNER_FIX': 4}
#Configuration['FITTING_DIR']=f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/Mass2.5e+12-i20d14.0-7.5pa35.0w0.0-0.0-No_Flare-ba12SNR8bm20.0-20.0ch4.0-Arms-No_Bar-rm0.0/"
#Configuration['FITTING_DIR']=f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/Mass5.0e+10-i42.0d15.0-12.0pa115w0.1-0.07-Flared-ba10SNR8.0bm10.0-10.0ch4.0-No_Arms-No_Bar-rm0.0/"
#Configuration['FITTING_DIR']= f"{homedir}/FAT_Main/Test_Sets/From_Bochum/Mass2.5e+11-i48.0d13.0-7.5pa115w0.07-0.15-No_Flare-ba15SNR1bm10.0-10.0ch4.0-No_Arms-No_Bar-rm0.0/"
#Configuration['FITTING_DIR']= f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/M_83_6.0Beams_3.0SNR/"
#Configuration['FITTING_DIR']= f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/Mass2.5e+12-i15d14.0-7.5pa35.0w0.0-0.0-No_Flare-ba12SNR8bm20.0-20.0ch4.0-Arms-No_Bar-rm0.0/"
#Configuration['FITTING_DIR']= f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/NGC_3198_36.9Beams_1.0SNR/"
#Configuration['FITTING_DIR']= f"{homedir}/FAT_Main/FAT_Testers/LVHIS-26_3/HPASS00018/"
Configuration['FITTING_DIR']='/home/peter/FAT_Main/FAT_Testers/Database-09-10-2020/M_83_6.0Beams_3.0SNR/'
#Configuration['FITTING_DIR']=f"{homedir}FAT_Main/Database/M_83_6.0Beams_3.0SNR/"
Fits_Files = {'ORIGINAL_CUBE': "Convolved_Cube.fits"}
Configuration['START_DIRECTORY']= f'{os.getcwd()}'
Configuration['CUBENAME']= 'Convolved_Cube'
Configuration['BASE_NAME']= 'Convolved_Cube_FAT'
Configuration['LOG_DIRECTORY']= f"{Configuration['FITTING_DIR']}Logs/"
Configuration['TIRIFIC'] = 'tirific'
Configuration['SHAKER_ITERATIONS'] = 2
#Fits_Files = {'ORIGINAL_CUBE': "Cube.fits"}
#Configuration['CUBENAME']= 'Cube'
#Configuration['BASE_NAME']= 'Cube_FAT'
Configuration['SUB_DIR']= 'Mass2.5e+12-i20d14.0-7.5pa35.0w0.0-0.0-No_Flare-ba12SNR8bm20.0-20.0ch4.0-Arms-No_Bar-rm0.0/'

#Configuration['CUBENAME']= 'Convolved_Cube'
#Configuration['BASE_NAME']= 'Convolved_Cube_FAT'
Fits_Files['NOISEMAP'] = f"{Configuration['BASE_NAME']}_noisemap.fits"
Fits_Files['FITTING_CUBE'] = f"{Configuration['CUBENAME']}_FAT.fits"
Fits_Files['OPTIMIZED_CUBE'] = f"{Configuration['CUBENAME']}_FAT_opt.fits"
Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"
Configuration['OUTPUTLOG'] = None

other_keys =  {'MINIMUM_WARP_SIZE': 3., # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc,set here.
               'MINIMUM_RINGS': 3,  # we need at least this amount of rings (Including 0 and 1/5 beam), set here
               'TOO_SMALL_GALAXY': 1., # if the number of beams across the major axis/2 is less than this we will not fit the galaxy, set here
               'DEBUG': True,
               'DISTANCE': 'Unset', # Distance to the galaxy, set from the catalogue at start of loop
               'ID_NR': 'Unset', # ID of the galaxy in the catalogue , set from the catalogue at start of loop
               'SOFIA_BASENAME': 'Unset', #Basename of pre-processed sofia products, only set when provided in catalogue at start of loop
               'BASENAME': 'Unset', #Basename for FAT products, typically {input_cube}_FAT, set at start of loop
               'LOG_DIR': 'Unset', #Directory to put log files from run, set at start of loop
               'USED_FITTING':'Fit_Tirific_OSC',
               'CURRENT_STAGE': 'initial', #Current stage of the fitting process, set at switiching stages
               'TIRIFIC_PID': 'Not Initialized', #Process ID of tirific that is running
               'ACCEPTED': False, #Whether a fit is accepted or not
               'FIXED_PARAMETERS': [['Z0','XPOS','YPOS','VSYS']],
               'MAX_SIZE_IN_BEAMS': 30, # The galaxy is not allowed to extend beyond this number of beams in radius, set in check_source
               'MIN_SIZE_IN_BEAMS': 0., # Minimum allowed radius in number of beams of the galaxy, set in check_source
               'SIZE_IN_BEAMS': 0, # The radius of the galaxy in number of beams, adapted after running Sofia
               'NO_RINGS': 0., # The number of rings in the fit,
               'LAST_RELIABLE_RINGS': [0.,0.], # Location of the rings where the SBR drops below the cutoff limits, adapted after every run. Should only be modified in check_size
               'LIMIT_MODIFIER': [1.], #Modifier for the cutoff limits based on the inclination , adapted after every run.
               'OLD_RINGS': [1.], # List to keep track of the ring sizes that have been fitted.
               'EXCLUDE_CENTRAL': False,
               'NO_POINTSOURCES': 0. , # Number of point sources, set in run_tirific

               'INNER_FIX': [4,4], #Number of rings that are fixed in the inner part for the INCL and PA, , adapted after every run in get_inner_fix in support_functions
               'WARP_SLOPE': [0.,0.], #Ring numbers from which outwards the warping should be fitted as a slope, set in get_warp_slope in modify_template
               'OUTER_SLOPE_START': 1, # Ring number from where the RC is fitted as a slope
               'RC_UNRELIABLE': 1, # Ring number from where the RC values are set flat. Should only be set in check_size

               'NOISE': 0. , #Noise of the input cube, set in main
               'BEAM_IN_PIXELS': [0.,0.,0.], #FWHM BMAJ, BMIN in pixels and total number of pixels in beam area, set in main
               'BEAM': [0.,0.,0.], #  FWHM BMAJ, BMIN in arcsec and BPA, set in main
               'NAXES': [0.,0.,0.], #  Size of the cube in pixels x,y,z arranged like sane people not python, set in main
               'NAXES_LIMITS': [[0.,0.],[0.,0.],[0.,0.]], #  Size of the cube in degree and km/s,  x,y,z arranged like sane people not python, set in main updated in cut_cubes
               'MAX_ERROR': {'SDIS': 4.,\
                             'PA' : 15.,\
                             'INCL': 15.}, #The maximum allowed errors for the parameters, set in main derived from cube
               'CHANNEL_WIDTH': 0., #Width of the channel in the cube in km/s, set in main derived from cube
               'PIXEL_SIZE': 0., #'Size of the pixels in degree'
               }

for key in other_keys:
    Configuration[key] = other_keys[key]

boundary_limit_keys = ['PA','INCL', 'SDIS', 'Z0','VSYS','XPOS','YPOS']
for key in boundary_limit_keys:
    Configuration[f"{key}_CURRENT_BOUNDARY"] = [[0.,0.],[0.,0.],[0.,0.]]

Configuration['RING_SIZE'] = 1.
Configuration['OPTIMIZED'] = False
Configuration['NCPU'] = 1.
Configuration['LOOPS'] = 3
Configuration['HANNING'] = False
Configuration['FIX_INCLINATION'] = [False,False]
Configuration['FIX_SDIS'] = [False,False]
Configuration['FIX_PA'] = [False,False]
Configuration['FIX_Z0'] = [True,True]
Configuration['FIX_SBR'] = [True,True]
Configuration['OUTER_RINGS_DOUBLED'] = False
Configuration['DISTANCE'] = 5.
Configuration['TIRIFIC_RUNNING'] = False
Configuration['TIMING'] = False
#Configuration['PA_CURRENT_BOUNDARY'] = [120.,280.]
cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}")

beamarea=(np.pi*abs(cube_hdr['BMAJ']*cube_hdr['BMIN']))/(4.*np.log(2.))
Configuration['CHANNEL_WIDTH'] = cube_hdr['CDELT3']/1000.
Configuration['BEAM_IN_PIXELS'] = [cube_hdr['BMAJ']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                   cube_hdr['BMIN']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                   beamarea/(abs(cube_hdr['CDELT1'])*abs(cube_hdr['CDELT2']))]
Configuration['BEAM'] = [float(cube_hdr['BMAJ']*3600.), float(cube_hdr['BMIN']*3600.),float(cube_hdr['BPA'])]

Configuration['NOISE'] = cube_hdr['FATNOISE']
Configuration['PIXEL_SIZE'] = np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])])
Configuration['NAXES'] = [cube_hdr['NAXIS1'],cube_hdr['NAXIS2'], cube_hdr['NAXIS3']]
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    coordinate_frame = WCS(cube_hdr)
    xlow,ylow,zlow = coordinate_frame.wcs_pix2world(1,1,1., 1.)
    xhigh,yhigh,zhigh = coordinate_frame.wcs_pix2world(*Configuration['NAXES'], 1.)
    Configuration['NAXES_LIMITS'] = [np.sort([xlow,xhigh]),np.sort([ylow,yhigh]),np.sort([zlow,zhigh])/1000.]


def Test_Regularise():
    #Variables = ['VROT','VROT_2','PA', 'PA_2','INCL','INCL_2','SDIS','SDIS_2','SBR','SBR_2','VSYS']
    Variables = ['SDIS']
    Tirific_Template = rf.tirific_template(f"{Configuration['FITTING_DIR']}/One_Step_Convergence/One_Step_Convergence.def")
    Configuration['NO_RINGS'] = int(Tirific_Template['NUR'])
    Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']-2
    Configuration['LIMIT_MODIFIER']= [3.469106943427586/2.5]
    radii = np.array(Tirific_Template['RADI'].split(),dtype=float)
    print(f"Size before checking {Configuration['SIZE_IN_BEAMS']}")
    accepted_size = mt.check_size(Configuration,Tirific_Template, fit_type = 'One_Step_Convergence', stage = 'after_os', debug=True,Fits_Files=Fits_Files)
    print(f"We are  modifying size {accepted_size}, {Configuration['SIZE_IN_BEAMS']}")

    if Configuration['NO_RINGS'] > 5:
        inner_slope = 0.
    else:
        inner_slope = NUR
    if Configuration['NO_RINGS'] > 15 and inner_slope > int(Configuration['NO_RINGS']*4./5.) :
        inner_slope = int(Configuration['NO_RINGS']*4./5.)
    Configuration['OUTER_SLOPE'] = inner_slope
    Variables_red = ['PA']
    for i,key in enumerate(Variables_red):
        sm_profile = mt.smooth_profile(Configuration,Tirific_Template, key,debug=True,no_apply=True)
        reg_prof = mt.regularise_profile(Configuration,Tirific_Template, key,min_error= 1.,debug = True, no_apply =True)
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii,np.array(Tirific_Template[key].split(),dtype=float),c='b',label='Original')
        plt.plot(radii,reg_prof[0],c='r',label='Regularised')
        plt.plot(radii,sm_profile[0],c='k',label='Smoothed')
        plt.legend()
        plt.savefig(f'Tests/{key}.png')
        plt.close()
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii,np.array(Tirific_Template[f'{key}_2'].split(),dtype=float),c='b',label='Original')
        plt.plot(radii,reg_prof[1],c='r',label='Regularised')
        plt.plot(radii,sm_profile[1],c='k',label='Smoothed')
        plt.legend()
        plt.savefig(f'Tests/{key}_2.png')
        plt.close()

def Test_Overview():
    wf.make_overview_plot(Configuration,Fits_Files, debug = True)

def Test_write_temp():
    Tirific_Template = rf.tirific_template(f'{homedir}/FAT_Main/FAT_Source_Dev/Support/template.def')
    key = 'VROT'
    Configuration = {'INNER_FIX': 4}
    Configuration['NO_RINGS'] = 10
    error = np.full([15,15],5.)
    print(error)
    Configuration['FITTING_DIR']= f'{homedir}/FAT_Main/FAT_Source_Dev/Installation_Check/'
    Tirific_Template.insert(key,f"# {key}_ERR",f"{' '.join([f'{x}' for x in error[0,:int(Configuration['NO_RINGS'])]])}")
    Tirific_Template.insert(f"{key}_2",f"# {key}_2_ERR",f"{' '.join([f'{x}' for x in error[1,:int(Configuration['NO_RINGS'])]])}")
    wf.tirific(Configuration,Tirific_Template, name = 'test.def', debug = False)

def Test_psutil_stats():
    import psutil as psu
    import subprocess
    import time
    startdir= os.getcwd()
    work_dir = '/home/peter/FAT_Main/FAT_Testers/Database-09-10-2020/M_83_6.0Beams_3.0SNR'

    deffile = 'Test.def'
    current_run = subprocess.Popen(['tirific',f"DEFFILE={deffile}","ACTION= 1"],\
                                   stdout = subprocess.PIPE, stderr = subprocess.PIPE,\
                                   cwd=work_dir,universal_newlines = True)

    Configuration['TIMING']=True
    time.sleep(0.1)
    current_proc = psu.Process(current_run.pid)
    for tir_out_line in current_run.stdout:
        if Configuration['TIMING']:

            cpu_duration,CPU,mem= get_usage_stats(Configuration,current_proc,debug=True)
            print(f'Attempting psutil statistics')
            print(f'clock_time = {cpu_duration} CPU%= {CPU} memory = {mem}')
            CPU,mem= sf.get_usage_statistics(Configuration,current_run.pid,debug=False)
            print(f'Attempting old statistics')
            print(CPU,mem)
            print(f'!!!!!!!!!!!yeah!!!!!!!'')
        time.sleep(1)

def get_usage_stats(Configuration,process, debug = False):

    Cpuduration=(process.cpu_times()[0]+process.cpu_times()[1])/60.
    memory_in_mb = (process.memory_info()[0]+process.memory_info()[1])/2**20. #psutilreturns bytes
    cpu_percent = process.cpu_percent(interval=1)
    return Cpuduration,cpu_percent,memory_in_mb

def Test_Ram():
    Configuration = {'FITTING_DIR': '/home/peter/FAT_Main/FAT_Source_Dev/Installation_Check/'}
    Configuration['OUTPUTLOG'] = None
    wf.plot_usage_stats(Configuration,debug = True)
def Test_One_Orientation():
    pa, inclination, SBR_initial, maj_extent,x_new,y_new,VROT_initial = rf.guess_orientation(Configuration,Fits_Files,debug=True)
    print(f"PA = {pa}, inclination = {inclination}, majextent = {maj_extent}")

def Test_Orientation():
    #Then read the input Catalogue
    Full_Catalogue = rf.catalogue( f'{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/Output_Summary.txt')
    with open(f"Inclinations.txt",'w') as file:
        file.write("#Guessing inclinations. \n")
    for current_galaxy_index in range(0,len(Full_Catalogue['DIRECTORYNAME'])):
        Configuration = {'RING_SIZE': 1.}
        Configuration['MAINDIR'] =  '/home/peter/FAT_Main/FAT_Testers/Database-09-10-2020/'
        Configuration['START_TIME'] = datetime.now()
        # First check the starttime
        Catalogue = {}
        for key in Full_Catalogue:
            if key != 'ENTRIES':
                Catalogue[key] = Full_Catalogue[key][current_galaxy_index]
            else:
                Catalogue[key] = Full_Catalogue[key]
        if 'BASENAME' in Catalogue['ENTRIES']:
            Configuration['BASE_NAME'] = Catalogue['BASENAME']
        else:
            Configuration['BASE_NAME'] = Catalogue['CUBENAME']+'_FAT'
        Fits_Files = {'ORIGINAL_CUBE': f"{Catalogue['CUBENAME']}.fits"}
        Fits_Files['FITTING_CUBE'] = f"{Catalogue['CUBENAME']}_FAT.fits"
        Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
        Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
        Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"
        if Catalogue['DIRECTORYNAME'] == './':
            Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/"
        else:
            Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/{Catalogue['DIRECTORYNAME']}/"
        if Configuration['FITTING_DIR'][-2:] == '//':
            Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'
        hdr = fits.getheader(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}")
        Configuration['NOISE'] = hdr['FATNOISE']
        Configuration['CHANNEL_WIDTH'] =  hdr['CDELT3']/1000.
        Configuration['OUTPUTLOG'] = None
        if os.path.exists(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['MOMENT0']}") \
            and os.path.exists(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['CHANNEL_MAP']}") \
            and os.path.exists(f"{Configuration['FITTING_DIR']}Sofia_Output/{Fits_Files['MOMENT1']}") :
            pa, inclination, sbr,majextent = rf.guess_orientation(Configuration, Fits_Files,debug =True)
        Vars_to_plot= ['INCL']
        if os.path.exists(f"{Configuration['FITTING_DIR']}ModelInput.def"):
            Input_Model = rf.load_tirific(f"{Configuration['FITTING_DIR']}ModelInput.def",Variables= Vars_to_plot,unpack=False,debug=True)
        else:
            Input_Model = []
        incl_in = Input_Model[0,0]
        with open(f"Inclinations.txt",'a') as file:
            file.write(f'For {incl_in} we guess {inclination[0]}+/-{inclination[1]} in the galaxy {Catalogue["DIRECTORYNAME"]} \n')

def Test_Parameterized():
    Configuration = {'INNER_FIX': 4}
    Configuration['FITTING_DIR']= f'{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/NGC_3198_6.0Beams_1.0SNR/'
    Configuration['OUTPUTLOG'] = None
    Configuration['RING_SIZE'] = 1.
    Configuration['NOISE'] = 0.003
    Configuration['MINIMUM_RINGS'] = 3
    Configuration['SIZE_IN_BEAMS'] = 3.5
    Configuration['CUBENAME']= 'Convolved_Cube'
    Configuration['BASE_NAME']= 'Convolved_Cube_FAT'

    Fits_Files = {'ORIGINAL_CUBE': "Convolved_Cube.fits"}
    Fits_Files['NOISEMAP'] = f"{Configuration['BASE_NAME']}_noisemap.fits"
    Fits_Files['FITTING_CUBE'] = f"{Configuration['CUBENAME']}_FAT.fits"
    Fits_Files['OPTIMIZED_CUBE'] = f"{Configuration['CUBENAME']}_FAT_opt.fits"
    Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
    Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
    Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
    Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
    Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"
    Configuration['FIX_INCLINATION'] = False
    Configuration['FIX_SDIS'] = False
    Configuration['FIX_PA'] = False
    Configuration['FIX_Z0'] = True
    hdr = fits.getheader(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}")
    Configuration['DISTANCE'] = 5.
    Configuration['TIRIFIC_RUNNING'] = False
    Configuration['TIMING'] = False
    # Initialize a def file.
    Tirific_Template = rf.tirific_template(f"{Configuration['FITTING_DIR']}/Centre_Convergence_In.def")
    Configuration['NO_RINGS'] = Tirific_Template['NUR']
    current_run = ''
    #runf.fit_warp(Configuration,Tirific_Template,current_run, fit_stage = 'None')
    accepted,current_run = runf.run_parameterized_tirific(Configuration,Tirific_Template,current_run,stage = 'paramiterized', fit_stage = 'Parameterized', debug= True)


def Test_Extract_PV():
    PA =53.53
    Center = [187.218,4.293,4221*1000.]
    cube = fits.open(f"{homedir}/Misc/Aditya_Stuff/PV_Scale/1440_data_nostokes.fits")
    PV = ff.extract_pv(cube,PA, \
                        center = Center, \
                        convert=1000.)

    fits.writeto(f"{homedir}/Misc/Aditya_Stuff/PV_Scale/test.fits",PV[0].data,PV[0].header,overwrite = True)

def Test_Set_Rings():
    Configuration = {'RING_SIZE': 1 ,
                     'MAX_SIZE_IN_BEAMS': 30, # The galaxy is not allowed to extend beyond this number of beams in radius, set in check_source
                     'MIN_SIZE_IN_BEAMS': 0., # Minimum allowed radius in number of beams of the galaxy, set in check_source
                     'SIZE_IN_BEAMS': 11.1975694623135514, # The radius of the galaxy in number of beams, adapted after running Sofia
                     'NO_RINGS': 0., # The number of rings in the fit
                     'OUTPUTLOG': None,
                     'BEAM': [20.,20.,0.]
                     }

    rad= sf.set_rings(Configuration,debug=True)


def Test_gflux_dH():
    Configuration = {'FITTING_DIR':f"{homedir}/FAT_Main/FAT_Testers/Database-09-10-2020/M_83_2.0Beams_1.0SNR/"}
    cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}Convolved_Cube_FAT.fits")
    beamarea=(np.pi*abs(cube_hdr['BMAJ']*cube_hdr['BMIN']))/(4.*np.log(2.))
    Configuration['CHANNEL_WIDTH'] = cube_hdr['CDELT3']/1000.
    Configuration['BEAM_IN_PIXELS'] = [cube_hdr['BMAJ']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                       cube_hdr['BMIN']/np.mean([abs(cube_hdr['CDELT1']),abs(cube_hdr['CDELT2'])]),\
                                       beamarea/(abs(cube_hdr['CDELT1'])*abs(cube_hdr['CDELT2']))]
    Configuration['NOISE'] = cube_hdr['FATNOISE']

    Totflux = rf.get_totflux(Configuration,f"/Finalmodel/Finalmodel_mom0.fits", debug=True)
    DHI = rf.get_DHI(Configuration,Model='One_Step_Convergence',debug=True)

    print(DHI,Totflux)

def Test_Tilts():
    Tirific_Template = rf.tirific_template(filename = f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def", debug = True)
    Configuration['NO_RINGS'] = int(Tirific_Template['NUR'])
    Configuration['INNER_FIX'] = [1.,1.]
    Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']-2
    tiltogram = sf.make_tiltogram(Configuration,Tirific_Template,debug=True)
    for i in [0,1]:
        print("Yes I like taco's")
        print(tiltogram[i])
        header = fits.PrimaryHDU(tiltogram[i])
        header.writeto(f'Tests/tilt_{i}.fits',overwrite =True)

    sf.check_tiltogram(Configuration,tiltogram)
    print(Configuration['INNER_FIX'])
def TestTir_Shaker():
    runf.tirshaker_call(Configuration,debug=Configuration['DEBUG'])

def Test_Check_Inclination():
    Tirific_Template = rf.tirific_template(filename = f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def", debug = True)
    print("These are the final original inclinations")
    print(Tirific_Template['INCL'],Tirific_Template['INCL_2'])
    print(Configuration['FITTING_DIR'])
    df.check_inclination(Configuration,Tirific_Template,Fits_Files, fit_type = 'One_Step_Convergence', debug =True)
def basic():
    print(f" fstring use double or triple quotes")
    print(f" Dictionary entries use single quotes")
    print(f"Below the basic doc structure")
basic.__doc__ =f'''
 NAME:

 PURPOSE:

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



if __name__ == '__main__':
    Test_psutil_stats()
