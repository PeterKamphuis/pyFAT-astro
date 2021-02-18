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

import pyFAT

import pyFAT.Support.modify_template as mt
import pyFAT.Support.read_functions as rf
import pyFAT.Support.run_functions as runf
import pyFAT.Support.write_functions as wf
import pyFAT.Support.support_functions as sf
import pyFAT.Support.fits_functions as ff

homedir = '/home/peter/'

def Test_Regularise():
    Variables = ['VROT','VROT_2','PA', 'PA_2','INCL','INCL_2','SDIS','SDIS_2','SBR','SBR_2','VSYS']
    dir = homedir+'/FAT_Main/FAT_Testers/Full_Database_pyFAT/'
    galaxy = 'Mass2.5e+12-i20d14.0-7.5pa35.0w0.0-0.0-No_Flare-ba12SNR8bm20.0-20.0ch4.0-Arms-No_Bar-rm0.0'
    profile = rf.load_tirific(dir+galaxy+'/Extent_Convergence/Extent_Convergence_Prev.def',Variables)
    radii = rf.load_tirific(dir+galaxy+'/Extent_Convergence/Extent_Convergence_Prev.def', ['RADI'])
    Tirific_Template = {}
    hdr = fits.getheader('Installation_Check/NGC_2903_FAT.fits')
    for i,key in enumerate(Variables):
        if key == 'SBR':
            format = '.2e'
        else:
            format = '.2f'
        Tirific_Template[key] = f"{' '.join([f'{x:{format}}' for x in profile[i]])}"
        if key == 'INCL':
            Inclination = np.array([(x+y)/2. for x,y in zip(profile[i],profile[i+1])],dtype=float)
    Tirific_Template['RADI'] = f"{' '.join([f'{x}' for x in radii[0]])}"
    print(Tirific_Template)
    Configuration = {'INNER_FIX': 4}
    Configuration['RING_SIZE'] = 1.
    Configuration['NO_RINGS'] = len(profile[0])
    Configuration['INNER_FIX']= sf.get_inner_fix(Configuration, Tirific_Template,debug=False)

    Configuration['BMMAJ'] = (radii[0][3]-radii[0][2])/Configuration['RING_SIZE']
    Configuration['OUTPUTLOG'] = None
    Configuration['OUTER_RINGS_DOUBLED'] = False
    Configuration['MINIMUM_RINGS'] = 3
    Configuration['LIMIT_MODIFIER'] = 1.
    Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']*Configuration['RING_SIZE']-2
    sf.set_limit_modifier(Configuration,Inclination, debug= False)
    NUR =  Configuration['NO_RINGS']
    if Configuration['NO_RINGS'] > 5:
        inner_slope =  int(round(sf.set_limits(mt.check_size(Configuration,Tirific_Template,hdr, debug = False,fix_rc = True),round(NUR/2.),NUR-1)))
    else:
        inner_slope = NUR
    if Configuration['NO_RINGS'] > 15 and inner_slope > int(Configuration['NO_RINGS']*4./5.) :
        inner_slope = int(Configuration['NO_RINGS']*4./5.)
    Configuration['OUTER_SLOPE'] = inner_slope
    Variables_red = ['VROT','PA','INCL','SDIS']
    for i,key in enumerate(Variables_red):
        sm_profile = mt.smooth_profile(Configuration,Tirific_Template, key,hdr ,debug=True,no_apply=True)
        reg_prof = mt.regularise_profile(Configuration,Tirific_Template, key ,hdr,min_error= 1.,debug = True, no_apply =True)
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii[0],profile[Variables.index(key)],c='b',label='Original')
        if key == 'VROT':
            print(reg_prof[0],sm_profile[0],profile[Variables.index(key)])
        plt.plot(radii[0],reg_prof[0],c='r',label='Regularised')
        plt.plot(radii[0],sm_profile[0],c='k',label='Smoothed')
        plt.legend()
        plt.savefig(f'Tests/{key}.png')
        plt.close()
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii[0],profile[Variables.index(f'{key}_2')],c='b',label='Original')
        plt.plot(radii[0],reg_prof[1],c='r',label='Regularised')
        plt.plot(radii[0],sm_profile[1],c='k',label='Smoothed')
        plt.legend()
        plt.savefig(f'Tests/{key}_2.png')
        plt.close()

def Test_Overview():
    Configuration = {'INNER_FIX': 4}
    Configuration['NO_RINGS'] = 10
    Configuration['OUTPUTLOG'] = None
    Configuration['RING_SIZE'] = 1.
    Configuration['NOISE'] = 0.003
    Configuration['MINIMUM_RINGS'] = 3
    Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']*Configuration['RING_SIZE']-2
    Configuration['CUBENAME']= 'NGC_2903'
    Configuration['BASE_NAME']= 'NGC_2903_FAT'
    Configuration['FITTING_DIR']= homedir+'/FAT_Main/FAT_Source_Dev/Installation_Check/'#Make a dictionary for the fitsfiles we use
    Fits_Files = {'ORIGINAL_CUBE': "NGC_2903.fits"}
    Fits_Files['NOISEMAP'] = f"{Configuration['BASE_NAME']}_noisemap.fits"
    Fits_Files['FITTING_CUBE'] = f"{Configuration['CUBENAME']}_FAT.fits"
    Fits_Files['OPTIMIZED_CUBE'] = f"{Configuration['CUBENAME']}_FAT_opt.fits"
    Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
    Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
    Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
    Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"
    Fits_Files['CHANNEL_MAP'] = f"{Configuration['BASE_NAME']}_chan.fits"
    Configuration['FIX_INCLINATION'] = True
    Configuration['FIX_SDIS'] = False
    Configuration['FIX_PA'] = False
    Configuration['FIX_Z0'] = True
    Configuration['DISTANCE'] = 5.
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

def Test_Ram():
    Configuration = {'FITTING_DIR': '/home/peter/FAT_Main/FAT_Source_Dev/Installation_Check/'}
    Configuration['OUTPUTLOG'] = None
    wf.plot_usage_stats(Configuration,debug = True)

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

def basic():
    print(f" fstring use double or triple quotes")
    print(f" Dictionary entries use single quotes")
    print(f"Below the basic doc structure")
basic.__doc__ =f'''
 NAME:

 PURPOSE:

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



if __name__ == '__main__':
    Test_gflux_dH()
