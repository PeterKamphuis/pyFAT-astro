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

sys.path.insert(1,'/home/peter/FAT_Main/FAT_Source_Dev/Support/')


import modify_template as mt
import read_functions as rf


def Test_Regularise():
    Variables = ['PA', 'PA_2','INCL','INCL_2','SDIS','SDIS_2']
    profile = rf.load_tirific('/home/peter/FAT_Main/FAT_Source_Dev/Installation_Check/Extent_Convergence/Extent_Convergence_Prev.def',Variables)
    radii = rf.load_tirific('/home/peter/FAT_Main/FAT_Source_Dev/Installation_Check/Extent_Convergence/Extent_Convergence_Prev.def', ['RADI'])
    Tirific_Template = {}
    hdr = fits.getheader('Installation_Check/NGC_2903_FAT.fits')
    for i,key in enumerate(Variables):
        if key == 'SBR':
            format = '.2e'
        else:
            format = '.2f'
        Tirific_Template[key] = f"{' '.join([f'{x:{format}}' for x in profile[i]])}"
    Tirific_Template['RADI'] = f"{' '.join([f'{x}' for x in radii[0]])}"
    print(Tirific_Template)
    Configuration = {'INNER_FIX': 4}
    Configuration['NO_RINGS'] = len(profile[0])
    Configuration['OUTPUTLOG'] = None
    Configuration['RING_SIZE'] = 1.
    Configuration['MINIMUM_RINGS'] = 3
    Configuration['SIZE_IN_BEAMS'] = Configuration['NO_RINGS']*Configuration['RING_SIZE']-2
    Variables_red = ['PA','INCL','SDIS']
    for i,key in enumerate(Variables_red):
        sm_profile = mt.smooth_profile(Configuration,Tirific_Template, key,hdr ,debug=True,no_apply=True)
        reg_prof = mt.regularise_profile(Configuration,Tirific_Template, key ,hdr,min_error= 1.,debug = True, no_apply =True)
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii[0],profile[Variables.index(key)],c='b')
        plt.plot(radii[0],reg_prof[0],c='r')
        plt.plot(radii[0],sm_profile[0],c='k')
        plt.savefig(f'Tests/{key}.png')
        plt.close()
        Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')
        plt.plot(radii[0],profile[Variables.index(f'{key}_2')],c='b')
        plt.plot(radii[0],reg_prof[1],c='r')
        plt.plot(radii[0],sm_profile[1],c='k')
        plt.savefig(f'Tests/{key}_2.png')
        plt.close()

if __name__ == '__main__':
    Test_Regularise()
