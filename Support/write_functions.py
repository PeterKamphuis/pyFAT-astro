#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to write text files to Disk

from support_functions import convert_type, print_log,convertRADEC,convertskyangle,set_limit_modifier
from modify_template import set_model_parameters, set_overall_parameters, set_fitting_parameters
import numpy as np

# create or append to the basic ifo file
def basicinfo(Configuration,initialize = False,first_fit = False, second_fit = False,
              RA=[0.,0.], DEC=[0.,0.], VSYS =[0.,0.],
              PA=[0,0], Inclination = [0,0], Max_Vrot = [0,0], Tot_Flux = [0,0], V_mask = [0,0],
              Distance = 0 , DHI = 0):
    if initialize:
        with open(Configuration['FITTING_DIR']+Configuration['BASE_NAME']+'basic_info.txt','w') as file:
            file.write(f'''#This file contains the basic parameters of the Galaxy
#The total flux is determined by the source finder in the initial estimates and the total emission in the masked cube of the final model in the fits
#D_HI is determined as the diameter where the major axis with given PA cuts the 10^20 column of the moment0 map
# {'RA':>25s} {'DEC':>25s} {'VSYS':>20s} {'PA':>20s} {'Inclination':>20s} {'Max VRot':>20s} {'V_mask':>20s} {'Tot FLux':>20s} {'D_HI':>20s} {'Distance':>20s} {'HI Mass':>20s} {'D_HI':>20s}
# {'hh:mm:ss':>25s} {'dd:mm:ss':>25s} {'km/s':>20s} {'deg':>20s} {'deg':>20s} {'km/s':>20s} {'km/s':>20s} {'Jy':>20s} {'arcsec':>20s} {'Mpc':>20s} {'M_solar':>20s} {'kpc':>20s}
# The initial input
''')
    with open(Configuration['FITTING_DIR']+Configuration['BASE_NAME']+'basic_info.txt','a') as file:
        if first_fit:
            file.write(f'''#After the center converged. \n''')
        elif second_fit:
            file.write(f'''#After the radii converged. \n''')
        RAhr,DEChr = convertRADEC(RA[0],DEC[0])
        RA_c = f'{RAhr}+/-{RA[1]*3600.:0.2f}'
        DEC_c = f'{DEChr}+/-{DEC[1]*3600.:0.2f}'
        VSYS_c = f'{VSYS[0]/1000.:.2f}+/-{VSYS[1]/1000.:.2f}'
        PA_c = f'{PA[0]:.2f}+/-{PA[1]:.2f}'
        INCL_c = f'{Inclination[0]:.2f}+/-{Inclination[1]:.2f}'
        MVROT_c = f'{Max_Vrot[0]/1000.:.2f}+/-{Max_Vrot[1]/1000.:.2f}'
        Vmask_c = f'{V_mask[0]/1000.:.2f}+/-{V_mask[1]/1000.:.2f}'
        DHI_a = f'{DHI:.2f}'
        Dist = f'{Distance:.2f}'
        HIMass  = f'{Tot_Flux[0]*2.36E5*Distance**2:.2e}'
        DHI_k = f'{convertskyangle(DHI,Distance):.2f}'
        Flux_c = f'{Tot_Flux[0]:.2f}+/-{Tot_Flux[1]:.2f}'
        file.write(f'''  {RA_c:>25s} {DEC_c:>25s} {VSYS_c:>20s} {PA_c:>20s} {INCL_c:>20s} {MVROT_c:>20s} {Vmask_c:>20s} {Flux_c:>20s} {DHI_a:>20s} {Dist:>20s} {HIMass:>20s} {DHI_k:>20s}
''')
basicinfo.__doc__ = '''

; NAME:
;       basicinfo(Configuration,initialize = True,
              RA=[0.,0.], DEC=[0.,0.], VSYS =[0.,0.],
              PA=[0,0], Inclination = [0,0], Max_Vrot = [0,0], Tot_Flux = [0,0], V_mask = [0,0],
              Distance = [0] , DHI = [0])
;
; PURPOSE:
;       write to the basic info file
;
; CATEGORY:
;       write
;
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''

# Function to write the first def file for a galaxy
def initial_def_file(Configuration, Fits_Files,Tirific_Template, Catalogue, cube_hdr,Initial_Parameters):
    #First we set some basic parameters that will hardly change
    if 'VSYS' in Initial_Parameters:
        Initial_Parameters['VSYS'] = [x/1000. for x in Initial_Parameters['VSYS']]
    set_overall_parameters(Configuration, Fits_Files,Tirific_Template,loops=7,outname = 'Cen_Conv',hdr=cube_hdr, flux = Initial_Parameters['FLUX'][0])
    # Then set the values for the various parameters of the model
    set_model_parameters(Configuration, Tirific_Template,Initial_Parameters, hdr=cube_hdr)
    # Finally we set how these parameters are fitted.
    set_limit_modifier(Configuration,Initial_Parameters['INCL'][0] )
    set_fitting_parameters(Configuration, Tirific_Template, hdr = cube_hdr,stage = 'initial',\
                           systemic =  Initial_Parameters['VSYS'], \
                           inclination = Initial_Parameters['INCL'],
                           pa = Initial_Parameters['PA'],
                           rotation =Initial_Parameters['VROT'],
                           ra = Initial_Parameters['RA'], dec = Initial_Parameters['DEC'])

    tirific(Configuration,Tirific_Template,name = 'Cen_Conv_In.def')
    if 'VSYS' in Initial_Parameters:
        Initial_Parameters['VSYS'] = [x*1000. for x in Initial_Parameters['VSYS']]

initial_def_file.__doc__ = '''

; NAME:
;        initial_def_file(Configuration, Fits_Files,Tirific_Template, Catalogue, cube_hdr,SBR_initial,pa,inclination):
;
; PURPOSE:
;       setup the first def file to be used in the fitting. As it is the first instance it has some special requirements.
;
; CATEGORY:
;       write
;
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''






def sofia(template,name):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")




sofia.__doc__ = '''

; NAME:
;       sofia
;
; PURPOSE:
;       write a sofia2 dictionary into file
;
; CATEGORY:
;       write
;
;
; INPUTS:
;       name = Name of the sofia file
;       template = a template dictionary that matches the sof read template
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''
def tirific(Configuration,Tirific_Template, name = 'tirific.def'):
    with open(Configuration['FITTING_DIR']+name, 'w') as file:
        for key in Tirific_Template:
            if key[0:5] == 'EMPTY':
                file.write('\n')
            else:
                file.write((f"{key} = {Tirific_Template[key]} \n"))

tirific.__doc__ = '''

; NAME:
;       sofia
;
; PURPOSE:
;       write a tirific.def file
;
; CATEGORY:
;       write
;
;
; INPUTS:
;       name = Name of the sofia file
;       template = a template dictionary that matches the sof read template
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''
