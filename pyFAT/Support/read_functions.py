# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used in FAT to read input files

from pyFAT.Support.support_functions import Proper_Dictionary,print_log,convertRADEC,set_limits, remove_inhomogeneities, \
                                obtain_border_pix, obtain_ratios, get_inclination_pa,get_vel_pa
from astropy.io import fits
from scipy import ndimage
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Ellipse
    import matplotlib.axes as maxes
import os
import time
import copy
import numpy as np
try:
    import importlib.resources as import_res
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as import_res

import shutil
import traceback


from pyFAT import Templates as templates
#Function to read a FAT input Catalogue
class BadCatalogueError(Exception):
    pass
class BadConfigurationError(Exception):
    pass
class NoConfigFile(Exception):
    pass
def catalogue(filename, debug = False):
    Catalogue = Proper_Dictionary({})
    tmpfile = open(filename,'r')
    #Define the exsiting catalogue input()
    input_columns = [x.strip().upper() for x in tmpfile.readline().split('|')]
    Catalogue['ENTRIES'] = ['ENTRIES']
    Catalogue['ENTRIES'].extend(input_columns)
    for key in input_columns:
        Catalogue[key] = []

    for line in tmpfile.readlines():
        input = [x.strip() for x  in line.split('|')]
        for i,key in enumerate(input_columns):
            if key == 'DISTANCE':
                Catalogue[key].append(float(input[i]))
            else:
                Catalogue[key].append(input[i])
    if 'NUMBER' in Catalogue['ENTRIES']:
        Catalogue['NUMBER'] = np.array(Catalogue['NUMBER'],dtype=int)

    return Catalogue
catalogue.__doc__ ='''
;+
; NAME:
;       catalogue(filename)
;
; PURPOSE:
;       Read the FAT input catalogue and write into the a dictionary
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the read file
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




#Function to read FAT configuration file into a dictionary
def config_file(input_parameters, start_dir, debug = False):
    No_File = True
    while No_File:
        try:
            if input_parameters.configfile == 'ChecK.ConfiG':
                from pyFAT import Installation_Check as IC
                tmpfile = import_res.open_text(IC,'FAT_INPUT.config')
            else:
                with open(input_parameters.configfile, 'r') as tmp:
                    tmpfile = tmp.readlines()
            No_File = False
        except:
            print(traceback.print_exc())
            input_parameters.configfile = input('''
                        You have provided a config file but it can't be found.
                        If you want to provide a config file please give the correct name.
                        Else press CTRL-C to abort.
                ''')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','WARP_OUTPUT','TWO_STEP']
    string_keys = ['OUTPUTLOG', 'OUTPUTCATALOGUE','MAINDIR','CATALOGUE']
    integer_keys = ['STARTGALAXY','ENDGALAXY','MAPS_OUTPUT','OPT_PIXELBEAM','FINISHAFTER']
    # Separate the keyword names
    for tmp in tmpfile:
        if tmp[0] != '#':
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
            add_key = tmp.split('=', 1)[0].strip().upper()
            if add_key in boolean_keys:
                invalid_input = True
                inp = tmp.split('=', 1)[1].strip()
                while invalid_input:
                    if inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes" or inp[0] == '1':
                        value = True
                        invalid_input = False
                    elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no" or inp[0] == '0':
                        value = False
                        invalid_input = False
                    else:
                        inp = input("The parameter {} in the configuration file  must be true/false or yes/no. Please give the correct value. \n".format(add_key))
                Configuration[add_key] = value
            elif add_key in string_keys:
                Configuration[add_key] = tmp.split('=', 1)[1].strip()
            elif add_key in integer_keys:
                Configuration[add_key] = int(tmp.split('=', 1)[1].strip())
            else:
                Configuration[add_key] = float(tmp.split('=', 1)[1].strip())
    #if we are checking the installation then the maindir, outputcatalogue and
    #Log go into the original_dir+installation_check.
    if input_parameters.installation_check:
        if Configuration['CATALOGUE'] != 'Installation_Check/FAT_Input_Catalogue.txt' or \
           Configuration['MAINDIR'] != 'Installation_Check/' or \
           Configuration['OUTPUTCATALOGUE'] != 'Installation_Check/Output_N2903.txt':
           raise BadConfigurationError('You can not modify the Installation Check input. It is solely for checking the Installation. Aborting')
        test_files = ['FAT_Input_Catalogue.txt','NGC_2903.fits']
        test_dir = start_dir+'/FAT_Installation_Check/'
        if not os.path.isdir(test_dir):
            print('We making a dir?)')
            os.mkdir(test_dir)
        else:
            for file in test_files:
                try:
                    os.remove(test_dir+file)
                except:
                    pass
        my_resources = import_res.files("pyFAT.Installation_Check")
        for file in test_files:
            data = (my_resources / file).read_bytes()
            print('Is it the opening')
            print(test_dir+file)
            with open(test_dir+file,'w+b') as tmp:
                tmp.write(data)
        Configuration['CATALOGUE'] = test_dir+'FAT_Input_Catalogue.txt'
        Configuration['MAINDIR'] = test_dir
        Configuration['OUTPUTCATALOGUE'] =test_dir+'Output_N2903.txt'
    #Make the input idiot safe
    if Configuration['MAINDIR'][-1] != '/':
        Configuration['MAINDIR'] = Configuration['MAINDIR']+'/'

    while not os.path.isdir(Configuration['MAINDIR']):
        Configuration['MAINDIR'] = input('''
                    Your main fitting directory ({}) does not exist.
                    Please provide the correct directory.
                    '''.format(Configuration['MAINDIR']))
    while not os.path.exists(Configuration['CATALOGUE']):
        Configuration['CATALOGUE'] = input('''
                    Your input catalogue ({}) does not exist.
                    Please provide the correct file name.
                    '''.format(Configuration['CATALOGUE']))
    #The output catalogue only needs to be in a valid directory as we create it
    output_catalogue_dir = Configuration['OUTPUTCATALOGUE'].split('/')
    if len(output_catalogue_dir) > 1:
        check_dir = '/'.join(output_catalogue_dir[:-1])
        while not os.path.isdir(check_dir):
            check_dir= input('''
                    The directory for your output catalogue ({}) does not exist.
                    Please provide the correct directory name.
                    '''.format(Configuration['OUTPUTCATALOGUE']))
            Configuration['OUTPUTCATALOGUE'] = check_dir+'/'+output_catalogue_dir[-1]


    required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','HANNING',\
                                   'STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT',\
                                   'RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR',\
                                    'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM',\
                                     'MAPS_OUTPUT','WARP_OUTPUT','TWO_STEP']

    for key in required_configuration_keys:
        if key not in Configuration:
            if key == 'TWO_STEP':
                Configuration[key] = False
            elif key == 'STARTGALAXY':
                Configuration[key] = 0
            elif key == 'FINISHAFTER':
                Configuration[key] = 2
            elif key == 'TESTING':
                Configuration[key] = 0
            elif key == 'START_POINT': #Previously calle allnew
                Configuration[key] = 1
            elif key == 'ENDGALAXY':
                Configuration[key] = -1
            elif key == 'NEW_OUTPUT':   # Called newresult in the gdl code
                Configuration[key] = True
            elif key == 'HANNING':
                Configuration[key] = False
            elif key == 'RING_SIZE': #Previosuly called RINGSPACING in
                Configuration[key] = 1.1
            elif key == 'FIX_INCLINATION': #Previosuly called fix_incl
                Configuration[key] = False
            elif key == 'FIX_PA':
                Configuration[key] = False
            elif key == 'FIX_SDIS':
                Configuration[key] = False
            elif key == 'FIX_Z0':
                Configuration[key] = True
            elif key == 'OPT_PIXELBEAM':
                Configuration[key] = 4
            elif key == 'MAPS_OUTPUT': # Previously called bookkeeping
                Configuration[key] = 3
            elif key == 'WARP_OUTPUT':
                Configuration[key] = False
            elif key == 'OUTPUTLOG':
                Configuration[key] = None
            else:
                raise BadConfigurationError('Something has gone wrong reading the required config key. This should never ever happen. Please file an issue on github')
    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['MAPS_OUTPUT'] == 5:
        Configuration['MAPS_OUTPUT'] = 4
    # We double the fix keys so we can  modify one while keeping the original as well
    fix_keys = ['FIX_PA','FIX_INCLINATION','FIX_SDIS','FIX_Z0']
    for key in fix_keys:
        Configuration[key] = [Configuration[key],Configuration[key]]
    return Configuration
config_file.__doc__ ='''
;+
; NAME:
;       config_file(input_parameters, start_dir)
;
; PURPOSE:
;       Read the FAT config file and write into the a dictionary
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the config file
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

def extract_vrot(Configuration, hdr,map ,angle,center, debug= False):
    if debug:
        print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
{'':8s} PA= {angle}
{'':8s} center= {center}
''',Configuration['OUTPUTLOG'], debug = True,screen = True)
    x1,x2,y1,y2 = obtain_border_pix(hdr,angle,center)
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    maj_resolution = np.sqrt(((x2-x1)/1000.)**2+((y2-y1)/1000.)**2)
    #maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(angle)))+abs(abs(y2-y1)/1000.*np.cos(np.radians(angle)))
    maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(angle)))+abs(abs(center[1])*np.cos(np.radians(angle))))
    #neg_index = np.where(maj_axis < 0.)[0]
    #pos_index = np.where(maj_axis > 0.)[0]
    neg_index = np.where(maj_profile > 0.)[0]
    pos_index = np.where(maj_profile < 0.)[0]
    if debug:
        print_log(f'''EXTRACT_VROT: The extracted major axis velocities
{'':8s} {maj_profile}
{'':8s} resolution = {maj_resolution}
''',Configuration['OUTPUTLOG'], debug = False)
    avg_profile = []
    neg_profile = []
    pos_profile = []
    diff = 0.
    for i in range(np.nanmin([neg_index.size,pos_index.size])):
        avg_profile.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],-1*maj_profile[pos_index[i]]]))
        neg_profile.append(maj_profile[neg_index[neg_index.size-i-1]])
        pos_profile.append(-1*maj_profile[pos_index[i]])
        #correct for beam smearing in the center
        beam_back = -1*int(hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution))

        if i > hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution) and i < -2.*beam_back:

        #if i*np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution < 1.5*hdr['BMAJ'] and i > int() :
            avg_profile[beam_back] = avg_profile[beam_back]+0.5*avg_profile[int(beam_back/2.)]+0.1*avg_profile[-1]
            neg_profile[beam_back] = neg_profile[beam_back]+0.5*neg_profile[int(beam_back/2.)]+0.1*neg_profile[-1]
            pos_profile[beam_back] = pos_profile[beam_back]+0.5*pos_profile[int(beam_back/2.)]+0.1*pos_profile[-1]


        #if neg_profile[-1]/pos_profile[-1] < 0.:
        #    avg_profile.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],-1*maj_profile[pos_index[i]]]))
        #elif neg_profile[-1] < 0.:
        #    avg_profile.append(np.nanmin([maj_profile[neg_index[neg_index.size-i-1]],-1*maj_profile[pos_index[i]]]))
        #else:
        #    avg_profile.append(np.nanmax([maj_profile[neg_index[neg_index.size-i-1]],-1*maj_profile[pos_index[i]]]))
    if debug:
        print_log(f'''EXTRACT_VROT: starting extraction of initial VROT.
{'':8s} negative= {neg_profile}
{'':8s} positive= {pos_profile}
{'':8s} avreage= {avg_profile}
''',Configuration['OUTPUTLOG'], debug = False)

    ring_size_req = hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution)
    profile = np.array(avg_profile[0::int(ring_size_req)],dtype=float)
    if debug:
        print_log(f'''EXTRACT_VROT: Constructing the final RC
{'':8s} initial RC= {profile}
{'':8s} at step width= {ring_size_req}
''',Configuration['OUTPUTLOG'], debug = False)
    try:
        profile[np.isnan(profile)] = profile[~np.isnan(profile)][-1]
    except IndexError:
        profile = []

    return profile

extract_vrot.__doc__ ='''
;+
; NAME:
;       extract_averaged(Configuration, hdr,map ,angle,center, debug= False):
;
; PURPOSE:
;       Extract a profile averaged over both sides.
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration =
;       Fits_Files =
;
;
; OPTIONAL INPUTS:
;       center=
;       pa =
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the config file
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

# Function to get the PA and inclination from the moment 0 for initial estimates

def guess_orientation(Configuration,Fits_Files, center = None, debug = False):
    #open the moment 0
    if debug:
        print_log(f'''GUESS_ORIENTATION: starting extraction of initial parameters.
''',Configuration['OUTPUTLOG'], debug = True)
    Image = fits.open(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files['MOMENT0'],\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = Image[0].data
    hdr = Image[0].header
    mom0 = copy.deepcopy(Image)
    Image.close()
    if not center:
        center = [hdr['NAXIS1']/2.,hdr['NAXIS2']/2]
    Image = fits.open(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files['CHANNEL_MAP'],\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    noise_map = np.sqrt(Image[0].data)*Configuration['NOISE']*Configuration['CHANNEL_WIDTH']
    SNR = np.nanmean(map[noise_map > 0.]/noise_map[noise_map > 0.])
    noise_hdr = Image[0].header

    median_noise_in_map = np.nanmedian(noise_map[noise_map > 0.])
    minimum_noise_in_map = np.nanmin(noise_map[noise_map > 0.])
    map[0.5*minimum_noise_in_map > noise_map] = 0.
    scale_factor = set_limits(SNR/3.*minimum_noise_in_map/median_noise_in_map, 0.05, 1.)
    if debug:
        print_log(f'''GUESS_ORIENTATION: We find SNR = {SNR} and a scale factor {scale_factor} and the noise median {median_noise_in_map}
{'':8s} minimum {minimum_noise_in_map}
''',Configuration['OUTPUTLOG'], debug = True)
    inclination, pa, maj_extent = get_inclination_pa(Configuration, map,mom0, hdr, center, cutoff = scale_factor* median_noise_in_map, debug = debug)

    maj_extent = maj_extent+(hdr['BMAJ']*0.2/scale_factor)
            #map[3*minimum_noise_in_map > noise_map] = 0.
    # From these estimates we also get an initial SBR
    x1,x2,y1,y2 = obtain_border_pix(hdr,pa[0],center)
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    maj_resolution = np.sqrt(((x2-x1)/1000.)**2+((y2-y1)/1000.)**2)
    #maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(pa[0])))+abs(abs(y2-y1)/1000.*np.cos(np.radians(pa[0])))
    maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa[0])))+abs(abs(center[1])*np.cos(np.radians(pa[0]))))

    # let's get an intensity weighted center for the extracted profile.

    center_of_profile = np.sum(maj_profile*maj_axis)/np.sum(maj_profile)
    neg_index = np.where(maj_axis < 0.)[0]
    pos_index = np.where(maj_axis > 0.)[0]
    avg_profile = []
    neg_profile = []
    pos_profile = []
    diff = 0.
    for i in range(np.nanmin([neg_index.size,pos_index.size])):
        avg_profile.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],maj_profile[pos_index[i]]]))
        neg_profile.append(maj_profile[neg_index[neg_index.size-i-1]])
        pos_profile.append(maj_profile[pos_index[i]])
        diff = diff+abs(avg_profile[i]-neg_profile[i])+abs(avg_profile[i]-pos_profile[i])
    diff = diff/np.nanmin([neg_index.size,pos_index.size])
    if debug:
        print_log(f'''GUESS_ORIENTATION:'Beam, center of profile, center
{'':8s}{hdr['BMAJ']*0.5/np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])} {center_of_profile} {center}
''',Configuration['OUTPUTLOG'], debug = False)
    # if the center of the profile is more than half a beam off from the Sofia center let's see which on provides a more symmetric profile
    if abs(center_of_profile) > hdr['BMAJ']*0.5/np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])]):
        if debug:
                print_log(f'''GUESS_ORIENTATION: The SoFiA center and that of the SBR profile are separated by more than half a bemea apart.
{'':8s}GUESS_ORIENTATION: Determining the more symmetric profile.
''',Configuration['OUTPUTLOG'], debug = False)
        neg_index = np.where(maj_axis < center_of_profile)[0]
        pos_index = np.where(maj_axis > center_of_profile)[0]
        avg_profile_new = []
        neg_profile_new = []
        pos_profile_new = []
        diff_new =0.
        for i in range(np.nanmin([neg_index.size,pos_index.size])):
            avg_profile_new.append(np.nanmean([maj_profile[neg_index[neg_index.size-i-1]],maj_profile[pos_index[i]]]))
            neg_profile_new.append(maj_profile[neg_index[neg_index.size-i-1]])
            pos_profile_new.append(maj_profile[pos_index[i]])
            diff_new = diff_new+abs(avg_profile_new[i]-neg_profile_new[i])+abs(avg_profile_new[i]-pos_profile_new[i])
        diff_new = diff_new/np.nanmin([neg_index.size,pos_index.size])
        if diff_new < diff:
            if debug:
                print_log(f'''GUESS_ORIENTATION: We are updating the center from {center[0]},{center[1]} to {center[0]+center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution},{center[1]+center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution}
''',Configuration['OUTPUTLOG'], debug = False)
            avg_profile = avg_profile_new
            center[0] = center[0]-center_of_profile/(2.*np.sin(np.radians(pa[0])))*maj_resolution
            center[1] = center[1]+center_of_profile/(2.*np.cos(np.radians(pa[0])))*maj_resolution
            maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa[0])))+abs(abs(center[1])*np.cos(np.radians(pa[0]))))


    ring_size_req = hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution)
    SBR_initial = avg_profile[0::int(ring_size_req)]/(np.pi*abs(hdr['BMAJ']*3600.*hdr['BMIN']*3600.)/(4.*np.log(2.))) # Jy*km/s
    SBR_initial =np.hstack((SBR_initial[0],SBR_initial,SBR_initial[-1]))

    SBR_initial[0:3] = SBR_initial[0:3] * (1.2 -float(inclination[0])/90.)


    #We need to know which is the approaching side and which is receding



    Image = fits.open(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files['MOMENT1'],\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = Image[0].data
    hdr = Image[0].header
    Image.close()
    map[3*minimum_noise_in_map > noise_map] = float('NaN')
    noise_map = []
    vel_pa = get_vel_pa(Configuration,map,center=center,debug=debug)

    x1,x2,y1,y2 = obtain_border_pix(hdr,pa[0],center)
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    #maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(pa[0])))+abs(abs(y2-y1)/1000.*np.cos(np.radians(pa[0])))
    maj_resolution = np.sqrt(((x2-x1)/1000.)**2+((y2-y1)/1000.)**2)
    maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa[0])))+abs(abs(center[1])*np.cos(np.radians(pa[0]))))
    loc_max = np.mean(maj_axis[np.where(maj_profile == np.nanmax(maj_profile))[0]])
    if loc_max > 0.:
            pa[0] = pa[0]+180
            print_log(f'''GUESS_ORIENTATION: We have modified the pa by 180 deg as we found the maximum velocity west of the center.
''' , Configuration['OUTPUTLOG'])
    if abs(pa[0]-vel_pa[0]) > 25. and ~np.isnan(vel_pa[0]):
        pa=vel_pa
    else:
        if ~np.isnan(vel_pa[0]):
            pa = [np.nansum([vel_pa[0]/vel_pa[1],pa[0]/pa[1]])/np.nansum([1./vel_pa[1],1./pa[1]]),2.*1./np.sqrt(np.nansum([1./vel_pa[1],1./pa[1]]))]
    #As python is utterly moronic the center goes in back wards to the map
    if debug:
        print_log(f'''GUESS_ORIENTATION: We found the following initial VSYS:
{'':8s}vsys = {map[int(round(center[1])),int(round(center[0]))]}, at {center}
''',Configuration['OUTPUTLOG'], debug = False)
    map = map  - map[int(round(center[1])),int(round(center[0]))]
    VROT_initial = extract_vrot(Configuration, hdr,map ,pa[0],center, debug= debug)
    min_RC_length= len(VROT_initial)
    if pa[1] < 10.:
        for x in [pa[0]-pa[1],pa[0]-pa[1]/2.,pa[0]+pa[1]/2.,pa[0]+pa[1]]:
            tmp  = extract_vrot(Configuration, hdr,map ,x,center, debug= debug)
            if len(tmp) > 0.:
                RC_length = len(tmp)
                if RC_length < min_RC_length:
                    min_RC_length = RC_length
                VROT_initial = np.vstack((VROT_initial[:min_RC_length],tmp[:min_RC_length]))
        VROT_initial = np.mean(VROT_initial,axis=0)
    map= []
    VROT_initial[0] = 0
    VROT_initial = VROT_initial/np.sin(np.radians(inclination[0]))
    if debug:
        print_log(f'''GUESS_ORIENTATION: We found the following initial rotation curve:
{'':8s}GUESS_ORIENTATION: RC = {VROT_initial}
''',Configuration['OUTPUTLOG'], debug = False)
    return np.array(pa),np.array(inclination),SBR_initial,maj_extent,center[0],center[1],VROT_initial
guess_orientation.__doc__ ='''
;+
; NAME:
;       guess_orientation(Configuration, Fits_Files,center = None)
;
; PURPOSE:
;       Read the moment map and noise map and estimate the the pa and inclination from axis ratios
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the config file
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      np.mean, astropy.io.fits, obtain_ratios, np.linspace
;
; EXAMPLE:
;
;
'''
# load the basic info file to get the Sofia FAT Initial_Estimates
def load_basicinfo(filename, Variables = ['RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask','Tot FLux','D_HI','Distance','HI_Mass' ,'D_HI' ], unpack = True, debug = False):
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
    del Var_inFile[0]
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


#Function for loading the variables of a tirific def file into a set of variables to be used
def load_template(Template,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True  ,debug = False):
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    if debug:
        print_log(f'''LOAD_TEMPLATE: We get the following number of rings from the template {Template['NUR']}
''',None, debug=debug)
    numrings = int(Template['NUR'])
    outputarray=np.zeros((numrings,len(Variables)),dtype=float)
    counter = 0
    for var in Variables:
        if debug:
            print_log(f'''LOAD_TEMPLATE: We are processing {var}.
{'':8s}LOAD_TEMPLATE: With the following values {Template[var]}''',None, debug=debug)
        tmp =  np.array(Template[var].rsplit(),dtype=float)
        outputarray[0:len(tmp),counter] = tmp[0:len(tmp)]
        counter +=1

    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray

#Function for loading the variables of a tirific def file into a set of variables to be used
def load_tirific(filename,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',
                 'Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True , debug = False ):
    if debug:
        print_log(f'''LOAD_TIRIFIC: Starting to extract the following paramaters:
{'':8s}{Variables}
''',None,screen=True, debug = True)
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    numrings = []
    while len(numrings) < 1:
        time.sleep(0.1)
        with open(filename, 'r') as tmp:
            numrings = [int(e.split('=')[1].strip()) for e in tmp.readlines() if e.split('=')[0].strip().upper() == 'NUR']



    #print(numrings)tmp
    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
    with open(filename, 'r') as tmp:
        unarranged = tmp.readlines()
    # Separate the keyword names
    for line in unarranged:
        var_concerned = str(line.split('=')[0].strip().upper())
        #if debug:
        #    print_log(f'''LOAD_TIRIFIC: extracting line
#{'':8s}{var_concerned}.
#''',None,screen=False, debug = True)
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
                varpos = np.where(var_concerned[2:].strip() == Variables)[0]
#                if debug:
#                    print_log(f'''LOAD_TIRIFIC: comparing {var_concerned[2:].strip()} to the variables.
#{'':8s}Found {varpos}.
#''',None,screen=True, debug = True)
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

# function to read the sofia catalogue
def sofia_catalogue(Configuration, Variables =['name','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra',\
                    'dec','v_app','f_sum','kin_pa','w50','err_f_sum','err_x','err_y','err_z'], header = None , debug = False):
    if debug:
        print_log(f'''SOFIA_CATLOGUE: Reading the source from the catalogue.
''',Configuration['OUTPUTLOG'],debug= debug)
    outlist = [[] for x in Variables]
    with open(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt') as sof_cat:
        for line in sof_cat.readlines():
            tmp =line.split()
            if line.strip() == '' or line.strip() == '#':
                pass
            elif tmp[0] == '#' and len(tmp) > 1:
                if tmp[1].strip().lower() == 'name':
                    # get the present columns
                    input_columns  = [x.strip() for x in tmp[1:]]
                    #determin their location in the line
                    column_locations = []
                    for col in input_columns:
                        column_locations.append(line.find(col)+len(col))
                    # check that we found all parameters
                    for value in Variables:
                        if value.lower() in input_columns:
                            continue
                        else:
                            log_statement=f'''READ_SOFIA_CATALOGUE: We cannot find the required column for {value} in the sofia catalogue.
{"":8s}READ_SOFIA_CATALOGUE: This can happen because a) you have tampered with the sofiainput.txt file in the Support directory,
{"":8s}READ_SOFIA_CATALOGUE: b) you are using an updated version of SoFiA2 where the names have changed and FAT is not yet updated.'
{"":8s}READ_SOFIA_CATALOGUE:    In this case please file a bug report at https://github.com/PeterKamphuis/FAT/issues/'
{"":8s}READ_SOFIA_CATALOGUE: c) You are using pre processed SoFiA output of your own and do not have all the output'
{"":8s}READ_SOFIA_CATALOGUE:    Required output is {','.join(Variables)})
'''
                            print_log(log_statement,Configuration['OUTPUTLOG'],debug= debug)
                            raise BadCatalogueError("READ_SOFIA_CATALOGUE: The required columns could not be found in the sofia catalogue.")
            else:
                for col in Variables:
                    if input_columns.index(col) == 0:
                        start = 0
                    else:
                        start = column_locations[input_columns.index(col)-1]
                    end = column_locations[input_columns.index(col)]
                    outlist[Variables.index(col)].append(line[start:end].strip())
    #if we have a header we convert w50 or w20
    if ('w50' in Variables or 'w20' in Variables or 'f_sum' in Variables or 'err_f_sum' in Variables) and header:
        if 'w50' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('w50')][i] = float(outlist[Variables.index('w50')][i])*header['CDELT3']
        if 'w20' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('w20')][i] = float(outlist[Variables.index('w20')][i])*header['CDELT3']
        if 'f_sum' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('f_sum')][i] = float(outlist[Variables.index('f_sum')][i])/Configuration['PIX_PER_BEAM']*header['CDELT3']/1000.
        if 'err_f_sum' in Variables:
            for i in range(len(outlist[0])):
                outlist[Variables.index('err_f_sum')][i] = float(outlist[Variables.index('err_f_sum')][i])/Configuration['PIX_PER_BEAM']*header['CDELT3']/1000.

    # we want to fit a specific source
    if len(outlist[0]) > 1 and 'f_sum' in Variables and header:
        if debug:
            print_log(f'''SOFIA_CATALOGUE: Multiple sources were found we will try to select the correct one.
''',Configuration['OUTPUTLOG'],debug= debug)
        found = False
        beam_edge=2.
        if Configuration['VEL_SMOOTH_EXTENDED']:
            vel_edge = 0.
        else:
            vel_edge = 1.
        while not found:
            many_sources  = copy.deepcopy(outlist)
            # We want to exclude any edge sources
            for i in range(len(many_sources[0])):
                edge = False
                diff = np.array([many_sources[Variables.index('x_min')][i],
                        abs(float(many_sources[Variables.index('x_max')][i])-header['NAXIS1']),
                        many_sources[Variables.index('y_min')][i],
                        abs(float(many_sources[Variables.index('y_max')][i])-header['NAXIS2'])
                        ],dtype=float)
                if debug:
                    print_log(f'''SOFIA_CATALOGUE: We find these differences and this edge size {beam_edge*header['BMAJ']/((abs(header['CDELT1'])+abs(header['CDELT2']))/2.)}
    {'':8s} diff  = {diff}
    ''',Configuration['OUTPUTLOG'],debug= debug,screen = True)
                if np.where(diff < beam_edge*header['BMAJ']/((abs(header['CDELT1'])+abs(header['CDELT2']))/2.))[0].size:
                    edge = True
                diff = np.array([many_sources[Variables.index('z_min')][i],
                        abs(float(many_sources[Variables.index('z_max')][i])-float(header['NAXIS3']))],dtype=float)
                if debug:
                    print_log(f'''SOFIA_CATALOGUE: And for velocity edge =  2
    {'':8s} diff  = {diff}
    ''',Configuration['OUTPUTLOG'],debug= debug,screen = True)
                if np.where(diff < vel_edge)[0].size:
                    edge = True
                if debug:
                    print_log(f'''SOFIA_CATALOGUE: Edge = {edge}
    ''',Configuration['OUTPUTLOG'],debug= debug,screen = True)
                if edge:
                    many_sources[Variables.index('f_sum')][i]=0.
            if np.nansum(many_sources[Variables.index('f_sum')]) == 0.:
                if beam_edge > 0.5:
                    beam_edge = beam_edge/2.
                elif vel_edge > 0.5:
                    vel_edge = vel_edge/2.
                else:
                    raise BadCatalogueError("The found sources are too close to the edges of the cube.")

            else:
                found = True
        if debug:
            print_log(f'''SOFIA_CATALOGUE: after checking edges we find these fluxes
{'':8s}{many_sources[Variables.index('f_sum')]}
''',Configuration['OUTPUTLOG'],debug= debug,screen = True)

        fluxes = np.array(many_sources[Variables.index('f_sum')],dtype= float)
        outlist = []
        #We want the source with the most total flux.
        index = np.where(np.nanmax(fluxes) == fluxes)[0][0]
        print_log(f'''SOFIA_CATALOGUE: We select the {index} source of this list.
''',Configuration['OUTPUTLOG'],debug= debug, screen =True)
        outlist = [x[index] for x in many_sources]
    else:
        outlist = [x[0] for x in outlist]
    if debug:
        print_log(f'''SOFIA_CATALOGUE: we found these values
{'':8s}{outlist}
''',Configuration['OUTPUTLOG'],debug= debug,screen = True)
    return outlist

sofia_catalogue.__doc__ ='''

;+
; NAME:
;       sofia_catalogue
;
; PURPOSE:
;       Read the sofia catalogue and extract several basic parameters from into a list.
;       In order to read
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the read file
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




def tirific_template(filename = '', debug = False):
    if filename == '':
        template = import_res.open_text(templates, 'template.def')
    elif filename == 'Installation_Check':
        from pyFAT import Installation_Check as IC
        template = import_res.open_text(IC, 'ModelInput.def')
    else:
        with open(filename, 'r') as tmp:
            template = tmp.readlines()
    result = Proper_Dictionary()
    counter = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip().upper())
        if key == '':
            result[f'EMPTY{counter}'] = line
            counter += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result
tirific_template.__doc__ ='''

;+
; NAME:
;       tirific_template
;
; PURPOSE:
;       Read a tirfic def file into a dictionary
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       filename = Name of the def file
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the read file
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

def sofia_template(debug = False):
    template = import_res.open_text(templates, 'sofia_template.par')
    result = {}
    counter = 0
    counter2 = 0
    # Separate the keyword names
    for line in template:
        key = str(line.split('=')[0].strip())
        if key == '':
            result[f'EMPTY{counter}'] = line
            counter += 1
        elif key[0] == '#':
            result[f'HASH{counter2}'] = line
            counter2 += 1
        else:
            result[key] = str(line.split('=')[1].strip())
    return result
sofia_template.__doc__ ='''

;+
; NAME:
;       sofia_template
;
; PURPOSE:
;       Read a sofia2 file into a dictionary
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       filename = Name of the sofia file
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the read file
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
