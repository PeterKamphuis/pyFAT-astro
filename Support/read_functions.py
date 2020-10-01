#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in FAT to read input files

from support_functions import Proper_Dictionary,print_log,convertRADEC,set_limits, remove_inhomogeneities
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
import copy
import numpy as np
import traceback
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
            tmpfile = open(input_parameters.configfile, 'r')
            No_File = False
        except:
            print(traceback.print_exc())
            input_parameters.configfile = input('''
                        You have provided a config file but it can't be found.
                        If you want to provide a config file please give the correct name.
                        Else press CTRL-C to abort.
                ''')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','WARP_OUTPUT']
    string_keys = ['OUTPUTLOG', 'OUTPUTCATALOGUE','MAINDIR','CATALOGUE']
    integer_keys = ['STARTGALAXY','ENDGALAXY','MAPS_OUTPUT','OPT_PIXELBEAM','FINISHAFTER']
    # Separate the keyword names
    for tmp in tmpfile.readlines():
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
        Configuration['CATALOGUE'] = start_dir+'/Installation_Check/FAT_Input_Catalogue.txt'
        Configuration['MAINDIR'] = start_dir+'/Installation_Check/'
        Configuration['OUTPUTCATALOGUE'] = start_dir+'/Installation_Check/Output_N2903.txt'
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


    required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','HANNING','STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT','RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR', 'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM', 'MAPS_OUTPUT','WARP_OUTPUT']

    for key in required_configuration_keys:
        if key not in Configuration:
            if key == 'STARTGALAXY':
                Configuration[key] = 0
            if key == 'FINISHAFTER':
                Configuration[key] = 2
            if key == 'TESTING':
                Configuration[key] = 0
            if key == 'START_POINT': #Previously calle allnew
                Configuration[key] = 1
            if key == 'ENDGALAXY':
                Configuration[key] = -1
            if key == 'NEW_OUTPUT':   # Called newresult in the gdl code
                Configuration[key] = True
            if key == 'HANNING':
                Configuration[key] = False
            if key == 'RING_SIZE': #Previosuly called RINGSPACING in
                Configuration[key] = 1.1
            if key == 'FIX_INCLINATION': #Previosuly called fix_incl
                Configuration[key] = False
            if key == 'FIX_PA':
                Configuration[key] = False
            if key == 'FIX_SDIS':
                Configuration[key] = False
            if key == 'FIX_Z0':
                Configuration[key] = True
            if key == 'OPT_PIXELBEAM':
                Configuration[key] = 4
            if key == 'MAPS_OUTPUT': # Previously called bookkeeping
                Configuration[key] = 3
            if key == 'WARP_OUTPUT':
                Configuration[key] = False
            if key == 'OUTPUTLOG':
                Configuration[key] = None
    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['MAPS_OUTPUT'] == 5:
        Configuration['MAPS_OUTPUT'] = 4
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
    noise_hdr = Image[0].header
    Image.close()
    median_noise_in_map = np.nanmedian(noise_map[noise_map > 0.])
    minimum_noise_in_map = np.nanmin(noise_map[noise_map > 0.])
    map[3*minimum_noise_in_map > noise_map] = 0.
    for i in [0,1]:
        # now we need to get profiles under many angles let's say 100
        #extract the profiles under a set of angles
        angles = np.linspace(0, 180, 180)

        ratios, maj_extent = obtain_ratios(map, hdr, center, angles)

        max_index = np.where(ratios == np.nanmax(ratios))[0]
        if max_index.size > 1:
            max_index =max_index[0]
        min_index = np.where(ratios == np.nanmin(ratios))[0]
        if min_index.size > 1:
            min_index =min_index[0]
        #get a 10% bracket

        tenp_max_index = np.where(ratios > np.nanmax(ratios)*0.9)[0]
        tenp_min_index = np.where(ratios < np.nanmin(ratios)*1.1)[0]
        if tenp_max_index.size <= 1:
            tenp_max_index= [max_index-2,max_index+2]
        if tenp_min_index.size <= 1:
            tenp_min_index= [min_index-2,min_index+2]
        pa = np.mean([angles[min_index],angles[max_index]-90.])
        pa_error = set_limits(np.mean([abs(angles[tenp_min_index[0]]-angles[min_index]),\
                            abs(angles[tenp_min_index[-1]]-angles[min_index]),\
                            abs(angles[tenp_max_index[0]]-angles[max_index]), \
                            abs(angles[tenp_min_index[-1]]-angles[max_index])]), \
                            0.5,15.)
        ratios[ratios < 0.204] = 0.204
        ratios[1./ratios < 0.204] = 1./0.204
        inclination = np.mean([np.degrees(np.arccos(np.sqrt((ratios[min_index]**2-0.2**2)/0.96))) \
                              ,np.degrees(np.arccos(np.sqrt(((1./ratios[max_index])**2-0.2**2)/0.96))) ])

        if i == 0:
            inclination_error = set_limits(np.nanmean([abs(np.degrees(np.arccos(np.sqrt((ratios[min_index]**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt((ratios[tenp_min_index[0]]**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt((ratios[min_index]**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt((ratios[tenp_min_index[-1]]**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt(((1./ratios[max_index])**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt(((1./ratios[tenp_max_index[0]])**2-0.2**2)/0.96)))),\
                                     abs(np.degrees(np.arccos(np.sqrt(((1./ratios[max_index])**2-0.2**2)/0.96))) - np.degrees(np.arccos(np.sqrt(((1./ratios[tenp_max_index[0]])**2-0.2**2)/0.96))))]), \
                                     2.5,25.)
            if not np.isfinite(inclination_error):
                inclination_error = 90./inclination                         
        if ratios[max_index]-ratios[min_index] < 0.4:
            inclination = float(inclination-(inclination*0.04/(ratios[max_index]-ratios[min_index])))
            if i == 0:
                inclination_error = float(inclination_error*0.4/(ratios[max_index]-ratios[min_index]))
        if maj_extent/hdr['BMAJ'] < 4:
            inclination = float(inclination+(inclination/10.*np.sqrt(4./(maj_extent/hdr['BMAJ']))))
            if i == 0:
                inclination_error = float(inclination_error*4./(maj_extent/hdr['BMAJ']))
        if i == 0 and inclination < 75.:
            mom0 = remove_inhomogeneities(Configuration,mom0,inclination=inclination, pa = pa, center = center,WCS_center = False, debug=debug)
            map = mom0[0].data
            #map[3*minimum_noise_in_map > noise_map] = 0.
    # From these estimates we also get an initial SBR
    x1,x2,y1,y2 = obtain_border_pix(hdr,pa,center)
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(pa)))+abs(abs(y2-y1)/1000.*np.cos(np.radians(pa)))
    maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa)))+abs(abs(center[1])*np.cos(np.radians(pa))))

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
                print_log(f'''GUESS_ORIENTATION: We are updating the center from {center[0]},{center[1]} to {center[0]+center_of_profile/(2.*np.sin(np.radians(pa)))*maj_resolution},{center[1]+center_of_profile/(2.*np.cos(np.radians(pa)))*maj_resolution}
''',Configuration['OUTPUTLOG'], debug = False)
            avg_profile = avg_profile_new
            center[0] = center[0]-center_of_profile/(2.*np.sin(np.radians(pa)))*maj_resolution
            center[1] = center[1]+center_of_profile/(2.*np.cos(np.radians(pa)))*maj_resolution
            maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa)))+abs(abs(center[1])*np.cos(np.radians(pa))))


    ring_size_req = Configuration['RING_SIZE']*hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])*maj_resolution)
    SBR_initial = avg_profile[0::int(ring_size_req)]/(np.pi*abs(hdr['BMAJ']*3600.*hdr['BMIN']*3600.)/(4.*np.log(2.))) # Jy*km/s
    SBR_initial =np.hstack((SBR_initial[0],SBR_initial,SBR_initial[-1]))
    #We need to know which is the approaching side and which is receding



    Image = fits.open(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files['MOMENT1'],\
            uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    map = Image[0].data/1000.
    hdr = Image[0].header
    Image.close()
    map[3*minimum_noise_in_map > noise_map] = float('NaN')
    noise_map = []

    x1,x2,y1,y2 = obtain_border_pix(hdr,pa,center)
    linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
    maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(pa)))+abs(abs(y2-y1)/1000.*np.cos(np.radians(pa)))
    maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
    maj_axis =  np.linspace(0,1000*maj_resolution,1000)- (abs((abs(center[0]))*np.sin(np.radians(pa)))+abs(abs(center[1])*np.cos(np.radians(pa))))
    map= []
    loc_max = np.mean(maj_axis[np.where(maj_profile == np.nanmax(maj_profile))[0]])
    if loc_max > 0.:
            pa = pa+180
            print_log(f'''GUESS_ORIENTATION: We have modified the pa by 180 deg as we found the maximum velocity west of the center.
''' , Configuration['OUTPUTLOG'])
    return np.array([pa,pa_error]),np.array([inclination,inclination_error]),SBR_initial,maj_extent,center[0],center[1]
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

    tmp = open(filename, 'r')

    numrings = [int(e.split('=')[1].strip()) for e in tmp.readlines() if e.split('=')[0].strip().upper() == 'NUR']
    tmp.seek(0)
    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
    unarranged = tmp.readlines()
    # Separate the keyword names
    for line in unarranged:
        var_concerned = str(line.split('=')[0].strip().upper())
        if debug:
            print_log(f'''LOAD_TIRIFIC: extracting line
{'':8s}{var_concerned}.
''',None,screen=True, debug = True)
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
                if debug:
                    print_log(f'''LOAD_TIRIFIC: comparing {var_concerned[2:].strip()} to the variables.
{'':8s}Found {varpos}.
''',None,screen=True, debug = True)
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

def obtain_ratios(map, hdr, center, angles, debug = False):
    ratios = []
    max_extent = 0.
    for angle in angles:
        #major axis

        x1,x2,y1,y2 = obtain_border_pix(hdr,angle,center)
        linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
        maj_resolution = abs((abs(x2-x1)/1000.)*np.sin(np.radians(angle)))+abs(abs(y2-y1)/1000.*np.cos(np.radians(angle)))
        maj_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
        maj_axis =  np.linspace(0,1000*maj_resolution,1000)
        tmp = np.where(maj_profile > 0.25*np.nanmax(maj_profile))[0]
        #gauss =fit_gaussian(maj_axis[tmp], maj_profile[tmp])
        #maj_gaus = gaussian_function(maj_profile, *gauss)
        #tmp = np.where(maj_gaus > 0.2*np.nanmax(maj_gaus))[0]
        width_maj = (tmp[-1]-tmp[0])*maj_resolution
        width_maj = np.sqrt(width_maj**2 - (hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])))**2)
        if width_maj > max_extent:
            max_extent = width_maj
        #minor axis
        if angle < 90:
            x1,x2,y1,y2 = obtain_border_pix(hdr,angle+90,center)
        else:
            x1,x2,y1,y2 = obtain_border_pix(hdr,angle-90,center)
        linex,liney = np.linspace(x1,x2,1000), np.linspace(y1,y2,1000)
        min_resolution =abs((abs(x2-x1)/1000.)*np.sin(np.radians(angle+90)))+abs(abs(y2-y1)/1000.*np.cos(np.radians(angle+90)))
        min_axis =  np.linspace(0,1000*min_resolution,1000)
        min_profile = ndimage.map_coordinates(map, np.vstack((liney,linex)),order=1)
        tmp = np.where(min_profile > 0.25*np.nanmax(min_profile))[0]
        #gauss =fit_gaussian(min_axis[tmp], maj_profile[tmp])
        #min_gaus = gaussian_function(min_profile,*gauss)
        #tmp = np.where(min_gaus > 0.2*np.nanmax(min_gaus))[0]
        width_min = (tmp[-1]-tmp[0])*min_resolution
        width_min = np.sqrt(width_min**2 - (hdr['BMAJ']/(np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])))**2)
        if width_min > max_extent:
            max_extent = width_min
        ratios.append(width_maj/width_min)
    #as the extend is at 25% let's take 2 time the sigma of that
    max_extent = (max_extent/(2.*np.sqrt(2*np.log(2))))*2.
    return np.array(ratios,dtype=float), max_extent*np.mean([abs(hdr['CDELT1']),abs(hdr['CDELT2'])])
obtain_ratios.__doc__ = '''
;+
; NAME:
;       obtain_ratios(map, hdr, center, angles)
;
; PURPOSE:
;       Obtain the ratio between an axis along the specified angle as well as one rotated 90 degrees to determine the pa and inclination.
;       Additionally keep track of the maximum width.
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
;          ratios =  ratios corresponding to the input angles
;         maj_width = the maximum extend of any profile analysed (in degree)
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      np.mean, ndimage.map_coordinates, np.linspace, np.vstack, np.cos, np.radians, np.sin
;
; EXAMPLE:
;
;
'''





def obtain_border_pix(hdr,angle,center, debug = False):

    if angle < 90.:
        x1 = center[0]-(hdr['NAXIS2']-center[1])*np.tan(np.radians(angle))
        x2 = center[0]+(center[1])*np.tan(np.radians(angle))
        if x1 < 0:
            x1 = 0
            y1 = center[1]+(center[0])*np.tan(np.radians(90-angle))
        else:
            y1 = hdr['NAXIS2']
        if x2 > hdr['NAXIS1']:
            x2 = hdr['NAXIS1']
            y2 = center[1]-(center[0])*np.tan(np.radians(90-angle))
        else:
            y2 = 0
    elif angle == 90:
        x1 = 0 ; y1 = center[1] ; x2 = hdr['NAXIS1'] ; y2 = center[1]
    else:
        x1 = center[0]-(center[1])*np.tan(np.radians(180.-angle))
        x2 = center[0]+(hdr['NAXIS2']-center[1])*np.tan(np.radians(180-angle))
        if x1 < 0:
            x1 = 0
            y1 = center[1]-(center[0])*np.tan(np.radians(angle-90))
        else:
            y1 = 0
        if x2 > hdr['NAXIS1']:
            x2 = hdr['NAXIS1']
            y2 = center[1]+(center[0])*np.tan(np.radians(angle-90))
        else:
            y2 = hdr['NAXIS2']

    return x1,x2,y1,y2

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
                if np.where(diff < 2)[0].size:
                    edge = True
                if debug:
                    print_log(f'''SOFIA_CATALOGUE: Edge = {edge}
    ''',Configuration['OUTPUTLOG'],debug= debug,screen = True)
                if edge:
                    many_sources[Variables.index('f_sum')][i]=0.
            if np.nansum(many_sources[Variables.index('f_sum')]) == 0.:
                    beam_edge = beam_edge/2.
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




def tirific_template(filename, debug = False):
    tmp = open(filename, 'r')
    result = Proper_Dictionary()
    counter = 0
    # Separate the keyword names
    for line in tmp.readlines():
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

def sofia_template(filename, debug = False):
    tmp = open(filename, 'r')
    result = {}
    counter = 0
    counter2 = 0
    # Separate the keyword names
    for line in tmp.readlines():
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
