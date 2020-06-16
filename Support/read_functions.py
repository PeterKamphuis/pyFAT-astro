#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in FAT to read input files
import os
from support_functions import Proper_Dictionary,print_log
import copy
import numpy as np
#Function to read a FAT input Catalogue
class BadCatalogueError(Exception):
    pass
def catalogue(filename):
    Catalogue = Proper_Dictionary({})
    tmpfile = open(filename,'r')
    #Define the exsiting catalogue input()
    input_columns = [x.strip().upper() for x in tmpfile.readline().split('|')]
    Catalogue['ENTRIES'] = ['COLUMNS']
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
    return Catalogue

#Function to read FAT configuration file into a dictionary
def config_file(input_parameters, start_dir):
    tmpfile = open(input_parameters.configfile, 'r')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','WARP_OUTPUT']
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


    required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','HANNING','STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT','RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR', 'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM', 'MAPS_OUTPUT','WARP_OUTPUT']

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

# function to read the sofia catalogue
def sofia_catalogue(Configuration, Variables =['name','x','x_min','x_max','y','y_min','y_max','z','z_min','z_max','ra','dec','v_app','f_sum','kin_pa'], header = None ):
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
                            print_log(log_statement,Configuration['OUTPUTLOG'])
                            raise BadCatalogueError("READ_SOFIA_CATALOGUE: The required columns could not be found in the sofia catalogue.")
            else:
                for col in Variables:
                    if input_columns.index(col) == 0:
                        start = 0
                    else:
                        start = column_locations[input_columns.index(col)-1]
                    end = column_locations[input_columns.index(col)]
                    outlist[Variables.index(col)].append(line[start:end].strip())
    # we want to fit a specific source
    if len(outlist[0]) > 1 and 'f_sum' in Variables and header:
        many_sources  = copy.deepcopy(outlist)
        # We want to exclude any edge sources
        for i in range(len(many_sources[0])):
            edge = False
            diff = np.array([many_sources[Variables.index('x_min')][i],
                    abs(float(many_sources[Variables.index('x_max')][i])-header['NAXIS1']),
                    many_sources[Variables.index('y_min')][i],
                    abs(float(many_sources[Variables.index('y_max')][i])-header['NAXIS2'])
                    ],dtype=float)

            index = np.where(diff < 2*header['BMAJ']/((abs(header['CDELT1'])+abs(header['CDELT2']))/2.))[0]
            if np.where(diff < 2*header['BMAJ']/((abs(header['CDELT1'])+abs(header['CDELT2']))/2.))[0].size:
                edge = True
            diff = np.array([many_sources[Variables.index('z_min')][i],
                    abs(float(many_sources[Variables.index('z_max')][i])-float(header['NAXIS3']))],dtype=float)
            if np.where(diff < 2)[0].size:
                edge = True
            if edge:
                many_sources[Variables.index('f_sum')][i]=0.


        fluxes = np.array(many_sources[Variables.index('f_sum')],dtype= float)
        outlist = []
        #We want the source with the most total flux.
        index = np.where(np.nanmax(fluxes) == fluxes)[0][0]
        outlist = [x[index] for x in many_sources]
    else:
        outlist = [x[0] for x in outlist]
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




def tirific_template(filename):
    tmp = open(filename, 'r')
    result = {}
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

def sofia_template(filename):
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
