#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in FAT to read input files
import os
from support_functions import Proper_Dictionary

#Function to read a FAT input Catalogue

def catalogue(filename):
    Catalogue = Proper_Dictionary({})
    tmpfile = open(filename,'r')
    #Define the exsiting catalogue input()
    input_columns = [x.strip().upper() for x in tmpfile.readline().split('|')]
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
def config_file(input_parameters):
    tmpfile = open(input_parameters.configfile, 'r')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','WARP_OUTPUT']
    string_keys = ['OUTPUTLOG', 'OUTPUTCATALOGUE','MAINDIR','CATALOGUE']
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
            else:
                Configuration[add_key] = float(tmp.split('=', 1)[1].strip())
    #if we are checking the installation then the maindir, outputcatalogue and
    #Log go into the original_dir+installation_check.
    if input_parameters.installation_check:
        Configuration['CATALOGUE'] = input_parameters.original_dir+'/Installation_Check/FAT_Input_Catalogue.txt'
        Configuration['MAINDIR'] = input_parameters.original_dir+'/Installation_Check/'
        Configuration['OUTPUTCATALOGUE'] =input_parameters.original_dir+'/Installation_Check/Output_N2903.txt'
    #Make the input idiot satified
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
    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['MAPS_OUTPUT'] == 5:
        Configuration['MAPS_OUTPUT'] = 4
    return Configuration
