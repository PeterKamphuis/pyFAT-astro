# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that written for developments purposes but are not included in the main part of the code.
# If you are too lazy to properly document please keep the functions in this file until you do document.

from pyFAT.Support.support_functions import print_log, convert_type,set_limits,rename_fit_products,\
                              set_ring_size,calc_rings,get_usage_statistics,get_inner_fix,convertskyangle,\
                              finish_current_run, remove_inhomogeneities,get_from_template,set_format, \
                              set_rings, convertRADEC, is_available
from pyFAT.Support.clean_functions import clean_before_sofia,clean_after_sofia
from pyFAT.Support.fits_functions import cut_cubes,extract_pv,make_moments
from pyFAT.Support.modify_template import write_new_to_template,smooth_profile,set_cflux,fix_sbr, \
                                          regularise_profile,set_fitting_parameters,check_size, \
                                          no_declining_vrot, set_new_size,set_errors,get_warp_slope

#from pyFAT.Support.run_functions import run_tirific

from pyFAT.Support.constants import H_0
from astropy.wcs import WCS
from astropy.io import fits
import pyFAT.Support.read_functions as rf
import pyFAT.Support.write_functions as wf
import pyFAT.Support.run_functions as runf
import os
import sys
import time
import copy
import subprocess
import numpy as np
import traceback
import warnings


'''
This is code snippet that belongs in main for when TWO_STEP = True However it is not developed at the moment as flat disks can be fitted by settinf FIX_INCLINATION and FIX_PA
multiline prints are excepted like  \'''
    if Configuration['START_POINT'] < 4:
        #We first fixe the variations
        Configuration['FIX_INCLINATION'][0] = True
        Configuration['FIX_SDIS'][0] = True
        Configuration['FIX_PA'][0] = True
        Configuration['FIX_Z0'][0] = True
        # setup the first def file to be used in the first loop
        wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                                cube_hdr,Initial_Parameters= Initial_Parameters,fit_stage='Centre_Convergence',debug=Configuration['DEBUG'])
        sf.print_log(f\'''The initial def file is written and we will now start fitting.
\''',Configuration['OUTPUTLOG'], screen =True, debug = Configuration['DEBUG'])
        Configuration['PREP_END_TIME'] = datetime.now()
        current_run = 'Not Initialized'
        # If we have no directory to put the output we create it
        if not os.path.isdir(Configuration['FITTING_DIR']+'Centre_Convergence'):
            os.mkdir(Configuration['FITTING_DIR']+'Centre_Convergence')
        #We skip the first fit atm
        #Configuration['CC_ACCEPTED'] = True
        #write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}Cen_Conv.def", Tirific_Template)
        #Upto here should be removed for real code
        while not Configuration['CC_ACCEPTED'] and Configuration['CC_LOOPS'] < 10:
            Configuration['CC_LOOPS'] = Configuration['CC_LOOPS']+1
            sf.print_log(f\'''We are starting loop {Configuration['CC_LOOPS']} of trying to converge the center.
\''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
            current_run = runf.central_converge(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,Initial_Parameters, debug = Configuration['DEBUG'])

        if Configuration['CC_ACCEPTED']:
            sf.print_log(f\'''The center has converged and we will adjust the smoothed profile and start to adjust the size of the galaxy.
\''',Configuration['OUTPUTLOG'],screen =True,debug=Configuration['DEBUG'])
        else:
            sf.print_log(f\''' We could not find a stable center for the the initial stages. We will now try while adapting the the size of the model.
\''',Configuration['OUTPUTLOG'],screen =True,debug=Configuration['DEBUG'])

        #Then we want to make a smoothed version that can be adapted
        #current_run = runf.fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,stage = 'after_cc', fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])
        #incl = rf.load_tirific(f"{Configuration['FITTING_DIR']}Centre_Convergence/Centre_Convergence.def",Variables = ['INCL'])
        #sf.print_log(f\'''BEFORE_CHECK_INCLINATION: CC_loops = {Configuration['CC_LOOPS']}
#{'':8s} Incl = {incl}
#{'':8s} Size in beams =  {Configuration['SIZE_IN_BEAMS']})
#\''',Configuration['OUTPUTLOG'],debug = Configuration['DEBUG'])
        #if float(incl[0][0]) < 40.:
            #If our fit stage is after cc we want to make sure we do an extra check on low inclinations or small Galaxies
        #    runf.check_inclination(Configuration,Tirific_Template,Fits_Files,fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])

        #if Configuration['OPTIMIZED']:
        #    runf.make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_stage = 'Centre_Convergence',debug=Configuration['DEBUG'])
    else:
        current_run = 'Not initialized'
        write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}Centre_Convergence/Centre_Convergence.def", Tirific_Template, \
                             Variables = ['VROT','Z0', 'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2',
                                             'RADI','INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2', 'CFLUX', 'CFLUX_2', 'NUR', 'CONDISP',\
                                             'BMIN','BMAJ','RMS','BPA','NCORES','INIMODE', 'VARY','VARINDX','MODERATE','DELEND','DELSTART'\
                                             ,'MINDELTA','PARMAX','PARMIN','DISTANCE','INSET'],debug=Configuration['DEBUG'])

    Configuration['CC_END_TIME'] = datetime.now()


    #If we only care about a centrally converged galaxy we stop here

    # if our current run is not broken then we want to stop it
    sf.finish_current_run(Configuration,current_run,debug= Configuration['DEBUG'])
    # write the new values to the basic info file
        #Write the info to the Basic info File

    wf.basicinfo(Configuration,first_fit = True, template=Tirific_Template,Fits_Files=Fits_Files)



    if Configuration['FINISHAFTER'] == 1:
        Configuration['FINAL_COMMENT'] = 'You have chosen to end the fitting after preprocessing and sofia.'
        cf.finish_galaxy(Configuration,maximum_directory_length,current_run =current_run, Fits_Files =Fits_Files,debug=Configuration['DEBUG'])
        continue
    # Now set our variations to the original values but only if the galaxy is large enough
    if Configuration['SIZE_IN_BEAMS'] > Configuration['MINIMUM_WARP_SIZE']:
        Configuration['FIX_INCLINATION'][0] = Original_Configuration['FIX_INCLINATION'][0]
        Configuration['FIX_SDIS'][0] = Original_Configuration['FIX_SDIS'][0]
        Configuration['FIX_PA'][0] = Original_Configuration['FIX_PA'][0]
        Configuration['FIX_Z0'][0] = Original_Configuration['FIX_Z0'][0]


    #Then we want to setup for the next fit.
    wf.initialize_def_file(Configuration, Fits_Files,Tirific_Template, \
                            cube_hdr,fit_stage='Extent_Convergence',debug=Configuration['DEBUG'])
    if not os.path.isdir(Configuration['FITTING_DIR']+'Extent_Convergence'):
        os.mkdir(Configuration['FITTING_DIR']+'Extent_Convergence')

    while not Configuration['EC_ACCEPTED'] and Configuration['EC_LOOPS'] < allowed_loops:
        Configuration['EC_LOOPS'] = Configuration['EC_LOOPS']+1
        sf.print_log(f\'''We are starting loop {Configuration['EC_LOOPS']} of trying to converge the extent.
\''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
        if Configuration['DEBUG']:
                sf.print_log(f\'''Settings for the variations will be.
{'':8s} INCLINATION: Fixed = {Original_Configuration['FIX_INCLINATION']}
{'':8s} PA: Fixed = {Original_Configuration['FIX_PA']}
{'':8s} SDIS: Fixed = {Original_Configuration['FIX_SDIS']}
{'':8s} Z0: Fixed = {Original_Configuration['FIX_Z0']}
\''',Configuration['OUTPUTLOG'],debug=Configuration['DEBUG'],screen =True)
        if Configuration['SIZE_IN_BEAMS'] > Configuration['MINIMUM_WARP_SIZE']:
            Configuration['FIX_INCLINATION'] = Original_Configuration['FIX_INCLINATION']
            Configuration['FIX_SDIS'] = Original_Configuration['FIX_SDIS']
            Configuration['FIX_PA'] = Original_Configuration['FIX_PA']
            Configuration['FIX_Z0'] = Original_Configuration['FIX_Z0']
        else:
            Configuration['FIX_INCLINATION'] = True
            Configuration['FIX_SDIS'] = True
            Configuration['FIX_PA'] = True
            Configuration['FIX_Z0'] = True
            flatten_the_curve(Configuration,Tirific_Template)
        current_run = runf.extent_converge(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,debug = Configuration['DEBUG'],allowed_loops = allowed_loops)


    if Configuration['EC_ACCEPTED']:
        sf.print_log(f\'''The extent has converged and we make a smoothed version.
\''',Configuration['OUTPUTLOG'],screen =True, debug = Configuration['DEBUG'])
        current_run = runf.fit_smoothed_check(Configuration, Fits_Files,Tirific_Template,current_run,cube_hdr,stage = 'after_ec', fit_stage = 'Extent_Convergence',debug = Configuration['DEBUG'])
        Configuration['FINAL_COMMENT'] = 'The galaxy has succesfully been fitted'
        if Configuration['OPTIMIZED']:
            runf.make_full_resolution(Configuration,Tirific_Template,Fits_Files,current_run = current_run,fit_stage = 'Extent_Convergence',debug=Configuration['DEBUG'])
    else:
        Configuration['FINAL_COMMENT'] = 'We could not converge on the extend of the galaxy'
        Configuration['MAPS_OUTPUT'] = 5
        cf.finish_galaxy(Configuration,maximum_directory_length, Fits_Files =Fits_Files,current_run =current_run,debug=Configuration['DEBUG'])
        continue
    Configuration['EC_END_TIME'] = datetime.now()
'''

def central_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr,Initial_Parameters, debug = False):
    #This function will run while it does not return 1, this means when it is called it is not in acceptance
    if debug:
        print_log('''CENTRAL_CONVERGE: Starting.
''',Configuration['OUTPUTLOG'],debug = debug)
    #First we actually run tirific
    accepted,current_run = run_tirific(Configuration,current_run,stage = 'run_cc', fit_stage = 'Centre_Convergence',debug=debug)
    #accepted,current_run = run_parameterized_tirific(Configuration,Tirific_Template,current_run,loops = 10, stage = 'run_cc', fit_stage = 'Centre_Convergence',debug=debug)

    accepted = check_central_convergence(Configuration,Tirific_Template,hdr,accepted, fit_stage = 'Centre_Convergence',debug=debug)
    set_cflux(Configuration,Tirific_Template,debug = debug)
    check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
    if accepted:
        Configuration['CC_ACCEPTED'] = True
    else:
        print_log('''CENTRAL_CONVERGE: Tirific ran the maximum amount of loops which means the fit is not accepted and we retry.
''',Configuration['OUTPUTLOG'],debug = debug)
        smoothed_vrot = smooth_profile(Configuration,Tirific_Template,'VROT',min_error = float(Configuration['CHANNEL_WIDTH']), debug=debug)
        smoothed_sbr = smooth_profile(Configuration,Tirific_Template,'SBR',min_error= np.max([float(Tirific_Template['CFLUX']),float(Tirific_Template['CFLUX_2'])]) ,debug=debug)
        #after smoothing the sbr we should check it
        check_sbr(Configuration,Tirific_Template,hdr,debug = debug)
        no_declining_vrot(Configuration,Tirific_Template,debug = debug)
        for key in ['XPOS','YPOS','VSYS','PA','INCL']:
            Initial_Parameters[key][0] = rf.load_template(Configuration,Tirific_Template,Variables = [key])[0][0]
        Initial_Parameters['VROT'] =  [np.nanmax(smoothed_vrot),set_limits(np.nanstd(smoothed_vrot[1:]),np.nanmax(smoothed_vrot)-np.nanmin(smoothed_vrot[1:]),np.nanmax(smoothed_vrot))]

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_cc',\
                               initial_estimates = Initial_Parameters, debug = debug)
        Tirific_Template['SIZE'] = f"{float(Tirific_Template['SIZE']) * 2.}"
        wf.tirific(Configuration,Tirific_Template,name = 'Centre_Convergence_In.def',debug=debug)

    return current_run



def check_inclination(Configuration,Tirific_Template,Fits_Files, fit_type = 'Undefined', debug =False):
    if debug:
        print_log(f'''CHECK_INCLINATION: estimating whether our inclination estimate is decent.
''',Configuration['OUTPUTLOG'],debug = True)

    current = rf.load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def",
                Variables= ['VROT','INCL','PA','XPOS','YPOS','INCL_ERR'])


    inclination = float(current[1][0])
    if debug:
        print_log(f'''CHECK_INCLINATION: This is the initial inclination {inclination}
''',Configuration['OUTPUTLOG'],debug = True)

    if Configuration['SIZE_IN_BEAMS'] < 5:
        incl_to_check = np.linspace(inclination-15.,inclination+15.,20)
    else:
        incl_to_check = np.linspace(10,50,20)
    # Let's make a directory to put things in
    tmp_stage = 'tmp_incl_check'
    try:
        os.mkdir(f"{Configuration['FITTING_DIR']}/{tmp_stage}")
    except:
        pass
    other_run = Configuration['TIRIFIC_RUNNING']
    Configuration['TIRIFIC_RUNNING'] = False
    #and a copy of the tirific template
    Check_Template = copy.deepcopy(Tirific_Template)
    write_new_to_template(Configuration, f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}.def", Check_Template,debug = debug)
    Check_Template['LOOPS'] = '0'
    Check_Template['INIMODE'] = '0'
    Check_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"
    Check_Template['RESTARTNAME'] = f"Logs/restart_{tmp_stage}.txt"
    out_keys = ['LOGNAME','OUTSET','TIRDEF']
    out_extensions = ['log','fits','def']
    incl_run= 'Not Initialized'
    for i,key in enumerate(out_keys):
        Check_Template[key] = f"{tmp_stage}/{tmp_stage}.{out_extensions[i]}"


    vobs = [x*np.sin(np.radians(y)) for x,y in zip(current[0][:],current[1][:])]
    if debug:
        print_log(f'''CHECK_INCLINATION: These are the values we get as input
{'':8s}Inclination = {current[1][:]}
{'':8s}Vrot = {current[0][:]}
{'':8s}Vobs = {vobs}
''',Configuration['OUTPUTLOG'])
    mom_chi = []
    model_mom0 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT0']}")
    #model_mom0 = remove_inhomogeneities(Configuration,model_mom0,inclination=float(current[1][0]), pa = float(current[2][0]), center = [current[3][0],current[4][0]],debug=debug)
    chan_map = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['CHANNEL_MAP']}")
    noisemap = np.sqrt(chan_map[0].data)*Configuration['NOISE']/np.nanmax(model_mom0[0].data)
    max_in_moment = np.nanmax(model_mom0[0].data)
    model_mom0[0].data = model_mom0[0].data
    #/max_in_moment
    for incl in incl_to_check:
        #print(f'We are doing this inclination {incl}')
        vrot = [x/np.sin(np.radians(incl)) for x in vobs]
        for key in ['INCL','INCL_2']:
            Check_Template[key]= f"{incl:.2f}"
        for key in ['VROT','VROT_2']:
            Check_Template[key]= f"{' '.join([f'{x:.2f}' for x in vrot])}"
        wf.tirific(Configuration,Check_Template,name = f'{tmp_stage}_In.def',debug=debug)
        accepted,incl_run = runf.run_tirific(Configuration,incl_run, stage = 'incl_check', fit_type=tmp_stage,debug=debug)
        make_moments(Configuration,Fits_Files,fit_type=tmp_stage,\
                     moments = [0], \
                     overwrite = True, vel_unit = 'm/s',debug=debug)
        #make_moments(filename = f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}.fits", basename = 'tmp_incl', directory = f"{Configuration['FITTING_DIR']}{tmp_stage}/",\
        #             moments = [0],level = 3.*Configuration['NOISE'], \
        #             overwrite = True, log= Configuration['OUTPUTLOG'], vel_unit = 'm/s',debug=debug)
        incl_mom0 = fits.open(f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}_mom0.fits")
        if debug:
            try:
                os.remove(f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
            except:
                pass
            os.rename(f"{Configuration['FITTING_DIR']}{tmp_stage}/{tmp_stage}_mom0.fits",f"{Configuration['FITTING_DIR']}{tmp_stage}/tmp_{incl:.1f}_mom0.fits")
        chi = np.nansum((model_mom0[0].data[noisemap > 0.]-incl_mom0[0].data[noisemap > 0.])**2/noisemap[noisemap > 0.]**2)
        mom_chi.append(abs(chi))
        incl_mom0.close()
    chan_map.close()
    model_mom0.close()
    finish_current_run(Configuration, incl_run)
    Configuration['TIRIFIC_RUNNING'] = other_run
    low= np.where(mom_chi == np.nanmin(mom_chi))[0]
    new_incl = float(incl_to_check[low])
    if debug:
        print_log(f'''CHECK_INCLINATION: This is the new inclination {new_incl} it was {current[1][0]}.
{'':8s}mom_chi = {mom_chi}
{'':8s}low = {low}
''',Configuration['OUTPUTLOG'])


    #exit()
    incl_err = np.mean(current[5])
    if incl_err < 5.:
        incl_err = 5.
    if not current[1][0]-incl_err < new_incl < current[1][0]+incl_err:
        if debug:
            print_log(f'''CHECK_INCLINATION: The inclination has changed, writing to file.
''',Configuration['OUTPUTLOG'])
        for key in ['INCL','INCL_2']:
                Check_Template[key]= f"{new_incl:.2f}"
        vrot = [x/np.sin(np.radians(new_incl)) for x in vobs]
        for key in ['VROT','VROT_2']:
                Check_Template[key]= f"{' '.join([f'{x:.2f}' for x in vrot])}"
        Tirific_Template = copy.deepcopy(Check_Template)
        Check_Template = []
        wf.tirific(Configuration,Tirific_Template,name = f"{fit_type}/{fit_type}.def",debug=debug)
check_inclination.__doc__ =f'''
 NAME:
    check_inclination

 PURPOSE:
    For low inclinations do a better check on the inclinations. Comparing the Chi^2 of the models over different inclinations

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:
    debug = False

    fit_type = 'Undefined'
    type of fitting

 OUTPUTS:
    A new file with modified inclination if it has decreased.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE: This should efinitely be further explored.
'''
#Function to read FAT configuration file into a dictionary
def config_file(input_parameters, start_dir, debug = False):
    No_File = True
    log_write_config = 'Empty'
    while No_File:
        try:
            if input_parameters.configfile == 'ChecK.ConfiG':
                from pyFAT import Installation_Check as IC
                with import_res.open_text(IC,'FAT_INPUT.config') as tmp:
                    tmpfile = tmp.readlines()
            elif input_parameters.configfile == 'No Default' and input_parameters.single_cube != 'CataloguE' :   #single_cube should be set, so can not be CataloguE
                import pyFAT
                with import_res.open_text(pyFAT,'FAT_INPUT.config') as tmp:
                    tmpfile = tmp.readlines()
            else:
                with open(input_parameters.configfile, 'r') as tmp:
                    tmpfile = tmp.readlines()
            No_File = False
        except:
            print(traceback.print_exc())
            input_parameters.configfile = input(f'''
                        You have provided a config file but it can't be found.
                        If you want to provide a config file please give the correct name.
                        Else press CTRL-C to abort.
                ''')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','FIX_SBR','FIX_VROT','WARP_OUTPUT']
    string_keys = ['OUTPUTLOG', 'OUTPUTCATALOGUE','MAINDIR','CATALOGUE']
    integer_keys = ['STARTGALAXY','ENDGALAXY','MAPS_OUTPUT','OPT_PIXELBEAM','FINISHAFTER','FITTING_TYPE']
    # Separate the keyword names
    for tmp in tmpfile:
        if tmp[0] != '#':
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
            add_key_in = tmp.split('=', 1)
            if len(add_key_in) > 1:
                add_key = add_key_in[0].strip().upper()
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
                            inp = input(f"The parameter {add_key} in the configuration file  must be true/false or yes/no. Please give the correct value. \n".format(add_key))
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
        if Configuration['CATALOGUE'] != f"Installation_Check/FAT_Input_Catalogue.txt" or \
           Configuration['MAIN_DIRECTORY'] != f"Installation_Check/" or \
           Configuration['OUTPUTCATALOGUE'] != f"Installation_Check/Output_N2903.txt":
           raise BadConfigurationError(f"You can not modify the Installation Check input. It is solely for checking the Installation. Aborting")
        test_files = ['FAT_Input_Catalogue.txt','NGC_2903.fits']
        test_dir = f"{start_dir}/FAT_Installation_Check/"
        if not os.path.isdir(test_dir):
            os.mkdir(test_dir)
        else:
            for file in test_files:
                try:
                    os.remove(test_dir+file)
                except:
                    pass
        my_resources = import_res.files('pyFAT.Installation_Check')
        for file in test_files:
            data = (my_resources / file).read_bytes()
            with open(test_dir+file,'w+b') as tmp:
                tmp.write(data)
        Configuration['CATALOGUE'] = f"{test_dir}FAT_Input_Catalogue.txt"
        Configuration['MAIN_DIRECTORY'] = test_dir
        Configuration['OUTPUTCATALOGUE'] =f"{test_dir}Output_N2903.txt"
    elif input_parameters.single_cube != 'CataloguE':
        from pyFAT.Support.write_functions import write_config
        file_location = input_parameters.single_cube.split('/')
        single_dir = f"{start_dir}{'/'.join(file_location[:-1])}/"
        Configuration['MAIN_DIRECTORY'] = single_dir
        Configuration['OUTPUTCATALOGUE'] = None
        Configuration['CATALOGUE'] = None
        if Configuration['OUTPUTLOG'] == "Logfileforthepergalaxyfit.txt":
            Configuration['OUTPUTLOG'] = f"Log_{os.path.splitext(file_location[-1])[0]}.txt"
        log_write_config = write_config(f"{single_dir}FAT_INPUT_{os.path.splitext(file_location[-1])[0]}.config", Configuration,debug=debug)

    #Make the input idiot safe
    if Configuration['MAIN_DIRECTORY'][-1] != '/':
        Configuration['MAIN_DIRECTORY'] = f"{Configuration['MAIN_DIRECTORY']}/"

    while not os.path.isdir(Configuration['MAIN_DIRECTORY']):
        Configuration['MAIN_DIRECTORY'] = input(f'''
                    Your main fitting directory ({Configuration['MAIN_DIRECTORY']}) does not exist.
                    Please provide the correct directory.
                    ''')
    if Configuration['CATALOGUE']:
        while not os.path.exists(Configuration['CATALOGUE']):
            Configuration['CATALOGUE'] = input(f'''
                        Your input catalogue ({Configuration['CATALOGUE']}) does not exist.
                        Please provide the correct file name.
                        ''')
    #The output catalogue only needs to be in a valid directory as we create it
    if Configuration['OUTPUTCATALOGUE']:
        output_catalogue_dir = Configuration['OUTPUTCATALOGUE'].split('/')
        if len(output_catalogue_dir) > 1:
            check_dir = '/'.join(output_catalogue_dir[:-1])
            while not os.path.isdir(check_dir):
                check_dir= input(f'''
                        The directory for your output catalogue ({Configuration['OUTPUTCATALOGUE']}) does not exist.
                        Please provide the correct directory name.
                        ''')
                Configuration['OUTPUTCATALOGUE'] = f"{check_dir}/{output_catalogue_dir[-1]}"


    required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','FIX_SBR','FIX_VROT','HANNING',\
                                   'STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT',\
                                   'RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR',\
                                    'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM',\
                                     'MAPS_OUTPUT','WARP_OUTPUT','FITTING_TYPE']

    for key in required_configuration_keys:
        if key not in Configuration:
            if key == 'STARTGALAXY':
                Configuration[key] = -1
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
            elif key == 'FITTING_TYPE':
                Configuration[key] = 'Fit_Tirific_OSC'
            elif key == 'FIX_INCLINATION': #Previosuly called fix_incl
                Configuration[key] = False
            elif key == 'FIX_PA':
                Configuration[key] = False
            elif key == 'FIX_SDIS':
                Configuration[key] = False
            elif key == 'FIX_SBR':
                Configuration[key] = False
            elif key == 'FIX_VROT':
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
                raise BadConfigurationError(f"Something has gone wrong reading the required config key. This should never ever happen. Please file an issue on github")
    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['OUTPUT_QUANTITY'] == 5:
        Configuration['OUTPUT_QUANTITY'] = 4
    # We double the fix keys so we can  modify one while keeping the original as well
    fix_keys = ['FIX_PA','FIX_INCLINATION','FIX_SDIS','FIX_Z0','FIX_SBR']
    for key in fix_keys:
        Configuration[key] = [Configuration[key],Configuration[key]]
    return Configuration, log_write_config
config_file.__doc__ =f'''
 NAME:
    config_file
 PURPOSE:
    Read the FAT config file and write into the a dictionary
 CATEGORY:
    read_functions

 INPUTS:
    input_parameters = input parameters for run
    start_dir = the directory where FAT is started from

 OPTIONAL INPUTS:
    debug = False
    fit_type = 'Undefined'

 OUTPUTS:
    Configuration = dictionary with the config file input

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def extent_converge(Configuration, Fits_Files,Tirific_Template,current_run,hdr, debug = False,allowed_loops = 10.):
    if debug:
        print_log(f'''EXTENT_CONVERGENCE: Starting with loop {Configuration['EC_LOOPS']} out of maximum {allowed_loops}.
''',Configuration['OUTPUTLOG'],debug = debug)
    fit_stage = 'Extent_Convergence'
    stage = 'run_ec'
    accepted,current_run = run_tirific(Configuration,current_run,stage = stage, fit_stage = fit_stage, debug= debug)
    write_new_to_template(Configuration,f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def" , Tirific_Template, debug = debug)
    accepted_size = check_size(Configuration,Tirific_Template, fit_stage = fit_stage, stage = stage,current_run = current_run, debug=debug,Fits_Files=Fits_Files)
    check_sbr(Configuration,Tirific_Template,hdr,stage = 'run_ec',debug=debug)
    if accepted and accepted_size:
            Configuration['EC_ACCEPTED'] = True
    else:
        if not accepted:
            print_log('''EXTENT_CONVERGE: Tirific ran the maximum amount of loops which means the fit is not accepted and we smooth and retry.
''',Configuration['OUTPUTLOG'],debug = debug)
        else:
            print_log(f'''EXTENT_CONVERGE: FAT adjusted the rings. Refitting with new settings after smoothing them.
''',Configuration['OUTPUTLOG'],debug = debug)
        Configuration['EC_ACCEPTED'] = False
        if Configuration['EC_LOOPS'] > allowed_loops:
            print_log(f'''EXTENT_CONVERGE: We have ran the convergence more than {allowed_loops} times aborting the fit.
''',Configuration['OUTPUTLOG'],debug = debug)
            return current_run

        Configuration['INNER_FIX'] = get_inner_fix(Configuration, Tirific_Template,debug=debug)
        set_cflux(Configuration,Tirific_Template,debug = debug)
        keys_to_smooth =['INCL','PA','SDIS','Z0','VROT']
        min_errors = [3.*np.mean(Configuration['LIMIT_MODIFIER']),2.,hdr['CDELT3']/(2000.*Configuration['LIMIT_MODIFIER']), \
                        convertskyangle(0.1,Configuration['DISTANCE'],physical= True)/Configuration['LIMIT_MODIFIER'],\
                        hdr['CDELT3']/(1000.*Configuration['LIMIT_MODIFIER'])]
        for j,key in enumerate(keys_to_smooth):
            smoothed = smooth_profile(Configuration,Tirific_Template,key,debug=debug,min_error=min_errors[j])

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'run_ec',\
                             debug = debug)

        #After smoothing the SBR we should check it again?
        check_sbr(Configuration,Tirific_Template,hdr,stage = 'run_ec',debug = debug)
        wf.tirific(Configuration,Tirific_Template,name = 'Extent_Convergence_In.def',debug = debug)
    return current_run

def fake_plastic_trees(Configuration,Tirific_Template,current_run,key,other_key,\
                        initial_profile,max_bend,step_bend,black_star = 0., \
                        counter = 0, total_fits = 1.,fit_stage= 'unknown',debug=False):
    prev_bend = -1*max_bend
    Chi = []
    profiles = []
    format = set_format(key)
    while prev_bend < max_bend:
        #first side 1
        bend = float(prev_bend+step_bend)
        new_profile = the_bends(Configuration,initial_profile[1],initial_profile[0],black_star,bend,debug=debug)
        diff = float(initial_profile[2][0]-new_profile[0])
        Tirific_Template[key] =  f"{' '.join([f'{x:{format}}' for x in new_profile])}"
        Tirific_Template[other_key] =  f"{' '.join([f'{x-diff:{format}}' for x in initial_profile[2]])}"
        wf.tirific(Configuration,Tirific_Template,name =  f"{fit_stage}_In.def",debug = False)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write(f"Restarting for bend = {bend} and RI = {Tirific_Template['RESTARTID']} \n")
        Chi_here = wait_for_tirific(Configuration,current_run, counter= counter, total_fits= total_fits)
        if ~np.isnan(Chi_here):
            Chi.append(Chi_here)
            profiles.append(rf.load_tirific(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def",Variables = [key]))
            prev_bend = copy.deepcopy(bend)
        else:
            prev_bend = copy.deepcopy(bend-step_bend)
        counter += 1
    return Chi,profiles,counter

def fit_warp(Configuration,Tirific_Template,current_run, fit_stage = 'None',debug=False):
    if debug:
        print_log(f'''FIT_WARP: Starting a new paremeterized run in fit_stage {fit_stage}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)

    Tirific_Template['LOOPS'] = 0
    incl_bend = planet_telex(Configuration, Tirific_Template,current_run,'INCL',fit_stage=fit_stage,debug=debug)
    incl_bend_2 = planet_telex(Configuration, Tirific_Template,current_run,'INCL_2',fit_stage=fit_stage,debug=debug)
    pa_bend = planet_telex(Configuration, Tirific_Template,current_run,'PA',max_bend = 50, fit_stage=fit_stage,debug=debug)
    pa_bend_2 = planet_telex(Configuration, Tirific_Template,current_run,'PA_2',max_bend = 50, fit_stage=fit_stage,debug=debug)
    if debug:
        print(f"PA = {Tirific_Template['PA']}" )
        print(f"PA_2 = {Tirific_Template['PA_2']}" )
        print(f"INCL = {Tirific_Template['INCL']}" )
        print(f"INCL_2 = {Tirific_Template['INCL_2']}" )

    return current_run,[Tirific_Template['INCL'],Tirific_Template['INCL_2'],Tirific_Template['PA'],Tirific_Template['PA_2']]


def planet_telex(Configuration, Tirific_Template,current_run, key, black_stars = 0.,max_bend =30 ,step_bend = 5.,fit_stage= 'None', debug=False):
    Chi = []
    thebends = []
    if key[-2:] == '_2':
        SBR_key = 'SBR_2'
        other_key = key[:-2]
        fitkey = other_key
    else:
        SBR_key = 'SBR'
        other_key = f"{key}_2"
        fitkey = key
    initial_profile = np.array(get_from_template(Configuration,Tirific_Template,['RADI',key,other_key]))
    if black_stars == 0.:
        black_stars = initial_profile[0][-1]/np.linspace(3,1.1,10)
        black_index = np.where(initial_profile[0] <= black_stars[0])[0][-1]

    # The start outermost bent of inclination
    SBRformat = set_format(SBR_key)
    format = set_format(key)
    SBR = np.array(get_from_template(Configuration,Tirific_Template,[SBR_key])[0],dtype=float)
    new_SBR = []
    for i,x in enumerate(SBR):
        if i < black_index:
            new_SBR.append(x)
        else:
            new_SBR.append(SBR[black_index])

    Tirific_Template[SBR_key] = f"{' '.join([f'{x:{SBRformat}}' for x in new_SBR ])}"
    Tirific_Template['LOOPS'] = 1.
    parameters = {key : [np.mean(initial_profile[1:2]), 10.]}
    set_fitting_parameters(Configuration, Tirific_Template,stage = 'parameterized',\
                            parameters_to_adjust = [fitkey],initial_estimates = parameters, debug = debug)
    Chi = []
    all_profs = []
    total_fits = len(black_stars)*2*max_bend/step_bend
    counter = 0.
    for black_star in black_stars:
        Chis,profiles,counter = fake_plastic_trees(Configuration,Tirific_Template,current_run,\
                                            key,other_key,initial_profile,max_bend,step_bend,\
                                            black_star= black_star,fit_stage= fit_stage, \
                                            counter = counter, total_fits = total_fits, debug=debug)
        Chi.append(Chis)
        all_profs.append(profiles)
    Chi = np.array(Chi,dtype=float)
    all_profs = np.array(all_profs,dtype=float)
    print(f"Chi = {Chi}")
    #print(np.where(Chi == np.min(Chi))[0][0],np.min(Chi))
    high_and_dry = np.where(Chi == np.min(Chi))
    print(f"This is the bend value we find {high_and_dry[0]}, {high_and_dry[1]}for {key}")
    new_profile = all_profs[high_and_dry][0][0]
    print(f" This is the profile {new_profile}")
    #Restore the brightness profiles
    Tirific_Template[SBR_key] = f"{' '.join([f'{x:{SBRformat}}' for x in SBR])}"
    diff = float(initial_profile[2][0]-new_profile[0])
    Tirific_Template[key] =  f"{' '.join([f'{x:{format}}' for x in new_profile])}"
    Tirific_Template[other_key] =  f"{' '.join([f'{x-diff:{format}}' for x in initial_profile[2]])}"

    return new_profile

def restore_bends(Configuration,Tirific_Template,bends,debug=False):
    if debug:
        print_log(f'''RESTORE_BENDS: setting the best fit bends to be output.
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    bend_keys = ['INCL','INCL_2','PA','PA_2']
    for bend,key in zip(bends,bend_keys):
        Tirific_Template[key] = bend


def run_parameterized_tirific(Configuration, Tirific_Template, current_run, loops = 5, stage = 'initial',fit_stage = 'Undefined_Stage', debug = False):
    if debug:
        print_log(f'''RUN_PARAMETIRIZED_TIRIFIC: Starting a new paremetirized run in stage {stage} and fit_stage {fit_stage}
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
    # First scale the inclination and the vrot

    Tirific_Template['LOOPS']  = 1
    wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = debug)
    if Configuration['TIRIFIC_RUNNING']:
        print_log('''RUN_TIRIFIC: We are using an initialized tirific
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
    else:
        print_log('''RUN_TIRIFIC: We are starting a new TiRiFiC
''',Configuration['OUTPUTLOG'], screen = True,debug = debug)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'w') as file:
            file.write("Initialized a new run")
        #current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_stage}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
        #                            stderr = subprocess.PIPE,cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        current_run = subprocess.Popen(["tirific",f"DEFFILE={fit_stage}_In.def","ACTION = 1"], stdout = subprocess.PIPE, \
                                    cwd=Configuration['FITTING_DIR'],universal_newlines = True)
        Configuration['TIRIFIC_RUNNING'] = True
        Configuration['TIRIFIC_PID'] = current_run.pid
    #rename_fit_products(Configuration,fit_stage = fit_stage, stage=stage, debug = debug)
    # Then if already running change restart file
    sys.stdout.flush()
    currentloop =1
    max_loop = 0
    counter = 0
    Chi_Current = 0.
    if Configuration['TIMING']:
        time.sleep(0.1)
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write(f"# Started Tirific at stage = {fit_stage}\n")
            CPU,mem = get_usage_statistics(Configuration,current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
    print(f"RUN_TIRIFIC: Starting loop 1")
    triggered = False

    for tir_out_line in current_run.stdout:
        tmp = re.split(r"[/: ]+",tir_out_line.strip())
        counter += 1
        if (counter % 50) == 0:
            if Configuration['TIMING']:
                with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
                    if tmp[0] == 'L' and not triggered:
                        if tmp[1] == '1':
                            file.write("# Started the actual fitting \n")
                            triggered = True
                    CPU,mem = get_usage_statistics(Configuration,current_run.pid)
                    file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
        if tmp[0] == 'L':
            if int(tmp[1]) != currentloop:
                print(f"RUN_TIRIFIC: Starting loop {tmp[1]} out of a maximum {max_loop}")
            currentloop  = int(tmp[1])
            if max_loop == 0:
                max_loop = int(tmp[2])
            try:
                Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
            except:
                #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                pass
        if tmp[0].strip() == 'Finished':
            Chi_Current = (float(tmp[tmp.index('C')+1]))
            break
        if tmp[0].strip() == 'Abort':
            break
    time.sleep(1.)
    # Now do the actual fitting loops
    for i in range(loops):
        print(f"RUN_Parameterized TIRIFIC: Starting loop {i+1}")
        resid= Tirific_Template['RESTARTID']
        resname =  Tirific_Template['RESTARTNAME']
        Tirific_Template = rf.tirific_template(f"{Configuration['FITTING_DIR']}{Tirific_Template['TIRDEF']}")
        if i == 0.:
            prev_bends = [Tirific_Template['INCL'],Tirific_Template['INCL_2'],Tirific_Template['PA'],Tirific_Template['PA_2']]
        Tirific_Template['RESTARTID'] = resid
        Tirific_Template['RESTARTNAME'] = resname
    # whatever has been fitted we will now attempt to modify by fitting a warped __version__
        current_run,bends = fit_warp(Configuration,Tirific_Template,current_run, fit_stage=fit_stage,debug=debug)
        Tirific_Template['LOOPS'] = 1
        wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = False)
        with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
            file.write("Restarting from previous run")
        for tir_out_line in current_run.stdout:
            tmp = re.split(r"[/: ]+",tir_out_line.strip())
            counter += 1
            if (counter % 50) == 0:
                if Configuration['TIMING']:
                    with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
                        if tmp[0] == 'L' and not triggered:
                            if tmp[1] == '1':
                                file.write("# Started the actual fitting \n")
                                triggered = True
                        CPU,mem = get_usage_statistics(Configuration,current_run.pid)
                        file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Mb \n")
            if tmp[0] == 'L':
                if int(tmp[1]) != currentloop:
                    print(f"RUN_TIRIFIC: Starting loop {tmp[1]} out of a maximum {max_loop}")
                currentloop  = int(tmp[1])
                if max_loop == 0:
                    max_loop = int(tmp[2])
                try:
                    Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
                except:
                    #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                    pass
            if tmp[0].strip() == 'Finished':
                Chi = float(tmp[tmp.index('C')+1])
                break
            if tmp[0].strip() == 'Abort':
                break
        comp_loop = i+1
        print(f"Chi = {Chi}, Chi_Current = {Chi_Current}, loop = {comp_loop}" )
        if Chi < Chi_Current:
            Chi_Current = Chi
            prev_bends = bends
        else:
            break
    restore_bends(Configuration,Tirific_Template,prev_bends,debug=debug)
    Tirific_Template['LOOPS'] = 0
    wf.tirific(Configuration,Tirific_Template,name = f"{fit_stage}_In.def",debug = False)
    with open(f"{Configuration['FITTING_DIR']}Logs/restart_{fit_stage}.txt",'a') as file:
        file.write("Restarting from previous run")
    if Configuration['TIMING']:
        with open(f"{Configuration['FITTING_DIR']}Logs/Usage_Statistics.txt",'a') as file:
            file.write("# Finished this run \n")
            CPU,mem = get_usage_statistics(Configuration,current_run.pid)
            file.write(f"{datetime.now()} CPU = {CPU} % Mem = {mem} Gb\n")
    print(f"RUN_TIRIFIC: Finished the current tirific run.")
    time.sleep(1.0)
    #The break off goes faster sometimes than the writing of the file so let's make sure it is present
    wait_counter = 0
    while not os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.def") and wait_counter < 1000.:
        time.sleep(0.1)
        wait_counter += 1

    if comp_loop != loops:
        return 1,current_run
    else:
        return 0,current_run

run_parameterized_tirific.__doc__= '''

; NAME:
;       tirific(Configuration)
;
; PURPOSE:
;       Fit a parametirized functions in tirific
;
; CATEGORY:
;       run function
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
;      set_limits, print_log
;
; EXAMPLE:
;
;
'''


def the_bends(Configuration,profile,radii,start,outer_diff,debug=False):
    profile = np.array(profile,dtype=float)
    mean = np.mean(profile[0])   # This is necessary to keep both sides conneccted and the same
    radii = np.array(radii,dtype=float)
    length = float((radii[-1] -start)/2.)
    amplitude =  float(outer_diff)
    center = float(radii[-1])
    mod_profile = np.arctan((radii-center)/(abs(length)))/(0.5*np.pi)
    mod_profile = mod_profile/np.nanmax(abs(mod_profile))
    add = mean-np.mean(mod_profile-mod_profile[0])*amplitude
    mod_profile = (mod_profile-mod_profile[0])*amplitude+mean
    return mod_profile

# Function to wait for tirific when it is too fast for stdout
def wait_for_tirific(Configuration,req_stamp,current_run,counter = 0. ,total_fits = 1.,debug =False):
    sys.stdout.flush()
    Chi = 0.
    attempt = 0.
    currentloop = 0.
    max_loop = 0
    time.sleep(1.0)
    for tir_out_line in current_run.stdout:
        tmp = re.split(r"[/: ]+",tir_out_line.strip())
        if tmp[0] == 'L':
            #if int(tmp[1]) != currentloop:
                #print(f"RUN_TIRIFIC: Starting loop {tmp[1]} out of a maximum {max_loop}")
            currentloop = int(tmp[1])
            if max_loop == 0:
                max_loop = int(tmp[2])
            try:
                Configuration['NO_POINTSOURCES'] = np.array([tmp[18],tmp[19]],dtype=float)
            except:
                #If this fails for some reason an old number suffices, if the code really crashed problems will occur elsewhere.
                pass
        if tmp[0].strip() == 'Finished':
            Chi = float(tmp[tmp.index('C')+1])
            break
        if tmp[0].strip() == 'Abort':
            break
        attempt += 1
        try:
            print(f"\rWarp Loop {counter/total_fits*100.:6.2f}% Complete.", end=" ",flush = True)
        except:
            pass
        if attempt > 1000.:
            Chi = float('NaN')
            break
    current_run.stdout.flush()
    return Chi

wait_for_tirific.__doc__ =f'''
  NAME:
    wait_for_tirific

  PURPOSE:
    Make the code wait for tirific to finish its fitting before it continues processing the fit.
    !!!!!!This is part of experimental code. Not functional currently!!!!!!

  CATEGORY:
     development_functions

  INPUTS:
     Configuration = Standard FAT configuration
     current_run = subproccess structure that is currently identified with tirific

  OPTIONAL INPUTS:
     debug = False

     counter
     Value for the current loop

     total_fits
     Value for the total number of fits to be done

  OUTPUTS:
     Chi
     Chi square calculated by tirific

  OPTIONAL OUTPUTS:

  PROCEDURES CALLED:
     Unspecified
 '''
