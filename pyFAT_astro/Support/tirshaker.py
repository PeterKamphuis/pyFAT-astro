# -*- coding: future_fstrings -*-
"""
Bootstrap errors in tirific, Written by G.I.G Josza.
"""
import matplotlib.ticker as ticker
import sys
import io
import random
import copy
import os
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.units as u
import astropy.constants as c
import scipy.special as spec
from scipy.optimize import least_squares
from matplotlib import colors as mplcolors
from scipy import stats

from pyFAT_astro.Support.write_functions import tirific as write_tirific
import pyFAT_astro.Support.support_functions as sf
import pyFAT_astro
from datetime import datetime
#
#tirshaker
def tirshaker(Configuration, Tirific_Template, outfilename = 'test_out.def', \
              outfileprefix = 'tirsh_', parameter_groups = None, rings = None,\
              block = None, variation = None, variation_type = None, iterations = 0, \
              random_seed = 0, mode = 'mad',fit_type='Error_Shaker'):
    Configuration['TIRSHAKER_TIME'][0] = datetime.now()
    if parameter_groups is None:
        parameter_groups = []
    if rings  is None:
        rings = []
    if block is None:
        block = []
    if  variation is None:
        variation = []
    if variation_type is None:
        variation_type = []
    # Initiate rng
    random.seed(random_seed)

    # Find the number of lines
    nur = int(Tirific_Template['NUR'])


    # Here we collect all parameter_groups as listed above and convert the lists into numbers
    # Then through the parameter groups and collect the parameter_groups
    allnumbers_in = []
    for j in range(len(parameter_groups)):
        numbers = []
        # Then go through the parameter_groups
        for k in range(len(parameter_groups[j])):
            #Once we checked all we can replace para with parameter_groups and use the input without =
            numbers.append([float(l) for l in Tirific_Template[parameter_groups[j][k]].split()])
            if parameter_groups[j][k] == 'CONDISP':
                pass
            else:
                while len(numbers[-1]) < nur:
                    numbers[-1].append(numbers[-1][-1])

        allnumbers_in.append(numbers)

    allnumbers_out = []
#Make sure some settings are blank
   
    Tirific_Template['OUTSET'] = ''
    Tirific_Template['PROGRESSLOG'] = ''
    Tirific_Template['TEXTLOG'] = ''
    Tirific_Template['TIRSMO'] = ''
    Tirific_Template['COOLGAL'] = ''
    Tirific_Template['TILT'] = ''
    Tirific_Template['BIGTILT'] = ''
    if nur < 15:
        Tirific_Template['INIMODE'] = 2
    else:
        Tirific_Template['INIMODE'] = 3

    Tirific_Template['LOGNAME'] = 'Error_Shaker.log'
    Tirific_Template['TIRDEF'] = 'Error_Shaker_Out.def'
    current_run='not set'
    original_running = copy.deepcopy(Configuration['TIRIFIC_RUNNING'])
    Configuration['TIRIFIC_RUNNING'] = False
    for i in range(iterations):
        Current_Template = copy.deepcopy(Tirific_Template)
        Current_Template['RESTARTID']= i
        # Provide some info where we are
        sf.print_log(f'''
        ******************************
        ******************************
        *** Tirshaker iteration {i:02d} ***
        ******************************
        ******************************
''',Configuration, case = ['screen'])

    #fortestin
    #    break

        # Find a name for a logfile and an output.def file
        #inname = 'blatirshin'
        #while os.path.exists(inname):
        #    inname = 'blatirshin_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        #logname = 'blatirshlog'
        #while os.path.exists(logname):
        #    logname = 'blatirshlog_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        #defname = 'blatirshdef'
        #while os.path.exists(defname):
        #    defname = 'blatirshdef_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        # Now random-generate the replacement for the numbers

        for j in range(len(parameter_groups)):
            if block[j]:
                variations = [variation[j]*random.uniform(-1.,1.)]*len(rings[j])
            else:
                variations = [variation[j]*random.uniform(-1.,1.) for klo in rings[j]]

            for k in range(len(parameter_groups[j])):
                current_list = [float(x) for x in Current_Template[parameter_groups[j][k]].split()]
                while len(current_list) < nur:
                    current_list.append(current_list[-1])
                for l in range(len(rings[j])):
                    if variation_type[j] == 'a':
                        current_list[rings[j][l]-1] = current_list[rings[j][l]-1]+variations[l]
                    else:
                        current_list[rings[j][l]-1] = current_list[rings[j][l]-1]*(1.+variations[l])
                format = sf.set_format(parameter_groups[j][k])
                Current_Template[parameter_groups[j][k]] = ' '.join([f'{x:{format}}' for x in current_list])



        #Current_Template['LOGNAME'] = logname
        #Current_Template['TIRDEF'] = defname

        write_tirific(Configuration,Current_Template, name =f'Error_Shaker/{fit_type}_In.def' )
        accepted,current_run = sf.run_tirific(Configuration, current_run, \
                                stage = 'shaker',fit_type = fit_type,\
                                max_ini_time= int(300*(int(Tirific_Template['INIMODE'])+1)))
        #os.system('tirific deffile= '+inname)

        # Read the values of the parameters requested

        allnumbers_out_part = []
        for j in range(len(parameter_groups)):
            numbers = []
            # Then go through the parameter_groups
            for k in range(len(parameter_groups[j])):
                numbers.append([float(x) for x in \
                    sf.load_tirific(Configuration,\
                        f"{Configuration['FITTING_DIR']}Error_Shaker/Error_Shaker_Out.def",\
                        Variables = [parameter_groups[j][k]],array=True)])
                if parameter_groups[j][k] == 'CONDISP':
                    pass
                else:
                    while len(numbers[-1]) < nur:
                        numbers[-1].append(numbers[-1][-1])

            allnumbers_out_part.append(numbers)
        allnumbers_out.append(allnumbers_out_part)

        # Now remove the files
        #os.remove(inname)
        #os.remove(logname)
        #os.remove(defname)

    #allnumbers_out = [[[[8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0], [36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0]], [[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], [22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0]], [[1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05], [1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05]]],
    #                 [[[9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0], [36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0]], [[2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], [22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0]], [[1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05], [1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05, 1e-05]]]]

    # Turn around the array, this could probably have been done above, but well...
    allparamsturned = []
    for j in range(len(allnumbers_in)):
        allparamsturned.append([])
        for k in range(len(allnumbers_in[j])):
            allparamsturned[j].append([])
            for l in range(len(allnumbers_in[j][k])):
                allparamsturned[j][k].append([])
                for m in range(len(allnumbers_out)):
                    allparamsturned[j][k][l].append(allnumbers_out[m][j][k][l])

    # Calculate mean and error
    allnumbers_final = []
    allnumbers_final_err = []
    #thenums = np.array(allnumbers_out)
    #print(len(allnumbers_in))
    for j in range(len(allnumbers_in)):
        allnumbers_final.append([])
        allnumbers_final_err.append([])
        for k in range(len(allnumbers_in[j])):
            allnumbers_final[j].append([])
            allnumbers_final_err[j].append([])
            for l in range(len(allnumbers_in[j][k])):
                # Attempt to use mad statistics for this
                if mode == 'mad':
#                    allparamsturned[j][k][l] = [1,2,3,1,2,3,1,2,3,1,2,3]
#                    meani = 0
#                    for ki in allparamsturned[j][k][l]:
#                        meani += float(ki)
#                    meani = meani/float(len(allparamsturned[j][k][l]))
#                    stdi = 0
#                    for ki in allparamsturned[j][k][l]:
#                        stdi += (meani-float(ki))*(meani-float(ki))
#                    stdi = np.sqrt(stdi/(len(allparamsturned[j][k][l])-1))
#                    print('{}+-{}'.format(meani, stdi))
                    median = np.median(np.array(allparamsturned[j][k][l]))
#                    mad = stats.median_absolute_deviation(np.array(allparamsturned[j][k][l]))
                    # Careful! This involves a scaling by 1.4826 by default as the default scale = 1.4826
                    madsigma = stats.median_abs_deviation(np.array(allparamsturned[j][k][l]))
                    average = np.average(np.array(allparamsturned[j][k][l]))
                    # Wow, np.std is the standard deviation using N and not N-1 in the denominator. So one has to use
                    std = np.sqrt(float(len(allparamsturned[j][k][l]))/float(len(allparamsturned[j][k][l])-1))*np.std(np.array(allparamsturned[j][k][l]))
                    allnumbers_final[j][k].append(stats.tmean(np.array(allparamsturned[j][k][l]), (median-3*madsigma, median+3*madsigma)))
                    allnumbers_final_err[j][k].append(stats.tstd(np.array(allparamsturned[j][k][l]), (median-3*madsigma, median+3*madsigma)))
                    #allnumbers_final[j][k].append(stats.tmean(np.array(allparamsturned[j][k][l])))
                    #allnumbers_final_err[j][k].append(stats.tstd(np.array(allparamsturned[j][k][l])))
                    sf.print_log('TIRSHAKER: Parameter: {:s} Ring: {:d} Pure average+-std: {:.3e}+-{:.3e} Median+-madsigma: {:.3e}+-{:.3e} Average+-sigma filtered: {:.3e}+-{:.3e} \n'.format(\
                                parameter_groups[j][k], l+1, average, std, median, madsigma, allnumbers_final[j][k][-1], allnumbers_final_err[j][k][-1])\
                                ,Configuration)
                else:
                    allnumbers_final[j][k].append(np.average(np.array(allparamsturned[j][k][l])))
                    allnumbers_final_err[j][k].append(np.sqrt(float(len(allparamsturned[j][k][l]))/float(len(allparamsturned[j][k][l])-1))*np.std(np.array(allparamsturned[j][k][l])))
    #print 'allnumbers_final'
    #print allnumbers_final
    #print 'allnumbers_final_err'
    #print allnumbers_final_err


    for j in range(len(parameter_groups)):
        for k in range(len(parameter_groups[j])):
            format = sf.set_format(parameter_groups[j][k])
            #Tirific_Template[parameter_groups[j][k]] = ' '.join([f'{x:{format}}' for x in allnumbers_final[j][k]])
            Tirific_Template.insert(f'{parameter_groups[j][k]}',f'# {parameter_groups[j][k]}_ERR',f"{' '.join([f'{x:{format}}' for x in allnumbers_final_err[j][k]])}")

            #Tirific_Template.insert(f'{parameter_groups[j][k]}',f'# {parameter_groups[j][k]}_ERR',' '.join([f'{x:{format}}' for x in allnumbers_final_err[j][k]]))

    # Put them into the output file
    # Write it to a copy of the file replacing the parameters


    sf.finish_current_run(Configuration, current_run)
    write_tirific(Configuration,Tirific_Template, name =f'Error_Shaker/{outfilename}' )
    Configuration['TIRIFIC_RUNNING'] = original_running
    Configuration['TIRSHAKER_TIME'][1] = datetime.now() 

tirshaker.__doc__ =f'''
 NAME:
    tirshaker

 PURPOSE:
    obtain errors through a FAT implemention of tirshaker developed by G.I.G. Jozsa.

 CATEGORY:
    run_functions

 INPUTS:
    Configuration = Standard FAT configuration
    outfileprefix (str)                    : Prefix to output parameters in outfilename
    parameter_groups (list of lists of str): List of parameters (parameter groups) that will be changed simultaneously
    rings (list of list of int)            : Ring numbers to be changed in parameter groups, starting at 1
    block (list of bool)                   : Change all rings by the same value (if True) or single rings (False). If the latter, the same ring for different parameters is changed by the same value
    variation (list of float)              : Amplitude of variation
    variation_type (list of str)           : Type of variation, 'a' for absolute, 'r' for relative
    iterations (int)                       : Number of re-runs
    random_seed (int)                      : Seed for random generator
    mode (str)                             : If 'mad' implements an outlier rejection.

 OPTIONAL INPUTS:


 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 This is a re-imaged version of the code developed by G.I.G. Jozsa found at  https://github.com/gigjozsa/tirshaker.git

 Takes a tirific def file filename and varies it iterations times
 and runs it as many times to then calculate the mean and the
 standard deviation of parameters that have been varied, both are
 put into the tirific deffile outfilename.

 With parameter_groups parameter groups are defined that are varied
 homogeneously. This is a list of list of parameter names including
 the '='-sign.  Caution! At this stage it is essential that in the
 .def file there is the '='-sign directly attached to the parameter
 names and that there is a space between parameter values and the
 '='.

 For each parameter group (list of parameter names), the values of
 the rings specified in rings (which is a list of integers, the
 list member corresponding to the parameter with the same index in
 parameter_groups) are varied. Block is a list of indicators if the
 values should be varied by ring (similar to the !-sign in the VARY
 parameter of tirific) or as a whole, again indices indicating the
 associated parameter group. Variation quantifies the maximum
 variation of the parameter (again indices matching), with
 variation type (indices matching) 'a' indicating an absolute
 variation, 'r' indicating a relative one. Parameters are varied by
 a uniform variation with maximal amplitude variation. So if a
 parameter is x and v the variation, the parameter gets changed by
 a number between -v and v in case the matching variation type is
 'a' and it gets changed by a number between -v*x and v*x if the
 matching variation type is 'r'. Tirific is started iterations
 times with varied input parameters but otherwise identical
 parameters (except for output files which are suppressed) and the
 results are recorded. For each varied parameter the mean and
 standard deviation is calculated and returned in the output .def
 file outfilename. In outfilename LOOPS is set to 0 and any output
 parameter is preceded by a prefix outfileprefix. Random_seed is
 the random seed to make the process deterministic. If mode ==
 'mad', the median and the MAD is calculated and, based on that,
 values beyond a 3-sigma (sigma estimated from MAD) bracket around
 the median are rejected before calculating mean and standard
 deviation.

'''
