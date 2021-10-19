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
from pyFAT_astro.Support.support_functions import run_tirific,print_log
# A cm in inch
cm_in_inch = 0.3937008

# A4widht_in_inch
A4widht_in_inch = 21.0*cm_in_inch

# A4height_in_inch
A4height_in_inch = 29.7*cm_in_inch

#tirshaker
def tirshaker_returnblanksatstart(string):
    """
    Helper to tirshaker
    """
    j = ''
    for i in string:
        if i != ' ':
            break
        j += ' '
    return j
#tirshaker
def tirshaker_convertostring(floatlist, parametername):
    """
    Helper to tirshaker
    """
#    print floatlist
    if 'XPOS' in parametername:
        return ' '.join(['{:+.7E}'.format(thing) for thing in floatlist])
    else:
        if 'YPOS' in parametername:
            return ' '.join(['{:+.7E}'.format(thing) for thing in floatlist])
        else:
            return ' '.join(['{:+.5E}'.format(thing) for thing in floatlist])
#tirshaker
def tirshaker(Configuration, filename = 'test.def', outfilename = 'test_out.def', \
              outfileprefix = 'tirsh_', parameter_groups = [], rings = [], block = []\
              , variation = [], variation_type = [], iterations = 0, random_seed = 0,\
               mode = 'mad',debug=False):
    # Initiate rng
    random.seed(random_seed)

    # Open the file
    lines = io.open(filename).read().split('\n')

    # Find the number of lines
    for line in lines:
        blankhere = tirshaker_returnblanksatstart(line)
        if line.find('NUR=',len(blankhere)) == len(blankhere):
            nur = int(line.split()[1])

    # Here we collect all parameter_groups as listed above and convert the lists into numbers
    # Then through the parameter groups and collect the parameter_groups
    allnumbers_in = []
    allblanks = []
    for j in range(len(parameter_groups)):
        numbers = []
        blanks = []
        # Then go through the parameter_groups
        for k in range(len(parameter_groups[j])):
            for line in lines:
                blankhere = tirshaker_returnblanksatstart(line)
                if line.find(parameter_groups[j][k],len(blankhere)) == len(blankhere):
                    blanks.append(blankhere)
                    numbers.append([float(l) for l in line.split()[1:]])

                    # Expand to number of rings by extrapolation unless CONDISP
                    if parameter_groups[j][k] == 'CONDISP=':
                       pass
                    else:
                        while len(numbers[-1]) < nur:
                            numbers[-1].append(numbers[-1][-1])
                    break
        allnumbers_in.append(numbers)
        allblanks.append(blanks)

    allnumbers_out = []
    for i in range(iterations):

        # Provide some info where we are
        print_log(f'''
        ******************************
        ******************************
        *** Tirshaker iteration {i:02d} ***
        ******************************
        ******************************
''',Configuration['OUTPUTLOG'])

    #fortestin
    #    break

        # Find a name for a logfile and an output.def file
        inname = 'blatirshin'
        while os.path.exists(inname):
            inname = 'blatirshin_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        logname = 'blatirshlog'
        while os.path.exists(logname):
            logname = 'blatirshlog_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        defname = 'blatirshdef'
        while os.path.exists(defname):
            defname = 'blatirshdef_'+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))

        replacement = copy.deepcopy(allnumbers_in)
        # Now random-generate the replacement for the numbers
        for j in range(len(parameter_groups)):
            if block[j]:
                variations = [variation[j]*random.uniform(-1.,1.)]*len(rings[j])
            else:
                variations = [variation[j]*random.uniform(-1.,1.) for klo in rings[j]]

            for k in range(len(parameter_groups[j])):
                for l in range(len(rings[j])):
                    #print(j,k,l,len(parameter_groups),len(parameter_groups[j]),len(rings[j]),len(variations))
                    if variation_type[j] == 'a':
                        replacement[j][k][rings[j][l]-1] = replacement[j][k][rings[j][l]-1]+variations[l]
                    else:
                        replacement[j][k][rings[j][l]-1] = replacement[j][k][rings[j][l]-1]*(1.+variations[l])

        # Write it to a copy of the file replacing the parameters
        #outstring = 'ACTION= 1\nPROMPT= 1\n'
        outstring = 'PROMPT= 1\n'
        #outstring=''
        # Create a firstfound array, all True
        firstfound = []
        for j in range(len(parameter_groups)):
            firstfound.append([])
            for k in parameter_groups[j]:
                firstfound[j].append(True)

        for linenum in range(len(lines)):
            anysearch = True
            for j in range(len(parameter_groups)):
                for k in range(len(parameter_groups[j])):
    #                print('pragr: {}'.format(parameter_groups[j][k]))
                    blankhere = tirshaker_returnblanksatstart(lines[linenum])
                    if lines[linenum].find(parameter_groups[j][k],len(blankhere)) == len(blankhere):
                        # Comment the line (not necessary, so we pass)
                        pass
                        # outstring = outstring+blankhere+'#'+lines[linenum][len(blankhere):]+'\n'

                        # Replace the first line
                        if firstfound[j][k]:
                            outstring = outstring+allblanks[j][k]+parameter_groups[j][k]+' '+tirshaker_convertostring(replacement[j][k],parameter_groups[j][k])+'\n'
                            firstfound[j][k] = False
                        anysearch = False
                        break
            if anysearch:

                # No output but the output.def file and the log
                blankhere = tirshaker_returnblanksatstart(lines[linenum])
                lbl = len(blankhere)
                if lines[linenum].find('LOGNAME=',lbl) == lbl:
                    outstring = outstring+'LOGNAME= '+logname+'\n'
                elif lines[linenum].find('OUTSET=',lbl) == lbl:
                    outstring = outstring+'OUTSET= '+'\n'
                elif lines[linenum].find('PROGRESSLOG=',lbl) == lbl:
                    outstring = outstring+'PROGRESSLOG= '+'\n'
                elif lines[linenum].find('TEXTLOG=',lbl) == lbl:
                    outstring = outstring+'TEXTLOG= '+'\n'
                elif lines[linenum].find('TIRDEF=',lbl) == lbl:
                    outstring = outstring+'TIRDEF= '+defname+'\n'
                elif lines[linenum].find('TIRSMO=',lbl) == lbl:
                    outstring = outstring+'TIRSMO= '+'\n'
                elif lines[linenum].find('COOLGAL=',lbl) == lbl:
                    outstring = outstring+'COOLGAL= '+'\n'
                elif lines[linenum].find('TILT=',lbl) == lbl:
                    outstring = outstring+'TILT= '+'\n'
                elif lines[linenum].find('BIGTILT=',lbl) == lbl:
                    outstring = outstring+'BIGTILT= '+'\n'
                elif lines[linenum].find('GR_DEVICE=',lbl) == lbl:
                    outstring = outstring+'GR_DEVICE= '+'\n'
                else:
                    outstring = outstring+lines[linenum]+'\n'

        # Dump it into a file and then run tirific
        thething = io.open(inname,'w')
        thething.write(outstring)
        thething.close()

        accepted,current_run = run_tirific(Configuration, 'Not_ZEd',deffile=inname, stage = 'shaker',fit_type = 'Error_Shaker', debug = debug)
        #os.system('tirific deffile= '+inname)

        # Read the values of the parameters requested
        thething = io.open(defname,'r')
        outlines = thething.read().split('\n')
        thething.close()
        allnumbers_out_part = []
        for j in range(len(parameter_groups)):
            numbers = []
            # Then go through the parameter_groups
            for k in range(len(parameter_groups[j])):
                for line in outlines:
                    blankhere = tirshaker_returnblanksatstart(line)
                    if line.find(parameter_groups[j][k],len(blankhere)) == len(blankhere):
                        blanks.append(blankhere)
                        numbers.append([float(l) for l in line.split()[1:]])
            allnumbers_out_part.append(numbers)
        allnumbers_out.append(allnumbers_out_part)

        # Now remove the files
        os.remove(inname)
        os.remove(logname)
        os.remove(defname)

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
                    madsigma = stats.median_absolute_deviation(np.array(allparamsturned[j][k][l]))
                    average = np.average(np.array(allparamsturned[j][k][l]))
                    # Wow, np.std is the standard deviation using N and not N-1 in the denominator. So one has to use
                    std = np.sqrt(float(len(allparamsturned[j][k][l]))/float(len(allparamsturned[j][k][l])-1))*np.std(np.array(allparamsturned[j][k][l]))
                    allnumbers_final[j][k].append(stats.tmean(np.array(allparamsturned[j][k][l]), (median-3*madsigma, median+3*madsigma)))
                    allnumbers_final_err[j][k].append(stats.tstd(np.array(allparamsturned[j][k][l]), (median-3*madsigma, median+3*madsigma)))
                    #allnumbers_final[j][k].append(stats.tmean(np.array(allparamsturned[j][k][l])))
                    #allnumbers_final_err[j][k].append(stats.tstd(np.array(allparamsturned[j][k][l])))
                    print_log('TIRSHAKER: Parameter: {:s} Ring: {:d} Pure average+-std: {:.3e}+-{:.3e} Median+-madsigma: {:.3e}+-{:.3e} Average+-sigma filtered: {:.3e}+-{:.3e} \n'.format(\
                                parameter_groups[j][k], l+1, average, std, median, madsigma, allnumbers_final[j][k][-1], allnumbers_final_err[j][k][-1])\
                                ,Configuration['OUTPUTLOG'],screen=False)
                else:
                    allnumbers_final[j][k].append(np.average(np.array(allparamsturned[j][k][l])))
                    allnumbers_final_err[j][k].append(np.sqrt(float(len(allparamsturned[j][k][l]))/float(len(allparamsturned[j][k][l])-1))*np.std(np.array(allparamsturned[j][k][l])))
    #print 'allnumbers_final'
    #print allnumbers_final
    #print 'allnumbers_final_err'
    #print allnumbers_final_err

    # Put them into the output file
    # Write it to a copy of the file replacing the parameters
    outstring = ''

    # Create a firstfound array, all True
    firstfound = []
    for j in range(len(parameter_groups)):
        firstfound.append([])
        for k in parameter_groups[j]:
            firstfound[j].append(True)

    for linenum in range(len(lines)):
        anysearch = True
        for j in range(len(parameter_groups)):
            for k in range(len(parameter_groups[j])):
    #                print('pragr: {}'.format(parameter_groups[j][k]))
                blankhere = tirshaker_returnblanksatstart(lines[linenum])
                if lines[linenum].find(parameter_groups[j][k],len(blankhere)) == len(blankhere):
                    # Comment the line (not necessary, so we pass)
                    outstring = outstring+blankhere+'#'+lines[linenum][len(blankhere):]+'\n'

                    # Replace the first line
                    if firstfound[j][k]:
                        outstring = outstring+allblanks[j][k]+parameter_groups[j][k]+' '+tirshaker_convertostring(allnumbers_final[j][k],parameter_groups[j][k])+'\n'
                        outstring = outstring+allblanks[j][k]+'ERR_'+parameter_groups[j][k]+' '+tirshaker_convertostring(allnumbers_final_err[j][k],parameter_groups[j][k])+'\n'
                        firstfound[j][k] = False
                    anysearch = False
                    break

        if anysearch:

            # No output but the output.def file and the log
            blankhere = tirshaker_returnblanksatstart(lines[linenum])
            lbl = len(blankhere)
            notfound = True
            for parametername in ['LOGNAME=', 'OUTSET=', 'PROGRESSLOG=', 'TEXTLOG=', 'TIRDEF=', 'TIRSMO=', 'COOLGAL=', 'TILT=', 'BIGTILT=', 'GR_DEVICE=']:
                if lines[linenum].find(parametername,lbl) == lbl:
                    try:
                        parameter = outfileprefix+lines[linenum].split()[1]
                        outstring = outstring+blankhere+parametername+' '+parameter+'\n'
                    except:
                        outstring = outstring+blankhere+parametername+'\n'
                    notfound = False
            if lines[linenum].find('LOOPS=',lbl) == lbl:
                outstring = outstring+'LOOPS= 0\n'
                notfound = False
            if notfound:
                outstring = outstring+lines[linenum]+'\n'

    # Dump it into a file

    thething = io.open(outfilename,'w')
    thething.write(outstring)
    thething.close()
tirshaker.__doc__ =f'''
 NAME:
    tirshaker

 PURPOSE:
    obtain mcmc errors through a FAT implemention of tirshaker developed by G.I.G. Jozsa.

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
    debug = False


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
