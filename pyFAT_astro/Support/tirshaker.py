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

def busyf(x,a,b1,b2,w,xe,xp,c,n):
    """
    The busy function, see Westmeier et al. 2014
    """
    return a/4.*(spec.erf(b1*(w+x-xe))+1)*(spec.erf(b2*(w-x+xe))+1)*(c*np.power(np.fabs(x-xp),n)+1)

def busyres(a, xes = None, ys = None):
    """
    Function fit for fitting with scipy, computing the residuals w.r.t. an observed spectrum, chi2 without weighting
    """
    return np.sum(np.power(ys - busyf(xes, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7]),2))

def read3dim(hdu):
    """
    Strip all non-3d info and return a newly created hdu.

    Input:
    hdu: astropy.io.fits HDU
    plane: number of plane, starting with 0
    """
    #hdu = fits.open(cube)[0]
    theshape = hdu.data.shape
    length = len(theshape)

    if length == 3:
        return hdu

    # We assume that the data are in that order
    if length > 3:
        newdata = hdu.data.reshape(theshape[-3],theshape[-2],theshape[-1])

    # Now rearrange header
    hdu.header['naxis'] = 3
    while length > 3:
        try:
            hdu.header.pop('NAXIS{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CDELT{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CRVAL{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CRPIX{:d}'.format(length))
        except:
            pass
        try:
            hdu.header.pop('CTYPE{:d}'.format(length))
        except:
            pass
        length -= 1

    # Re-generate the thing
    return(fits.PrimaryHDU(data=newdata, header=hdu.header))


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
    thenums = np.array(allnumbers_out)
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
                    print_log('TIRSHAKER: Parameter: {:s} Ring: {:d} Pure average+-std: {:.3e}+-{:.3e} Median+-madsigma: {:.3e}+-{:.3e} Average+-sigma filtered: {:.3e}+-{:.3e}'.format(\
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


def replacekeys(instr = '', outstr = '', replace = {}):
    """
    Take string instr, scan for keywords in replace to replace the line with the value of the keys, append if key does not exist. Blanks at the beginning don't count.
    """

    lines = instr.split('\n')
    for linenum in range(len(lines)):
        blankhere = tirshaker_returnblanksatstart(lines[linenum])
        lbl = len(blankhere)
        notfound = True
        for key in replace.keys():
            foundkeys = []
            if lines[linenum].find(key,lbl) == lbl:
                outstr = outstr+key+' '+replace[key]+'\n'
                notfound = False
                foundkeys.append(key)
        for key in foundkeys:
            replace.pop(key)
        if notfound:
            outstr = outstr+lines[linenum]+'\n'
    for key in replace.keys():
        outstr = outstr+key+' '+replace[key]+'\n'

    return outstr

def getkey(instr = '', key = ''):
    """
    find key in instr and return a list of value (separated by spaces in instr), including the key. Blanks don't count. First match counts.
    """
    lines = instr.split('\n')
    for linenum in range(len(lines)):
        blankhere = tirshaker_returnblanksatstart(lines[linenum])
        if lines[linenum].find(key,len(blankhere)) == len(blankhere):
            return lines[linenum].split()
    return [key,'0.0']

def findnfilename(prefix = ''):
    """
    Generate a string which does not represent any existing file name by appending a random string to prefix
    """
    while os.path.exists(prefix):
        prefix = prefix+''.join(random.choice('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWxyz') for i in range(6))
    return prefix

def gethirad(indef = '', dr = 1., random_seed = 1):
    """
    Spit out HI radius and velocity at that radius from file resulting from tirshaker run, where dr is the precision in arcsec.
    """

    asolarmass = 1/(1.248683E24*8.01325e-21) # One solar mass per square parsec in Jy km/s / sqare arcsec for redshift 0

    random.seed(random_seed)

    # Open the file
    if indef == '':
        return
    dinstream = io.open(indef)
    dinstr = dinstream.read()
    dinstream.close

    # Find maxrad and create many radii
    maxrad = float(getkey(instr = dinstr, key = 'RADI=')[-1])
    newradi = np.arange(0., maxrad+dr, dr)

    # Get three surface brightness lists
    sbr   = [float(i) for i in getkey(instr = dinstr, key = 'SBR=')[1:]]
    sberr = [float(i) for i in getkey(instr = dinstr, key = 'ERR_SBR=')[1:]]
    vroterr = getkey(instr = dinstr, key = 'ERR_VROT=')[1:]

    sbrno = ' '.join(['{:.5E}'.format(sbr[i])          for i in range(len(sbr))])
    sbrpl = ' '.join(['{:.5E}'.format(sbr[i]+sberr[i]) for i in range(len(sbr))])
    sbrmi = ' '.join(['{:.5E}'.format(sbr[i]-sberr[i]) for i in range(len(sbr))])
#    print(sbr, '\n', sberr, '\n', sbrno, '\n', sbrpl, '\n', sbrmi)

    hirad = []
    hivrot = []
    dhivrot = []
#    for i in [sbrmi, sbrno, sbrpl]:
    for i in [sbrmi, sbrno, sbrpl]:
        # Find names for new files
        outname = findnfilename('blatirshout_')
        inname  = findnfilename('blatirshin_')
        logname = findnfilename('blatirslog_')

        # Write new deffiles
        replace = {
            'LOGNAME='    : logname,
            'OUTSET='     : '',
            'PROGRESSLOG=': '',
            'TEXTLOG='    : '',
            'TIRSMO='     : '',
            'COOLGAL='    : '',
            'TILT='       : '',
            'BIGTILT='    : '',
            'GR_DEVICE='  : '',
            'TIRDEF='     : outname,
            'TIRNR='      : str(len(newradi)),
            'TIRRAD='     : ' '.join([str(round(i,1)) for i in newradi]),
            'INCL='       : ' '.join(vroterr),
            'SBR='        : i
        }
        instr = replacekeys(instr = dinstr, outstr = '', replace = replace)

        # Put it into inname
        instream = io.open(inname, 'w')
        instream.write(instr)
        instream.close()

        # Run tirific
        os.system('tirific deffile= '+inname)

        # Open, scan through lines and determine HI radius
        doutstream = io.open(outname)
        doutstr = doutstream.read()
        doutstream.close()

        radi = [float(i) for i in getkey(instr = doutstr, key = 'RADI=')[1:]]
        vrot = [float(i) for i in getkey(instr = doutstr, key = 'VROT=')[1:]]
        dvrot = [float(i) for i in getkey(instr = doutstr, key = 'INCL=')[1:]]
        sbr = [float(i) for i in getkey(instr = doutstr, key = 'SBR=')[1:]]

        hiradhere = -1.
        for j in range(1,len(radi)):
            if sbr[j-1] > asolarmass and sbr[j] < asolarmass:
                hiradhere = (radi[j]+radi[j-1])/2
                hivrothere = (vrot[j]+vrot[j-1])/2
                dhivrothere = (dvrot[j]+dvrot[j-1])/2

        if hiradhere < 0.:
            if sbr[0] > asolarmass:
                hiradhere = radi[-1]
                hivrothere = vrot[-1]
                dhivrothere = vrot[-1]
            else:
                hiradhere = 0.
                hivrothere = 0.
                dhivrothere = 0.

        hirad.append(hiradhere)
        hivrot.append(hivrothere)
        dhivrot.append(dhivrothere)
        os.remove(inname)
        os.remove(logname)
        os.remove(outname)
    hiradius = hirad[1]
    dhiradius = (hirad[2]-hirad[0])/2.
    hirvelocity = hivrot[1]
    dhirvelocity = dhivrot[1]

    distance = float(getkey(instr = dinstr, key = 'DISTANCE=')[-1])*1000000.
    arcsectokpc = distance*np.pi/(180.*3600.*1000.)
    hiradiuskpc = hiradius*arcsectokpc
    dhiradiuskpc = dhiradius*arcsectokpc
    print('HI Radius: {:.2f}+-{:.2f} arcsec = {:.2f}+-{:.2f} kpc'.format(hiradius, dhiradius, hiradiuskpc, dhiradiuskpc))
    print('VROT(HI Radius): {:.2f}+-{:.2f} km/s'.format(hirvelocity, dhirvelocity))
    print('HI Diameter: {:.2f}+-{:.2f} arcsec = {:.2f}+-{:.2f} kpc'.format(hiradius*2., dhiradius*2., hiradiuskpc*2., dhiradiuskpc*2.))
    hiradm = (hiradiuskpc*1000.*u.pc).to(u.m)
    dhiradm = (dhiradiuskpc*1000.*u.pc).to(u.m)
    vrotm = hirvelocity*1000.*u.m/u.s
    dvrotm = dhirvelocity*1000.*u.m/u.s
    mdyn = (hiradm*np.power(vrotm,2.)/c.G).to(u.solMass)
    dmdyn = (np.sqrt(np.power(dvrotm*2*vrotm*hiradm/c.G,2)+np.power(dhiradm*vrotm*vrotm/c.G,2))).to(u.solMass)
    print('Mdyn(HI Radius): {:.2E}+-{:.2E} solar masses'.format(mdyn.value, dmdyn.value))
    return

def getspectrandhiwidth(indef = '', figname = None, maskfile = None, distance = None, plotbusy = True, widfrac = 0.2, modelcolour = 'DeepPink'):
    '''
Only for high S/N: extracts and plots spectra of the original file
    as provided in the input .def file, then of model. Input .def file should be same format as a tirshaker output (i.e. contain ERR_INCL and ERR_VROT and ERR_SBR lines). If maskfile is
    given, will use it as a mask. Maskfile should be same size as .def
    input file and be 1 where there is emission and 0 elsewhere. Then determines wwidfrac*100 two ways for the data, the model, and model+error and model-error, once finding the first channels above widfrac*100% of the peak on either side of the peak and then linearly interpolating to find the positon of the intersection with the 20% line. Then by fitting a busy function to the profile. Widfrac is the fraction of the peak that we are supposed to trace. Typically 0.5 (W50) or 0.2 (W20)
    '''
    if indef == '':
        return
    dinstream = io.open(indef)
    dinstr = dinstream.read()
    dinstream.close

    if distance == None:
        distance = 0.1
    else:
        distance = distance.to(u.pc).value/1000000.

    # Get inclination and its errors, sbr and its errors, vrot and its errors, sbr and its errors
    pargroups = {'INCL': ['INCL=', 'ERR_INCL='], 'VROT' : ['VROT=','ERR_VROT='], 'SBR': ['SBR=', 'ERR_SBR='], 'SDIS': ['SDIS=', 'ERR_SDIS='], 'CONDISP': ['CONDISP=', 'ERR_CONDISP=']}

    # Filling dicts with values from file
    for parname in pargroups.keys():
        # Parameter
        pargroups[parname].append([float(i) for i in getkey(instr = dinstr, key = pargroups[parname][0])[1:]])
        # Error
        pargroups[parname].append([float(i) for i in getkey(instr = dinstr, key = pargroups[parname][1])[1:]])
        # Min
        pargroups[parname].append([pargroups[parname][2][i]-pargroups[parname][3][i] for i in range(len(pargroups[parname][2]))])
        # Max
        pargroups[parname].append([pargroups[parname][2][i]+pargroups[parname][3][i] for i in range(len(pargroups[parname][2]))])
        # Now in letters, 6: normal, 7: min, 8: max
        for j in [2, 4, 5]:
            pargroups[parname].append(' '.join(['{:.5E}'.format(i) for i in pargroups[parname][j]]))

    # Make models
    models = {6: None, 7: None, 8: None}

    for i in [6,7,8]:
        models[i] = findnfilename('bloutfits_')
        inname  = findnfilename('blatirshin_')
        logname = findnfilename('blatirslog_')

        replace = {
            'LOGNAME='    : logname,
            'OUTSET='     : models[i],
            'PROGRESSLOG=': '',
            'TEXTLOG='    : '',
            'TIRSMO='     : '',
            'COOLGAL='    : '',
            'TILT='       : '',
            'BIGTILT='    : '',
            'GR_DEVICE='  : '',
            'TIRDEF='     : '',
        }
        for parname in pargroups.keys():
            replace[parname+'='] = pargroups[parname][i]
            replace[parname+'_2='] = pargroups[parname][i]
            replace[parname+'_3='] = pargroups[parname][i]
            replace[parname+'_4='] = pargroups[parname][i]
            replace[parname+'_5='] = pargroups[parname][i]
            replace[parname+'_6='] = pargroups[parname][i]
        instr = replacekeys(instr = dinstr, outstr = '', replace = replace)

        # Put it into inname
        instream = io.open(inname, 'w')
        instream.write(instr)
        instream.close()

        # Run tirific
        os.system('tirific deffile= '+inname)
        # Remove files
        os.remove(inname)
        os.remove(logname)

    # Now get spectra and inclinations, i.e. read triplets [velocity, flux density in mJy, inclination]
    spectra = {'orig': [getkey(instr = dinstr, key = 'INSET=')[-1], np.array(pargroups['INCL'][2]).mean()], 'model': [models[6], np.array(pargroups['INCL'][2]).mean()], 'model_low': [models[7], np.array(pargroups['INCL'][4]).mean()], 'model_hi': [models[8], np.array(pargroups['INCL'][5]).mean()]}

    for key in spectra.keys():
        filename = spectra[key][0]
        fset = fits.open(filename)
        hdu  = read3dim(fset[0])

        inclination =spectra[key][1]

        # inclination = np.array([float(i) for i in getkey(instr = dinstr, key = 'INCL=')[1:]]).mean()

        # Open model
        #model = fits.open(outfitname)
        #hdu_model = read3dim(model[0])

        # Open maskfile or create a mask which is no mask
        if maskfile != None:
            fset_mask = fits.open(maskfile)
            hdu_mask = read3dim(fset_mask[0])
            data = hdu.data*hdu_mask.data
            #data_model = hdu_model.data*hdu_mask.data
        else:
            data = hdu.data
            #data_model = hdu_model.data

        # Get intel
        crpix3 = int(hdu.header['CRPIX3'])
        crval3 = float(hdu.header['CRVAL3'])/1000.
        cdelt3 = float(hdu.header['CDELT3'])/1000.
        naxis3 = int(hdu.header['NAXIS3'])

        bmaj = float(hdu.header['BMAJ'])
        bmin = float(hdu.header['BMIN'])
        pixx = np.fabs(float(hdu.header['CDELT1']))
        pixy = np.fabs(float(hdu.header['CDELT2']))

        # Collapse on first axis (velocity)
        spectrum = 1000.*np.sum(data, axis = (1,2))*pixx*pixy/(1.1330900354567984*bmaj*bmin)
        # spectrum_model = 1000.*np.sum(data_model, axis = (1,2))*pixx*pixy/(1.1330900354567984*bmaj*bmin)
        # Calculate a velocity axis
        firstvalvel = (1.-crpix3)*cdelt3+crval3

        # This creates naxis3+1 samples for some extremely strange reason not possible with arange
        vaxis  = np.array([firstvalvel+i*cdelt3 for i in range(naxis3)])

        peak = spectrum.max()
        # This is just by linear interpolation between two channels
        chanmax_20 = -np.argmax(spectrum[-1:0:-1] > (widfrac*peak))-1
        chanmax_20_aft = chanmax_20+1
        vmax_20 = vaxis[chanmax_20]+cdelt3*(spectrum[chanmax_20]-widfrac*peak)/((spectrum[chanmax_20]-widfrac*peak)+(widfrac*peak-(spectrum[chanmax_20_aft])))
        chanmin_20 = np.argmax(spectrum > (widfrac*peak))
        chanmin_20_bef = chanmin_20-1
        vmin_20 = vaxis[chanmin_20]-cdelt3*(spectrum[chanmin_20]-widfrac*peak)/((spectrum[chanmin_20]-widfrac*peak)+(widfrac*peak-(spectrum[chanmin_20_bef])))
        #vmin_20 = vaxis[chanmin_20_bef]
        w20 = (vmax_20-vmin_20)/np.sin(np.pi*inclination/180.)
        w20uc = (vmax_20-vmin_20)

        vmax_50 = vaxis[-np.argmax(spectrum[-1:0:-1] > (0.5*peak))-1]
        vmin_50 = vaxis[np.argmax(spectrum > (0.5*peak))]
        fwhm = vmax_50-vmin_50
        #peak_model = spectrum_model.max()
        #vmax_20_model = vaxis[-np.argmax(spectrum_model[-1:0:-1] > (widfrac*peak_model))-1]
        #vmin_20_model = vaxis[np.argmax(spectrum_model > (widfrac*peak_model))]

        #/np.sin(np.pi*inclination/180.)
        xpos = (vmax_50+vmin_50)/2.

        # Now get a busy fit
        reo = least_squares(busyres, np.array([4*peak, np.sqrt(np.pi)/(np.sqrt(8.)*fwhm/2.3548), np.sqrt(np.pi)/(np.sqrt(8.)*fwhm/2.3548), 1., xpos, xpos, 0.0001, 2.]), kwargs = {'xes' : vaxis, 'ys' : spectrum}, max_nfev = 1000000)['x']

        # Then find vmax and vmin again
        delta = 0.5*(vaxis[-1]-vaxis[0])/abs(vaxis[-1]-vaxis[0])
        vaxisbus = np.arange(vaxis[0],vaxis[-1]+delta,delta)
        spectrumbus = busyf(vaxisbus, reo[0], reo[1], reo[2], reo[3], reo[4], reo[5], reo[6], reo[7])
        maxibus = np.max(spectrumbus)
        vmaxbus_20 = vaxisbus[-np.argmax(spectrumbus[-1:0:-1] > (widfrac*maxibus))-1]
        vminbus_20 = vaxisbus[np.argmax(spectrumbus > (widfrac*maxibus))]
        w20bus = (vmaxbus_20-vminbus_20)/np.sin(np.pi*inclination/180.)
        w20ucbus = (vmaxbus_20-vminbus_20)

        #w20_model = (vmax_20-vmin_20)/np.sin(np.pi*inclination/180.)

        totflux = spectrum.sum()*cdelt3/1000.
        HImass = totflux*2.36E5*distance*distance

        spectrump = np.append(spectrum[:1],spectrum)
        delta = vaxis[1]-vaxis[0]
        vaxisp = np.append(vaxis[:1]-delta/2., vaxis+delta/2.)

        delta = 0.5
        vaxispbus = np.arange(vaxisp[0],vaxisp[-1]+delta,delta)
        spectrumpbus = busyf(vaxispbus, reo[0], reo[1], reo[2], reo[3], reo[4], reo[5], reo[6], reo[7])
        #print(reo)
        #spectrumpbus_in = busyf(vaxispbus, 4*peak, np.sqrt(np.pi)/(np.sqrt(8.)*fwhm/2.3548), np.sqrt(np.pi)/(np.sqrt(8.)*fwhm/2.3548), 0., xpos, xpos, 0., 2.)

        # Now put this into the dict
        spectra[key] = [vaxis, spectrum, vaxisp, spectrump, peak, vmax_20, vmin_20, w20, totflux, HImass, vaxispbus, spectrumpbus, vmaxbus_20, vminbus_20, w20bus, inclination, w20uc, w20ucbus]

        #print('w20: {:.2f} km/s'.format(w20))
        #print('w20_model: {:.2f} km/s'.format(w20_model))

        #totflux_model = spectrum_model.sum()
        fset.close()
        if maskfile != None:
            fset_mask.close()

    # Models not needed any more
    for key in [6,7,8]:
        os.remove(models[key])

    #####
    # Plotting from here onwards
    #####
    #spectrump = np.append(spectrum[:1],spectrum)
    #spectrump_model = np.append(spectrum_model[:1],spectrum_model)
    #delta = vaxis[1]-vaxis[0]
    #vaxisp = np.append(vaxis[:1]-delta/2., vaxis+delta/2.)

    fig = plt.figure(figsize=(A4height_in_inch, A4widht_in_inch))
#    plt.rc('text', usetex=True)
    ax = fig.subplots()
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
    ax.tick_params(which='both', direction = 'in', width = 2, labelsize='xx-large', top = True, bottom = True, left = True, right = True)
    ax.tick_params(which='major', length = 10)
    ax.tick_params(which='minor', length = 5)
    ax.set_xlabel(r'$v$($\,\mathrm{{km}}\,\mathrm{{s}}^{{-1}}$)')
    ax.set_ylabel(r'$S$($\,\mathrm{{mJy}}$)')
    ax.yaxis.label.set_size('xx-large')
    ax.xaxis.label.set_size('xx-large')
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2)
    ax.axhline(y = 0., linewidth=2, color='k', linestyle= 'dashed')

    # Possible for key are: 'orig', 'model', 'model_low', 'model_hi'
    # spectra[key] = [vaxis, spectrum, vaxisp, spectrump, peak, vmax_20, vmin_20, w20_prep, totflux, HImass]
    thingtoplot = ['model', 'orig']
    lws = [3,2]
    colors = [modelcolour, 'DarkBlue']
    for i in range(len(thingtoplot)):
#        ax.plot(spec[10], spec[12], color = 'black', lw = lws[i], linestyle = 'dotted')
        spec = spectra[thingtoplot[i]]
        ax.axhline(y = widfrac*spec[4], linewidth=lws[i], color=colors[i], linestyle= 'dashed')
        if plotbusy:
            ax.axvline(x = spec[12], linewidth=lws[i], color=colors[i], linestyle= 'dashed')
            ax.axvline(x = spec[13], linewidth=lws[i], color=colors[i], linestyle= 'dashed')
            hungl = mplcolors.to_rgba(colors[i])
            newcolor = (hungl[0], hungl[1], hungl[2], 0.6)
            ax.plot(spec[10], spec[11], color = newcolor, lw = lws[i], linestyle = 'dotted')
        else:
            ax.axvline(x = spec[5], linewidth=lws[i], color=colors[i], linestyle= 'dashed')
            ax.axvline(x = spec[6], linewidth=lws[i], color=colors[i], linestyle= 'dashed')
        ax.step(spec[2], spec[3], color = colors[i], lw = lws[i])
        ax.set_xbound(lower=spec[2][0], upper=spec[2][-1])

    fig.tight_layout()
    #ax.plot(showgouse, calculated, 'g-')
    #ax.plot(showgouse, fitted, 'r-')
    #ax.axvline(x = average+threshold*stdev, linewidth=2, color='k')
    #ax.set_xlim(min(bin_edges), max(bin_edges))
    #ax.set_title(title)
    #    grofo = mpl.colors.to_rgba('DarkBlue')
    #    ourcolor = (grofo[0], grofo[1], grofo[2], 0.5)
    #    for i in dir(ax):
    #        if '' in i:
    #            print(i)

    # Finally analysis
    # Possible for key are: 'orig', 'model', 'model_low', 'model_hi'
    # spectra[key] = [vaxis, spectrum, vaxisp, spectrump, peak, vmax_20, vmin_20, w20, totflux, HImass]

    print('The peak flux of the original is:        {:.1f} mJy'.format(spectra['orig'][4]))
    print('The peak flux of the model is:           {:.1f}+{:.1f}-{:.1f} mJy'.format(spectra['model'][4], spectra['model_hi'][4]-spectra['model'][4], spectra['model'][4]-spectra['model_low'][4]))
    print('Just using the channel values:')
    print('The w{:d} of the original is:              {:.1f} km/s'.format(int(widfrac*100.),spectra['orig'][16]))
    print('The w{:d} of the model is:                 {:.1f}+{:.1f}-{:.1f} km/s'.format(int(widfrac*100.),spectra['model'][16], spectra['model_hi'][16]-spectra['model'][16], spectra['model'][16]-spectra['model_low'][16]))
    print('The w{:d}c of the original is:              {:.1f} km/s'.format(int(widfrac*100.),spectra['orig'][7]))
    print('The w{:d}c of the model is:                 {:.1f}+{:.1f}-{:.1f} km/s'.format(int(widfrac*100.),spectra['model'][7], spectra['model_hi'][7]-spectra['model'][7], spectra['model'][7]-spectra['model_low'][7]))
    print('Using a fit with the busy function first:')
    print('The w{:d} of the original is:              {:.1f} km/s'.format(int(widfrac*100.), spectra['orig'][17]))
    print('The w{:d} of the model is:                 {:.1f}+{:.1f}-{:.1f} km/s'.format(int(widfrac*100.), spectra['model'][17], spectra['model_hi'][17]-spectra['model'][17], spectra['model'][17]-spectra['model_low'][17]))
    print('The w{:d}c of the original is:              {:.1f} km/s'.format(int(widfrac*100.), spectra['orig'][14]))
    print('The w{:d}c of the model is:                 {:.1f}+{:.1f}-{:.1f} km/s'.format(int(widfrac*100.), spectra['model'][14], spectra['model_hi'][14]-spectra['model'][14], spectra['model'][14]-spectra['model_low'][14]))
    print('The channelwidth is: {:.1f} km/s'.format(spec[2][1]-spec[2][0]))
    print('The channelwidth / sin(i) is: {:.1f} km/s'.format((spec[2][1]-spec[2][0])/np.sin(np.pi*spectra['model'][15]/180.)))
    print('The total flux of the original is:       {:.1f} Jy km/s'.format(spectra['orig'][8]))
    print('The total flux of the model is:          {:.1f}+{:.1f}-{:.1f} Jy km/s'.format(spectra['model'][8], spectra['model_hi'][8]-spectra['model'][8], spectra['model'][8]-spectra['model_low'][8]))
    print('The HI mass of the original is:          {:.3E} solar masses'.format(spectra['orig'][9]))
    print('The HI mass of the model is:             {:.3E}+{:.2E}-{:.2E} solar masses'.format(spectra['model'][9], spectra['model_hi'][9]-spectra['model'][9], spectra['model'][9]-spectra['model_low'][9]))


    if figname != None:
        fig.savefig(figname)
    return

def replaceforplot_lines(lines = None, gr_parms = None, gr_device = None):
    """
    helper working on a list of lines, see replaceforplot
    """

    # Create a block containing the correct output lines
    block = ['GR_PARMS= {:s}'.format(gr_parms)]
    grnr = 0
    for parameter in gr_parms.split(' ')[1:]:
        grnr += 1
        block += ['GR_COL_{:d}= 1'.format(grnr)]
        block += ['GR_LINES_{:d}= 1'.format(grnr)]

        for line in lines:
            blankhere = tirshaker_returnblanksatstart(line)
            lbl = len(blankhere)
            if line.find('ERR_'+parameter+'=',lbl) == lbl:
                block += ['GR_ERRB_{:d}= 1'.format(grnr)]
                block += ['GR_ERRV_{:d}='.format(grnr)+line[lbl+len('ERR_'+parameter+'='):]]


    # Then put it in after GR_PARMS
    outputlist = []
    for line in lines:
        blankhere = tirshaker_returnblanksatstart(line)
        lbl = len(blankhere)
        if line.find('GR_PARMS=',lbl) == lbl:
            outputlist += ['#'+line]
            outputlist.extend(block)
        elif line.find('GR_DEVICE=',lbl) == lbl:
            if gr_device == None:
                outputlist += [line]
            else:
                outputlist += ['#'+line]
                outputlist += ['GR_DEVICE= '+gr_device+'\n']
        else:
            outputlist += [line]
    return outputlist

def replaceforplot_lines_smooth(lines = None, gr_parms = None, parms = None, gr_device = None):
    """
    helper working on a list of lines, see replaceforplot
    """

    # Create a block containing the correct output lines
    block = ['GR_PARMS= {:s}'.format(gr_parms)]
    grnr = 0
    outputlist = []
    for line in lines:
        blankhere = tirshaker_returnblanksatstart(line)
        lbl = len(blankhere)
        notfound = True
        for parameter in parms.split(' '):
            if line.find(parameter+'=',lbl) == lbl:
                notfound = False
                outputlist += ['#'+line]
                ourline = [float(i) for i in line[lbl:].split(' ')[1:]]
                hanfloats = [(3.*ourline[0]+ourline[1])/4]
                [hanfloats.append((ourline[i-1]+2.*ourline[i]+ourline[i+1])/4) for i in range(1,len(ourline)-1)]
                hanfloats.append((3.*ourline[-1]+ourline[-2])/4)
                outstr = ''
                for hanfloat in hanfloats:
                    if 'XPOS' in parameter:
                        outstr += '{:+.7E} '.format(hanfloat)
                    elif 'YPOS' in parameter:
                         outstr += '{:+.7E} '.format(hanfloat)
                    else:
                         outstr += '{:+.5E} '.format(hanfloat)
                outputlist += [parameter+'= '+ outstr]
        if notfound:
            outputlist += [line]

    for parameter in gr_parms.split(' ')[1:]:
        grnr += 1
        block += ['GR_COL_{:d}= 1'.format(grnr)]
        block += ['GR_LINES_{:d}= 1'.format(grnr)]

    # Then put it in after GR_PARMS
    outputlist2 = []
    for line in outputlist:
        blankhere = tirshaker_returnblanksatstart(line)
        lbl = len(blankhere)
        if line.find('GR_PARMS=',lbl) == lbl:
            outputlist2 += ['#'+line]
            outputlist2.extend(block)
        elif line.find('GR_DEVICE=',lbl) == lbl:
            if gr_device == None:
                outputlist2 += [line]
            else:
                outputlist2 += ['#'+line]
                outputlist2 += ['GR_DEVICE= '+gr_device+'\n']
        else:
            outputlist2 += [line]
    return outputlist2

def replaceforplot(indef = None, outdef = None, parms = None, gr_parms = 'RADI VROT', smooth = False, gr_device = None):
    """
    Add lines to the .def file preparing it to plot error bars based on tirshaker errors. If smooth = True perform a Hanning smoothing for the parameters mentioned. Replace GR_DEVICE with the string given for gr_device unless gr_device == None.
    """

    # Open the file
    thestream = io.open(indef)
    lines = thestream.read().split('\n')
    thestream.close()

    # Add lines to the lines depending on finding the errors
    if smooth:
        newlines = replaceforplot_lines_smooth(lines = lines, gr_parms = gr_parms, parms = parms, gr_device = gr_device)
    else:
        newlines = replaceforplot_lines(lines = lines, gr_parms = gr_parms, gr_device = gr_device)

    outstring = ''
    for line in newlines:
        outstring = outstring+line+'\n'

    # Dump it into a file
    thething = io.open(outdef,'w')
    thething.write(outstring)
    thething.close()
