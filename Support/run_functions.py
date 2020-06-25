#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to run external programs

class BadSourceError(Exception):
    pass
class SofiaRunError(Exception):
    pass
class TirificRunError(Exception):
    pass

from support_functions import print_log, convert_type,set_limits,sbr_limits
from clean_functions import clean_before_sofia,clean_after_sofia
from fits_functions import cut_cubes,extract_pv
from astropy.wcs import WCS
from astropy.io import fits

import read_functions as rf
import write_functions as wf
import os
import subprocess
import numpy as np
import traceback



def central_converge(Configuration, Fits_Files,Tirific_Template, Catalogue):
    return 1

def check_source(Configuration, Fits_Files, Catalogue, header):


    name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa, \
        w50,err_f_sum, err_x,err_y,err_z= rf.sofia_catalogue(Configuration,header=header)

    x_min,x_max,y_min,y_max,z_min,z_max = convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int')
    x,y,z,ra,dec,v_app,f_sum,kin_pa,f_sum_err , err_x,err_y,err_z= convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa,err_f_sum, err_x,err_y,err_z])
    #How does sofia 2 deal with the fully flagged channels?
    # Need to check for that here if NaNs are included
    if f_sum < 0.:
            print_log(f'''CHECK_SOURCE: This galaxy has negative total flux. That will not work. Aborting.
    ''',Configuration['OUTPUTLOG'])
            raise BadSourceError('We found an initial negative total flux.')
    galaxy_box = [[z_min,z_max],[y_min,y_max],[x_min,x_max]]

    # If the provided distance  = -1 we assume a Hubble follow
    if Catalogue['DISTANCE'] == -1:
        Catalogue['DISTANCE'] == v_app/(1000.*H_0)
    if Catalogue['DISTANCE'] < 0.5:
        Catalogue['DISTANCE'] == 0.5


    #Check whether the cube is very large, if so cut it down

    new_box = cut_cubes(Configuration, Fits_Files, galaxy_box, header)
    #update our pixel values to match the new sizes
    for i in range(len(new_box)):
        shift = new_box[i,0]
        if i == 0:
            z -= shift; z_min -= shift; z_max -= shift
        elif i == 1:
            y -= shift; y_min -= shift; y_max -= shift
        elif i == 2:
            x -= shift; x_min -= shift; x_max -= shift

    Cube = fits.open(Configuration['FITTING_DIR']+Fits_Files['FITTING_CUBE'],uint = False, do_not_scale_image_data=True,ignore_blank = True, output_verify= 'ignore')
    data = Cube[0].data
    header = Cube[0].header

    Configuration['CHANNEL_WIDTH'] = header['CDELT3']/1000.
    cube_wcs = WCS(header)
    # convert the boundaries to real coordinates
    ralow,declow,vellow = cube_wcs.wcs_pix2world(x_min,y_min,z_min,0.)
    rahigh,dechigh,velhigh = cube_wcs.wcs_pix2world(x_max,y_max,z_max,0.)
    DECboun = np.sort([float(declow),float(dechigh)])
    RAboun = np.sort([float(ralow),float(rahigh)])
    VELboun = np.sort([float(vellow),float(velhigh)])
    # We write the results of the cut cube to the log
    print_log(f'''CHECK_SOURCE: The source finder found the following center in pixels.
{"":8s}CHECK_SOURCE: RA center = {x} with boundaries {x_min}, {x_max}
{"":8s}CHECK_SOURCE: DEC center = {y} with boundaries {y_min}, {y_max}
{"":8s}CHECK_SOURCE: V_sys center = {z} with boundaries {z_min}, {z_max}
{"":8s}CHECK_SOURCE: This should correspond to the WCS coordinates:
{"":8s}CHECK_SOURCE: RA center = {ra} with boundaries {','.join(convert_type(RAboun,type='str'))}
{"":8s}CHECK_SOURCE: DEC center = {dec} with boundaries {','.join(convert_type(DECboun,type='str'))}
{"":8s}CHECK_SOURCE: V_sys center = {v_app/1000.:.2f} with boundaries {','.join(convert_type(VELboun/1000.,type='str'))}
''', Configuration['OUTPUTLOG'])

    #get the maximum amount of rings
    ringbuffer = set_limits(round(3./Configuration['RING_SIZE']),2,6)
    Configuration['MAX_RINGS'] = int(round(np.sqrt(((x_max-x_min)/2.)**2+((y_max-y_min)/2.)**2) \
                /((header['BMAJ']*Configuration['RING_SIZE'])/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])]))\
                +ringbuffer))

    # Determine whether the centre is blanked or not
    Central_Flux = data[int(round(z)),int(round(y)),int(round(x))]
    print_log(f'''CHECK_SOURCE: In the center we find {Central_Flux} Jy/beam at the location:
{"":8s}CHECK_SOURCE: x,y,z = {int(round(x))}, {int(round(y))}, {int(round(z))}.
''',Configuration['OUTPUTLOG'])

    if not np.isfinite(Central_Flux):
        Configuration['EXCLUDE_CENTRAL'] = True
        print_log(f'''CHECK_SOURCE: The flux in the central part is blanked. We exclude the central rings.
''',Configuration['OUTPUTLOG'])
    else:
        Configuration['EXCLUDE_CENTRAL'] = False

    #Check that the source is bright enough
    Max_SNR = np.nanmax(data)/Configuration['NOISE']
    if Max_SNR < 2.5:
        log_statement = f'''CHECK_SOURCE: The maximum Signal to Noise in this cube is {Max_SNR} that is not enough for a fit.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'])
        raise BadSourceError(log_statement)
    #Let's get the initial estimates of the PA and inclination from the axis ratios

    try:
        pa, inclination, SBR_initial, maj_extent = rf.guess_orientation(Configuration,Fits_Files, center = [x,y])
    except:
        print_log(f'''CHECK_SOURCE: We could not establish proper initial estimates from the moment maps.
{"":8s} {traceback.print_exc()}
''', Configuration['OUTPUTLOG'], screen = True)
        raise BadSourceError("No initial estimates. Likely the source is too faint.")
    if abs(pa[0]-kin_pa) < 25:
        pa[0] = (pa[0]/pa[1]+kin_pa/10.)/(1./pa[1]+1./10.)
    else:
        if np.isfinite(kin_pa):
            pa[0] = kin_pa

    # Let's see how many ring we need. #The +1 accounts for the inner 1./5. of the beam.
    Configuration['NO_RINGS'] = int(round(maj_extent/(header['BMAJ']*Configuration['RING_SIZE'])+1.))

    print_log(f'''CHECK_SOURCE: From the original Configuration and SoFiA we find:
{"":8s}CHECK_SOURCE: The Maximum amount of rings possible {Configuration['MAX_RINGS']}
{"":8s}CHECK_SOURCE: We start with {Configuration['NO_RINGS']} Rings in the model.
{"":8s}CHECK_SOURCE: SoFiA found a PA of {kin_pa:.2f} and we use a PA = {pa[0]:.2f} +/- {pa[1]:.2f}
{"":8s}CHECK_SOURCE: We start with an inclination of {inclination[0]:.2f} +/- {inclination[1]:.2f}{"":8s}
{"":8s}CHECK_SOURCE: SoFiA found a W50 of {w50/1000.:.2f} km/s
{"":8s}CHECK_SOURCE: We will now check whether the size of the galaxy is sufficient.
''',Configuration['OUTPUTLOG'])

    initial_ringsize = Configuration['RING_SIZE']
    while Configuration['RING_SIZE'] > 0.5 and  Configuration['NO_RINGS'] <= 4.:
        previous_ringsize = Configuration['RING_SIZE']
        Configuration['RING_SIZE'] = set_limits(Configuration['RING_SIZE']/1.5,0.5,float('NaN'))
        if  Configuration['NO_RINGS'] <= 2.:
            Configuration['NO_RINGS'] = int(round(Configuration['NO_RINGS']*previous_ringsize/Configuration['RING_SIZE']))
        else:
            Configuration['NO_RINGS'] = int(round(Configuration['NO_RINGS']*previous_ringsize/Configuration['RING_SIZE']-0.66))
        print_log(f'''CHECK_SOURCE: Because we had less than four rings we have reduced the ring size from {previous_ringsize} to {Configuration['RING_SIZE']}
''',Configuration['OUTPUTLOG'])

    if Configuration['NO_RINGS'] < 3:
        print_log(f'''CHECK_SOURCE: With a ring size of {Configuration['RING_SIZE']} we still only find {Configuration['NO_RINGS']}.
{"":8s}CHECK_SOURCE: This is not enough for a fit.
''',Configuration['OUTPUTLOG'])
        raise BadSourceError('The extracted source is too small to fit.')

    if initial_ringsize != Configuration['RING_SIZE']:
        Configuration['MAX_RINGS'] = Configuration['MAX_RINGS']*initial_ringsize/Configuration['RING_SIZE']
    #For galaxies smaller than a certain size we do not fit a warp
    if Configuration['NO_RINGS']*Configuration['RING_SIZE'] < Configuration['MINIMUM_WARP_SIZE']/2.:
        print_log(f'''CHECK_SOURCE: We force the fit of a flat disk because the estimated size of {Configuration['RING_SIZE']*Configuration['NO_RINGS']*2.} is less than {Configuration['MINIMUM_WARP_SIZE']}.
''',Configuration['OUTPUTLOG'])
        Configuration['FIX_INCLINATION'] = True
        Configuration['FIX_PA'] = True
        Configuration['FIX_SDIS'] = True
    # Make sure we did not get over extended
    if Configuration['NO_RINGS'] > Configuration['MAX_RINGS']:
        print_log(f'''CHECK_SOURCE: As we had more rings (# Rings = {Configuration['NO_RINGS']}) than the maximum ({Configuration['MAX_RINGS']}) we limit the amount of rings to the maximum.
''',Configuration['OUTPUTLOG'])
        Configuration['NO_RINGS'] = Configuration['MAX_RINGS']

    if Configuration['NO_RINGS'] > 20 and Configuration['MAX_RINGS'] > 25:
        Configuration['OUTER_RINGS_DOUBLED'] = True
        print_log(f'''CHECK_SOURCE: This is a large galaxy (# Rings = {Configuration['NO_RINGS']}) Therefore we use twice the ring size in the outer parts.
''',Configuration['OUTPUTLOG'])



    if Configuration['HANNING']:
        vres = header['CDELT3']*2.
    else:
        vres = header['CDELT3']

    if inclination[0] < 40:
        max_vrot=w50/2./np.sin(np.radians(abs(inclination[0]+5.)))
    else:
        max_vrot=w50/2./np.sin(np.radians(abs(inclination[0])))
    max_vrot_dev=set_limits((VELboun[1]-VELboun[0])/4./np.sin(np.radians(inclination[0])),4.*vres,20*vres)
    #Write the info to the Basic info File
    wf.basicinfo(Configuration,initialize = True,
              RA=[ra,abs(err_x*header['CDELT1'])],
              DEC=[dec,abs(err_y*header['CDELT2'])],
              VSYS =[v_app,abs(err_z*header['CDELT3'])],
              PA=pa, Inclination = inclination, Max_Vrot = [max_vrot,max_vrot_dev], Tot_Flux = [f_sum,f_sum_err], V_mask = [VELboun[1]-VELboun[0],vres],
              Distance = Configuration['DISTANCE'] , DHI = maj_extent*3600.)

    if Configuration['MAPS_OUTPUT'] < 4:
        # extract a PV-Diagram
        print(z_max,z_min)
        PV =extract_pv(Cube, pa[0], center=[ra,dec,v_app], convert = 1000.,
                       finalsize = [int(round(maj_extent/np.mean([abs(header['CDELT1']),abs(header['CDELT2'])])*1.25+header['NAXIS1']*0.2)),
                                    int(round(z_max-z_min)+10.)])
        if not os.path.isdir(Configuration['FITTING_DIR']+'/PV-Diagrams'):
                os.mkdir(Configuration['FITTING_DIR']+'/PV-Diagrams')
        fits.writeto(Configuration['FITTING_DIR']+'/PV-Diagrams/'+Configuration['BASE_NAME']+'_sofia_xv.fits',PV[0].data,PV[0].header)
    Cube.close()
    Initial_Parameters = {}
    Initial_Parameters['RA'] = [ra,abs(err_x*header['CDELT1'])]
    Initial_Parameters['DEC'] = [dec,abs(err_y*header['CDELT2'])]
    Initial_Parameters['VSYS'] =[v_app,abs(err_z*header['CDELT3'])]
    Initial_Parameters['SBR'] = SBR_initial
    Initial_Parameters['VROT'] = [max_vrot/1000.,max_vrot_dev/1000.]
    Initial_Parameters['PA'] = pa
    Initial_Parameters['INCL'] = inclination

    return Initial_Parameters

check_source.__doc__='''

; NAME:
;       check_source(Configuration, Fits_Files, Tirific_Template, Catalogue)
;
; PURPOSE:
;       Check that the source found by sofia can be fitted and the cube is suitable for fitting
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
;      set_limits, print_log
;
; EXAMPLE:
;
;
'''


def extent_converge(Configuration, Fits_Files,Tirific_Template, Catalogue):
    return 1

def sofia(Configuration, Fits_Files, hdr, supportdir):

    sofia_template = rf.sofia_template(supportdir+'/sofia_template.par')
    if not os.path.isdir(Configuration['FITTING_DIR']+'/Sofia_Output'):
        os.mkdir(Configuration['FITTING_DIR']+'/Sofia_Output')
    os.chdir(Configuration['FITTING_DIR'])
    threshold = 5.
    counter = 3
    sofia_template['input.data'] = Fits_Files['FITTING_CUBE']
    beam_in_pixels = int(round(hdr['BMAJ']/((abs(hdr['CDELT1'])+abs(hdr['CDELT2']))/2.)))
    spatial_kernels = [0,beam_in_pixels,beam_in_pixels*2]
    if (hdr['NAXIS1']+hdr['NAXIS2'])/(2.*beam_in_pixels)  > 30:
        spatial_kernels.append(beam_in_pixels*3)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube is more than 30 beams across.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'])

    velocity_kernels = [0,3,6,12]
    if hdr['NAXIS3'] > 52:
        velocity_kernels.append(16)
        log_statement=f'''RUN_SOFIA: Adding an extra kernel scale as the cube has more than 52 channels.
'''
        print_log(log_statement, Configuration['OUTPUTLOG'])
    sofia_template['scfind.kernelsXY'] = ','.join([str(x) for x in spatial_kernels])
    sofia_template['scfind.kernelsZ'] = ','.join([str(x) for x in velocity_kernels])
    sofia_ok = False
    while not sofia_ok:
        clean_before_sofia(Configuration)
        sofia_template['scfind.threshold'] = str(threshold)
        wf.sofia(sofia_template,'sofia_input.par')
        print_log("RUN_SOFIA: Running SoFiA.",Configuration['OUTPUTLOG'],screen = True)
        sfrun = subprocess.Popen(["sofia2",'sofia_input.par'], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
        sofia_run, sofia_warnings_are_annoying = sfrun.communicate()
        print_log(sofia_run.decode("utf-8"), Configuration['OUTPUTLOG'])

        if sfrun.returncode == 8:
            if threshold > 3.:
                log_statement = f'''RUN_SOFIA: We did not find a source at a threshold of {threshold}
{"":8s} RUN_SOFIA: Lowering the threshold and trying again."
'''
                print_log(log_statement,Configuration['OUTPUTLOG'])
                threshold -= 1
            else:
                clean_after_sofia(Configuration)
                log_statement = f'''RUN_SOFIA: We did not find a source above a threshold of {threshold}.
{"":8s}RUN_SOFIA: We cannot lower the threshold lower as the risk of fitting noise becomes too high.
{"":8s}Continuing to the next galaxy.
'''
                print_log(log_statement,Configuration['OUTPUTLOG'])
                raise SofiaRunError("RUN_SOFIA:Sofia cannot find a source above a threshold of 3.")
        elif sfrun.returncode == 0:
            sofia_ok = True
        else:
            print_log(sofia_warnings_are_annoying.decode("utf-8"), Configuration['OUTPUTLOG'])
            raise SofiaRunError("RUN_SOFIA:Sofia did not execute properly. See log for details")

    #Move sofia output to the desired Directory
    clean_after_sofia(Configuration)

    os.chdir(Configuration['START_DIR'])


sofia.__doc__ ='''
;+
; NAME:
;       SOFIA
;
; PURPOSE:
;       Run SoFiA to create a mask Moment maps and a catalogue
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;        sofia(Configuration, Fits_Files)
;       RUN_SOFIA,allnew,new_dir,currentfitcube,catcatalogname,supportdirchecked,pixfwhm,header,errormessage,VSYSpix,RApix,DECpix,Totflux,log=log
;
;
; INPUTS:
;      Configuration, Fits_files hdr and supportdir
; OPTIONAL INPUTS:
;       LOG = name of the tracing log
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;       FILE_TEST(),READ_FITS(),READ_TEMPLATE,SPAWN
;
'''
