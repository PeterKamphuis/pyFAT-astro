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
    return SBR_initial,pa,inclination
'''







    
                                ;break if we do not want to use tirific
     IF finishafter EQ 0. then begin
        if setfinishafter NE 1 then begin
           comment = 'You have chosen to skip the fitting process after all preparations for the fit'
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
        ENDIF
        if optimized then begin
           tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
           currentfitcube=tmp[0]
        endif
                                ;       bookkeeping=0
        goto,finishthisgalaxy
     ENDIF
                                ;build up moment 0 axis
     ;buildaxii,headermap,xaxmom0,yaxmom0
                                ;some counters for keeping track
     prevmodification=0.
     overwrite=0.
     lastsbr=100
     counter=0.
     countsbr=0.
     restart_counter = 0
     continue_tirific = 'initialized'
     curr_run_id = 'stopped'
     case 1 of
       norings[0] LT 5: inimode = 1
       norings[0] LT 15: inimode = 2
       else: inimode = 3
     endcase
     ringmodifier=0
     plus2beamshift=0
     shiftcentercounter=0.
                                ;CD to the correct directory
     CD, new_dir, CURRENT=old_dir
     newrings=fix(norings[0])
                                ;Calculate the beam in pixels  and get an estimate for the peak sbr
     beaminpixels=fix(catmajbeam[i]/((ABS(pixelsizeRA)+ABS(pixelsizeDEC))/2.*3600.))
     centralarea=moment0map[fix(RApix[0]-beaminpixels):fix(RApix[0]+beaminpixels),fix(DECpix[0]-beaminpixels):fix(DECpix[0]+beaminpixels)]
     cenav=TOTAL(centralarea[WHERE(FINITE(centralarea))])/n_elements(centralarea[WHERE(FINITE(centralarea))])
     peaksbr=cenav*channelwidth/(catmajbeam[i]*catminbeam[i])
                                ; at high inclinations rings overlap
                                ; so the peaksbr needs to be reduced
                                ; according to the inclination and the
                                ; physical size of the galaxy
     IF catinc[i] GT 70. then begin
        IF convertskyanglefunction(DHI,catDistance[i]) GT 1. then $
           peaksbr=peaksbr/(convertskyanglefunction(DHI,catDistance[i])*SIN(catinc[i]*!DtoR)*0.5) else $
              peaksbr=peaksbr/(convertskyanglefunction(3.*catmajbeam[i],catDistance[i])*SIN(catinc[i]*!DtoR)*0.5)
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We estimate the peaksbr to be "+string(peaksbr)
        close,66
     ENDIF

                                ;Calculate an initial surface
                                ;brightness profile based on the
                                ;central brightness and an exponential
     velaverage=catmaxrot[i]
     signvelaverage=catmaxrot[i]
     stringsbr='SBR= '+string(peaksbr)+' '+string(peaksbr)
     scalelength=convertskyanglefunction(7.5,catDistance[i],/PHYSICAL)
     for j=1,norings[0]-2 do begin
        exposbr=peaksbr*exp(-(rings[j]-rings[0])/(scalelength))
        stringsbr=stringsbr+string(exposbr)+' '
     endfor
     tmppos=where('SBR' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringsbr
                                ;setting the values for the 2nd disk
     tmppos=where('SBR_2' EQ tirificfirstvars)
     tmp=str_sep(strtrim(strcompress(stringsbr),2),'=')
     tirificfirst[tmppos]='SBR_2='+tmp[1]
     tmp2=str_sep(strtrim(strcompress(tmp[1]),2),' ')


                                ;but that does mean we'd like
                                ;to change SBR a lot so we set quite
                                ;large starting steps
     string1=''
     string2=''
     string3=''
     string4=''
     string5=''
     string6=''
     string7=''

                                ;The minimum is based on the cutoff
                                ;values and therefore we need to set
                                ;every ring in the tirific file
     IF norings[0]*ring_spacing LT 7 then begin
        for j=norings[0],3,-1 do begin
           string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
           string2=string2+' 1'
           string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
           string4=string4+' 7.5E-5'
           string5=string5+' 2E-6'
           string6=string6+' 5E-6'
           string7=string7+' 3'
        endfor
     ENDIF else begin
        for j=norings[0],4,-1 do begin
           string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+', SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
           string2=string2+' 1 1'
           string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
           string4=string4+' 7.5E-5 7.5E-5'
           string5=string5+' 2E-6 2E-6'
           string6=string6+' 5E-6 5E-6'
           string7=string7+' 3 3'
        endfor
     ENDELSE
     string1 = STRMID(string1,1,STRLEN(string1)-1)
     SBRinput1=[string1,string2,string3,string4,string5,string6,string7]
     SBRinput2=[' SBR 1 2 3 SBR_2 1 2 3',$
                strtrim(strcompress(string(peaksbr,format='(E12.5)'))),'0','1E-5','1E-6','5E-5','1E-6','3','70','70']
                                ;Now setting some general values
     tmppos=where('INSET' EQ tirificfirstvars)
     tirificfirst[tmppos]='INSET='+strtrim(strcompress(string(currentfitcube+'.fits')))
     tmppos=where('BMIN' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMIN='+strtrim(strcompress(string(catminbeam[i])))
     tmppos=where('BMAJ' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMAJ='+strtrim(strcompress(string(catmajbeam[i])))
     tmppos=where('BPA' EQ tirificfirstvars)
     tirificfirst[tmppos]='BPA='+strtrim(strcompress(string(catbpa[i])))
     tmppos=where('SDIS' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS=  8.'
     tmppos=where('SDIS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS_2=  8.'
     tmppos=where('CONDISP' EQ tirificfirstvars)
     IF vresolution EQ 0. then begin
        tirificfirst[tmppos]='CONDISP=  '+string(channelwidth*1.2/(2*SQRT(2*ALOG(2.))))
     ENDIF ELSE BEGIN
        tirificfirst[tmppos]='CONDISP=  '+string((channelwidth*(1.+vresolution))/(2*SQRT(2*ALOG(2.))))
     ENDELSE
     tmppos=where('RMS' EQ tirificfirstvars)
     tirificfirst[tmppos]='RMS='+strtrim(strcompress(string(catnoise[i])))
     tmppos=where('DISTANCE' EQ tirificfirstvars)
     tirificfirst[tmppos]='DISTANCE='+strtrim(strcompress(string(catDistance[i])))
     prevrings=norings[0]
                                ; defining some flags that we will use
                                ; for checks on what is going on in
                                ; the fitting
                                ;the First estimate of Z0 is 0.2 kpc
                                ; if 0.2 kpc is less
                                ;than 1/4th of the beam we want to
                                ;make it a quarter of the beam at high inclinations
     If catdistance[i] EQ 1. then begin
        IF catinc[i] GT 80. then inpzval=MAX([(catmajbeam[i]*ring_spacing*norings[0])/150.,catmajbeam[i]/4.]) else  inpzval=(catmajbeam[i]*ring_spacing*norings[0])/150.
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string((inpzval))))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string((inpzval))))
     ENDIF ELSE BEGIN
        IF catinc[i] GT 80. then inpzval=MAX([convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL),catmajbeam[i]/4.]) else inpzval=convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL)
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string(inpzval)))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string(inpzval)))
     ENDELSE
     constring=0.
     forcedring=0.
     lastcutrings=0.
     lastaddrings=0.
     secondtime=0.
     sofiarings=norings[0]
     stringvelocities='VROT= 0. '+string(ABS(double(catmaxrot[i])))
     tmppos=where('CFLUX' EQ tirificfirstvars)
     IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX= 1e-5'
     tmppos=where('CFLUX_2' EQ tirificfirstvars)
     IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX_2= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX_2= 1e-5'
                                ;If we change the amount of rings we need to come back here
     sbrshift:
     case finishafter of
        2.1:begin
           tmpring=norings[0]-10.
           rings=dblarr(norings[0])
           rings[0:9]=(findgen(10))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
           rings[10:norings[0]-1]=(findgen(fix(tmpring)))*catmajbeam[i]*ring_spacing*2+catmajbeam[i]/5.+11.*catmajbeam[i]*ring_spacing
        end
       else:rings=(findgen(norings[0]))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
     endcase


     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
        close,66
     ENDIF ELSE begin
        print,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
     endelse
     fluxadjust=0.
                                ;let's write the number of rings to the tirific file
     tmppos=where('NUR' EQ tirificfirstvars)
     tirificfirst[tmppos]='NUR='+strtrim(strcompress(string(norings[0])),1)
                                ;and cflux


                                ;Now we write the radii of our rings
     stringring='RADI=0.0 '
     for j=0,norings[0]-2 do begin
        stringring=stringring+string(rings[j])+' '
     endfor
     tmppos=where('RADI' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringring
                                ;Using the parameters from sofia or a
                                ;previous fit
     tmppos=where('VROT' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringvelocities
     tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
     tmppos=where('VROT_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VROT_2='+tmp[1]
     tmppos=where('INCL' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
     tmppos=where('PA' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA='+strtrim(strcompress(string(catPA[i])))
     tmppos=where('XPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS='+strtrim(strcompress(string(catvsys[i])))
     tmppos=where('INCL_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
     tmppos=where('PA_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA_2='+strtrim(strcompress(string(catPA[i])))
     tmppos=where('XPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS_2='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS_2='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS_2='+strtrim(strcompress(string(catvsys[i])))

                                ; Set the tirific fitting parameters for the inclination
                                ;If the inclination is above 75 the
                                ;initial estimates must be good so we
                                ;want less change
     case 1 of
        catinc[i] LT 40: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.2),string(0.1),string(5.0),string(0.1),'5','70','70']
           ;We put this inside might not work then take it out again
           IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
           IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
        end
        catinc[i] LT 75: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(5.),string(0.1),string(1.0),string(0.1),'3','70','70']

           IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
           IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
        end
        else:INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                         string(catinc[i]+catincdev[i]),string(catinc[i]-catincdev[i]),string(1.0),string(0.01),string(1.0),string(0.001),'3','70','70']
     endcase
                                ; have to ensure that the parameters are within the limits
     IF INCLinput1[1] GT 90. then INCLinput1[1]='90.'
     IF INCLinput1[1] LT 30.then INCLinput1[1]='30.'

                                ;Set the input for the PA

     PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string((catPA[i])+catPAdev[i]+10.),string((catPA[i])-catPAdev[i]-10),string(1.),string(0.1),string(1.),string(0.1),'3','70','70']

                                ;If the estimate PA isn't close
                                ;to 180 or 360 we want to block
                                ;changes that can flip the rotation curve
     IF PAfixboun NE 'n' then begin
        IF PAfixboun EQ 'w' then begin
           IF catPA[i]+catPAdev[i]+10. GT 360 then PAinput1[1]='360'
           IF catPA[i]-catPAdev[i]-10. LT 180 then PAinput1[2]='180'
        ENDIF
        IF PAfixboun EQ 'e' then begin
           IF catPA[i]+catPAdev[i]+10. GT 180. then PAinput1[1]='180'
           IF catPA[i]-catPAdev[i]-10. LT 0. then PAinput1[2]='0'
        ENDIF
     ENDIF
                                ;Also the minimum and maximum of Z0 which have to be based on physical
                                ;values
     Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.1,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(1.0,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),'3','70','70']
                                ;if the inclination is high we want to
                                ;make sure that the scale height can
                                ;be at least half a beam or 1 kpc
                                ;whichever is larger because otherwise
                                ;it might be pushing the inclination
                                ;down.

     IF catinc[i] GT 80 then begin
        Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
     ENDIF
     IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([catmajbeam[i]*norings[0]*ring_spacing/66.,catmajbeam[i]/2.]))
        ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*ring_spacing*norings[0]/66.)
        Z0input1[2]='0.'
        Z0input1[3]=string(-1*catmajbeam[i]*ring_spacing*norings[0]/1.5E5)
        Z0input1[4]=string(catmajbeam[i]*ring_spacing*norings[0]/1.5E6)
        Z0input1[5]=string(catmajbeam[i]*ring_spacing*norings[0]/1500)
        Z0input1[6]=string(catmajbeam[i]*ring_spacing*norings[0]/1.5E6)
     ENDIF





                                ;Now The rotation

     VROTmax=catmaxrot[i]+catmaxrotdev[i]
     VROTmin=channelwidth
                                ;ensure reasonable vrotmax dependent
                                ;on inclination as more unsure at low inclination
     IF VROTmax LT 80. then VROTmax=80.
     IF VROTmax GT 600. then VROTmax=600.

     string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
     string2=string(VROTmax)
     string3=string(VROTmin)
     string4=string(channelwidth)
                                ;At low inclinations we want to be
                                ;very accepting of the rotational
                                ;values as they are very unsure.
     case 1 of
        catinc[i] LT 30:begin
           string4=string(2*channelwidth)
           string5=string(0.02*channelwidth)
           string6=string(0.2*channelwidth)
           string7=string(0.02*channelwidth)
        end
        else:begin
           string4=string(channelwidth)
           string5=string(0.01*channelwidth)
           string6=string(0.1*channelwidth)
           string7=string(0.01*channelwidth)
        end
     ENDCASE
     string8='3'
     string9='70'
                                ;If the model is large enough we fit
                                ;only a slope to the outer quarter of
                                ;the model at low inclinations we want
                                ;this slope to be at least half of the
                                ;rings and be flat
     IF centralexclude then begin
        IF catinc[i] LT 40 then begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDIF ELSE BEGIN
        IF catinc[i] LT 40 then begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDELSE
     IF norings[0] GT 4 then begin
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
     ENDIF ELSE BEGIN
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
     ENDELSE
     vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                  strtrim(strcompress(string(catvsys[i]+100.)),1),strtrim(strcompress(string(catvsys[i]-100.)),1),'10','0.01','0.5','0.01','3','70','70']

                                ;IF we have a small cube we will only accept very small changes to the
                                ;central position
     IF norings[0] LT 4 then begin
        maxxpos=strcompress(string(RAdeg+catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-catmajbeam[i]/3600.))
        ;We'll be satified with a change less than a tenth of a pixel
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), +$
                      maxxpos,minxpos,strtrim(string(pixelsizeRA)),strtrim(string(pixelsizeRA/20.)),+$
                      strtrim(string(pixelsizeRA/10.)),strtrim(string(pixelsizeRA/20.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                     maxypos,minypos,strtrim(string(pixelsizeDEC)),strtrim(string(pixelsizeDEC/20.)),+$
                      strtrim(string(pixelsizeDEC/10.)),strtrim(string(pixelsizeDEC/20.)),'3','70','70']
     ENDIF ELSE BEGIN
        maxxpos=strcompress(string(RAdeg+3.*catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-3.*catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+3.*catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-3.*catmajbeam[i]/3600.))
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), +$
                      maxxpos,minxpos,strtrim(string(pixelsizeRA*2.)),strtrim(string(pixelsizeRA/10.)),+$
                      strtrim(string(pixelsizeRA/2.)),strtrim(string(pixelsizeRA/10.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                     maxypos,minypos,strtrim(string(pixelsizeDEC*2.)),strtrim(string(pixelsizeDEC/10.)),+$
                      strtrim(string(pixelsizeDEC/2.)),strtrim(string(pixelsizeDEC/10.)),'3','70','70']
     ENDELSE
     tryone=0.
                                ;   IF we have a low inclination we first want to fit the PA by itself
     PAest=0.
     INCLest=0.
     INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                 ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                 string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(2),string(0.1),string(0.1),'3']
     IF (catinc[i]  LT 50 OR mismatch EQ 1 OR ceil(norings[0]*COS(catinc[i]*!DtoR))*ring_spacing LE 5) AND counter EQ 0. then begin
                                ;If we have a small number of beams
                                ;acros the minor axis the inclination
                                ;is very unsure and we first want to
                                ;fit the inclination. If the KinPA and VelPA do not match we want to fit the PA
        PAinincl=catPA[i]
        INCLinincl=catinc[i]
        maxrotinincl=catmaxrot[i]
        fixstring='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        IF ceil(norings[0]*COS(catinc[i]*!DtoR))*ring_spacing LE 5 OR catinc[i] LT 30. OR mismatch EQ 1 then begin

                                ; Avoid runaway rotation
           IF mismatch EQ 1 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Because of a difference between the kinematic PA and the morphological PA we will scramble the inclination."
                 printf,66,linenumber()+"Original inclination ="+string(catinc[i])+" New = "+string(catinc[i]-5.)
                 close,66
              ENDIF
              catinc[i]=catinc[i]-5.
              catincdev[i]=catincdev[i]+10.
           ENDIF ELSE BEGIN
              catinc[i]=catinc[i]-3.
              catincdev[i]=catincdev[i]+5.
           ENDELSE
           IF catinc[i] LT 5 then catinc[i]=5.

           tmppos=where('INCL' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
           tmppos=where('INCL_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
           TMP_VEL=W50/2./SIN(ABS(INCLinincl)*!pi/180.)
           IF TMP_VEL GT 300. then TMP_VEL= 300.
           tmppos=where('VROT' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT= 0. '+strtrim(strcompress(string(tmp_vel)))
           tmppos=where('VROT_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT_2= 0. '+strtrim(strcompress(string(tmp_vel)))
           INCLest=1
           ;INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
           ;             ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
           ;             '90.',string(catinc[i]-catincdev[i]-5),string(1.),string(0.1),string(0.5),string(0.05),'5','70','70']
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                        ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                        '90.',string(catinc[i]-catincdev[i]-5),string(5.),string(0.5),string(0.1),'5']
           ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
           ;           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
           ;           string((catPA[i])+catPAdev[i]+40.),string((catPA[i])-catPAdev[i]-40),string(6),string(0.1),string(1.0),string(0.01),'3','70','70']
           IF mismatch EQ 1 THEN BEGIN

             PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                      string((catPA[i])+catPAdev[i]+40.),string((catPA[i])-catPAdev[i]-40),string(10),string(1.),string(0.5),'3']
           ENDIF ELSE BEGIN
             PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                      string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(1.),string(0.25),'3']
           ENDELSE

           IF INCLinput1[2] LT 2. then INCLinput1[2]='2'

           IF PAfixboun NE 'n' then begin
              IF PAfixboun EQ 'w' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
                 IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
              ENDIF
              IF PAfixboun EQ 'e' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
                 IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
              ENDIF
           ENDIF

           ;VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
          ;                ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
          ;                ,string(VROTmax),string(0.),$
          ;                string(channelwidth),string(0.1*channelwidth),string(0.5*channelwidth),string(0.01*channelwidth),'3','70','70',fixstring]
           VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
                          ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
                          ,string(VROTmax),string(0.), string(4*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth),'5',fixstring]
           Writefittingvariables,tirificfirst,inclinput1,VROTinputINCL,painput1,sbrinput1,sbrinput2,INIMODE=inimode

           againINCLestimate:
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because there are only a few beams across the minor axis or we have low inclination we first adjust the inclination "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           restart_counter++
           tmppos=where('RESTARTID' EQ tirificfirstvars)
           tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))

           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificfirst)-1 do begin
              printf,1,tirificfirst[index]
           endfor
           close,1

           IF testing GE 1 then goto,testing1INCL
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           print,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()


           run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name='1stfit.', $
                  store_name='1stfitold.', log=log, loops=loops,nopoints=nopoints, $
                  AC=AC1,toymodels=toymodels,run_unit=run_unit
           IF bookkeeping EQ 5 THEN goto,finishthisgalaxy

           ;If not accepted we try again
           ;IF AC1 EQ 0. and INCLest LT 1 then begin
            ;  IF size(log,/TYPE) EQ 7 then begin
            ;     openu,66,log,/APPEND
            ;     printf,66,linenumber()+"The INCL estimate was not accepted try again. Tries: "+string(INCLest)
            ;     close,66
            ;  ENDIF
            ;  VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2','SBR','SBR_2']
            ;  firstfitvalues=0.
            ;  writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
            ;  catPA[i]=firstfitvalues[0,0]
            ;  catinc[i]=firstfitvalues[0,2]
            ;  catmaxrot[i]=firstfitvalues[1,4]
            ;  IF catinc[i] LT 5 then catinc[i]=5.
            ;  PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(1.),string(PAinput1[5]*(1+INCLest*0.25)),'5']
            ;  ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;  ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;  ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
            ;  IF PAfixboun NE 'n' then begin
            ;     IF PAfixboun EQ 'w' then begin
            ;        IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
            ;        IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
            ;     ENDIF
            ;     IF PAfixboun EQ 'e' then begin
            ;        IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
            ;        IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
            ;     ENDIF
            ;  ENDIF
            ;  INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;            '90.',string(catinc[i]-catincdev[i]-5),string(5.),string(0.5),string(0.1*(1.+INCLest*0.5)),'5']
            ;  VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
            ;              ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
            ;              ,string(VROTmax),string(0.), string(4*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth*(1.+INCLest*0.5)),'5',fixstring]
              ;INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
              ;            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
              ;            '90.',string(catinc[i]-catincdev[i]-20),string(1.),string(0.5),string(0.1),string(0.01),'3','70','70']
            ;  Writefittingvariables,tirificfirst,painput1,inclinput1,vrotinputINCL,INIMODE=inimode


            ;  INCLest=INCLest+1.
            ;  goto,againINCLestimate
           ;ENDIF
           ;ELSE BEGIN

              ;When accepted update our values
            ;  VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2','SBR','SBR_2']
            ;  firstfitvalues=0.
          ;    writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted

          ;    IF size(log,/TYPE) EQ 7 then begin
          ;       openu,66,log,/APPEND
          ;       printf,66,linenumber()+"Because of low inclination (<30) or small number of beams across minor axis (< 5):"
          ;       printf,66,linenumber()+"We have adjusted the PA from "+string(PAinincl)+" to "+string(firstfitvalues[0,0])+'+/-'+string(newPA[1])
          ;       printf,66,linenumber()+"The inclination from "+string(INCLinincl)+" to "+string(firstfitvalues[0,2])+'+/-'+string(newinclination[1])
          ;       printf,66,linenumber()+"The max rotation from "+string(maxrotinincl)+" to "+string(firstfitvalues[n_elements(firstfitvalues[*,0])-2,4])+'+/-'+string(catmaxrotdev[1])
          ;       printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT."
          ;       close,66
          ;    ENDIF
          ;    catPA[i]=firstfitvalues[0,0]
          ;    catinc[i]=firstfitvalues[0,2]
          ;    IF catinc[i] LT 5 then catinc[i]=5.
          ;    catmaxrot[i]=firstfitvalues[1,4]


           ;ENDELSE
        ENDIF
        testing1INCL:


        ;We only want to check the PA if we have not checked the INCL
        IF INCLest EQ 0 THEN BEGIN
          ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
          ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
          ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(2.5),string(0.1),string(0.05),string(0.01),'3','70','70']
          PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.5),string(0.1),'3']
          IF PAfixboun NE 'n' then begin
             IF PAfixboun EQ 'w' then begin
                IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
             ENDIF
             IF PAfixboun EQ 'e' then begin
                IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
             ENDIF
          ENDIF
          ;VROTinputPA=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
          ;             ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
          ;             ,string(VROTmax),string(0.),$
          ;             string(channelwidth),string(0.01*channelwidth),string(channelwidth),string(0.01*channelwidth),'3','70','70',fixstring]
          VROTinputPA=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
                       ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
                       ,string(VROTmax),string(0.), string(3*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth),'5',fixstring]
          Writefittingvariables,tirificfirst,painput1,sbrinput1,sbrinput2,VROTinputPA,INIMODE=inimode


          againPAestimate:
          restart_counter++
          tmppos=where('RESTARTID' EQ tirificfirstvars)
          tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))
          openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
          for index=0,n_elements(tirificfirst)-1 do begin
             printf,1,tirificfirst[index]
          endfor
          close,1

          IF testing GE 1 then goto,testing1PA
          IF size(log,/TYPE) EQ 7 then begin
             openu,66,log,/APPEND
             printf,66,linenumber()+"Starting tirific the PA estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
             close,66
          ENDIF
          print,linenumber()+"Starting tirific the PA estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()

          run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name='1stfit.', $
                    store_name='1stfitold.', log=log,loops=loops,nopoints=nopoints,AC=AC1, $
                    toymodels=toymodels,run_unit=run_unit
          IF bookkeeping EQ 5 THEN goto,finishthisgalaxy


'''

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
