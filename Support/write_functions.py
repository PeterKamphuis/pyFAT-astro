#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to write text files to Disk

from support_functions import convert_type, print_log,convertRADEC,convertskyangle,set_limit_modifier
from modify_template import set_model_parameters, set_overall_parameters, set_fitting_parameters
import numpy as np

# create or append to the basic ifo file
def basicinfo(Configuration,initialize = False,first_fit = False, second_fit = False,
              RA=[0.,0.], DEC=[0.,0.], VSYS =[0.,0.],
              PA=[0,0], Inclination = [0,0], Max_Vrot = [0,0], Tot_Flux = [0,0], V_mask = [0,0],
              Distance = 0 , DHI = 0):
    if initialize:
        with open(Configuration['FITTING_DIR']+Configuration['BASE_NAME']+'basic_info.txt','w') as file:
            file.write(f'''#This file contains the basic parameters of the Galaxy
#The total flux is determined by the source finder in the initial estimates and the total emission in the masked cube of the final model in the fits
#D_HI is determined as the diameter where the major axis with given PA cuts the 10^20 column of the moment0 map
# {'RA':>25s} {'DEC':>25s} {'VSYS':>20s} {'PA':>20s} {'Inclination':>20s} {'Max VRot':>20s} {'V_mask':>20s} {'Tot FLux':>20s} {'D_HI':>20s} {'Distance':>20s} {'HI Mass':>20s} {'D_HI':>20s}
# {'hh:mm:ss':>25s} {'dd:mm:ss':>25s} {'km/s':>20s} {'deg':>20s} {'deg':>20s} {'km/s':>20s} {'km/s':>20s} {'Jy':>20s} {'arcsec':>20s} {'Mpc':>20s} {'M_solar':>20s} {'kpc':>20s}
# The initial input
''')
    with open(Configuration['FITTING_DIR']+Configuration['BASE_NAME']+'basic_info.txt','a') as file:
        if first_fit:
            file.write(f'''#After the center converged. \n''')
        elif second_fit:
            file.write(f'''#After the radii converged. \n''')
        RAhr,DEChr = convertRADEC(RA[0],DEC[0])
        RA_c = f'{RAhr}+/-{RA[1]*3600.:0.2f}'
        DEC_c = f'{DEChr}+/-{DEC[1]*3600.:0.2f}'
        VSYS_c = f'{VSYS[0]/1000.:.2f}+/-{VSYS[1]/1000.:.2f}'
        PA_c = f'{PA[0]:.2f}+/-{PA[1]:.2f}'
        INCL_c = f'{Inclination[0]:.2f}+/-{Inclination[1]:.2f}'
        MVROT_c = f'{Max_Vrot[0]/1000.:.2f}+/-{Max_Vrot[1]/1000.:.2f}'
        Vmask_c = f'{V_mask[0]/1000.:.2f}+/-{V_mask[1]/1000.:.2f}'
        DHI_a = f'{DHI:.2f}'
        Dist = f'{Distance:.2f}'
        HIMass  = f'{Tot_Flux[0]*2.36E5*Distance**2:.2e}'
        DHI_k = f'{convertskyangle(DHI,Distance):.2f}'
        Flux_c = f'{Tot_Flux[0]:.2f}+/-{Tot_Flux[1]:.2f}'
        file.write(f'''  {RA_c:>25s} {DEC_c:>25s} {VSYS_c:>20s} {PA_c:>20s} {INCL_c:>20s} {MVROT_c:>20s} {Vmask_c:>20s} {Flux_c:>20s} {DHI_a:>20s} {Dist:>20s} {HIMass:>20s} {DHI_k:>20s}
''')
basicinfo.__doc__ = '''

; NAME:
;       basicinfo(Configuration,initialize = True,
              RA=[0.,0.], DEC=[0.,0.], VSYS =[0.,0.],
              PA=[0,0], Inclination = [0,0], Max_Vrot = [0,0], Tot_Flux = [0,0], V_mask = [0,0],
              Distance = [0] , DHI = [0])
;
; PURPOSE:
;       write to the basic info file
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
;      split, strip, open
;
; EXAMPLE:
;
;
'''

# Function to write the first def file for a galaxy
def initial_def_file(Configuration, Fits_Files,Tirific_Template, Catalogue, cube_hdr,Initial_Parameters):



        #First we set some basic parameters that will hardly changed
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template,loops=7,outname = 'Cen_Conv',hdr=cube_hdr)
        # Then set the values for the various parameters
        set_model_parameters(Configuration, Tirific_Template,Initial_Parameters, hdr=cube_hdr)
        # Finally we set how these parameters are fitted.
        set_limit_modifier(Configuration,Initial_Parameters['INCL'][0] )
        set_fitting_parameters(Configuration, Tirific_Template,stage = 'initial',systemic = Initial_Parameters['VSYS'][0]/1000. )

        for key in Tirific_Template:
            print(f"{key} = {Tirific_Template[key]}")
        print("Not yet implemented")

'''






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


initial_def_file.__doc__ = '''

; NAME:
;        initial_def_file(Configuration, Fits_Files,Tirific_Template, Catalogue, cube_hdr,SBR_initial,pa,inclination):
;
; PURPOSE:
;       setup the first def file to be used in the fitting. As it is the first instance it has some special requirements.
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
;      split, strip, open
;
; EXAMPLE:
;
;
'''






def sofia(template,name):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")




sofia.__doc__ = '''

; NAME:
;       sofia
;
; PURPOSE:
;       write a sofia2 dictionary into file
;
; CATEGORY:
;       write
;
;
; INPUTS:
;       name = Name of the sofia file
;       template = a template dictionary that matches the sof read template
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
;      split, strip, open
;
; EXAMPLE:
;
;
'''
