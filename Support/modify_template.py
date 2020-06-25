#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to Modify the Tirific_Template


class InitializeError(Exception):
    pass

from support_functions import set_rings,convertskyangle
import numpy as np




#Function
def set_model_parameters(Configuration, Tirific_Template,Model_Values, hdr = None):
    parameters_to_set = ['RADI','VROT','Z0','SBR','INCL','PA','XPOS','YPOS','VSYS','SDIS']
    if 'VSYS' in Model_Values:
        vsys =Model_Values['VSYS'][0]/1000.
    else:
        vsys=100.
    for key in parameters_to_set:
        if key in Model_Values:
            # if 2 long we have a value and error
            if len(Model_Values[key]) == 2:
                Tirific_Template[key]= f"{Model_Values[key][0]}"
            else:
                Tirific_Template[key]= f"{' '.join([str(x) for x in Model_Values[key][:int(Tirific_Template['NUR'])]])}"
            if key != 'RADI':
                key_write = f"{key}_2"
                if f"{key}_2" in Model_Values:
                    key = f"{key}_2"
                if len(Model_Values[key]) == 2:
                    Tirific_Template[key_write]= f"{Model_Values[key][0]}"
                else:
                    Tirific_Template[key_write]= f"{' '.join([str(x) for x in Model_Values[key][:int(Tirific_Template['NUR'])]])}"
        else:
            if key == 'RADI':
                if hdr:
                    rad = set_rings(Configuration,hdr)
                    Tirific_Template['RADI']= f"{' '.join([str(x) for x in rad])}"
                    Tirific_Template['NUR']=str(len(rad))
                    Configuration['NO_RINGS'] = len(rad)-2
                else:
                    raise InitializeError('We cannot guess the radi without a header')
            elif key == 'Z0':
                if hdr:
                    if Model_Values['INCL'][0] > 80:
                        Tirific_Template['Z0'] = f"{np.max([convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True),hdr['BMAJ']/4.*3600.])}"

                    else:
                        Tirific_Template['Z0'] = f"{convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True)}"

                else:
                    Tirific_Template['Z0'] = f"{convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True)}"

                Tirific_Template['Z0_2'] = Tirific_Template['Z0']
            elif key == 'SDIS':
                Tirific_Template['SDIS'] = '8.'
                Tirific_Template['SDIS_2'] = '8.'

set_model_parameters.__doc__ = '''

    ; NAME:
    ;      set_model_parameters(Configuration, Tirific_Template,Model_Values, hdr = None):
    ;
    ; PURPOSE:
    ;      Set the model values parameters in the tirific file that are singular values that apply to all of the fitting such as names and other issues
    ;
    ; CATEGORY:
    ;       modify_template
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

def set_fitting_parameters(Configuration, Tirific_Template,stage = 'initial',systemic = 100.):

    sbr_input = set_sbr_fitting(Configuration, Tirific_Template,stage = stage, systemic = systemic)






    print("Not yet implemented")

def set_overall_parameters(Configuration, Fits_Files,Tirific_Template, loops = 0, outname = 'random_fit',hdr= None):

            if Configuration['OPTIMIZED']:
                Tirific_Template['INSET'] = f"{Fits_Files['OPTIMIZED_CUBE']}"
            else:
                Tirific_Template['INSET'] = f"{Fits_Files['FITTING_CUBE']}"

            if Configuration['NO_RINGS'] < 3:
                Tirific_Template['INIMODE'] = '1'
            elif Configuration['NO_RINGS'] < 13:
                Tirific_Template['INIMODE'] = '2'
            else:
                Tirific_Template['INIMODE'] = '3'

            Tirific_Template['NUR'] = f"{Configuration['NO_RINGS']+2}"

            #this could be fancier
            if Configuration['NO_RINGS'] < 3:
                Tirific_Template['NCORES'] = '2'
            elif Configuration['NO_RINGS'] < 6:
                Tirific_Template['NCORES'] = '3'
            elif Configuration['NO_RINGS'] < 12:
                Tirific_Template['NCORES'] = '4'
            else:
                Tirific_Template['NCORES'] = '6'

            Tirific_Template['LOOPS'] = f"{int(loops)}"
            Tirific_Template['DISTANCE'] = f"{Configuration['DISTANCE']}"
            out_keys = ['LOGNAME','OUTSET', 'GR_DEVICE','TIRDEF']
            out_extensions = ['log','fits', 'ps/vcps','def']
            for i,key in enumerate(out_keys):
                Tirific_Template[key] = f"{outname}.{out_extensions[i]}"
            #some things we only set if a header is provided
            if hdr:
                Tirific_Template['BMAJ'] = f"{hdr['BMAJ']*3600}"
                Tirific_Template['BMIN'] = f"{hdr['BMIN']*3600}"
                Tirific_Template['RMS'] = f"{hdr['FATNOISE']}"
                try:
                    Tirific_Template['BPA'] = f"{hdr['BPA']}"
                except:
                    Tirific_Template['BPA'] = '0'
                if Configuration['HANNING']:
                    instrumental_vres = (hdr['CDELT3']/1000.*2)/(2.*np.sqrt(2.*np.log(2)))
                else:
                    instrumental_vres = (hdr['CDELT3']/1000.*1.2)/(2.*np.sqrt(2.*np.log(2)))
                Tirific_Template['CONDISP'] = f"{instrumental_vres}"

set_overall_parameters.__doc__ = '''

    ; NAME:
    ;      set_overall_parameters(Configuration, Fits_Files,Tirific_Template, loops = 0, outname = 'random_fit')
    ;
    ; PURPOSE:
    ;      Set the parameters in the tirific file that are singular values that apply to all of the fitting such as names and other issues
    ;
    ; CATEGORY:
    ;       modify_template
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
def set_sbr_fitting(Configuration,Tirific_Template,hdr,systemic = 100., stage = 'initial'):
    radii,sbr_ring_limits = sbr_limits(Configuration,hdr,systemic = systemic)
    sbr_input = {}



    sbr_input['PARMAX'] = [1 for x in range(3,len(radii))]
    sbr_input['PARMIN'] = [sbr_ring_limits[x]/2. for x in range(3,len(radii))]
    sbr_input['MODERATE'] = [5 for x in range(3,len(radii))] #How many steps from del start to del end
    sbr_input['DELSTART'] = [7.5e-5 for x in range(3,len(radii))] # Starting step
    sbr_input['DELEND'] = [2e-6 for x in range(3,len(radii))]
    sbr_input['MINDELTA'] = [2e-6 for x in range(3,len(radii))]
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
