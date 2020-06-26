#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used to Modify the Tirific_Template


class InitializeError(Exception):
    pass

from support_functions import set_rings,convertskyangle,sbr_limits,set_limits
import numpy as np

def set_fitting_parameters(Configuration, Tirific_Template, \
                           parameters_to_adjust  = ['NO_ADJUSTMENT'],
                           hdr = None,stage = 'initial',systemic = [100.,2], \
                           inclination = [60.,2.], pa = [90,1], \
                           rotation = [100.,5.],ra = [180,1e-4], dec= [0,1e-4]):
    fitting_settings = {}
    fitting_keys = ['VARY','VARINDX','MODERATE','DELEND','DELSTART','MINDELTA','PARMAX','PARMIN']
    if stage == 'initial':
        fitting_settings['SBR'] = set_sbr_fitting(Configuration, hdr = hdr,stage = stage, systemic = systemic[0])
        fitting_settings['INCL'] = set_incl_fitting(Configuration, hdr = hdr,stage = stage, inclination = inclination)
        fitting_settings['PA'] = set_pa_fitting(Configuration, hdr = hdr,stage = stage, pa = pa )
        fitting_settings['Z0'] = set_z0_fitting(Configuration, hdr = hdr,stage = stage, z0 = [0.2,2.] )
        fitting_settings['VROT'] = set_vrot_fitting(Configuration, hdr = hdr,stage = stage, rotation = rotation )
        fitting_settings['XPOS'],fitting_settings['YPOS'] = set_spatial_fitting(Configuration, hdr = hdr,stage = stage, ra = ra, dec = dec )
        fitting_settings['VSYS'] = set_vsys_fitting(Configuration, hdr = hdr,stage = stage, systemic = systemic )
        if inclination[0] < 30.:
            parameters_to_adjust = ['XPOS','YPOS','VSYS','PA','VROT','SBR','INCL']
        elif inclination[0] < 50.:
            parameters_to_adjust = ['XPOS','YPOS','VSYS','PA','VROT','INCL','SBR']
        elif inclination[0] > 75.:
            parameters_to_adjust = ['XPOS','YPOS','VSYS','PA','VROT','SBR','INCL','Z0']
        else:
            parameters_to_adjust = ['XPOS','YPOS','VSYS','PA','INCL','VROT','SBR']
    # Reset the fitting parameters
    for fit_key in fitting_keys:
        Tirific_Template[fit_key]= ''
    #write the new parameters
    for key in parameters_to_adjust:
        if key in fitting_settings:
            for fit_key in fitting_keys:
                if  fit_key in fitting_settings[key]:
                    if fit_key == 'VARY':
                        if len(Tirific_Template[fit_key]) == 0:
                            Tirific_Template[fit_key] = ', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])
                        else:
                            Tirific_Template[fit_key] = f"{Tirific_Template[fit_key]}, {', '.join([f'{x:<10s}' for x in fitting_settings[key][fit_key]])}"
                    else:
                        if fit_key == 'VARINDX':
                            format = '<10s'
                        else:
                            if key == 'SBR':
                                format = '<10.2e'
                            else:
                                format = '<10.2f'
                        Tirific_Template[fit_key] = f"{Tirific_Template[fit_key]} {' '.join([f'{x:{format}}' for x in fitting_settings[key][fit_key]])}"

set_fitting_parameters.__doc__ = '''

    ; NAME:
    ;      set_fitting_parameters(Configuration, Tirific_Template, \
                               parameters_to_adjust  = ['NO_ADJUSTMENT'],
                               hdr = None,stage = 'initial',systemic = [100.,2], \
                               inclination = [60.,2.], pa = [90,1], \
                               rotation = [100.,5.],ra = [180,1e-4], dec= [0,1e-4]):
    ;
    ; PURPOSE:
    ;      Set the parameters that control the fitting in the Tirific template
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


def set_incl_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', inclination = [60,5.]):
    NUR = Configuration['NO_RINGS']+2
    incl_limits = [set_limits(inclination[0]-inclination[1]-10.,5.,60.),
                    set_limits(inclination[0]+inclination[1]+10.,30.,90.)]
    incl_input= {}
    if stage == 'initial':
        incl_input['VARY'] =  np.array([f"INCL 0:{NUR} INCL_2 0:{NUR}"])
        incl_input['PARMAX'] = np.array([incl_limits[1]])
        incl_input['PARMIN'] = np.array([incl_limits[0]])
        incl_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        incl_input['DELSTART'] = np.array([inclination[1]]) # Starting step
        incl_input['DELEND'] = np.array([0.1*inclination[1]]) #Ending step
        incl_input['MINDELTA'] = np.array([0.1*inclination[1]]) #saturation criterum when /SIZE SIZE should be 10 troughout the code

    return incl_input
set_incl_fitting.__doc__ = '''

    ; NAME:
    ;      set_incl_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', incl_limits = [5.,90.]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the incl if required
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
            if key == 'SBR':
                format = '.2e'
            else:
                format = '.2f'
            if len(Model_Values[key]) == 2:
                Tirific_Template[key]= f"{Model_Values[key][0]:{format}}"
            else:
                Tirific_Template[key]= f"{' '.join([f'{x:{format}}' for x in Model_Values[key][:int(Configuration['NO_RINGS']+2)]])}"
            if key != 'RADI':
                key_write = f"{key}_2"
                if f"{key}_2" in Model_Values:
                    key = f"{key}_2"
                if len(Model_Values[key]) == 2:
                    Tirific_Template[key_write]= f"{Model_Values[key][0]:{format}}"
                else:
                    Tirific_Template[key_write]= f"{' '.join([f'{x:{format}}' for x in Model_Values[key][:int(Configuration['NO_RINGS']+2)]])}"


        else:
            if key == 'RADI':
                if hdr:
                    rad = set_rings(Configuration,hdr)
                    Tirific_Template['RADI']= f"{' '.join([f'{x:.2f}' for x in rad])}"
                    Tirific_Template['NUR']=str(len(rad))
                    Configuration['NO_RINGS'] = len(rad)-2
                else:
                    raise InitializeError('We cannot guess the radi without a header')
            elif key == 'Z0':
                if hdr:
                    if Model_Values['INCL'][0] > 80:
                        Tirific_Template['Z0'] = f"{np.max([convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True),hdr['BMAJ']/4.*3600.]):.3f}"

                    else:
                        Tirific_Template['Z0'] = f"{convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True):.3f}"

                else:
                    Tirific_Template['Z0'] = f"{convertskyangle(0.2,distance=Configuration['DISTANCE'],physical= True):.3f}"

                Tirific_Template['Z0_2'] = Tirific_Template['Z0']
            elif key == 'SDIS':
                Tirific_Template['SDIS'] = '8.'
                Tirific_Template['SDIS_2'] = '8.'
    #make sure that the first value in VROT = 0
    vrot = Tirific_Template['VROT'].split()
    if float(vrot[0]) != 0.:
        #These are +1 becuase +2-1
        Tirific_Template['VROT']=f" 0. {' '.join([f'{x}' for x in vrot[:int(Configuration['NO_RINGS']+1)]])}"
        Tirific_Template['VROT_2']=f" 0. {' '.join([f'{x}' for x in vrot[:int(Configuration['NO_RINGS']+1)]])}"

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

def set_overall_parameters(Configuration, Fits_Files,Tirific_Template, loops = 0, outname = 'random_fit',hdr= None, flux = None):

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
                Tirific_Template['BMAJ'] = f"{hdr['BMAJ']*3600:.2f}"
                Tirific_Template['BMIN'] = f"{hdr['BMIN']*3600:.2f}"
                Tirific_Template['RMS'] = f"{hdr['FATNOISE']:.2e}"
                try:
                    Tirific_Template['BPA'] = f"{hdr['BPA']:.2f}"
                except:
                    Tirific_Template['BPA'] = '0'
                if Configuration['HANNING']:
                    instrumental_vres = (hdr['CDELT3']/1000.*2)/(2.*np.sqrt(2.*np.log(2)))
                else:
                    instrumental_vres = (hdr['CDELT3']/1000.*1.2)/(2.*np.sqrt(2.*np.log(2)))
                Tirific_Template['CONDISP'] = f"{instrumental_vres:.2f}"
            if flux:
                Tirific_Template['CFLUX'] = f"{set_limits(flux/1.5e7,1e-6,1e-3):.2e}"
                Tirific_Template['CFLUX_2'] = f"{set_limits(flux/1.5e7,1e-6,1e-3):.2e}"

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


def set_pa_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', pa = [90,2.]):
    NUR = Configuration['NO_RINGS']+2
    pa_input = {}
    pa_limits = [0.,0.]
    if pa[0] < 180:
        pa_limits[1] = set_limits(pa[0]+pa[1]+10,0.,180.)
        pa_limits[0] = set_limits(pa[0]-pa[1]-10,0.,180.)
    else:
        pa_limits[1] = set_limits(pa[0]+pa[1]+10,180.,360.)
        pa_limits[0] = set_limits(pa[0]-pa[1]-10,180.,360.)
    if stage == 'initial':
        pa_input['VARY'] =  np.array([f"INCL 0:{NUR} pa_2 0:{NUR}"])
        pa_input['PARMAX'] = np.array([pa_limits[1]])
        pa_input['PARMIN'] = np.array([pa_limits[0]])
        pa_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        pa_input['DELSTART'] = np.array([0.5*pa[1]]) # Starting step
        pa_input['DELEND'] = np.array([0.05*pa[1]]) #Ending step
        pa_input['MINDELTA'] = np.array([0.1*pa[1]]) #saturation criterum when /SIZE SIZE should be 10 troughout the code

    return pa_input
set_pa_fitting.__doc__ = '''

    ; NAME:
    ;      set_pa_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', pa_limits = [5.,90.]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the incl if required
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

def set_sbr_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial'):
    if hdr:
        radii,sbr_ring_limits = sbr_limits(Configuration,hdr,systemic = systemic)
    else:
        sbr_ring_limits = np.zeros(Configuration['NO_RINGS']+2.)
        radii = np.zeros(Configuration['NO_RINGS']+2.)
    sbr_input = {}

    if Configuration['NO_RINGS']*Configuration['RING_SIZE'] < 7:
        sbr_input['VARY'] =  np.array([f"SBR {x} SBR_2 {x}" for x in range(len(radii),3,-1)])
        sbr_input['PARMAX'] = np.array([1 for x in range(len(radii),3,-1)])
        sbr_input['PARMIN'] = np.array([sbr_ring_limits[x]/2. if x <= (3./4.)*len(radii) else 0 for x in range(len(radii),3,-1)])
        sbr_input['MODERATE'] = np.array([5 for x in range(len(radii),3,-1)]) #How many steps from del start to del end
        sbr_input['DELSTART'] = np.array([7.5e-5 for x in range(len(radii),3,-1)]) # Starting step
        sbr_input['DELEND'] = np.array([2e-6 for x in range(len(radii),3,-1)]) #Ending step
        sbr_input['MINDELTA'] = np.array([2e-6 for x in range(len(radii),3,-1)]) #saturation criterum when /SIZE SIZE should be 10 troughout the code
    else:
        sbr_input['VARY'] =  np.array([[f"SBR {x}",f"SBR_2 {x}"] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2)
        sbr_input['PARMAX'] = np.array([[1,1] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2)
        sbr_input['PARMIN'] = np.array([[sbr_ring_limits[x]/2.,sbr_ring_limits[x]/2.] if x <= (3./4.)*len(radii) else [0.,0.] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2)
        sbr_input['MODERATE'] = np.array([[5,5] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2) #How many steps from del start to del end
        sbr_input['DELSTART'] = np.array([[7.5e-5,7.5e-5] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2) # Starting step
        sbr_input['DELEND'] = np.array([[2e-6,2e-6] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2)
        sbr_input['MINDELTA'] = np.array([[2e-6,2e-6] for x in range(len(radii),3,-1)]).reshape((len(radii)-3)*2)

    sbr_input['VARY'] = np.concatenate((sbr_input['VARY'],['SBR 1 2 3 SBR_2 1 2 3']),axis=0)
    sbr_input['PARMAX'] = np.concatenate((sbr_input['PARMAX'],[1]))
    sbr_input['PARMIN'] = np.concatenate((sbr_input['PARMIN'],[0]))
    sbr_input['MODERATE'] = np.concatenate((sbr_input['MODERATE'],[5]))
    sbr_input['DELSTART'] = np.concatenate((sbr_input['DELSTART'],[1e-5]))
    sbr_input['DELEND'] = np.concatenate((sbr_input['DELEND'],[1e-6]))
    sbr_input['MINDELTA'] = np.concatenate((sbr_input['MINDELTA'],[1e-6]))
    return sbr_input
set_sbr_fitting.__doc__ = '''

    ; NAME:
    ;      set_sbr_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial'):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the sbr if required
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
def set_spatial_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', ra = [180,1e-4], dec= [0,1e-4]):
    NUR = Configuration['NO_RINGS']+2
    xpos_input = {}
    ypos_input = {}
    if hdr:
        rabox_limits = np.sort([hdr['CRVAL1']-(hdr['CRPIX1']*hdr['CDELT1']),hdr['CRVAL1']+((hdr['NAXIS1']-hdr['CRPIX1'])*hdr['CDELT1']) ])
        decbox_limits = np.sort([hdr['CRVAL2']-(hdr['CRPIX2']*hdr['CDELT2']),hdr['CRVAL2']+((hdr['NAXIS2']-hdr['CRPIX2'])*hdr['CDELT2']) ])
    else:
        rabox_limits = [0,360]
        decbox_limits = [-90,90]
    xpos_limits = [0.,0.]
    ypos_limits = [0.,0.]
    if Configuration['NO_RINGS']*Configuration['RING_SIZE'] < 2:
        xpos_limits[1] = set_limits(ra[0]+3.*ra[1],*rabox_limits)
        xpos_limits[0] = set_limits(ra[0]-3.*ra[1],*rabox_limits)
        ypos_limits[1] = set_limits(dec[0]+3.*dec[1],*decbox_limits)
        ypos_limits[0] = set_limits(dec[0]-3.*dec[1],*decbox_limits)
    else:
        xpos_limits[1] = set_limits(ra[0]+np.max([3.*ra[1],3.*float(hdr['BMAJ'])]),*rabox_limits)
        xpos_limits[0] = set_limits(ra[0]-np.max([3.*ra[1],3.*float(hdr['BMAJ'])]),*rabox_limits)
        ypos_limits[1] = set_limits(dec[0]+np.max([3.*dec[1],3.*float(hdr['BMAJ'])]),*decbox_limits)
        ypos_limits[0] = set_limits(dec[0]-np.max([3.*dec[1],3.*float(hdr['BMAJ'])]),*decbox_limits)

    if stage == 'initial':
        xpos_input['VARY'] =  np.array([f"XPOS 0:{NUR} XPOS_2 0:{NUR}"])
        xpos_input['PARMAX'] = np.array([xpos_limits[1]])
        xpos_input['PARMIN'] = np.array([xpos_limits[0]])
        xpos_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        xpos_input['DELSTART'] = np.array([hdr['CDELT1']]) # Starting step
        xpos_input['DELEND'] = np.array([0.05*float(hdr['CDELT1'])]) #Ending step
        xpos_input['MINDELTA'] = np.array([0.1*float(hdr['CDELT1'])]) #saturation criterum when /SIZE SIZE should be 10 troughout the code
        ypos_input['VARY'] =  np.array([f"YPOS 0:{NUR} YPOS_2 0:{NUR}"])
        ypos_input['PARMAX'] = np.array([ypos_limits[1]])
        ypos_input['PARMIN'] = np.array([ypos_limits[0]])
        ypos_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        ypos_input['DELSTART'] = np.array([hdr['CDELT2']]) # Starting step
        ypos_input['DELEND'] = np.array([0.05*float(hdr['CDELT2'])]) #Ending step
        ypos_input['MINDELTA'] = np.array([0.1*float(hdr['CDELT2'])])

    return xpos_input, ypos_input
set_spatial_fitting.__doc__ = '''

    ; NAME:
    ;      set_spatial_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', ra = [180,1e-4.], dec= [0,1e-4]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the xpos and ypos if required
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
def set_vrot_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', rotation = [100,5.]):
    NUR = Configuration['NO_RINGS']+2
    vrot_input = {}
    vrot_limits = [set_limits(rotation[0]-rotation[1]-10,hdr['CDELT3']/1000.,360.), \
                   set_limits(rotation[0]+rotation[1]+10,80.,600.)]
    if stage == 'initial':
        inner_slope = int(round(set_limits(NUR*(4.-Configuration['LIMIT_MODIFIER'][0])/4.,round(NUR/2.),NUR-2)))
        vrot_input['VARY'] =  np.array([f"!VROT {NUR}:2 VROT_2 {NUR}:2"])

        vrot_input['PARMAX'] = np.array([vrot_limits[1]])
        vrot_input['PARMIN'] = np.array([vrot_limits[0]])
        vrot_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        vrot_input['DELSTART'] = np.array([hdr['CDELT3']/1000.*Configuration['LIMIT_MODIFIER'][0]]) # Starting step
        #These were lower in the original fat
        vrot_input['DELEND'] = np.array([0.1*hdr['CDELT3']/1000.*Configuration['LIMIT_MODIFIER'][0]]) #Ending step
        vrot_input['MINDELTA'] = np.array([0.05*hdr['CDELT3']/1000.*Configuration['LIMIT_MODIFIER'][0]]) #saturation criterum when /SIZE SIZE should be 10 troughout the code
        #if there is not values in the center we connect the inner ring to the next ring
        forvarindex = ''
        if Configuration['EXCLUDE_CENTRAL']:
            forvarindex = forvarindex+'VROT 2 VROT_2 2'
        if Configuration['NO_RINGS'] > 3:
            inner_slope = int(round(set_limits(NUR*(4.-Configuration['LIMIT_MODIFIER'][0])/4.,round(NUR/2.),NUR-2)))
            forvarindex = forvarindex+f"VROT {NUR-1}:{inner_slope} VROT_2 {NUR-1}:{inner_slope}"
        vrot_input['VARINDX'] = np.array([forvarindex])

    return vrot_input
set_vrot_fitting.__doc__ = '''

    ; NAME:
    ;      set_vrot_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', rotation = [100,5.]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the incl if required
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

def set_vsys_fitting(Configuration,hdr = None,systemic = [100.,2], stage = 'initial'):
    NUR = Configuration['NO_RINGS']+2
    vsys_input = {}
    vsys_limits = [set_limits(systemic[0]-3.*systemic[1],hdr['CDELT3']/1000.,200000.), \
                      set_limits(systemic[0]+3.*systemic[1],hdr['CDELT3']/1000.,200000.)]
    if stage == 'initial':
        vsys_input['VARY'] =  np.array([f"VSYS 0:{NUR} VSYS_2 0:{NUR}"])
        vsys_input['PARMAX'] = np.array([vsys_limits[1]])
        vsys_input['PARMIN'] = np.array([vsys_limits[0]])
        vsys_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        vsys_input['DELSTART'] = np.array([hdr['CDELT3']/1000.]) # Starting step
        vsys_input['DELEND'] = np.array([0.1*hdr['CDELT3']/1000.]) #Ending step
        vsys_input['MINDELTA'] = np.array([0.05*hdr['CDELT3']/1000.]) #saturation criterum when /SIZE SIZE should be 10 troughout the code

    return vsys_input


set_vsys_fitting.__doc__ = '''

    ; NAME:
    ;      set_vsys_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', z0 = [0.2,2]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the Z0 if required
    ;      Z0 input is in kpc
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



def set_z0_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', z0 = [0.2,2.]):
    NUR = Configuration['NO_RINGS']+2
    z0_input = {}
    z0_limits = [set_limits(z0[0]-z0[1]-10.,0.05,0.5),
                    set_limits(np.sum(z0),0.5,5)]
    z0_limits = convertskyangle(z0,Configuration['DISTANCE'], physical = True)
    #if our limit is smaller than 1.5* beam make it 1.5*beam # This originally only applied at high inclination.
    if hdr:
        if z0_limits[1] < 1.5*float(hdr['BMAJ']):
            z0_limits[1] = 1.5*float(hdr['BMAJ'])
    if stage == 'initial':
        z0_input['VARY'] =  np.array([f"Z0 0:{NUR} Z0_2 0:{NUR}"])
        z0_input['PARMAX'] = np.array([z0_limits[1]])
        z0_input['PARMIN'] = np.array([z0_limits[0]])
        z0_input['MODERATE'] = np.array([5]) #How many steps from del start to del end
        z0_input['DELSTART'] = np.array([convertskyangle(0.1,Configuration['DISTANCE'], physical = True)]) # Starting step
        z0_input['DELEND'] = np.array([convertskyangle(0.05,Configuration['DISTANCE'], physical = True)]) #Ending step
        z0_input['MINDELTA'] = np.array([convertskyangle(0.05,Configuration['DISTANCE'], physical = True)]) #saturation criterum when /SIZE SIZE should be 10 troughout the code

    return z0_input
set_z0_fitting.__doc__ = '''

    ; NAME:
    ;      set_z0_fitting(Configuration,hdr = None,systemic = 100., stage = 'initial', z0 = [0.2,2]):
    ;
    ; PURPOSE:
    ;      the fitting parameter for the Z0 if required
    ;      Z0 input is in kpc
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
