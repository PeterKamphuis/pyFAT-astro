# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to write text files to Disk

from pyFAT_astro.Support.log_functions import print_log,extract_date
from pyFAT_astro.Support.modify_template import set_model_parameters, set_overall_parameters,\
                                set_fitting_parameters,get_warp_slope, update_disk_angles
from make_moments.functions import extract_pv                                
from pyFAT_astro.Support.read_functions import load_basicinfo

import pyFAT_astro.Support.support_functions as sf

import copy
import numpy as np
import os
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Ellipse
    import matplotlib.axes as maxes
    import matplotlib.font_manager as mpl_fm
from astropy.io import fits
from astropy.wcs import WCS
# create or append to the basic ifo file
def basicinfo(Configuration,initialize = False,stage='TiRiFiC' ,
              RA= None, DEC= None,VSYS = None, PA= None,Inclination = None, \
              Max_Vrot = None,Tot_Flux = None, V_mask = None, \
              Distance = float('NaN') , DHI = float('NaN'), template = None):
    if RA is None:
        RA=[float('NaN'),float('NaN')]
    if DEC is None:
        DEC=[float('NaN'),float('NaN')]
    if VSYS is None:
        VSYS =[float('NaN'),float('NaN')]
    if PA is None:
        PA=[float('NaN'),float('NaN')]
    if Inclination is None:
        Inclination = [float('NaN'),float('NaN')]
    if Max_Vrot is None:
        Max_Vrot = [float('NaN'),float('NaN')]
    if Tot_Flux is None:
        Tot_Flux = [float('NaN'),float('NaN')]
    if V_mask is None:
        if np.sum(Configuration['VROT_CURRENT_BOUNDARY'])== 0.:
            V_mask = [float('NaN'),float('NaN')]
        else:
            V_mask = [2.*np.max(Configuration['VROT_CURRENT_BOUNDARY']),\
                        2* Configuration['CHANNEL_WIDTH']]
    if template is None:
        template = {}

    if initialize:
        with open(f"{Configuration['FITTING_DIR']}{Configuration['BASE_NAME']}-Basic_Info.txt",'w') as file:
            file.write(f'''#This file contains the basic parameters of the Galaxy
#The total flux is determined by the source finder in the initial estimates and the total emission in the masked cube of the final model in the fits
#D_HI is determined as the diameter where the major axis with given PA cuts the 10^20 column of the moment0 map
# {'RA':>25s} {'DEC':>25s} {'VSYS':>20s} {'PA':>20s} {'Inclination':>20s} {'Max VRot':>20s} {'V_mask':>20s} {'Tot FLux':>20s} {'D_HI':>20s} {'Distance':>20s} {'HI Mass':>20s} {'D_HI_kpc':>20s}
# {'hh:mm:ss':>25s} {'dd:mm:ss':>25s} {'km/s':>20s} {'deg':>20s} {'deg':>20s} {'km/s':>20s} {'km/s':>20s} {'Jy':>20s} {'arcsec':>20s} {'Mpc':>20s} {'M_solar':>20s} {'kpc':>20s}
# The initial input
''')
    else:
        Vars_to_Set =  ['XPOS','YPOS','VSYS','VROT','INCL','PA','SDIS','SBR','SBR_2','Z0']
        FAT_Model = sf.load_tirific(Configuration,template,\
            Variables= Vars_to_Set,array=True )
        RA=[FAT_Model[Vars_to_Set.index('XPOS'),0],abs(Configuration['BEAM'][0]/(3600.*2.))]
        DEC=[FAT_Model[Vars_to_Set.index('YPOS'),0],abs(Configuration['BEAM'][0]/(3600.*2.))]
        VSYS =np.array([FAT_Model[Vars_to_Set.index('VSYS'),0],Configuration['CHANNEL_WIDTH']],dtype=float)
        PA=[FAT_Model[Vars_to_Set.index('PA'),0], 3.]
        Inclination = [FAT_Model[Vars_to_Set.index('INCL'),0], 3.]
        Max_Vrot = [np.max(FAT_Model[Vars_to_Set.index('VROT'),:]),np.max(FAT_Model[Vars_to_Set.index('VROT'),:])-np.min(FAT_Model[Vars_to_Set.index('VROT'),1:]) ]
        Distance = Configuration['DISTANCE']

    with open(f"{Configuration['FITTING_DIR']}{Configuration['BASE_NAME']}-Basic_Info.txt",'a') as file:
        if not initialize:
            file.write(f'''#These are the values from {stage}. \n''')
        RAhr,DEChr = sf.convertRADEC(Configuration,RA[0],DEC[0] )
        RA_c = f'{RAhr}+/-{RA[1]*3600.:0.2f}'
        DEC_c = f'{DEChr}+/-{DEC[1]*3600.:0.2f}'
        VSYS_c = f'{VSYS[0]:.2f}+/-{VSYS[1]:.2f}'
        PA_c = f'{PA[0]:.2f}+/-{PA[1]:.2f}'
        INCL_c = f'{Inclination[0]:.2f}+/-{Inclination[1]:.2f}'
        MVROT_c = f'{Max_Vrot[0]:.2f}+/-{Max_Vrot[1]:.2f}'
        Vmask_c = f'{float(V_mask[0])/1000.:.2f}+/-{float(V_mask[1])/1000.:.2f}'
        DHI_a = f'{DHI[0]:.2f}+/-{DHI[1]:.2f}'
        Dist = f'{Distance:.2f}'
        HIMass  = f'{Tot_Flux[0]*2.36E5*Distance**2:.2e}'
        DHI_k = f'{sf.convertskyangle(Configuration,DHI[0]):.2f}'
        Flux_c = f'{Tot_Flux[0]:.2f}+/-{Tot_Flux[1]:.2f}'
        file.write(f'''  {RA_c:>25s} {DEC_c:>25s} {VSYS_c:>20s} {PA_c:>20s} {INCL_c:>20s} {MVROT_c:>20s} {Vmask_c:>20s} {Flux_c:>20s} {DHI_a:>20s} {Dist:>20s} {HIMass:>20s} {DHI_k:>20s}
''')

basicinfo.__doc__ =f'''
 NAME:
    basicinfo

 PURPOSE:
    Create or update the file with the basic information of the fit

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


    initialize = False
    If true a new file with header but and first guess values will be created.
    These parameters need to be provided through the keywords below.

    stage = 'Tirific'

    RA=[float('NaN'),float('NaN')], DEC=[float('NaN'),float('NaN')]
    RA and DEC in degrees

    VSYS =[float('NaN'),float('NaN')]
    Systemic velocity in m/s

    PA=[float('NaN'),float('NaN')], Inclination = [float('NaN'),float('NaN')]
    Position angle and inclination in degrees

    Max_Vrot = [float('NaN'),float('NaN')]
    The maximum rotational velocity

    Tot_Flux = [float('NaN'),float('NaN')]
    The total flux in the source in Jy

    V_mask = [float('NaN'),float('NaN')],
    extend of the mask in the cube

    Distance = float('NaN')
    Distance to the source in Mpc

    DHI = float('NaN')
    Diameter of the 1 M_solar/pc^2 iso density surface ellipse. in arcsec

    Standard FAT dictionary with filenames

    template = ['EMPTY']
    Tirific template to obtain values from

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def create_tirific_run_cube(Configuration,Fits_Files):
    # if we have an optimized size cube this is only used in the run and we need t o change the ctyp 
    if Configuration['OPTIMIZED']:
        optimized_cube = f'{Configuration["FITTING_DIR"]}/{Fits_Files["OPTIMIZED_CUBE"]}'
        cube = fits.open(optimized_cube)
        cube[0].header['CTYPE3'] = 'VELO'
        fits.writeto(optimized_cube,cube[0].data,cube[0].header,overwrite=True)
    source = f'{Configuration["FITTING_DIR"]}/{Fits_Files["FITTING_CUBE"]}'
    stripped_file_name = os.path.splitext(Fits_Files["FITTING_CUBE"])[0]
    target = f'{Configuration["FITTING_DIR"]}/{stripped_file_name}_tirific.fits'
    os.system(f'''cp {source} {target}''')
    cube = fits.open(target)
    cube[0].header['CTYPE3'] = 'VELO'
    fits.writeto(target,cube[0].data,cube[0].header,overwrite=True)
    Fits_Files['TIR_RUN_CUBE'] = target    
    return Fits_Files

# Function to write the first def file for a galaxy
def initialize_def_file(Configuration, Fits_Files,Tirific_Template,\
        Initial_Parameters = None, fit_type = 'Undefined' ):

    if Initial_Parameters is None:
        Initial_Parameters = {}


    #First we set some basic parameters that will hardly change
    if fit_type == 'Centre_Convergence':

        #if 'VSYS' in Initial_Parameters:
        #    Initial_Parameters['VSYS'] = [x/1000. for x in Initial_Parameters['VSYS']]
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template,\
            fit_type=fit_type, flux = Initial_Parameters['FLUX'][0] )
        # Then set the values for the various parameters of the model

        set_model_parameters(Configuration, Tirific_Template,Initial_Parameters )

        sf.set_limit_modifier(Configuration,Tirific_Template  )

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'initial',
                                initial_estimates = Initial_Parameters )
        #if 'VSYS' in Initial_Parameters:
        #    Initial_Parameters['VSYS'] = [x*1000. for x in Initial_Parameters['VSYS']]
    elif fit_type in ['Extent_Convergence','Fit_Tirific_OSC']:
        #if 'VSYS' in Initial_Parameters and fit_type == 'Fit_Tirific_OSC':
        #    Initial_Parameters['VSYS'] = [x/1000. for x in Initial_Parameters['VSYS']]
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,\
            fit_type=fit_type)
        Vars_to_Set =  ['XPOS','YPOS','VSYS','VROT','INCL','PA','SDIS','SBR','SBR_2','Z0']
        if fit_type == 'Fit_Tirific_OSC':
            set_model_parameters(Configuration, Tirific_Template,Initial_Parameters,stage='initialize_def_file' )

        # Finally we set how these parameters are fitted.
        sf.set_limit_modifier(Configuration,Tirific_Template )
        #get_inner_fix(Configuration,Tirific_Template )
        get_warp_slope(Configuration,Tirific_Template )
        #if Configuration['OUTER_RINGS_DOUBLED']:
        Initial_Parameters['XPOS'][1]= sf.set_limits(Initial_Parameters['XPOS'][1],Configuration['BEAM'][0]/3600.,Configuration['BEAM'][0]/3600.*5)
        Initial_Parameters['YPOS'][1]= sf.set_limits(Initial_Parameters['YPOS'][1],Configuration['BEAM'][0]/3600.,Configuration['BEAM'][0]/3600.*5)
        Initial_Parameters['VSYS'][1]= sf.set_limits(Initial_Parameters['VSYS'][1],Configuration['CHANNEL_WIDTH']/2.,Configuration['CHANNEL_WIDTH']*5)
        if Initial_Parameters['INCL'][0] >35:
            Initial_Parameters['PA'][1]= sf.set_limits(Initial_Parameters['PA'][1],1.,15)
        Initial_Parameters['INCL'][1]= sf.set_limits(Initial_Parameters['INCL'][1],3.,15)


        set_fitting_parameters(Configuration, Tirific_Template,stage = 'initialize_os',
                               initial_estimates=Initial_Parameters )

    tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def' )

initialize_def_file.__doc__ =f'''
 NAME:
    initialize_def_file

 PURPOSE:
    setup the first def file to be used in the fitting. As it is the first instance it has some special requirements.

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Standard FAT dictionary with filenames
    Tirific_Template = Standard FAT Tirific Template


 OPTIONAL INPUTS:


    Initial_Parameters = ['EMPTY']
    The initial parameters to be used to set up the def file.

    fit_type = 'Undefined'
    type of fitting being done

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:"/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
    Unspecified

 NOTE:
'''

def beam_artist(ax,hdr,im_wcs):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ghxloc, ghyloc = im_wcs.wcs_pix2world(float(xmin+(xmax-xmin)/18.), float(ymin+(ymax-ymin)/18.), 1.)
    localoc = [float(ghxloc),float(ghyloc) ]
    widthb = hdr['BMIN']
    heightb = hdr['BMAJ']
    try:
        angleb  = hdr['BPA']
    except KeyError:
        angleb = 0.
    #either the location or the beam has to be transformed 
    beam = Ellipse(xy=localoc, width=widthb, height=heightb, angle=angleb, transform=ax.get_transform('world'),
           edgecolor='k', lw=1, facecolor='none', hatch='/////',zorder=15)
    return beam

beam_artist.__doc__ =f'''
 NAME:
    ist

 PURPOSE:
    create a beam patch
        
 CATEGORY:
    write_functions

 INPUTS:
    ax = is the axes object where the beam is intended to go
    hdr = the image header
    im_wcs = WCS frame of the image

 OPTIONAL INPUTS:


 OUTPUTS:
    IST
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:

 NOTE:
'''

def create_plot_stats(Configuration,loads,labels):
    combined_time =  np.sort(np.array(loads['Tirific']['Time']+loads['FAT']['Time'],dtype=float))

    combined_loads ={'Tirific':{'CPU':np.interp(combined_time,np.array(loads['Tirific']['Time'],dtype=float),np.array(loads['Tirific']['CPU'],dtype=float)),\
                                'MEM':np.interp(combined_time,np.array(loads['Tirific']['Time'],dtype=float),np.array(loads['Tirific']['MEM'],dtype=float))},\
                    'FAT':{'CPU':np.interp(combined_time,np.array(loads['FAT']['Time'],dtype=float),np.array(loads['FAT']['CPU'],dtype=float)),\
                                                'MEM':np.interp(combined_time,np.array(loads['FAT']['Time'],dtype=float),np.array(loads['FAT']['MEM'],dtype=float))}

    }
    comb_list= labels['Tirific']['Time']+labels['FAT']['Time']
    comb_label = labels['Tirific']['label']+labels['FAT']['label']

    labels_times=np.array([x for x, _ in sorted(zip(comb_list, comb_label))],dtype=float)
    labels_comb = np.array([x for _, x in sorted(zip(comb_list, comb_label))],dtype=str)
    if Configuration['MULTIPROCESSING']:
        non_split_ct = copy.deepcopy(combined_time)
        combined_time =[[],[]]
        non_split_cl = copy.deepcopy(combined_loads)
        combined_loads = [[],[]]
        non_split_lt = copy.deepcopy(labels_times)
        labels_times =[[],[]]
        non_split_lc = copy.deepcopy(labels_comb)
        labels_comb =[[],[]]
        split_time_lab = non_split_lt[np.where(non_split_lc == 'Pausing FAT')[0]] 
        split_time = non_split_lt[np.where(non_split_lc == 'Pausing FAT')[0]+1] 
     
        for i in  [0,1]:
            if i == 0:
                indxs = np.where(non_split_ct <= split_time)[0]
                indxslab =  np.where(non_split_lt <= split_time_lab)[0]
              
            else:
                indxs = np.where(non_split_ct > split_time)[0]  
                indxslab =  np.where(non_split_lt > split_time_lab)[0]
          
            combined_time[i] = non_split_ct[indxs]
            labels_comb[i] = non_split_lc[indxslab]
            labels_times[i] = non_split_lt[indxslab ]
            combined_loads[i] = {'Tirific':{'CPU': non_split_cl['Tirific']['CPU'][indxs],\
                                            'MEM': non_split_cl['Tirific']['MEM'][indxs] },\
                                'FAT': {'CPU': non_split_cl['FAT']['CPU'][indxs],\
                                        'MEM': non_split_cl['FAT']['MEM'][indxs] }}
           
    else:
        combined_time = [combined_time ]
        combined_loads = [combined_loads]
        labels_times = [labels_times]
        labels_comb = [labels_comb ]

    return combined_time, combined_loads, labels_times, labels_comb

create_plot_stats.__doc__ =f'''
 NAME:
    create_plot_stats(

 PURPOSE:
    create the lists and dictionaries that are required to plot a single frame in the usage plot
    from the over loads and times dictionaries
        
 CATEGORY:
    write_functions

 INPUTS:
    Configuration = is the standard FAT configuration
    loads = is the overall loads dictionary
    labels = the overall labels


 OPTIONAL INPUTS:


 OUTPUTS:
    combined_time = list with all times for individual axis object
    combined_loads = dictionary with all loads  for individual axis object
    labels_times = list with all label times  for individual axis object
    labels_comb = liast with all labels  for individual axis object
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:

 NOTE:
'''

def make_overview_plot(Configuration,Fits_Files ):
    fit_type = Configuration['USED_FITTING']
    print_log(f'''MAKE_OVERVIEW_PLOT: We are starting the overview plot.
''',Configuration,case=['debug_start'])

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cube_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
        moment0_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom0.fits")
        moment1_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom1.fits")
        moment2_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom2.fits")
        cube = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}")
        moment0 = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT0']}")
        moment1 = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT1']}")
        moment2 = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['MOMENT2']}")
        channels_map = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['CHANNEL_MAP']}")
        im_wcs = WCS(moment0[0].header).celestial
   

    # Open the model info
    print_log(f'''MAKE_OVERVIEW_PLOT: Reading the variables from the final model
''',Configuration,case=['debug_add'])
    Vars_to_plot_short= ['RADI','XPOS','YPOS','VSYS','VROT','INCL','PA','SDIS',\
                    'SBR','Z0']
    Vars_to_plot=copy.deepcopy(Vars_to_plot_short)
    for x in Vars_to_plot_short:
        if x != 'RADI':
            Vars_to_plot.append(f'# {x}_ERR')
        if x not in ['XPOS','YPOS','VSYS']:
            Vars_to_plot.append(f'{x}_2')
            Vars_to_plot.append(f'# {x}_2_ERR')
    FAT_Model = sf.load_tirific(Configuration,\
        f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def",\
        Variables= Vars_to_plot,array=True ,brightness_check=True)
    Extra_Model_File = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Iteration_{Configuration['ITERATIONS']}.def"

    if os.path.exists(Extra_Model_File):
        Extra_Model = sf.load_tirific(Configuration,Extra_Model_File,\
            Variables= Vars_to_plot,array=True ,brightness_check=True)
        
    else:
        Extra_Model = []
    print_log(f'''MAKE_OVERVIEW_PLOT: We find the following model values.
{'':8s}{[f"{x} = {FAT_Model[i,:]}" for i,x in enumerate(Vars_to_plot)]}
''',Configuration,case=['debug_add'])

    if os.path.exists(f"{Configuration['FITTING_DIR']}ModelInput.def"):
        Input_Model = sf.load_tirific(Configuration,\
            f"{Configuration['FITTING_DIR']}ModelInput.def",\
            Variables= Vars_to_plot,array=True,brightness_check=False )
    else:
        Input_Model = []
    sof_basic_ra,sof_basic_dec, sof_basic_vsys,sof_basic_maxrot,sof_basic_pa,sof_basic_inclination,sof_basic_extent = load_basicinfo(Configuration,
        f"{Configuration['FITTING_DIR']}{Configuration['BASE_NAME']}-Basic_Info.txt",Variables=['RA','DEC','VSYS','Max VRot','PA','Inclination','D_HI'])


    #Let's start plotting
    ysize = 23.2/2.
    xsize = 0.7*ysize
    Overview = plt.figure(2, figsize=(xsize, ysize), dpi=300, facecolor='w', edgecolor='k')

    size_factor= ysize/11.6
    size_ratio = ysize/xsize

    #stupid pythonic layout for grid spec, which means it is yx instead of xy like for normal human beings
    gs = Overview.add_gridspec(int(20*size_ratio),20)
    try:
        mpl_fm.fontManager.addfont(Configuration['FONT_FILE'])
        font_name = mpl_fm.FontProperties(fname=Configuration['FONT_FILE']).get_name()
    except FileNotFoundError:
        font_name = 'DejaVu Sans'
    labelfont = {'family': font_name,
             'weight': 'normal',
             'size': 8*size_factor}
    plt.rc('font', **labelfont)
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
# First some textual information

    if Configuration['SUB_DIR'] == './':
        name = Configuration['BASE_NAME']
    else:
        name = Configuration['SUB_DIR']
    #plt.title(f"Overview for {name}")
    ax_text = Overview.add_subplot(gs[0:1,0:21],frameon = False)
    plt.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False,
        right = False,
        left= False,
        labelleft = False)
    ax_text.text(0.5,1.0,f'''Overview for {name}''',rotation=0, va='top',ha='center', color='black',transform = ax_text.transAxes,
      bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=14*size_factor)
    ax_text.text(0.5,0.25,f'''The ring size used in the model is {Configuration['RING_SIZE']:.2f} x BMAJ, with BMAJ = {Configuration['BEAM'][0]:.1f} arcsec. We assumed a distance  of {Configuration['DISTANCE']:.1f} Mpc.'''
      ,rotation=0, va='top',ha='center', color='black',transform = ax_text.transAxes,
      bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)


#-----------------------------------------------------------------Moment 0 ------------------------------------------------------
    ax_moment0 = Overview.add_subplot(gs[2:8,0:6], projection=im_wcs)
    # we need contour levels and
    min_color = 0.
    max_color = np.nanmax(moment0[0].data)*0.8
    moment0_plot = ax_moment0.imshow(moment0[0].data, origin='lower', alpha=1,\
        vmin = min_color, vmax = max_color,cmap='hot_r',transform=ax_moment0.get_transform(im_wcs) )
    moment0_plot.set_label('Intensity Map')
    plt.ylabel('DEC')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA')

    median_noise_in_map = np.sqrt(np.nanmedian(channels_map[0].data[channels_map[0].data > 0.]))*Configuration['NOISE']*Configuration['CHANNEL_WIDTH']
    mindism0 = median_noise_in_map
    mindism0 = median_noise_in_map
    #"We find this {} as the minimum of the moment0 map".format(mindism0))
    if mindism0 < 0.:
        mindism0  =abs(mindism0)*2.
    if mindism0 == 0:
        mindism0 = np.max(moment0[0].data)/64.
    #print("We find this {} as the minimum of the moment0 map".format(mindism0))
    maxdism0 = np.max(moment0[0].data) * 0.8
    if mindism0 > maxdism0:
        mindism0 = 0.1*maxdism0
    if maxdism0 < 16*mindism0:
        momlevel = np.array([1,4,8,12],dtype=float)* mindism0
    elif maxdism0 < 32*mindism0:
        momlevel = np.array([1,4,8,16,24,32],dtype=float)* mindism0
    elif maxdism0 < 256*mindism0:
        momlevel = np.array([1,4,8,32,64,128],dtype=float) * mindism0
    else:
        momlevel = np.array([0,1,2,3,4,5,6,7],dtype=float) * (maxdism0-mindism0)/7+ mindism0

    #print("We find this {} as the minimum of the moment0 map".format(mindism0))
    momlevel = np.array([x for x in momlevel if x < np.max(moment0[0].data)*0.95],dtype=float)
    if momlevel.size == 0:
        momlevel=0.5*mindism0
    ax_moment0.contour(moment0[0].data, transform=ax_moment0.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.5*size_factor , zorder =4)
    ax_moment0.contour(moment0[0].data, transform=ax_moment0.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=1.2*size_factor)
    ax_moment0.contour(moment0_mod[0].data, transform=ax_moment0.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2*size_factor , zorder =7)
    ax_moment0.contour(moment0_mod[0].data, transform=ax_moment0.get_transform(im_wcs),
              levels=momlevel, colors='r',zorder=8, linewidths=0.9*size_factor)
  
    square_plot(ax_moment0)
    ax_moment0.grid()    
    beam = beam_artist(ax_moment0,moment0[0].header,im_wcs)
    ax_moment0.add_patch(beam)
    #center_x,center_y = im_wcs.wcs_world2pix(FAT_Model[Vars_to_plot.index('XPOS'),0],\
    #                    FAT_Model[Vars_to_plot.index('YPOS'),0], 1.)

    center_x = float(FAT_Model[Vars_to_plot.index('XPOS'),0])
    center_y = float(FAT_Model[Vars_to_plot.index('YPOS'),0])
   
    ax_moment0.text(center_x,center_y,'X', size= 5.,va='center',ha='center', \
                    color='white',zorder=17, transform=ax_moment0.get_transform('world'))
    #ax_moment0.text(center_x,center_y,'X', size= 7.,va='center',ha='center', \
     #               color='black',zorder=16, transform=ax_moment0.get_transform('world'))


    #                ,transform = ax_moment0.transAxes,zorder=7)
    #
    # colorbar
    divider = make_axes_locatable(ax_moment0)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment0_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    cbar.ax.set_title(f"{moment0[0].header['BUNIT']}", y= 0.2*size_factor**2)

    '''the channel width should be 1 as the map is already in Jy/beam*km/s'''
    column_levels = sf.columndensity(Configuration,momlevel*1000.,\
        systemic = FAT_Model[Vars_to_plot.index('VSYS'),0], channel_width=1)

    if 1e21 < np.min(column_levels):
        fact= 1e21
        str_fact = r'$\times 10^{21} {\rm cm}^{-2}$'
    elif 1e20 < np.min(column_levels) < 1e21:
        fact= 1e20
        str_fact = r'$\times 10^{20} {\rm cm}^{-2}$'
    elif 1e19 < np.min(column_levels) < 1e20:
        fact= 1e19
        str_fact = r'$\times 10^{19} {\rm cm}^{-2}$'
    elif 1e18 < np.min(column_levels) < 1e19:
        fact= 1e18
        str_fact = r'$\times 10^{18} {\rm cm}^{-2}$'
    else:
        fact= 1e17
        str_fact = r'$\times 10^{17} {\rm cm}^{-2}$'
    if len(momlevel) < 4:
        info_string = f"The contours are at {', '.join(['{:.1f}'.format(x/fact) for x in column_levels])} {str_fact}."
    else:
        info_string = f"The contours are at {', '.join(['{:.1f}'.format(x/fact) for x in column_levels[0:4]])}"
        counter = 4
        while counter < len(column_levels):
            info_string = info_string+f"\n {', '.join(['{:.1f}'.format(x/fact) for x in column_levels[counter:counter+7]])}"
            counter += 7
        info_string = info_string+f" {str_fact}."
    #info_string = f"The contours are at {column_levels}."
    ax_moment0.text(-0.1,-0.2,info_string, va='top',ha='left', color='black',transform = ax_moment0.transAxes,
              bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7)
    # No further need for the moment maps
    moment0.close()
    moment0_mod.close()
#-----------------------------------------------------------------Velocity Field------------------------------------------------------
    ax_moment1 = Overview.add_subplot(gs[10:16, 0:6], projection=im_wcs)
    ax_moment1.set_label('Velocity Field')
    #Comp_ax1.set_facecolor('black')
    # we need contour levels and
    inclination_correction = sf.set_limits(FAT_Model[Vars_to_plot.index('INCL'),0]+12.5,20.,90.)
    velocity_width= sf.set_limits(1.25*np.nanmax(FAT_Model[Vars_to_plot.index('VROT'),:])*np.sin(np.radians(inclination_correction)),30.,700.)

    max_color= FAT_Model[Vars_to_plot.index('VSYS'),0]+velocity_width
    min_color= FAT_Model[Vars_to_plot.index('VSYS'),0]-velocity_width

    moment1_plot = ax_moment1.imshow(moment1[0].data, cmap='rainbow', origin='lower', alpha=1, vmin = min_color, vmax = max_color )
    plt.ylabel('DEC')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA')

    # contours
    velocity_step=sf.set_limits(int((int(max_color-min_color)*0.9)/20.),1.,30.)
    integer_array = np.linspace(0,20,21)-10
    momlevel = [FAT_Model[Vars_to_plot.index('VSYS'),0]+x*velocity_step for x in integer_array if min_color < FAT_Model[Vars_to_plot.index('VSYS'),0]+x*velocity_step < max_color]
    ax_moment1.contour(moment1[0].data, transform=ax_moment1.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.5 *size_factor, zorder =4)
    ax_moment1.contour(moment1[0].data, transform=ax_moment1.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=0.9*size_factor)
    ax_moment1.contour(moment1_mod[0].data, transform=ax_moment1.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2*size_factor , zorder =7)
    #ax_moment1.contour(moment1_mod[0].data, transform=ax_moment1.get_transform(im_wcs),
    #          levels=momlevel, colors='r',zorder=8, linewidths=0.9)
   
    square_plot(ax_moment1)
    beam = beam_artist(ax_moment1,moment1[0].header,im_wcs)
    ax_moment1.add_patch(beam)
    ax_moment1.grid()
    # colorbar
    divider = make_axes_locatable(ax_moment1)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment1_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    cbar.ax.set_title(f"{moment1[0].header['BUNIT']}", y= 0.2*size_factor**2)

    column_levels = ', '.join(["{:.1f}".format(x) for x in momlevel])
    info_string= f'''The contours start at {float(momlevel[0]):.1f} km/s'''
    if len(momlevel) > 1:
        info_string=f'''{info_string} \n and increase with {float(momlevel[1])-float(momlevel[0]):.1f} km/s.'''
    else:
        info_string=f'''{info_string}.'''


    ax_moment1.text(-0.1,-0.2,info_string, va='top',ha='left', color='black',transform = ax_moment1.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7)
    # No further need for the moment maps
    moment1.close()
    moment1_mod.close()
#-----------------------------------------------------------------Moment 2------------------------------------------------------
    ax_moment2 = Overview.add_subplot(gs[10:16:, 8:14], projection=im_wcs)
    ax_moment2.set_label('Moment2')
    #Comp_ax1.set_facecolor('black')
    # we need contour levels and
    max_color= sf.set_limits(np.nanmax(moment2[0].data),15,50)
    min_color= 0.

    moment2_plot = ax_moment2.imshow(moment2[0].data, cmap='rainbow' ,origin='lower', alpha=1, vmin = min_color, vmax = max_color )
    plt.ylabel('DEC')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA')

    # contours

    momlevel = np.linspace(min_color,max_color*0.8,5)

    ax_moment2.contour(moment2[0].data, transform=ax_moment2.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.5*size_factor , zorder =4)
    ax_moment2.contour(moment2[0].data, transform=ax_moment2.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=1.2*size_factor)
    ax_moment2.contour(moment2_mod[0].data, transform=ax_moment2.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2*size_factor , zorder =7)

   
    square_plot(ax_moment2)
    beam = beam_artist(ax_moment2,moment2[0].header,im_wcs)
    ax_moment2.add_patch(beam)
    ax_moment2.grid()

    # colorbar
    divider = make_axes_locatable(ax_moment2)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment2_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    #cbar.set_title(label=f"{moment2[0].header['BUNIT']}")
    cbar.ax.set_title(f"{moment2[0].header['BUNIT']}", y= 0.2*size_factor**2)

    column_levels = ', '.join(["{:.1f}".format(x) for x in momlevel])
    if len(momlevel) < 4:
        info_string = f"The contours are at {column_levels} km/s."
    else:
        info_string = f"The contours are at {', '.join(['{:.1f}'.format(x) for x in momlevel[0:4]])}"
        counter = 4
        while counter < len(momlevel):
            info_string = info_string+f"\n {', '.join(['{:.1f}'.format(x) for x in momlevel[counter:counter+7]])}"
            counter += 7
            info_string = info_string+" km/s."

    ax_moment2.text(-0.1,-0.2,info_string, va='top',ha='left', color='black',transform = ax_moment2.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7)
    # No further need for the moment maps
    moment2.close()
    moment2_mod.close()



#__________________------------------------------------------------------------PV Diagram

    extract_angle = np.mean(FAT_Model[Vars_to_plot.index('PA'),0:round(len(FAT_Model[Vars_to_plot.index('PA'),:])/2.)])

    messages = extract_pv(cube = cube,\
                overwrite = False,PA=extract_angle,\
                center=  [float(FAT_Model[Vars_to_plot.index('XPOS'),0]),\
                        float(FAT_Model[Vars_to_plot.index('YPOS'),0]),\
                        float(FAT_Model[Vars_to_plot.index('VSYS'),0]*1000.)],\
                map_velocity_unit= 'km/s',log = True,silent = True,\
                output_directory = f"{Configuration['FITTING_DIR']}Finalmodel",\
                output_name =f"{Configuration['BASE_NAME']}_final_xv.fits")
    print_log(messages,Configuration,case=["verbose"])
    PV = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/{Configuration['BASE_NAME']}_final_xv.fits")
   
    messages = extract_pv(cube = cube_mod,\
                overwrite = False,PA=extract_angle,\
                center=  [float(FAT_Model[Vars_to_plot.index('XPOS'),0]),\
                        float(FAT_Model[Vars_to_plot.index('YPOS'),0]),\
                        float(FAT_Model[Vars_to_plot.index('VSYS'),0]*1000.)],\
                cube_velocity_unit='m/s',\
                map_velocity_unit= 'km/s',log = True,silent = True,\
                output_directory = f"{Configuration['FITTING_DIR']}/Finalmodel/",\
                output_name =f"Finalmodel_xv.fits")
    print_log(messages,Configuration,case=["verbose"])
    PV_model = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_xv.fits")


    ax_PV = plot_PV(Configuration,image=PV, model = PV_model, figure = Overview, \
        location = gs[2:8,8:14],size_factor = size_factor,
        tirific_model = f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def")
    PV.close()
    PV_model.close()

# ADD a legend for the plots

    ax_legend = Overview.add_subplot(gs[18:27,0:2],frameon = False)
    plt.tick_params(
            axis='both',          # changes apply to the x-axis
            which='both',      # both major and minor ticks are affected
            direction = 'in',
            bottom=False,      # ticks along the bottom edge are off
            top=False,         # ticks along the top edge are off
            labelbottom=False,
            right = False,
            left= False,
            labelleft = False)
    legend = ['Appr.','Rec.','Correct Appr.','Correct Rec.','Unsm. Appr.','Unsm. Rec.']
    plt.plot([0,1],[0,1],'k',marker='o',label=legend[0])
    plt.plot([0,1],[0,1],'r',marker='o',label=legend[1])
    if len(Input_Model) > 0:
        plt.plot([0,1],[0,1],'b',marker='o',label=legend[2])
        plt.plot([0,1],[0,1],'yellow',marker='o',label=legend[3])
    if len(Extra_Model) > 0:
        plt.plot([0,1],[0,1],'k',marker='o',label=legend[4],alpha=0.2)
        plt.plot([0,1],[0,1],'r',marker='o',label=legend[5],alpha=0.2)
    plt.plot([0,1],[0,1],'white',marker='o',linewidth=5)
    plt.ylim(0.4,0.6)
    plt.xlim(0.4,0.6)

    chartBox = ax_legend.get_position()
    ax_legend.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
    ax_legend.legend(loc='upper left', bbox_to_anchor = (-1.0,1.03),shadow=True, ncol=1)


# ------------------------------Rotation curves------------------------------------

    labelfont= {'family':font_name,
            'weight':'normal',
            'size':10*size_factor}
    ax_RC = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[19:22,3:9],\
                            Overview,'VROT',Input_Model = Input_Model,initial_extent= sof_basic_extent[0], \
                            initial = sof_basic_maxrot[0],Extra_Model = Extra_Model)
    ymin =np.nanmin([FAT_Model[Vars_to_plot.index('VROT'),1:],FAT_Model[Vars_to_plot.index('VROT_2'),1:]])
    ymax =np.nanmax([FAT_Model[Vars_to_plot.index('VROT'),1:],FAT_Model[Vars_to_plot.index('VROT_2'),1:]])

    if len(Extra_Model) > 0:
        ymin2 =np.nanmin([Extra_Model[Vars_to_plot.index('VROT'),1:],Extra_Model[Vars_to_plot.index('VROT_2'),1:]])
        ymax2 =np.nanmax([Extra_Model[Vars_to_plot.index('VROT'),1:],Extra_Model[Vars_to_plot.index('VROT_2'),1:]])
    else:
        ymin2 = ymin
        ymax2 = ymax
  
    if len(Input_Model) > 0:
        ymin3 =np.nanmin([Input_Model[Vars_to_plot.index('VROT'),1:],Input_Model[Vars_to_plot.index('VROT_2'),1:]])
        ymax3 =np.nanmax([Input_Model[Vars_to_plot.index('VROT'),1:],Input_Model[Vars_to_plot.index('VROT_2'),1:]])
    else:
        ymin3=ymin
        ymax3= ymax
    ymin = np.nanmin([ymin,ymin2,ymin3])
    ymax = np.nanmax([ymax,ymax2,ymax3])

    buffer = np.nanmean([ymin,ymax])/20.
    ymin= ymin-buffer
    ymax = ymax+buffer
    if not np.isnan(ymin) and not np.isnan(ymax):
        ax_RC.set_ylim(ymin,ymax)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=False)
    plt.ylabel('RC (km s$^{-1}$)',**labelfont)
    arcmin,arcmax = ax_RC.get_xlim()
    sec_ax = ax_RC.twiny()
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.figure.canvas.draw()
    sec_ax.set_xlabel('Radius (kpc)',va='bottom',**labelfont)

# ------------------------------Inclination------------------------------------
    ax_INCL = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[22:25,3:9],\
                              Overview,'INCL',Input_Model = Input_Model,initial_extent= sof_basic_extent[0], \
                              initial =sof_basic_inclination[0],Extra_Model = Extra_Model  )

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=False)
    if 'INCL' in Configuration['FIXED_PARAMETERS'][0]:
        ax_INCL.text(1.01,0.5,'Forced Flat', rotation =-90,va='center',ha='left', color='black',transform = ax_INCL.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)
    plt.ylabel('Incl ($^{\circ}$)',**labelfont)
    arcmin,arcmax = ax_INCL.get_xlim()
    sec_ax = ax_INCL.twiny()
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.set_xticklabels([])
    sec_ax.figure.canvas.draw()

# ------------------------------PA------------------------------------
    ax_PA = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[25:28,3:9],\
                            Overview,'PA',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                            initial = sof_basic_pa[0],Extra_Model = Extra_Model )

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=True)
    if 'PA' in Configuration['FIXED_PARAMETERS'][0]:
        ax_PA.text(1.01,0.5,'Forced Flat', va='center',ha='left', color='black',rotation = -90, transform = ax_PA.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)
    plt.xlabel('Radius (arcsec)',**labelfont)
    plt.ylabel('PA ($^{\circ}$)',**labelfont)
    arcmin,arcmax = ax_PA.get_xlim()
    sec_ax = ax_PA.twiny()
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.set_xticklabels([])
    sec_ax.figure.canvas.draw()

# ------------------------------SDIS------------------------------------
    ax_SDIS = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[19:22,12:18],\
                              Overview,'SDIS',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                              Extra_Model = Extra_Model,initial = 8. )
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=False,
        labeltop= True)
    if 'SDIS' in Configuration['FIXED_PARAMETERS'][0]:
        ax_SDIS.text(1.01,0.5,'Forced Flat',rotation=-90, va='center',ha='left', color='black',transform = ax_SDIS.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)

    plt.ylabel('Disp (km s$^{-1}$)',**labelfont)
    arcmin,arcmax = ax_SDIS.get_xlim()
    sec_ax = ax_SDIS.twiny()
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.figure.canvas.draw()
    sec_ax.set_xlabel('Radius (kpc)',va='bottom',**labelfont)


# ------------------------------Scale height------------------------------------
    ax_Z0 = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[22:25,12:18],\
                            Overview,'Z0',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                            Extra_Model = Extra_Model,initial = sf.convertskyangle(Configuration,0.2,physical=True) )

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=False,
        labeltop=False)
    if 'Z0' in Configuration['FIXED_PARAMETERS'][0]:
        ax_Z0.text(1.25,0.5,'Forced Flat',rotation=-90, va='center',ha='left', color='black',transform = ax_Z0.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)
    plt.ylabel('Z0 (arcsec)',**labelfont)
    arcmin,arcmax = ax_Z0.get_ylim()
    sec_ax = ax_Z0.twinx()
    sec_ax.set_ylim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.figure.canvas.draw()
    sec_ax.set_ylabel('Z0 (kpc)',rotation=-90,va='bottom',**labelfont)
    arcmin,arcmax = ax_Z0.get_xlim()
    sec_ax = ax_Z0.twiny()
    sec_ax.set_xticklabels([])
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.figure.canvas.draw()


# ------------------------------SBR------------------------------------
    ax_SBR = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[25:28,12:18],\
                             Overview,'SBR',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                             Extra_Model = Extra_Model )

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labeltop=False)


    plt.ylabel('SBR \n (Jy km s$^{-1}$ arcsec$^{-2}$)',**labelfont)
    plt.xlabel('Radius (arcsec)',**labelfont)

    jymin,jymax = ax_SBR.get_ylim()
    sec_ax = ax_SBR.twinx()
    sec_ax.set_ylim(sf.columndensity(Configuration,jymin*1000.,arcsquare = True)/1e20,sf.columndensity(Configuration,jymax*1000.,arcsquare = True)/1e20)
    sec_ax.set_ylabel('Col. Dens. \n (x10$^{20}$ cm$^{-2}$)',rotation=-90,va='bottom',**labelfont)
    sec_ax.figure.canvas.draw()
    if 'SBR' in Configuration['FIXED_PARAMETERS'][0]:
        ax_SBR.text(1.25,0.5,'Forced Gaussian',rotation=-90, va='center',ha='left', color='black',transform = ax_SBR.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10*size_factor)
    sec_ax = ax_SBR.twiny()
    arcmin,arcmax = ax_SBR.get_xlim()
    sec_ax.set_xticklabels([])
    sec_ax.set_xlim(sf.convertskyangle(Configuration,arcmin),sf.convertskyangle(Configuration,arcmax))
    sec_ax.figure.canvas.draw()

#
#----------------------------------------------Distance vs VSYS -----------------------------------------
    ax_VSYS = Overview.add_subplot(gs[2:6,16:20])
    plt.xlabel('Sys. Vel. (km s$^{-1}$)',**labelfont)
    plt.ylabel('Distance (Mpc)',**labelfont)

    plt.errorbar(float(FAT_Model[Vars_to_plot.index('VSYS'),0]),float(Configuration['DISTANCE'])\
                ,xerr=float(FAT_Model[Vars_to_plot.index('# VSYS_ERR'),0]), c='k',zorder= 3,fmt="o")
    #plt.scatter(float(FAT_Model[0,Vars_to_plot.index('VSYS')]),float(Configuration['DISTANCE']),c='k',zorder= 3)
    if len(Extra_Model) > 0:
        plt.scatter(float(Extra_Model[Vars_to_plot.index('VSYS'),0]),float(Configuration['DISTANCE']),c='r',alpha = 0.5,zorder=1)
    if len(Input_Model) > 0:
        plt.scatter(float(Input_Model[Vars_to_plot.index('VSYS'),0]),float(Configuration['DISTANCE']),c='b',zorder= 2)
    plt.scatter(sof_basic_vsys[0],float(Configuration['DISTANCE']),marker='x',alpha=0.5, c = 'k')
    xmin,xmax = ax_VSYS.get_xlim()
    ymin,ymax = ax_VSYS.get_ylim()

    #plt.scatter(sof_basic_vsys[0],float(Configuration['DISTANCE']),marker='x',alpha=0.5, c = 'k')
    vall= np.linspace(0,15000,100)
    Distanc=vall/70.
    plt.plot(vall,Distanc,'k--',alpha=0.5)
    ax_VSYS.set_xlim(xmin, xmax)
    left = float(FAT_Model[Vars_to_plot.index('VSYS'),0])-cube[0].header['CDELT3']/2000.
    bottom = ymin-5
    width =  cube[0].header['CDELT3']/1000.
    height = ymax-ymin+50
#plt.fill([sof_basic_vsys[0]-Cube[0].header['CDELT3']/2000.,sof_basic_vsys[0]+Cube[0].header['CDELT3']/2000.],[ymin-50,ymax+50],c='k',alpha=0.5)
    why = plt.Rectangle((left,bottom), width,height,facecolor='black',alpha=0.4 , zorder=1)
    ax_VSYS.add_patch(why)
    ax_VSYS.set_ylim(ymin, ymax)
#----------------------------------------------RA vs DEC -----------------------------------------
    ax_RAD = Overview.add_subplot(gs[8:12,16:20])

    plt.errorbar(float(FAT_Model[Vars_to_plot.index('XPOS'),0]),\
                 float(FAT_Model[Vars_to_plot.index('YPOS'),0]),
                xerr=float(FAT_Model[Vars_to_plot.index('# XPOS_ERR'),0]),\
                yerr=float(FAT_Model[Vars_to_plot.index('# YPOS_ERR'),0]), c='k',zorder= 3,fmt="o")

    plt.scatter(float(FAT_Model[Vars_to_plot.index('XPOS'),0]),float(FAT_Model[Vars_to_plot.index('YPOS'),0]),c='k',zorder=3,label = 'Final')
    if len(Extra_Model) > 0:
        lab = 'Unsmoothed'
        alpha =0.5
        plt.scatter(float(Extra_Model[Vars_to_plot.index('XPOS'),0]),float(Extra_Model[Vars_to_plot.index('YPOS'),0]),c='r',zorder=1,alpha=alpha,label=lab)
    if len(Input_Model) > 0:
        plt.scatter(float(Input_Model[Vars_to_plot.index('XPOS'),0]),float(Input_Model[Vars_to_plot.index('YPOS'),0]),c='b',zorder =2,label='Input')
    plt.scatter(sof_basic_ra[0],sof_basic_dec[0],marker='x',alpha=0.5, c = 'k',label='Initial')
    mod_ell = Ellipse(xy=[float(FAT_Model[Vars_to_plot.index('XPOS'),0]),float(FAT_Model[Vars_to_plot.index('YPOS'),0])], width=cube[0].header['BMAJ'] , height=cube[0].header['BMAJ'], angle=0,
               edgecolor='none', alpha=0.4, lw=4, facecolor='k', hatch = None, zorder=-1)
    ax_RAD.add_patch(mod_ell)
    ax_RAD.legend(loc='upper left', bbox_to_anchor=(0.0, -0.3), shadow=True, ncol=1)
    plt.xlabel('RA ($^{\circ}$)',**labelfont)
    plt.ylabel('DEC ($^{\circ}$)',**labelfont)
    cube_mod.close()
    cube.close()
    channels_map.close()

    #
    plt.savefig(f"{Configuration['FITTING_DIR']}Overview.png", bbox_inches='tight')
    #plt.savefig(f"Overview_Test.png", bbox_inches='tight')
    plt.close()

make_overview_plot.__doc__ =f'''
 NAME:
    make_overview_plot(Configuration,Fits_Files ):

 PURPOSE:
    Create a plot that shows the various stages of the fitting and output in a handy overview

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Locations of the Fits files used by FAT

 OPTIONAL INPUTS:


 OUTPUTS:
    The overview plot

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def plot_parameters(Configuration,Vars_to_plot,FAT_Model,location,Figure,\
        parameter, Input_Model = None,legend = None,initial = None, \
        initial_extent=  None, Extra_Model = None):
    if Input_Model is None:
        Input_Model = []
    if Extra_Model is None:
        Extra_Model = []
    if legend is None:
        legend = ['Empty','Empty','Empty','Empty']
    print_log(f'''PLOT_PARAMETERS: We are starting to plot {parameter}
''', Configuration,case=['debug_start'])
    ax = Figure.add_subplot(location)
    try:
        yerr = FAT_Model[Vars_to_plot.index(f'# {parameter}_ERR'),:]
    except ValueError:
        yerr =np.zeros(len(FAT_Model[Vars_to_plot.index('RADI'),:]))
    print_log(f'''PLOT_PARAMETERS: We found these errors {yerr}
''', Configuration,case=['debug_add'])

    ax.errorbar(FAT_Model[Vars_to_plot.index('RADI'),:],FAT_Model[Vars_to_plot.index(f'{parameter}'),:],yerr= yerr, c ='k', label=f'{legend[0]}',zorder=3)
    ax.plot(FAT_Model[Vars_to_plot.index('RADI'),:],FAT_Model[Vars_to_plot.index(f'{parameter}'),:],'ko', ms = 3.,zorder=3)
    if np.sum(FAT_Model[Vars_to_plot.index(f'{parameter}_2'),:]) != 0.:
        diff = np.sum(abs(FAT_Model[Vars_to_plot.index(f'{parameter}_2'),:]-FAT_Model[Vars_to_plot.index(f'{parameter}'),:]))
        if diff != 0.:
            try:
                yerr = FAT_Model[Vars_to_plot.index(f'# {parameter}_2_ERR'),:]
            except ValueError:
                yerr =np.zeros(len(FAT_Model[Vars_to_plot.index('RADI'),:]))
            ax.errorbar(FAT_Model[Vars_to_plot.index('RADI'),:],FAT_Model[Vars_to_plot.index(f'{parameter}_2'),:],yerr= yerr, c ='r', label=f'{legend[1]}',zorder=3)
            ax.plot(FAT_Model[Vars_to_plot.index('RADI'),:],FAT_Model[Vars_to_plot.index(f'{parameter}_2'),:],'ro', ms = 3.,zorder=3)

    if len(Input_Model) > 0:
        if np.sum(Input_Model[Vars_to_plot.index(f'{parameter}'),:]) != 0.:
            last_index = int(np.where(Input_Model[Vars_to_plot.index(f'{parameter}'),:] != 0.)[0][-1])
            try:
                yerr = Input_Model[Vars_to_plot.index(f'# {parameter}_ERR'),:last_index]
            except ValueError:
                yerr =np.zeros(len(Input_Model[Vars_to_plot.index('RADI'),:last_index]))
            ax.errorbar(Input_Model[Vars_to_plot.index('RADI'),:last_index],Input_Model[Vars_to_plot.index(f'{parameter}'),:last_index],yerr= yerr, c ='b',linestyle='-', label=f'{legend[2]}',zorder=2)
            ax.plot(Input_Model[Vars_to_plot.index('RADI'),:last_index],Input_Model[Vars_to_plot.index(f'{parameter}'),:last_index],'bo',linestyle='-', ms = 3.,zorder=2)
        if np.sum(Input_Model[Vars_to_plot.index(f'{parameter}_2'),:]) != 0.:
            diff = np.sum(abs(Input_Model[Vars_to_plot.index(f'{parameter}_2'),:]-Input_Model[Vars_to_plot.index(f'{parameter}'),:]))
            if diff != 0.:
                last_index = int(np.where(Input_Model[Vars_to_plot.index(f'{parameter}_2'),:] != 0.)[0][-1])
                try: ###### Keep swapping indices
                    yerr = Input_Model[Vars_to_plot.index(f'# {parameter}_2_ERR'),:last_index]
                except ValueError:
                    yerr =np.zeros(len(Input_Model[Vars_to_plot.index('RADI'),:last_index]))

                ax.errorbar(Input_Model[Vars_to_plot.index('RADI'),:last_index],Input_Model[Vars_to_plot.index(f'{parameter}_2'),:last_index],yerr= yerr, c ='yellow', label=f'{legend[3]}',zorder=2)
                ax.plot(Input_Model[Vars_to_plot.index('RADI'),:last_index,],Input_Model[Vars_to_plot.index(f'{parameter}_2'),:last_index],'yellow',zorder=2,marker ='o',linestyle='-' , ms = 3.)
    ymin,ymax = ax.get_ylim()
    if len(Extra_Model) > 0:
        ax.plot(Extra_Model[Vars_to_plot.index('RADI'),:],Extra_Model[Vars_to_plot.index(f'{parameter}'),:],'ko',linestyle='-', ms = 3., alpha=0.2,zorder=1)
        if np.sum(Extra_Model[Vars_to_plot.index(f'{parameter}_2'),:]) != 0.:
            diff = np.sum(abs(Extra_Model[Vars_to_plot.index(f'{parameter}_2'),:]-Extra_Model[Vars_to_plot.index(f'{parameter}'),:]))
            if diff != 0.:
                ax.plot(Extra_Model[Vars_to_plot.index('RADI'),:],Extra_Model[Vars_to_plot.index(f'{parameter}_2'),:],'r', alpha=0.2,zorder=1,marker ='o',linestyle='-' , ms = 3.)
    xmin,xmax = ax.get_xlim()
    if initial:
        ax.plot([xmin-5,xmax+5],[float(initial),float(initial)],c='k',alpha=0.4, linestyle ='--')
    if initial_extent:
        ax.plot([float(initial_extent)/2.,float(initial_extent)/2.],[ymin-5,ymax+5],c='k',alpha=0.4, linestyle ='--')
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)
    return ax

plot_parameters.__doc__ =f'''
     NAME:
        plot_parameters

     PURPOSE:
        Sub plot a specific parameter from the fit

     CATEGORY:
        write_functions

     INPUTS:
        Configuration = Standard FAT configuration
        Vars_to_plot =List with order of the variables in the FAT_Model
        FAT_Model = Array with the values for all parameters
        location = location in the figure where to put this subplot
        Figure = the figure in which we are plotting
        parameter = the parameter to plot

     OPTIONAL INPUTS:


        Input_Model = []
        A model to compare to, should be ordered the same as the FAT_Model

        legend = ['Empty','Empty','Empty','Empty']
        legend names for fat_model and Input model

        initial = 'No Value'
        Initial guess of the parameter

        Extra_Model = []
        A third model to plot weakly behing the main model

     OUTPUTS:
        ax object from matplot lib

     OPTIONAL OUTPUTS:

     PROCEDURES CALLED:
        Unspecified

     NOTE:
'''

def plot_PV(Configuration,image=None, model = None, figure = None, \
    location = [0.1,0.1,0.8,0.8], tirific_model=None, size_factor = 1.):
    if image == None:
        print('plot_PV will not work witout an image.')
        return 'Empty'
    else:
        try:
            hdr = image[0].header
            data = image[0].data
        except:
            hdr = image.header
            data = image.data

    ratio=hdr['NAXIS2']/hdr['NAXIS1']
    # Then we want to plot our PV-Diagram
    if figure == None:
        figure = plt.figure(2, figsize=(8, 8), dpi=300, facecolor='w', edgecolor='k')

    ax = figure.add_subplot(location)
    #Comp_ax2.set_title('PV-Diagram')

    maxint= np.nanmax(data)*0.85
    minint= np.nanmin(data)/3.
    PV_plot = ax.imshow(data,  cmap='hot_r', origin='lower', alpha=1, \
                            vmin=minint, vmax=maxint,aspect='auto')
    xaxis = [hdr['CRVAL1'] + (i - hdr['CRPIX1'] + 1) \
        * (hdr['CDELT1']) for i in range(hdr['NAXIS1'])]
    yaxis = [hdr['CRVAL2'] + (i - hdr['CRPIX2'] + 1) * (hdr['CDELT2']) for i in
         range(hdr['NAXIS2'])]
    #something is going wrong here, Fixed as usual astropy was fucking thing up
    step = int(abs((xaxis[-1]-xaxis[0])/7))
    ticks = np.array([-3*step,-2*step,-1*step,0.,step,2*step,3*step],dtype=float)+hdr['CRVAL1']
    pix_ticks =(ticks - hdr['CRVAL1'])/hdr['CDELT1']+(hdr['CRPIX1']-1)
    ax.set_xticks(pix_ticks)
    ax.set_xticklabels([f'{int(i):d}' for i in ticks],size=5,ha='center')
    ax.set_yticks(range(len(yaxis))[0:-1:int(len(yaxis) / 5)])
    ax.set_yticklabels(['{:.1f}'.format(i) for i in yaxis[0:-1:int(len(yaxis) / 5)]])
    ax.grid()
   
    #Add some contours
    neg_cont = np.array([-3,-1.5],dtype=float)*Configuration['NOISE']
    if  np.nanmax(data) * 0.95 < 96*Configuration['NOISE']:
        pos_cont =  np.array([1.5,3.,6,12,24,48,96],dtype=float)*Configuration['NOISE']
    else:
        pos_cont =  np.array([0,1,2,3,4,5,6,7],dtype=float)*(np.nanmax(data) * 0.95-3.*Configuration['NOISE'])/7. +3.*Configuration['NOISE']
    pos_cont = np.array([x for x in pos_cont if x <= np.nanmax(data) * 0.95],dtype=float)
    print_log(f'''PV_PLOT: postive {pos_cont}, negative {neg_cont}, noise {Configuration['NOISE']}
''',Configuration,case=['debug_add'] )
    if pos_cont.size == 0:
        pos_cont = 0.5 * np.nanmax(data) * 0.95

#ax_PV.contour(PV[0].data, levels=pos_cont, colors='k',transform=ax_PV.get_transform(xv_proj))
#ax_PV.contour(PV[0].data, levels=neg_cont, colors='grey',linestyles='--',transform=ax_PV.get_transform(xv_proj))
#ax_PV.contour(PV_model[0].data, levels=pos_cont, colors='b',transform=ax_PV.get_transform(xv_model_proj),linewidths=1.)
    ax.contour(data, levels=pos_cont,linewidths=1.*size_factor, colors='k')
    ax.contour(data, levels=neg_cont,linewidths=1.*size_factor, colors='grey',linestyles='--')
    if model != None:
        try:
            model_data = model[0].data
        except:
            model_data = model.data
        ax.contour(model_data, levels=pos_cont, colors='b',linewidths=1.*size_factor)

    momlevel = np.hstack((neg_cont,pos_cont))
    column_levels = ', '.join(["{:.1f}".format(x*1000.) for x in momlevel])
    if len(momlevel) < 4:
        info_string = f"The contours are at {column_levels} mJy/beam"
    else:
        info_string = f"The contours are at {', '.join(['{:.1f}'.format(x*1000.) for x in momlevel[0:4]])}"
        counter = 5
        while counter < len(momlevel):
            info_string = info_string+f"\n {', '.join(['{:.1f}'.format(x*1000.) for x in momlevel[counter:counter+7]])}"
            counter += 7
    info_string = info_string+" mJy/beam."

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(PV_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([minint, maxint])
    cbar.ax.set_title(f"{hdr['BUNIT']}", y= 0.2*size_factor**2)

    ax.text(-0.1,-0.2,info_string, va='top',ha='left', color='black',transform = ax.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7)
    if tirific_model != None:
        parameters = sf.load_tirific(Configuration,tirific_model,\
            Variables= ['RADI','VROT','VROT_2','INCL','INCL_2','VSYS','VSYS_2']\
                ,array=True ,brightness_check=True)
        if np.sum(parameters[2]) == 0.:
            parameters[2] =   parameters[1]
            parameters[4] =   parameters[3]
            parameters[6] =   parameters[5]
        if 'DRVAL1' in hdr:
            center = [float(x) for x in hdr['DRVAL1'].split('+')]
            model_center = [float(x[0]) for x in sf.load_tirific(Configuration,tirific_model,\
            Variables= ['XPOS','XPOS_2','YPOS','YPOS_2'],array=True ) ]
        if 'DRVAL2' in hdr:
            vsys = float(hdr['DRVAL2'])
            model_vsys = [float(x[0]) for x in sf.load_tirific(Configuration,tirific_model,\
            Variables= ['VSYS','VSYS_2'],array=True ) ]

    
     
                            
        plotrc  = np.array([-1.*x*np.sin(np.radians(y))+parameters[5][0] for x,y \
            in zip(parameters[1],parameters[3])],dtype=float)
        plotrc2  = np.array([x*np.sin(np.radians(y))+parameters[6][0] for x,y in\
            zip(parameters[2],parameters[4])],dtype=float)
        # As python is a piece of shit that only does things half we need to convert to pixels

        plotrc= (plotrc - hdr['CRVAL2'])/hdr['CDELT2']+(hdr['CRPIX2']-1)
        radius = (-1.*parameters[0]-hdr['CRVAL1'])/hdr['CDELT1']+(hdr['CRPIX1']-1)
        plotrc2= (plotrc2 - hdr['CRVAL2'])/hdr['CDELT2']+(hdr['CRPIX2']-1)
        radius2 = (parameters[0]-hdr['CRVAL1'])/hdr['CDELT1']+(hdr['CRPIX1']-1)
        ax.plot(radius,plotrc,'o',c='r')
        ax.plot(radius2,plotrc2,'o',c='r')
    #cf.plot_fits(filename, Comp_ax2, cmap='hot_r', aspect=ratio, cbar ='horizontal')
    ax.set_xlabel("Offset (arcsec)")
    ax.set_ylabel("Velocity (km s$^{-1}$)")
    return ax
plot_PV.__doc__ =f'''
     NAME:
        plot_PV

     PURPOSE:
        Plot the PV with th RC overlaid.

     CATEGORY:
        write_functions

     INPUTS:
        Configuration = Standard FAT configuration
        image = fits object of the PV

     OPTIONAL INPUTS:
        model = fits object of model to overplot
        figure = figure to plot the ax in if not defined it will be an 8 x 8 figure \
        location = [0.1,0.1,0.8,0.8]
            location of the ax
        tirific_model=None
            Model to obtain the RC from
        size_factor = 1.
            scaling factor for text and line widths



     OUTPUTS:
        ax object from matplot lib

     OPTIONAL OUTPUTS:

     PROCEDURES CALLED:
        Unspecified

     NOTE:
'''

def plot_individual_ax(Configuration,ax,combined_time,combined_loads):
    ax.plot(combined_time,combined_loads['Tirific']['MEM'],'b-',lw=0.5)
    ax.plot(combined_time,combined_loads['FAT']['MEM'],'b--',lw=0.5)
    ax.set_ylim(0,np.max([combined_loads['Tirific']['MEM'],combined_loads['FAT']['MEM']]) \
                      +np.max([combined_loads['Tirific']['MEM'],combined_loads['FAT']['MEM']])/10.)
   
   
    ax2 = ax.twinx()
    ax2.plot(combined_time,combined_loads['Tirific']['CPU'],'r-',lw=0.5)
    ax2.plot(combined_time,combined_loads['FAT']['CPU'],'r--',lw=0.5)
    return ax,ax2
plot_individual_ax.__doc__ =f'''
     NAME:
        plot_individual_ax

     PURPOSE:
        setup a double plot with 2 axes 

     CATEGORY:
        write_functions

     INPUTS:
        Configuration = Standard FAT configuration
        ax = the axis object to hold the double plot
        combined_time = list containing the x axis
        combined_load = dictionary that that contains the Memory (blue) and CPU (red) loads split
                         according to tirific (solid line)  and FAT (dashed line)

     OPTIONAL INPUTS:


     OUTPUTS:
        ax = axis object with the memory plot
        ax2 = axis object with the CPU plot

     OPTIONAL OUTPUTS:

     PROCEDURES CALLED:
        Unspecified

     NOTE:
'''

def plot_usage_stats(Configuration ):
    with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt") as file:
        lines = file.readlines()

    
    try:
        mpl_fm.fontManager.addfont(Configuration['FONT_FILE'])
        font_name = mpl_fm.FontProperties(fname=Configuration['FONT_FILE']).get_name()
    except FileNotFoundError:
        font_name = 'DejaVu Sans'
    labelfont = {'family': font_name,
             'weight': 'normal',
             'size': 4}
    loads,labels,maxCPU, maxMEM = read_statistics(Configuration)
   
    # Below thanks to P. Serra
    # Make single-PID figures and total figure
    #print(loads['Tirific']['Time'],loads['FAT']['Time'])
    if len(loads['Tirific']['Time']) > 0.:
        combined_time, combined_loads, labels_times, labels_comb = \
                create_plot_stats(Configuration,loads,labels)
          
        if Configuration['MULTIPROCESSING']:
            fig, ax = plt.subplots(1,2,sharey=True,figsize = (8,6))
            ax[0].spines['right'].set_visible(False)
            ax[1].spines['left'].set_visible(False)
            ax[0].yaxis.tick_left()
            ax[0].tick_params(labelright='off')
            ax[1].yaxis.tick_right()
          
           
            
        else:
            fig, ax = plt.subplots(figsize = (8,6))
           
            ax = [ax]
        fig.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.3, top = 0.7)
        left_bottom = 0.1
        tot_time = 0.
        for timeone in combined_time:
            tot_time += timeone[-1]-timeone[0]
        normalize = tot_time/0.8
        lengths = [np.max([(x[-1]-x[0])/normalize,0.15]) for x in combined_time]
        lengths = np.array(lengths,dtype=float) * 0.8/np.sum(lengths) 
        for i,axplot in enumerate(ax):
            normalize = (labels_times[i][-1]-labels_times[i][0])/lengths[i]
            ax_MEM,ax_CPU = plot_individual_ax(Configuration,axplot,combined_time[i],combined_loads[i])
            ax_MEM,ax_CPU,left_bottom = set_proper_edges(Configuration, i, ax, ax_MEM,\
                                         ax_CPU,lengths[i],left_bottom)
            ax_CPU.set_ylim(0.,maxCPU*1.05)
            ax_MEM.set_ylim(0.,maxMEM*1.05)
            ax_MEM.set_xlim(combined_time[i][0],combined_time[i][-1])
            ax_CPU.set_xlim(combined_time[i][0],combined_time[i][-1])
                
            ax_CPU = set_timing_labels(Configuration,ax_CPU, labels_times[i],labels_comb[i],labels,labelfont,normalize)
        fig.text(0.5,0.25,'time (min)', va='center',ha='center',rotation= 0, color='black',
                          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7)
        fig.savefig(f"{Configuration['LOG_DIRECTORY']}ram_cpu.pdf")
        plt.close()

plot_usage_stats.__doc__ =f'''
 NAME:
    plot_usage_stats

 PURPOSE:
    Plot the RAM and CPU usage over time for the fit

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:


 OUTPUTS:
    ram_cpu.pdf in the Logs directory

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def read_statistics(Configuration):
    labels = {'FAT': {'label':[], 'Time':[]}, 'Tirific':{'label':[], 'Time':[]}}
    loads = {'FAT':{'CPU':[],'MEM':[],'Time':[]},'Tirific':{'CPU':[],'MEM':[],'Time': []}}
    with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt") as file:
        lines = file.readlines()
    current_stage = 'Not_Found'
    #current_module = 'Unknown'
    startdate = 0
    maxCPU = 0.
    maxMEM= 0.
    paused_time = 0.
    gap_time= 0.
    for line in lines:
        line = line.strip()
        tmp = line.split(' ')
        if line[0] == '#':
            date = extract_date(f"{tmp[-2]} {tmp[-1]}")
        else:
            date = extract_date(f"{tmp[0]} {tmp[1]}")
        if startdate == 0:
            startdate = date
        diff = date - startdate
        time = diff.total_seconds()/60.-gap_time
    
        if line[0] == '#':
            if tmp[1] == 'TIRIFIC:':

                if tmp[2].lower() == 'initializing':
                    tmp2 = line.split('=')[1].split()
                    current_stage = tmp2[0].strip()
                    labels['Tirific']['label'].append(f'Initializing {current_stage}')
                    labels['Tirific']['Time'].append(time)
                elif tmp[2].lower() == 'finished':
                    labels['Tirific']['label'].append(f'Ended {current_stage}')
                    labels['Tirific']['Time'].append(time)
                    #current_stage = 'No Tirific'
                elif tmp[2].lower() == 'started':
                    labels['Tirific']['label'].append(f'Started {current_stage}')
                    labels['Tirific']['Time'].append(time)
            elif tmp[1] == 'MP_FITTING_LOOP:':
                labels['Tirific']['label'].append(f'Started {tmp[1]}')  
                if Configuration['MULTIPROCESSING']:
                    gap_time = time-paused_time
                labels['Tirific']['Time'].append(time-gap_time)
            else:
                if tmp[2].lower() == 'pause':
                    labels['FAT']['label'].append(f'Pausing FAT')
                    labels['FAT']['Time'].append(time)
                    if Configuration['MULTIPROCESSING']:
                        '''If we are pausing due to splitting sofia and fit we subtract the pause time'''
                        paused_time = time
                else:
                    labels['FAT']['label'].append(f'Starting {tmp[1]}')
                    labels['FAT']['Time'].append(time)
        else:
            if tmp[-1].lower() == 'tirific':
                loads['Tirific']['Time'].append(time)
                loads['Tirific']['CPU'].append(tmp[4])
               
                loads['Tirific']['MEM'].append(tmp[8])
            else:
                loads['FAT']['Time'].append(time)
                loads['FAT']['CPU'].append(tmp[4])
                loads['FAT']['MEM'].append(tmp[8])
            if float(tmp[4]) > maxCPU:
                    maxCPU = float(tmp[4])
            if float(tmp[8]) > maxMEM:
                    maxMEM = float(tmp[8])
    return loads,labels,maxCPU,maxMEM
read_statistics.__doc__ =f'''
 NAME:
    read_statistics

 PURPOSE:
    Read the file Usage_Statistics.txt in the FAT log directory and transform the input into two 
    dictionaries. In Multiprocessing the pool time is subtracted
  

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    

 OPTIONAL INPUTS:


 OUTPUTS:
    loads = dictionary with the times and Memory and CPU loads, split in system and tirific
    labels = labels indicating specifics of what was happening at times
    maxCPU = the max CPU load encountered
    maxMEM = the maximum memory load encountered

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def reorder_output_catalogue(Configuration,Full_Catalogue):
    #Read in the output catalogue
    with open(Configuration['OUTPUT_CATALOGUE']) as file:
        lines = file.readlines()
    #Determine column sizes
    header = lines[0]
    output = lines[1:]
    #IDs cannot have spaces
    outputIDs = []
    for line in output:
        outputIDs.append(line.split()[0].strip())
    #Sort based on the inpu IDs
    with open(Configuration['OUTPUT_CATALOGUE'],'w') as file:
        file.write(header)
        for galaxy in Full_Catalogue['DIRECTORYNAME']:
            try:
                index_no = np.where(galaxy == \
                    np.array(outputIDs))[0][0]
                file.write(output[index_no])
            except IndexError:
                pass

reorder_output_catalogue.__doc__ =f'''
 NAME:
    reorder_output_catalogue

 PURPOSE:
    When running in multiprocessing mode the output catalogue can be
    in a different order as the input. This function makes sure the are sorted

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Full_Catalogue = Full input catalogue

 OPTIONAL INPUTS:


 OUTPUTS:
    reorder output catalogue

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_proper_edges(Configuration,i,ax,ax_MEM,ax_CPU,length,left_bottom):
    if i+1 == len(ax):
        if i > 0:
            ax_MEM.spines['left'].set_visible(False) 
            ax_MEM.tick_params(axis='y', length = 0)    
        ax_CPU.spines['left'].set_visible(False)
        ax_CPU.set_ylabel('CPUs (%)',color='r')
        ax_CPU.tick_params(axis='y', labelcolor='r')
    else:
        '''This currently only works with a single break'''
        ax_CPU.spines['right'].set_visible(False)
        ax_CPU.axis('off')    
        ax_CPU.tick_params(labelleft='off')
    if i == 0:
        ax_MEM.set_ylabel('RAM (Mb) ', color='b')
        ax_MEM.tick_params(axis='y', labelcolor='b')
    ax_MEM.set_position([left_bottom,0.3,length,0.4])
    ax_CPU.set_position([left_bottom,0.3,length,0.4])
    left_bottom += length+0.01  

    d =0.025
    dx = d/(6*length)
    kwargs = dict(transform=ax_MEM.transAxes, color='k', clip_on=False)
    if i != 0:
        ax_MEM.plot((-dx, dx), (1-d, 1+d), **kwargs)
        ax_MEM.plot((-dx, dx), (-d, +d), **kwargs)
    
    if i+1 != len(ax):
        ax_MEM.plot((1-dx, 1+dx), (-d, +d), **kwargs)
        ax_MEM.plot((1-dx, 1+dx), (1-d, 1+d), **kwargs)   
    return ax_MEM,ax_CPU,left_bottom
set_proper_edges.__doc__ =f'''
 NAME:
    set_proper_edges

 PURPOSE:
    Set the proper edges form the timing plot  
 
 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    i = index of current ax
    ax = all axes
    ax_MEM = the memory axis
    ax_CPU = the CPU axis
    length = size of the plot in fig coordinate
    left_bottom = location of the right side of the plot in fig coordinate

 OPTIONAL INPUTS:


 OUTPUTS:
    ax_MEM = updated memory axis object
    ax_CPU = updated CPU axis object
    left_bottom = right hand side for next axis object

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def set_timing_labels(Configuration,ax, labels_times,labels_comb,labels,labelfont,normalize):
    labelfont['size'] = 6
    miny, maxy = ax.get_ylim()
    last_label = -100

    label_sep = 0.015*normalize
    #label_sep = 0.001
    color, linest = '0.5', '--'
   
    prev_label = ''

    last_label_top = -100.
    last_label_bottom = -100.
        #for label,time in zip(labels['Tirific']['label'],labels['Tirific']['Time']):
    for label,time in zip(labels_comb,labels_times):
        if label in labels['Tirific']['label']:
            offset = 20.
            xoffset = 0.005*normalize
            vertical_start = maxy
            va= 'bottom'
            ha= 'left'
            color = 'k'
            linest = '-'
        else:
            offset=-50
            xoffset = -0.005*normalize
            vertical_start = miny
            va= 'top'
            ha='right'
            color = '0.5'
            linest = '--'


        if (prev_label == 'Initializing tmp_incl_check' or prev_label == 'Ended tmp_incl_check'):
            if (label != 'Initializing tmp_incl_check' and label != 'Ended tmp_incl_check') or \
                    time == labels_times[-1]    :

                if time != labels_times[-1]:
                    ax.axvline(x=prev_time, linestyle=linest, color=color, linewidth=0.05)
                    if label in labels['Tirific']['label']:
                        last_label = last_label_top = max(prev_time,last_label_top+label_sep)
                    else:
                        last_label = last_label_bottom = max(prev_time,last_label_bottom+label_sep)
                    ax.text(last_label,vertical_start+offset,prev_label, va=va,ha=ha,rotation= 60, color='black',
                            bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                    ax.plot([prev_time,last_label+xoffset],[vertical_start,vertical_start+offset],linest,color=color,linewidth=0.05,clip_on=False)
                ax.axvline(x=time, linestyle=linest, color=color, linewidth=0.05)
                if label in labels['Tirific']['label']:
                    last_label = last_label_top = max(time,last_label_top+label_sep)
                else:
                    last_label = last_label_bottom = max(time,last_label_bottom+label_sep)
                ax.text(last_label,vertical_start+offset,label, va=va,ha=ha,rotation= 60, color='black',
                        bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                ax.plot([time,last_label+xoffset],[vertical_start,vertical_start+offset],linest,color=color,linewidth=0.05,clip_on=False)
            else:
                prev_time = time
        elif (prev_label == 'Initializing Error_Shaker' or prev_label == 'Ended Error_Shaker' or prev_label == 'Started Error_Shaker'):
            if (label != 'Initializing Error_Shaker' and label != 'Ended Error_Shaker' and  label != 'Started Error_Shaker') or \
                    time == labels_times[-1]:

                    #ax2.axvline(x=prev_time, linestyle=linest, color=color, linewidth=0.05)
                    #last_label = max(prev_time,last_label+label_sep)
                    #ax2.text(last_label,ax2maxy+20.,prev_label, va='bottom',ha='left',rotation= 60, color='black',
                    #      bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                    #ax2.plot([prev_time,last_label+0.1],[ax2maxy,ax2maxy+20.],linest,color=color,linewidth=0.05,clip_on=False)
                ax.axvline(x=time, linestyle=linest, color=color, linewidth=0.05)
                if label in labels['Tirific']['label']:
                    last_label = last_label_top = max(time,last_label_top+label_sep)
                else:
                    last_label = last_label_bottom = max(time,last_label_bottom+label_sep)
                ax.text(last_label,vertical_start+offset,label, va=va,ha=ha,rotation= 60, color='black',
                          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                ax.plot([time,last_label+xoffset],[vertical_start,vertical_start+offset],linest,color=color,linewidth=0.05,clip_on=False)
            else:
                prev_label = label
                prev_time = time
        else:
            ax.axvline(x=time, linestyle=linest, color=color, linewidth=0.05)
          
            if label in labels['Tirific']['label']:
                last_label = last_label_top = max(time,last_label_top+label_sep)
            else:
                last_label = last_label_bottom= max(time,last_label_bottom+label_sep)
           
            #exit()
            ax.text(last_label,vertical_start+offset,label,va=va,ha=ha,rotation= 60, color='black',
                      bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                #This should be the line to the label
            ax.plot([time,last_label+xoffset],[vertical_start,vertical_start+offset],linest,color=color,linewidth=0.05,clip_on=False)
        prev_label = label

    ax.set_ylim(miny,maxy)    
    return ax

set_timing_labels.__doc__ =f'''
NAME:
set_timing_labels

PURPOSE:
Set the timing labels for the various accournces

CATEGORY:
write_functions

INPUTS:
Configuration = Standard FAT configuration
ax = the axis object to attach the labels to
labels_times = list of label times
labels_comb = list of labels
labels = dictionary with all labels split into upper and lower 
labelfont = label font settings
normalize = conversion factor from image coordinate to plot

OPTIONAL INPUTS:


OUTPUTS:
ax = label axis object

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
Unspecified

NOTE:
'''

def sofia(template,name):
    with open(name,'w') as file:
        for key in template:
            if key[0] == 'E' or key [0] == 'H':
                file.write(template[key])
            else:
                file.write(f"{key} = {template[key]}\n")

sofia.__doc__ =f'''
NAME:
sofia
PURPOSE:
write a sofia2 dictionary into file
CATEGORY:
write_functions

INPUTS:
template = sofia template
name = name of the file to write to

OPTIONAL INPUTS:

OUTPUTS:

OPTIONAL OUTPUTS:

PROCEDURES CALLED:
Unspecified

NOTE:
'''

def square_plot(ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    if xmax > ymax:
        diff = int(xmax-ymax)/2.
        ax.set_ylim(ymin-diff,ymax+diff)
        ymin, ymax = ax.get_ylim()
    else:
        diff = int(ymax-xmax)/2.
        ax.set_xlim(xmin-diff,xmax+diff)
        xmin, xmax = ax.get_xlim()
square_plot.__doc__ =f'''
 NAME:
    square_plot

 PURPOSE:
    square the axes object
        
 CATEGORY:
    write_functions

 INPUTS:
    ax = is the axes object to be squared
        
 OPTIONAL INPUTS:


 OUTPUTS:
    square axes
 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:

 NOTE:
'''

def tirific(Configuration,Tirific_Template, name = 'tirific.def',\
                full_name = False  ):
    #IF we're writing we bump up the restart_ID and adjust the AZ1P angles to the current warping
    update_disk_angles(Configuration,Tirific_Template )
    try:
        Tirific_Template['RESTARTID'] = str(int(Tirific_Template['RESTARTID'])+1)
    except ValueError:
        Tirific_Template['RESTARTID'] = 0 
    except KeyError:
        Tirific_Template['RESTARTID'] = 0
    if full_name:
        file_name = name
    else:
        file_name = f'{Configuration["FITTING_DIR"]}{name}'
    with open(file_name, 'w') as file:
        for key in Tirific_Template:
            if key[0:5] == 'EMPTY':
                file.write('\n')
            else:
                file.write((f"{key}= {Tirific_Template[key]} \n"))
tirific.__doc__ =f'''
 NAME:
    tirific

 PURPOSE:
    Write a tirific template to file

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Tirific_Template = Standard FAT Tirific Template

 OPTIONAL INPUTS:


    name = 'tirific.def'
    name of the file to write to

 OUTPUTS:
    Tirific def file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 '''
