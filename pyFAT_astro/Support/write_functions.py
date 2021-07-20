# -*- coding: future_fstrings -*-
# This module contains a set of functions and classes that are used to write text files to Disk

from pyFAT_astro.Support.support_functions import print_log,convertRADEC,convertskyangle,set_limit_modifier,columndensity,set_limits,get_inner_fix,linenumber
from pyFAT_astro.Support.modify_template import set_model_parameters, set_overall_parameters, set_fitting_parameters,get_warp_slope, update_disk_angles
from pyFAT_astro.Support.fits_functions import extract_pv
from pyFAT_astro.Support.read_functions import load_tirific,load_basicinfo, load_template
import numpy as np
import warnings
import datetime
import os
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    from matplotlib.patches import Ellipse
    import matplotlib.axes as maxes
from astropy.io import fits
from astropy.wcs import WCS
# create or append to the basic ifo file
def basicinfo(Configuration,initialize = False,stage='TiRiFiC', debug = False,
              RA=[float('NaN'),float('NaN')], DEC=[float('NaN'),float('NaN')],
              VSYS =[float('NaN'),float('NaN')], PA=[float('NaN'),float('NaN')],
              Inclination = [float('NaN'),float('NaN')], Max_Vrot = [float('NaN'),float('NaN')],
              Tot_Flux = [float('NaN'),float('NaN')], V_mask = [float('NaN'),float('NaN')],
              Distance = float('NaN') , DHI = float('NaN'), template = ['EMPTY']):

    try:
        if template[0] == 'EMPTY':
            template = {}
    except KeyError:
        pass
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
        FAT_Model = load_template(Configuration,template,Variables= Vars_to_Set,unpack=False, debug=debug)
        RA=[FAT_Model[0,Vars_to_Set.index('XPOS')],abs(Configuration['BEAM'][0]/(3600.*2.))]
        DEC=[FAT_Model[0,Vars_to_Set.index('YPOS')],abs(Configuration['BEAM'][0]/(3600.*2.))]
        VSYS =np.array([FAT_Model[0,Vars_to_Set.index('VSYS')],Configuration['CHANNEL_WIDTH']],dtype=float)
        PA=[FAT_Model[0,Vars_to_Set.index('PA')], 3.]
        Inclination = [FAT_Model[0,Vars_to_Set.index('INCL')], 3.]
        Max_Vrot = [np.max(FAT_Model[:,Vars_to_Set.index('VROT')]),np.max(FAT_Model[:,Vars_to_Set.index('VROT')])-np.min(FAT_Model[1:,Vars_to_Set.index('VROT')]) ]
        Distance = Configuration['DISTANCE']

    with open(f"{Configuration['FITTING_DIR']}{Configuration['BASE_NAME']}-Basic_Info.txt",'a') as file:
        if not initialize:
            file.write(f'''#These are the values from {stage}. \n''')
        RAhr,DEChr = convertRADEC(Configuration,RA[0],DEC[0],debug=debug)
        RA_c = f'{RAhr}+/-{RA[1]*3600.:0.2f}'
        DEC_c = f'{DEChr}+/-{DEC[1]*3600.:0.2f}'
        VSYS_c = f'{VSYS[0]:.2f}+/-{VSYS[1]:.2f}'
        PA_c = f'{PA[0]:.2f}+/-{PA[1]:.2f}'
        INCL_c = f'{Inclination[0]:.2f}+/-{Inclination[1]:.2f}'
        MVROT_c = f'{Max_Vrot[0]:.2f}+/-{Max_Vrot[1]:.2f}'
        Vmask_c = f'{V_mask[0]/1000.:.2f}+/-{V_mask[1]/1000.:.2f}'
        DHI_a = f'{DHI[0]:.2f}+/-{DHI[1]:.2f}'
        Dist = f'{Distance:.2f}'
        HIMass  = f'{Tot_Flux[0]*2.36E5*Distance**2:.2e}'
        DHI_k = f'{convertskyangle(Configuration,DHI[0],Distance):.2f}'
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
    debug = False

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

# Function to write the first def file for a galaxy
def initialize_def_file(Configuration, Fits_Files,Tirific_Template,Initial_Parameters = ['EMPTY'], fit_type = 'Undefined',debug = False):
    try:
        if Initial_Parameters[0] == 'EMPTY':
            Initial_Parameters = {}
    except KeyError:
        pass

    #First we set some basic parameters that will hardly change
    if fit_type == 'Centre_Convergence':

        if 'VSYS' in Initial_Parameters:
            Initial_Parameters['VSYS'] = [x/1000. for x in Initial_Parameters['VSYS']]
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template,fit_type=fit_type, flux = Initial_Parameters['FLUX'][0], debug=debug)
        # Then set the values for the various parameters of the model

        set_model_parameters(Configuration, Tirific_Template,Initial_Parameters, debug=debug)

        set_limit_modifier(Configuration,Initial_Parameters['INCL'][0], debug=debug )

        set_fitting_parameters(Configuration, Tirific_Template,stage = 'initial',
                                initial_estimates = Initial_Parameters, debug=debug)
        if 'VSYS' in Initial_Parameters:
            Initial_Parameters['VSYS'] = [x*1000. for x in Initial_Parameters['VSYS']]
    elif fit_type in ['Extent_Convergence','Fit_Tirific_OSC']:
        if 'VSYS' in Initial_Parameters and fit_type == 'Fit_Tirific_OSC':
            Initial_Parameters['VSYS'] = [x/1000. for x in Initial_Parameters['VSYS']]
        set_overall_parameters(Configuration, Fits_Files,Tirific_Template ,fit_type=fit_type, debug=debug,stage='initialize_ec')
        Vars_to_Set =  ['XPOS','YPOS','VSYS','VROT','INCL','PA','SDIS','SBR','SBR_2','Z0']
        if fit_type == 'Fit_Tirific_OSC':
            set_model_parameters(Configuration, Tirific_Template,Initial_Parameters,stage='initialize_def_file', debug=debug)
        FAT_Model = load_template(Configuration,Tirific_Template,Variables= Vars_to_Set,unpack=False, debug=debug)

        # Finally we set how these parameters are fitted.
        set_limit_modifier(Configuration,FAT_Model[0,Vars_to_Set.index('INCL')], debug=debug)
        get_inner_fix(Configuration,Tirific_Template, debug=debug)
        get_warp_slope(Configuration,Tirific_Template, debug=debug)

        parameters = {'VSYS': [FAT_Model[0,Vars_to_Set.index('VSYS')], Configuration['CHANNEL_WIDTH']], \
                      'XPOS': [FAT_Model[0,Vars_to_Set.index('XPOS')], Configuration['BEAM'][0]/3600.] ,
                      'YPOS': [FAT_Model[0,Vars_to_Set.index('YPOS')], Configuration['BEAM'][0]/3600.],
                      'INCL': [FAT_Model[0,Vars_to_Set.index('INCL')], 3.],
                      'PA':  [FAT_Model[0,Vars_to_Set.index('PA')], 3.],
                      'VROT':[np.mean(FAT_Model[:,Vars_to_Set.index('VROT')]),np.max(FAT_Model[:,Vars_to_Set.index('VROT')])-np.min(FAT_Model[1:,Vars_to_Set.index('VROT')]) ]  }
        set_fitting_parameters(Configuration, Tirific_Template,stage = 'initialize_os',
                               initial_estimates=parameters, debug=debug)

    tirific(Configuration,Tirific_Template,name = f'{fit_type}_In.def', debug=debug)

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
    debug = False

    Initial_Parameters = ['EMPTY']
    The initial parameters to be used to set up the def file.

    fit_type = 'Undefined'
    type of fitting being done

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def make_overview_plot(Configuration,Fits_Files, debug = False):
    fit_type = Configuration['USED_FITTING']
    if debug:
        print_log(f'''MAKE_OVERVIEW_PLOT: We are starting the overview plot.
''',Configuration['OUTPUTLOG'],debug =True )

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        cube_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel.fits")
        moment0_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom0.fits")
        moment1_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom1.fits")
        moment2_mod = fits.open(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_mom2.fits")
        cube = fits.open(f"{Configuration['FITTING_DIR']}{Fits_Files['FITTING_CUBE']}")
        moment0 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT0']}")
        moment1 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT1']}")
        moment2 = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['MOMENT2']}")
        channels_map = fits.open(f"{Configuration['FITTING_DIR']}/Sofia_Output/{Fits_Files['CHANNEL_MAP']}")
        im_wcs = WCS(moment0[0].header)

    # Open the model info
    Vars_to_plot= ['RADI','XPOS','YPOS','VSYS','VROT','VROT_ERR','VROT_2','VROT_2_ERR','INCL','INCL_ERR','INCL_2',
                    'INCL_2_ERR','PA','PA_ERR','PA_2','PA_2_ERR','SDIS','SDIS_ERR','SDIS_2','SDIS_2_ERR','SBR',
                    'SBR_2','Z0','Z0_2','Z0_ERR','Z0_2_ERR']
    FAT_Model = load_tirific(Configuration,f"{Configuration['FITTING_DIR']}Finalmodel/Finalmodel.def",Variables= Vars_to_plot,unpack=False,debug=debug)
    Extra_Model_File = f"{Configuration['FITTING_DIR']}{fit_type}/{fit_type}_Iteration_{Configuration['ITERATIONS']}.def"

    if os.path.exists(Extra_Model_File):
        Extra_Model = load_tirific(Configuration,Extra_Model_File,Variables= Vars_to_plot,unpack=False,debug=debug)
    else:
        Extra_Model = []
    if debug:
        print_log(f'''MAKE_OVERVIEW_PLOT: We find the following model values.
{'':8s}{[f"{x} = {FAT_Model[:,i]}" for i,x in enumerate(Vars_to_plot)]}
''',Configuration['OUTPUTLOG'])

    if os.path.exists(f"{Configuration['FITTING_DIR']}ModelInput.def"):
        Input_Model = load_tirific(Configuration,f"{Configuration['FITTING_DIR']}ModelInput.def",Variables= Vars_to_plot,unpack=False,debug=debug)
    else:
        Input_Model = []
    sof_basic_ra,sof_basic_dec, sof_basic_vsys,sof_basic_maxrot,sof_basic_pa,sof_basic_inclination,sof_basic_extent = load_basicinfo(Configuration,
        f"{Configuration['FITTING_DIR']}{Configuration['BASE_NAME']}-Basic_Info.txt",Variables=['RA','DEC','VSYS','Max VRot','PA','Inclination','D_HI'])


    #Let's start plotting

    Overview = plt.figure(2, figsize=(8.2, 11.6), dpi=300, facecolor='w', edgecolor='k')

    size_ratio = 11.6/8.2
    #stupid pythonic layout for grid spec, which means it is yx instead of xy like for normal human beings
    gs = Overview.add_gridspec(int(20*size_ratio),20)
    labelfont = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 8}
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
      bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=14)
    ax_text.text(0.5,0.25,f'''The ring size used in the model is {Configuration['RING_SIZE']:.2f} x BMAJ, with BMAJ = {Configuration['BEAM'][0]:.1f} arcsec. We assumed a distance  of {Configuration['DISTANCE']:.1f} Mpc.'''
      ,rotation=0, va='top',ha='center', color='black',transform = ax_text.transAxes,
      bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)


#-----------------------------------------------------------------Moment 0 ------------------------------------------------------
    ax_moment0 = Overview.add_subplot(gs[2:8,0:6], projection=im_wcs)
    ax_moment0.set_label('Intensity Map')
    #Comp_ax1.set_facecolor('black')
    # we need contour levels and
    min_color = 0.
    max_color = np.nanmax(moment0[0].data)*0.8
    moment0_plot = ax_moment0.imshow(moment0[0].data, origin='lower', alpha=1, vmin = min_color, vmax = max_color,cmap='hot_r' )
    moment0_plot.set_label('Intensity Map')
    plt.ylabel('DEC (J2000)')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA J2000')

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
    else:
        momlevel = np.array([1,4,8,32,64,128],dtype=float) * mindism0
    #print("We find this {} as the minimum of the moment0 map".format(mindism0))
    momlevel = np.array([x for x in momlevel if x < np.max(moment0[0].data)*0.95],dtype=float)
    if momlevel.size == 0:
        momlevel=0.5*mindism0
    ax_moment0.contour(moment0[0].data, transform=ax_moment0.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.5 , zorder =4)
    ax_moment0.contour(moment0[0].data, transform=ax_moment0.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=1.2)
    ax_moment0.contour(moment0_mod[0].data, transform=ax_moment0.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2 , zorder =7)
    ax_moment0.contour(moment0_mod[0].data, transform=ax_moment0.get_transform(im_wcs),
              levels=momlevel, colors='r',zorder=8, linewidths=0.9)
    xmin, xmax = ax_moment0.get_xlim()
    ymin, ymax = ax_moment0.get_ylim()
    if xmax > ymax:
        diff = int(xmax-ymax)/2.
        ax_moment0.set_ylim(ymin-diff,ymax+diff)
        ymin, ymax = ax_moment0.get_ylim()
    else:
        diff = int(ymax-xmax)/2.
        ax_moment0.set_xlim(xmin-diff,xmax+diff)
        xmin, xmax = ax_moment0.get_xlim()
    ghxloc, ghyloc = im_wcs.wcs_pix2world(float(xmin+(xmax-xmin)/18.), float(ymin+(ymax-ymin)/18.), 1.)
    localoc = [float(ghxloc),float(ghyloc) ]
    widthb = moment0[0].header['BMIN']
    heightb = moment0[0].header['BMAJ']
    try:
        angleb  = moment0[0].header['BPA']
    except:
        angleb = 0.
    beam = Ellipse(xy=localoc, width=widthb, height=heightb, angle=angleb, transform=ax_moment0.get_transform('fk4'),
           edgecolor='k', lw=1, facecolor='none', hatch='/////',zorder=15)
    ax_moment0.add_patch(beam)
    # colorbar
    divider = make_axes_locatable(ax_moment0)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment0_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    cbar.ax.set_title(f"{moment0[0].header['BUNIT']}", y= 0.2)


    column_levels = columndensity(Configuration,momlevel*1000.,systemic = FAT_Model[0,Vars_to_plot.index('VSYS')])

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
    inclination_correction = set_limits(FAT_Model[0,Vars_to_plot.index('INCL')]+12.5,20.,90.)
    velocity_width= set_limits(1.25*np.nanmax(FAT_Model[:,Vars_to_plot.index('VROT')])*np.sin(np.radians(inclination_correction)),30.,700.)

    max_color= FAT_Model[0,Vars_to_plot.index('VSYS')]+velocity_width
    min_color= FAT_Model[0,Vars_to_plot.index('VSYS')]-velocity_width

    moment1_plot = ax_moment1.imshow(moment1[0].data, cmap='rainbow', origin='lower', alpha=1, vmin = min_color, vmax = max_color )
    plt.ylabel('DEC (J2000)')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA J2000')

    # contours
    velocity_step=set_limits(int((int(max_color-min_color)*0.9)/20.),1.,30.)
    integer_array = np.linspace(0,20,21)-10
    momlevel = [FAT_Model[0,Vars_to_plot.index('VSYS')]+x*velocity_step for x in integer_array if min_color < FAT_Model[0,Vars_to_plot.index('VSYS')]+x*velocity_step < max_color]
    ax_moment1.contour(moment1[0].data, transform=ax_moment1.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2 , zorder =4)
    ax_moment1.contour(moment1[0].data, transform=ax_moment1.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=0.9)
    ax_moment1.contour(moment1_mod[0].data, transform=ax_moment1.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2 , zorder =7)
    #ax_moment1.contour(moment1_mod[0].data, transform=ax_moment1.get_transform(im_wcs),
    #          levels=momlevel, colors='r',zorder=8, linewidths=0.9)
    xmin, xmax = ax_moment1.get_xlim()
    ymin, ymax = ax_moment1.get_ylim()
    if xmax > ymax:
        diff = int(xmax-ymax)/2.
        ax_moment1.set_ylim(ymin-diff,ymax+diff)
        ymin, ymax = ax_moment1.get_ylim()
        xmin, xmax = ax_moment1.get_xlim()
    else:
        diff = int(ymax-xmax)/2.
        ax_moment1.set_xlim(xmin-diff,xmax+diff)
    ghxloc, ghyloc = im_wcs.wcs_pix2world(float(xmin+(xmax-xmin)/18.), float(ymin+(ymax-ymin)/18.), 1.)
    localoc = [float(ghxloc),float(ghyloc) ]
    widthb = moment1[0].header['BMIN']
    heightb = moment1[0].header['BMAJ']
    try:
        angleb  = moment1[0].header['BPA']
    except:
        angleb = 0.
    beam = Ellipse(xy=localoc, width=widthb, height=heightb, angle=angleb, transform=ax_moment1.get_transform('fk4'),
           edgecolor='k', lw=1, facecolor='none', hatch='/////',zorder=15)
    ax_moment1.add_patch(beam)
    # colorbar
    divider = make_axes_locatable(ax_moment1)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment1_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    cbar.ax.set_title(f"{moment1[0].header['BUNIT']}", y= 0.2)

    column_levels = ', '.join(["{:.1f}".format(x) for x in momlevel])
    info_string= f'''The contours start at {float(momlevel[0]):.1f} km/s'''
    if len(momlevel) > 1:
        info_string=f'''{info_string} and increase with {float(momlevel[1])-float(momlevel[0]):.1f} km/s.'''
    else:
        info_string=f'''{info_string}.'''


    '''
    if len(column_levels) < 4:
        info_string = f"The contours are at {column_levels} km/s."
    else:
        info_string = f"The contours are at {', '.join(['{:.1f}'.format(x) for x in momlevel[0:4]])}"
        counter = 4
        while counter < len(momlevel):
            info_string = info_string+f"\n {', '.join(['{:.1f}'.format(x) for x in momlevel[counter:counter+7]])}"
            counter += 7
        info_string = info_string+" km/s."
    '''
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
    max_color= set_limits(np.nanmax(moment2[0].data),15,50)
    min_color= 0.

    moment2_plot = ax_moment2.imshow(moment2[0].data, cmap='rainbow' ,origin='lower', alpha=1, vmin = min_color, vmax = max_color )
    plt.ylabel('DEC (J2000)')
    #Stupid python suddenly finds its own labels
    plt.xlabel('RA J2000')

    # contours

    momlevel = np.linspace(min_color,max_color*0.8,5)

    ax_moment2.contour(moment2[0].data, transform=ax_moment2.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.5 , zorder =4)
    ax_moment2.contour(moment2[0].data, transform=ax_moment2.get_transform(im_wcs),
              levels=momlevel, colors='k',zorder=6, linewidths=1.2)
    ax_moment2.contour(moment2_mod[0].data, transform=ax_moment2.get_transform(im_wcs),
               levels=momlevel, colors='white',linewidths=1.2 , zorder =7)
    #ax_moment2.contour(moment2_mod[0].data, transform=ax_moment1.get_transform(im_wcs),
    #          levels=momlevel, colors='yellow',zorder=8, linewidths=0.9)
    xmin, xmax = ax_moment2.get_xlim()
    ymin, ymax = ax_moment2.get_ylim()
    if xmax > ymax:
        diff = int(xmax-ymax)/2.
        ax_moment2.set_ylim(ymin-diff,ymax+diff)
        ymin, ymax = ax_moment2.get_ylim()
    else:
        diff = int(ymax-xmax)/2.
        ax_moment2.set_xlim(xmin-diff,xmax+diff)
        xmin, xmax = ax_moment2.get_xlim()

    ghxloc, ghyloc = im_wcs.wcs_pix2world(float(xmin+(xmax-xmin)/18.), float(ymin+(ymax-ymin)/18.), 1.)
    localoc = [float(ghxloc),float(ghyloc) ]
    widthb = moment2[0].header['BMIN']
    heightb = moment2[0].header['BMAJ']
    try:
        angleb  = moment2[0].header['BPA']
    except:
        angleb = 0.

    beam = Ellipse(xy=localoc, width=widthb, height=heightb, angle=angleb, transform = ax_moment2.get_transform('fk4'),
           edgecolor='k', lw=1, facecolor='none', hatch='/////',zorder=15)
    ax_moment2.add_patch(beam)


    # colorbar
    divider = make_axes_locatable(ax_moment2)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(moment2_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([min_color, max_color])
    #cbar.set_title(label=f"{moment2[0].header['BUNIT']}")
    cbar.ax.set_title(f"{moment2[0].header['BUNIT']}", y= 0.2)

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

    extract_angle = np.mean(FAT_Model[0:round(len(FAT_Model[:,Vars_to_plot.index('PA')])/2.),Vars_to_plot.index('PA')])
    PV = extract_pv(Configuration,cube,extract_angle, \
                    center = [float(FAT_Model[0,Vars_to_plot.index('XPOS')]),float(FAT_Model[0,Vars_to_plot.index('YPOS')]),float(FAT_Model[0,Vars_to_plot.index('VSYS')]*1000.)], \
                    convert=1000.)
    if not os.path.exists(f"{Configuration['FITTING_DIR']}/Finalmodel/{Configuration['BASE_NAME']}_final_xv.fits"):
        fits.writeto(f"{Configuration['FITTING_DIR']}/Finalmodel/{Configuration['BASE_NAME']}_final_xv.fits",PV[0].data,PV[0].header)
    PV_model = extract_pv(Configuration,cube_mod,extract_angle, \
                    center = [float(FAT_Model[0,Vars_to_plot.index('XPOS')]),float(FAT_Model[0,Vars_to_plot.index('YPOS')]),float(FAT_Model[0,Vars_to_plot.index('VSYS')]*1000.)], \
                    convert=1000.)
    if not os.path.exists(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_xv.fits"):
        fits.writeto(f"{Configuration['FITTING_DIR']}/Finalmodel/Finalmodel_xv.fits",PV_model[0].data,PV_model[0].header)
    ratio=PV[0].header['NAXIS2']/PV[0].header['NAXIS1']
    # Then we want to plot our PV-Diagram
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        xv_proj = WCS(PV[0].header)
        xv_model_proj = WCS(PV_model[0].header)
    ax_PV = Overview.add_subplot(gs[2:8,8:14], projection=xv_proj)
    #Comp_ax2.set_title('PV-Diagram')

    maxint= np.nanmax(PV[0].data)*0.85
    minint= np.nanmin(PV[0].data)/3.


    PV_plot = ax_PV.imshow(PV[0].data,  cmap='hot_r', origin='lower', alpha=1, vmin=minint, vmax=maxint,aspect='auto')
    xaxis = [PV[0].header['CRVAL1'] + (i - PV[0].header['CRPIX1'] + 1) * (PV[0].header['CDELT1']) for i in
             range(PV[0].header['NAXIS1'])]
    yaxis = [PV[0].header['CRVAL2'] + (i - PV[0].header['CRPIX2'] + 1) * (PV[0].header['CDELT2']) for i in
             range(PV[0].header['NAXIS2'])]
    plt.gca().set_xticks(range(len(xaxis))[0:-1:int(len(xaxis) / 5)])
    plt.gca().set_yticks(range(len(yaxis))[0:-1:int(len(yaxis) / 5)])
    plt.gca().set_xticklabels(['{:10.0f}'.format(i) for i in xaxis[0:-1:int(len(xaxis) / 5)]])
    plt.gca().set_yticklabels(['{:10.1f}'.format(i) for i in yaxis[0:-1:int(len(yaxis) / 5)]])

    #Add some contours
    neg_cont = np.array([-3,-1.5],dtype=float)*Configuration['NOISE']
    pos_cont =  np.array([1.5,3.,6,12,24,48,96],dtype=float)*Configuration['NOISE']
    pos_cont = np.array([x for x in pos_cont if x < np.nanmax(PV[0].data) * 0.95],dtype=float)
    if debug:
            print_log(f'''MAKE_OVERVIEW_PLOT: postive {pos_cont}, negative {neg_cont}, noise {Configuration['NOISE']}
    ''',Configuration['OUTPUTLOG'],debug =True )
    if pos_cont.size == 0:
        pos_cont = 0.5 * np.nanmax(PV[0].data) * 0.95
    try:
        ax_PV.contour(PV[0].data, levels=pos_cont, colors='k',transform=ax_PV.get_transform(xv_proj))
        ax_PV.contour(PV[0].data, levels=neg_cont, colors='grey',linestyles='--',transform=ax_PV.get_transform(xv_proj))
        ax_PV.contour(PV_model[0].data, levels=pos_cont, colors='b',transform=ax_PV.get_transform(xv_model_proj),linewidths=1.)
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
    except:
        info_string = f"Something went wrong plotting the contours."
    divider = make_axes_locatable(ax_PV)
    cax = divider.append_axes("top", size="5%", pad=0.05, axes_class=maxes.Axes)
    cbar = plt.colorbar(PV_plot, cax=cax, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    cbar.set_ticks([minint, maxint])
    cbar.ax.set_title(f"{PV[0].header['BUNIT']}", y= 0.2)

    ax_PV.text(-0.1,-0.2,info_string, va='top',ha='left', color='black',transform = ax_PV.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7)

    #cf.plot_fits(filename, Comp_ax2, cmap='hot_r', aspect=ratio, cbar ='horizontal')
    ax_PV.set_xlabel("Offset (arcsec)")
    ax_PV.set_ylabel("Velocity (km s$^{-1}$)")
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
    labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':10}
    ax_RC = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[19:22,3:9],\
                            Overview,'VROT',Input_Model = Input_Model,initial_extent= sof_basic_extent[0], \
                            initial = sof_basic_maxrot[0],Extra_Model = Extra_Model,\
                            debug=debug)
    ymin =np.min([FAT_Model[1:,Vars_to_plot.index('VROT')],FAT_Model[1:,Vars_to_plot.index('VROT_2')]])
    if len(Extra_Model) > 0:
        ymin2 =np.min([Extra_Model[1:,Vars_to_plot.index('VROT')],Extra_Model[1:,Vars_to_plot.index('VROT_2')]])
        ymax2 =np.max([Extra_Model[1:,Vars_to_plot.index('VROT')],Extra_Model[1:,Vars_to_plot.index('VROT_2')]])
    else:
        ymin2 = ymin
        ymax2 = ymax
    ymax =np.max([FAT_Model[1:,Vars_to_plot.index('VROT')],FAT_Model[1:,Vars_to_plot.index('VROT_2')]])

    if len(Input_Model) > 0:
        ymin3 =np.min([Input_Model[1:,Vars_to_plot.index('VROT')],Input_Model[1:,Vars_to_plot.index('VROT_2')]])
        ymax3 =np.max([Input_Model[1:,Vars_to_plot.index('VROT')],Input_Model[1:,Vars_to_plot.index('VROT_2')]])
    else:
        ymin3=ymin
        ymax3= ymax
    ymin = np.min([ymin,ymin2,ymin3])
    ymax = np.max([ymax,ymax2,ymax3])

    buffer = np.mean([ymin,ymax])/20.
    ymin= ymin-buffer
    ymax = ymax+buffer

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
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.figure.canvas.draw()
    sec_ax.set_xlabel('Radius (kpc)',va='bottom',**labelfont)

# ------------------------------Inclination------------------------------------
    ax_INCL = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[22:25,3:9],\
                              Overview,'INCL',Input_Model = Input_Model,initial_extent= sof_basic_extent[0], \
                              initial =sof_basic_inclination[0],Extra_Model = Extra_Model ,debug=debug)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=False)
    if 'INCL' in Configuration['FIXED_PARAMETERS'][0]:
        ax_INCL.text(1.01,0.5,'Forced Flat', rotation =-90,va='center',ha='left', color='black',transform = ax_INCL.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)
    plt.ylabel('Incl ($^{\circ}$)',**labelfont)
    arcmin,arcmax = ax_INCL.get_xlim()
    sec_ax = ax_INCL.twiny()
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.set_xticklabels([])
    sec_ax.figure.canvas.draw()

# ------------------------------PA------------------------------------
    ax_PA = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[25:28,3:9],\
                            Overview,'PA',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                            initial = sof_basic_pa[0],Extra_Model = Extra_Model,debug=debug)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=True)
    if 'PA' in Configuration['FIXED_PARAMETERS'][0]:
        ax_PA.text(1.01,0.5,'Forced Flat', va='center',ha='left', color='black',rotation = -90, transform = ax_PA.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)
    plt.xlabel('Radius (arcsec)',**labelfont)
    plt.ylabel('PA ($^{\circ}$)',**labelfont)
    arcmin,arcmax = ax_PA.get_xlim()
    sec_ax = ax_PA.twiny()
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.set_xticklabels([])
    sec_ax.figure.canvas.draw()

# ------------------------------SDIS------------------------------------
    ax_SDIS = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[19:22,12:18],\
                              Overview,'SDIS',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                              Extra_Model = Extra_Model,initial = 8.,debug=debug)
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
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)

    plt.ylabel('Disp (km s$^{-1}$)',**labelfont)
    arcmin,arcmax = ax_SDIS.get_xlim()
    sec_ax = ax_SDIS.twiny()
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.figure.canvas.draw()
    sec_ax.set_xlabel('Radius (kpc)',va='bottom',**labelfont)


# ------------------------------Scale height------------------------------------
    ax_Z0 = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[22:25,12:18],\
                            Overview,'Z0',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                            Extra_Model = Extra_Model,initial = convertskyangle(Configuration,0.2,Configuration['DISTANCE'],physical=True),debug=debug)

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
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)
    plt.ylabel('Z0 (arcsec)',**labelfont)
    arcmin,arcmax = ax_Z0.get_ylim()
    sec_ax = ax_Z0.twinx()
    sec_ax.set_ylim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.figure.canvas.draw()
    sec_ax.set_ylabel('Z0 (kpc)',rotation=-90,va='bottom',**labelfont)
    arcmin,arcmax = ax_Z0.get_xlim()
    sec_ax = ax_Z0.twiny()
    sec_ax.set_xticklabels([])
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.figure.canvas.draw()


# ------------------------------SBR------------------------------------
    ax_SBR = plot_parameters(Configuration,Vars_to_plot,FAT_Model,gs[25:28,12:18],\
                             Overview,'SBR',Input_Model = Input_Model,initial_extent= sof_basic_extent[0],\
                             Extra_Model = Extra_Model,debug=debug)

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
    sec_ax.set_ylim(columndensity(Configuration,jymin*1000.,arcsquare = True)/1e20,columndensity(Configuration,jymax*1000.,arcsquare = True)/1e20)
    sec_ax.set_ylabel('Col. Dens. \n (x10$^{20}$ cm$^{-2}$)',rotation=-90,va='bottom',**labelfont)
    sec_ax.figure.canvas.draw()
    if 'SBR' in Configuration['FIXED_PARAMETERS'][0]:
        ax_SBR.text(1.25,0.5,'Forced Gaussian',rotation=-90, va='center',ha='left', color='black',transform = ax_SBR.transAxes,
          bbox=dict(facecolor='white',edgecolor='white',pad=0.,alpha=0.),zorder=7,fontsize=10)
    sec_ax = ax_SBR.twiny()
    arcmin,arcmax = ax_SBR.get_xlim()
    sec_ax.set_xticklabels([])
    sec_ax.set_xlim(convertskyangle(Configuration,arcmin,Configuration['DISTANCE']),convertskyangle(Configuration,arcmax,Configuration['DISTANCE']))
    sec_ax.figure.canvas.draw()


#----------------------------------------------Distance vs VSYS -----------------------------------------
    ax_VSYS = Overview.add_subplot(gs[2:6,16:20])
    plt.xlabel('Sys. Vel. (km s$^{-1}$)',**labelfont)
    plt.ylabel('Distance (Mpc)',**labelfont)
    plt.scatter(float(FAT_Model[0,Vars_to_plot.index('VSYS')]),float(Configuration['DISTANCE']),c='k',zorder= 3)
    if len(Extra_Model) > 0:
        plt.scatter(float(Extra_Model[0,Vars_to_plot.index('VSYS')]),float(Configuration['DISTANCE']),c='r',alpha = 0.5,zorder=1)
    if len(Input_Model) > 0:
        plt.scatter(float(Input_Model[0,Vars_to_plot.index('VSYS')]),float(Configuration['DISTANCE']),c='b',zorder= 2)
    plt.scatter(sof_basic_vsys[0],float(Configuration['DISTANCE']),marker='x',alpha=0.5, c = 'k')
    xmin,xmax = ax_VSYS.get_xlim()
    ymin,ymax = ax_VSYS.get_ylim()

    #plt.scatter(sof_basic_vsys[0],float(Configuration['DISTANCE']),marker='x',alpha=0.5, c = 'k')
    vall= np.linspace(0,15000,100)
    Distanc=vall/70.
    plt.plot(vall,Distanc,'k--',alpha=0.5)
    ax_VSYS.set_xlim(xmin, xmax)
    left = float(FAT_Model[0,Vars_to_plot.index('VSYS')])-cube[0].header['CDELT3']/2000.
    bottom = ymin-5
    width =  cube[0].header['CDELT3']/1000.
    height = ymax-ymin+50
#plt.fill([sof_basic_vsys[0]-Cube[0].header['CDELT3']/2000.,sof_basic_vsys[0]+Cube[0].header['CDELT3']/2000.],[ymin-50,ymax+50],c='k',alpha=0.5)
    why = plt.Rectangle((left,bottom), width,height,facecolor='black',alpha=0.4 , zorder=1)
    ax_VSYS.add_patch(why)
    ax_VSYS.set_ylim(ymin, ymax)
#----------------------------------------------RA vs DEC -----------------------------------------
    ax_RAD = Overview.add_subplot(gs[8:12,16:20])
    plt.scatter(float(FAT_Model[0,Vars_to_plot.index('XPOS')]),float(FAT_Model[0,Vars_to_plot.index('YPOS')]),c='k',zorder=3,label = 'Final')
    if len(Extra_Model) > 0:
        lab = 'Unsmoothed'
        alpha =0.5
        plt.scatter(float(Extra_Model[0,Vars_to_plot.index('XPOS')]),float(Extra_Model[0,Vars_to_plot.index('YPOS')]),c='r',zorder=1,alpha=alpha,label=lab)
    if len(Input_Model) > 0:
        plt.scatter(float(Input_Model[0,Vars_to_plot.index('XPOS')]),float(Input_Model[0,Vars_to_plot.index('YPOS')]),c='b',zorder =2,label='Input')
    plt.scatter(sof_basic_ra[0],sof_basic_dec[0],marker='x',alpha=0.5, c = 'k',label='Initial')
    mod_ell = Ellipse(xy=[float(FAT_Model[0,Vars_to_plot.index('XPOS')]),float(FAT_Model[0,Vars_to_plot.index('YPOS')])], width=cube[0].header['BMAJ'] , height=cube[0].header['BMAJ'], angle=0,
               edgecolor='none', alpha=0.4, lw=4, facecolor='k', hatch=' ',zorder=-1)
    ax_RAD.add_patch(mod_ell)
    ax_RAD.legend(loc='upper left', bbox_to_anchor=(0.0, -0.3), shadow=True, ncol=1)
    plt.xlabel('RA ($^{\circ}$)',**labelfont)
    plt.ylabel('DEC ($^{\circ}$)',**labelfont)
    cube_mod.close()
    cube.close()
    channels_map.close()

    plt.savefig(f"{Configuration['FITTING_DIR']}Overview.png", bbox_inches='tight')
    plt.close()

make_overview_plot.__doc__ =f'''
 NAME:
    make_overview_plot(Configuration,Fits_Files, debug = False):

 PURPOSE:
    Create a plot that shows the various stages of the fitting and output in a handy overview

 CATEGORY:
    write_functions

 INPUTS:
    Configuration = Standard FAT configuration
    Fits_Files = Locations of the Fits files used by FAT

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    The overview plot

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def plot_parameters(Configuration,Vars_to_plot,FAT_Model,location,Figure,parameter,\
                    Input_Model = [],legend = ['Empty','Empty','Empty','Empty'],
                    initial = None, initial_extent=  None, Extra_Model = [], debug = False):
    if debug:
        print_log(f'''PLOT_PARAMETERS: We are starting to plot {parameter}
''', Configuration['OUTPUTLOG'], debug = True)
    ax = Figure.add_subplot(location)
    try:
        yerr = FAT_Model[:,Vars_to_plot.index(f'{parameter}_ERR')]
    except:
        yerr =np.zeros(len(FAT_Model[:,Vars_to_plot.index('RADI')]))
    if debug:
        print_log(f'''PLOT_PARAMETERS: We found these errors {yerr}
''', Configuration['OUTPUTLOG'])

    ax.errorbar(FAT_Model[:,Vars_to_plot.index('RADI')],FAT_Model[:,Vars_to_plot.index(f'{parameter}')],yerr= yerr, c ='k', label=f'{legend[0]}',zorder=3)
    ax.plot(FAT_Model[:,Vars_to_plot.index('RADI')],FAT_Model[:,Vars_to_plot.index(f'{parameter}')],'ko', ms = 3.,zorder=3)
    if np.sum(FAT_Model[:,Vars_to_plot.index(f'{parameter}_2')]) != 0.:
        diff = np.sum(abs(FAT_Model[:,Vars_to_plot.index(f'{parameter}_2')]-FAT_Model[:,Vars_to_plot.index(f'{parameter}')]))
        if diff != 0.:
            try:
                yerr = FAT_Model[:,Vars_to_plot.index(f'{parameter}_2_ERR')]
            except:
                yerr =np.zeros(len(FAT_Model[:,Vars_to_plot.index('RADI')]))
            ax.errorbar(FAT_Model[:,Vars_to_plot.index('RADI')],FAT_Model[:,Vars_to_plot.index(f'{parameter}_2')],yerr= yerr, c ='r', label=f'{legend[1]}',zorder=3)
            ax.plot(FAT_Model[:,Vars_to_plot.index('RADI')],FAT_Model[:,Vars_to_plot.index(f'{parameter}_2')],'ro', ms = 3.,zorder=3)

    if len(Input_Model) > 0:
        if np.sum(Input_Model[:,Vars_to_plot.index(f'{parameter}')]) != 0.:
            last_index = int(np.where(Input_Model[:,Vars_to_plot.index(f'{parameter}')] != 0.)[0][-1])
            try:
                yerr = Input_Model[:last_index,Vars_to_plot.index(f'# {parameter}_ERR')]
            except:
                yerr =np.zeros(len(Input_Model[:last_index,Vars_to_plot.index('RADI')]))
            ax.errorbar(Input_Model[:last_index,Vars_to_plot.index('RADI')],Input_Model[:last_index,Vars_to_plot.index(f'{parameter}')],yerr= yerr, c ='b',linestyle='-', label=f'{legend[2]}',zorder=2)
            ax.plot(Input_Model[:last_index,Vars_to_plot.index('RADI')],Input_Model[:last_index,Vars_to_plot.index(f'{parameter}')],'bo',linestyle='-', ms = 3.,zorder=2)
        if np.sum(Input_Model[:,Vars_to_plot.index(f'{parameter}_2')]) != 0.:
            diff = np.sum(abs(Input_Model[:,Vars_to_plot.index(f'{parameter}_2')]-Input_Model[:,Vars_to_plot.index(f'{parameter}')]))
            if diff != 0.:
                last_index = int(np.where(Input_Model[:,Vars_to_plot.index(f'{parameter}_2')] != 0.)[0][-1])
                try:
                    yerr = Input_Model[:last_index,Vars_to_plot.index(f'# {parameter}_2_ERR')]
                except:
                    yerr =np.zeros(len(Input_Model[:last_index,Vars_to_plot.index('RADI')]))

                ax.errorbar(Input_Model[:last_index,Vars_to_plot.index('RADI')],Input_Model[:last_index,Vars_to_plot.index(f'{parameter}_2')],yerr= yerr, c ='yellow', label=f'{legend[3]}',zorder=2)
                ax.plot(Input_Model[:last_index,Vars_to_plot.index('RADI')],Input_Model[:last_index,Vars_to_plot.index(f'{parameter}_2')],'yellow',zorder=2,marker ='o',linestyle='-' , ms = 3.)
    ymin,ymax = ax.get_ylim()
    if len(Extra_Model) > 0:
        ax.plot(Extra_Model[:,Vars_to_plot.index('RADI')],Extra_Model[:,Vars_to_plot.index(f'{parameter}')],'ko',linestyle='-', ms = 3., alpha=0.2,zorder=1)
        if np.sum(Extra_Model[:,Vars_to_plot.index(f'{parameter}_2')]) != 0.:
            diff = np.sum(abs(Extra_Model[:,Vars_to_plot.index(f'{parameter}_2')]-Extra_Model[:,Vars_to_plot.index(f'{parameter}')]))
            if diff != 0.:
                ax.plot(Extra_Model[:,Vars_to_plot.index('RADI')],Extra_Model[:,Vars_to_plot.index(f'{parameter}_2')],'r', alpha=0.2,zorder=1,marker ='o',linestyle='-' , ms = 3.)
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
        debug = False

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

def plot_usage_stats(Configuration,debug = False):
    with open(f"{Configuration['LOG_DIRECTORY']}Usage_Statistics.txt") as file:
        lines = file.readlines()
    labels = []
    label_times = []
    times = []
    CPU = []
    mem = []
    labelfont = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 4}
    current_stage = 'Not_Found'
    startdate = 0
    for line in lines:
        line = line.strip()
        if line[0] == '#':
            tmp = line.split('=')
            if len(tmp) == 2:
                current_stage = tmp[1].strip()
                labels.append(f'Initializing {current_stage}')
                label_times.append('No Time')
            else:
                tmp = line.split(' ')
                if tmp[1].lower() == 'finished':
                    labels.append(f'Ended {current_stage}')
                    label_times.append('No Time')
                elif tmp[3].lower() == 'actual':
                    labels.append(f'Started {current_stage}')
                    label_times.append('No Time')
                    times.append('No Time')
                    CPU.append(CPU[-1])
                    mem.append(mem[-1])
        else:
            tmp = line.split(' ')
            tmp2 = tmp[0].split('-')
            if len(tmp2) == 3:
                try:
                    date =  datetime.datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S.%f')
                except ValueError:
                    date =  datetime.datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S')

                if startdate == 0:
                    try:
                        startdate  =  datetime.datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S.%f')
                    except ValueError:
                        startdate  =  datetime.datetime.strptime(f"{tmp[0]} {tmp[1]}", '%Y-%m-%d %H:%M:%S')
                diff = date - startdate
                if len(times) > 0 :
                    if times[-1] == 'No Time':
                        times[-1] = diff.total_seconds()/60.
                times.append(diff.total_seconds()/60.)
                if label_times[-1] == 'No Time':
                    label_times[-1] = diff.total_seconds()/60.
                CPU.append(tmp[4])
                mem.append(tmp[8])
    # Below thanks to P. Serra
    # Make single-PID figures and total figure
    if len(mem) > 0.:
        fig, ax1 = plt.subplots(figsize = (8,4))
        fig.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.15, top = 0.7)
        times = np.array(times, dtype = float)
        mem = np.array(mem, dtype = float)
        CPU = np.array(CPU, dtype = float)

        ax1.plot(times,mem,'b-',lw=0.5)
        ax1.set_ylim(0,np.max(mem)+np.max(mem)/10.)
        ax1.set_ylabel('RAM (Mb) ', color='b')
        ax1.tick_params(axis='y', labelcolor='b')
        ax1.set_xlabel('time (min)', color='k')
        ax2 = ax1.twinx()
        ax2.plot(times,CPU,'r-',lw=0.5)
        ax2.set_ylabel('CPUs (%)',color='r')
        ax2miny,ax2maxy = ax2.get_ylim()
        ax2.tick_params(axis='y', labelcolor='r')
        last_label = -100
        label_sep = label_times[-1]/40.
        color, linest = '0.5', '--'
        labelfont = {'family': 'Times New Roman',
                 'weight': 'normal',
                 'size': 6.5}
        prev_label = ''
        for label,time in zip(labels,label_times):
            if color == '0.5':
                color = 'k'
            elif color == 'k':
                color = '0.5'
            if linest == '--':
                linest = '-'
            elif linest == '-':
                linest = '--'
            if (prev_label == 'Initializing tmp_incl_check' or prev_label == 'Ended tmp_incl_check'):
                if (label != 'Initializing tmp_incl_check' and label != 'Ended tmp_incl_check'):
                    ax2.axvline(x=prev_time, linestyle=linest, color=color, linewidth=0.05)
                    last_label = max(prev_time,last_label+label_sep)
                    ax2.text(last_label,ax2maxy+20.,prev_label, va='bottom',ha='left',rotation= 60, color='black',
                          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                    ax2.plot([prev_time,last_label+0.1],[ax2maxy,ax2maxy+20.],'k'+linest,color=color,linewidth=0.05,clip_on=False)
                    ax2.axvline(x=time, linestyle=linest, color=color, linewidth=0.05)
                    last_label = max(time,last_label+label_sep)
                    ax2.text(last_label,ax2maxy+20.,label, va='bottom',ha='left',rotation= 60, color='black',
                          bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                    ax2.plot([time,last_label+0.1],[ax2maxy,ax2maxy+20.],'k'+linest,color=color,linewidth=0.05,clip_on=False)
                else:
                    prev_time = time
            else:
                ax2.axvline(x=time, linestyle=linest, color=color, linewidth=0.05)
                last_label = max(time,last_label+label_sep)
                ax2.text(last_label,ax2maxy+20.,label, va='bottom',ha='left',rotation= 60, color='black',
                      bbox=dict(facecolor='white',edgecolor='white',pad= 0.,alpha=0.),zorder=7,fontdict = labelfont)
                ax2.plot([time,last_label+0.1],[ax2maxy,ax2maxy+20.],'k'+linest,color=color,linewidth=0.05,clip_on=False)
            prev_label = label
        #This is beyond stupid again, but hey it is python so needed to make things work.
        ax2.set_ylim([ax2miny,ax2maxy])
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
    debug = False

 OUTPUTS:
    ram_cpu.pdf in the Logs directory

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

def tirific(Configuration,Tirific_Template, name = 'tirific.def', debug = False):
    #IF we're writing we bump up the restart_ID and adjust the AZ1P angles to the current warping
    update_disk_angles(Configuration,Tirific_Template, debug= debug)
    Tirific_Template['RESTARTID'] = str(int(Tirific_Template['RESTARTID'])+1)
    with open(Configuration['FITTING_DIR']+name, 'w') as file:
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
    debug = False

    name = 'tirific.def'
    name of the file to write to

 OUTPUTS:
    Tirific def file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
 '''

def write_config(file,Configuration,debug = False):
    if debug:
        print_log(f'''WRITE_CONFIG: writing the configuration to {file}
''',Configuration['OUTPUTLOG'], screen = True)
    # Separate the keyword names
    with open(file,'w') as tmp:
        for key in Configuration:
            tmp.write(f'{key} = {Configuration[key]} \n')

write_config.__doc__ =f'''
 NAME:
    write_config

 PURPOSE:
    Write a config file to the fitting directory.

 CATEGORY:
    write_functions

 INPUTS:
    file = name of the file to write to
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    A FAT config file.

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''
