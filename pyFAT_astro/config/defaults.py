# -*- coding: future_fstrings -*-

import omegaconf
import os
import pyFAT_astro
from dataclasses import dataclass, field
from datetime import datetime
from omegaconf import MISSING
#from multiprocessing import cpu_count
import psutil
from typing import List, Optional


@dataclass
class Fitting:

    # Catalogue ID of the first galaxy to be fitted. -1 Means start at the beginning
    catalogue_start_id: str = '-1'
    # the last galaxy to be fitted, if set to -1 the whole catalogue will be fitted
    catalogue_end_id: str = '-1'
    fitting_stages: List = field(default_factory=lambda: [
                                 'Create_FAT_Cube', 'Run_Sofia', 'Fit_Tirific_OSC'])
    # Possible stages are
    # Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory
    # Sofia_Catalogue: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted
    # Run_Sofia: Run Sofia on the FAT cube and  process the output
    # Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet
    # Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth
    # Tirshaker: Bootstrap errors for the final model. This can take a long time
    ring_size: float = 1.1  # The size of the rings in number of beams
    fixed_parameters:  List = field(default_factory=lambda: [
                                    'Z0', 'XPOS', 'YPOS', 'VSYS'])  # Options are INCL, PA, SDIS, SBR
    opt_pixel_beam: int = 4

    distance: float = -1.  # Distance to the galaxy, set from the catalogue at start of loop in case of batch fitting


@dataclass
class Input:
    main_directory: str = f'{os.getcwd()}'
    # 'Options are independent, sinusoidal, hanning
    channel_dependency: str = 'independent'
    catalogue: Optional[str] = None
    tirific: str = "tirific"  # Command to call tirific
    sofia2: str = "sofia"  # Command to call sofia 2
    sofia_basename: Optional[str] = None
    # Directory of the existing Sofia output. Only used if the input catalogue is sofia or pre-processed sofia is somewhere
    sofia_dir: Optional[str] = None


@dataclass
class Output:
    # Name of the log dir
    log_directory: str = f'Logs/{datetime.now().strftime("%d-%m-%Y")}'
    log_file: str = 'Log.txt'  # Name of the log file
    catalogue: str = 'pyFAT_results.txt'
    new_output: bool = True  # Create all output anew
    # How much output you would like to maintain for each galaxy. 0 just organize the output and keep all (This will also happen when a fit is unsuccesful, this can be a lot of files); 1 remove optimized files, log files and input files; 2  remove optimized files, log files, input files, ps files and unsmoothed files; 3 (Default) remove optimized files, log files, input files, ps files, unsmoothed files and all model fits files except the final model; 4 keep only the def files and remove all other output. 5 indicates a failed fit clean up. >6 is the same as 0. Residuals are created for all cases where the fits files are maintained.
    output_quantity: int = 3
    # If you want FAT to output a warp radius, tiltograms and warp radius set warp_output (Default = n)
    warp_output: bool = False
    # The font to be used in the plots
    font_file: str = "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"
    debug: bool = False
    timing: bool = False
    verbose_log: bool = False
    verbose_screen: bool = False

@dataclass
class Advanced:
    start_directory: str = f'{os.getcwd()}'
    # The maximum number of iterations that FAT tries to calls trific bfeore it call it quits
    max_iterations: int = 15
    loops: int = 7  # The number of full loops set for tirific in a  single iteration, this is increased every time the maximum is reach upto max 15 loops
    # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc,set here.
    minimum_warp_size: float = 3.
    # we need at least this amount of rings (Including 0 and 1/5 beam), set here
    minimum_rings: int = 3
    too_small_galaxy: float = 1.  # if the number of beams across the major axis/2 is less than this we will not fit the galaxy, set here
    unreliable_size: float = 2.  # If the final diameter is smaller than this the fit is considered unreliable
    unreliable_inclination: float = 10.  # If the final inclination is below this the fit is considered unreliable
    shaker_iterations: int = 20
    multiprocessing: bool = True
    number_of_disks: int = 2
    # The whole mechanism beheind the sbr_limits (which weigh the smoothing) is emperically determined. They are multiplied with this factor to allow optimization
    limit_modifier_factor: float = 1.05
    #We do not want to use too many cores per galaxy.
    per_galaxy_ncpu: int = 4
    catalogue_split_character: str = '|'
    # Allow for the user to set the boundaries in the fitting
    pa_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    incl_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    sdis_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    z0_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    vsys_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    xpos_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    ypos_input_boundary:  List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    vrot_input_boundary: List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    sbr_input_boundary: List = field(
        default_factory=lambda: [[0., 0.], [0., 0.], [0., 0.]])
    # The brightest pixels need to have a SNR above this value
    source_max_snr: float = 2.5
    # The fraction of pixels in the source that need to be above max_snr
    source_max_fraction: float = 0.075
    # The mean SNR required in the source
    source_mean_snr: float = 0.75
    sofia_threshold: int = 5
    #option to bypass all initial checks on the intial sofia source and simply force the TRM fitting
    force_fit: bool = False
    # Add the channel dependency, minimum inclination,
    debug_function: List = field(default_factory=lambda: ['ALL'] )

    # If we are using a sofia catalogue and want tor overwrite the create FAT Cubes
    sofia_overwrite: bool = False

@dataclass
class defaults:
    try:
        ncpu: int = len(psutil.Process().cpu_affinity())
    except AttributeError:
        ncpu: int = psutil.cpu_count()

    print_examples: bool = False
    installation_check: bool = False
    cube_name: Optional[str] = None
    configuration_file: Optional[str] = None
    input: Input = Input()
    output: Output = Output()
    fitting: Fitting = Fitting()
    advanced: Advanced = Advanced()
