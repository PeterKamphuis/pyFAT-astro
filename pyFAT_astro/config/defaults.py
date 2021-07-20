# -*- coding: future_fstrings -*-

from dataclasses import dataclass,field
import omegaconf
from omegaconf import MISSING
from typing import List,Optional
import os
from datetime import datetime

@dataclass
class Fitting:
    catalogue_start_id: int = -1 #Catalogue ID of the first galaxy to be fitted. -1 Means start at the beginning
    catalogue_end_id: int= -1 #the last galaxy to be fitted, if set to -1 the whole catalogue will be fitted
    fitting_stages: List = field(default_factory=lambda: ['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC'])
    # Possible stages are
    # Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory
    # Catalogue_Sofia: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted
    # Run_Sofia: Run Sofia on the FAT cube and  process the output
    # Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet
    # Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth
    ring_size: float = 1.1 # The size of the rings in number of beams
    fixed_parameters:  List = field(default_factory=lambda: ['Z0','XPOS','YPOS','VSYS']) #Options are INCL, PA, SDIS, SBR
    opt_pixel_beam: int=4
    ncpu: int = 6
    max_iterations: int=15
    distance: float = -1. # Distance to the galaxy, set from the catalogue at start of loop in case of batch fitting

@dataclass
class Input:
    main_directory: str =f'{os.getcwd()}'
    hanning_smoothed: bool = False # is the input cube hanning smoothed
    catalogue: Optional[str] = None
    tirific: str = "tirific" #Command to call tirific
    sofia2: str = "sofia2"   #Command to call sofia 2
    sofia_basename: Optional[str] = None
    sofia_dir: Optional[str] = None #Directory of the existing Sofia output. Only used if the input catalogue is sofia or pre-processed sofia is somewhere

@dataclass
class Output:
    log_directory: str = f'Logs/{datetime.now().strftime("%d-%m-%Y")}' # Name of the log dir
    log_file: str = 'Log.txt' #Name of the log file
    catalogue: str = 'FAT_results.txt'
    new_output: bool = True # Create all output anew
    # How much output you would like to maintain for each galaxy. 0 just organize the output and keep all (This will also happen when a fit is unsuccesful, this can be a lot of files); 1 remove optimized files, log files and input files; 2  remove optimized files, log files, input files, ps files and unsmoothed files; 3 (Default) remove optimized files, log files, input files, ps files, unsmoothed files and all model fits files except the final model; 4 keep only the def files and remove all other output. 5 indicates a failed fit clean up. >6 is the same as 0. Residuals are created for all cases where the fits files are maintained.
    output_quantity: int = 3
    warp_output: bool = False #If you want FAT to output a warp radius, tiltograms and warp radius set warp_output (Default = n)
    debug: bool = False
    timing: bool = False

@dataclass
class Advanced:
    start_directory: str =f'{os.getcwd()}'
    max_iterations: int=15 #The maximum number of iterations that FAT tries to calls trific bfeore it call it quits
    loops: int =10 #The number of full loops set for tirific in a  single iteration
    minimum_warp_size: float = 3. # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc,set here.
    minimum_rings: int = 3  # we need at least this amount of rings (Including 0 and 1/5 beam), set here
    too_small_galaxy: float = 1. # if the number of beams across the major axis/2 is less than this we will not fit the galaxy, set here

@dataclass
class defaults:
    print_examples: bool=False
    installation_check: bool=False
    cube_name: Optional[str] = None
    configuration_file: Optional[str] = None
    input: Input = Input()
    output: Output = Output()
    fitting: Fitting = Fitting()
    advanced: Advanced = Advanced()
