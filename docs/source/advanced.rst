Setting your fit preferences through a yaml file.
=================================

Introduction
--------

In comparision to FAT pyFAT delivers much more control over the fitting procedure. This is made possible by the new modular setup and the use of omegaconf for the input and default settings.
In this page we explain how the fitting in FAT can be adapted by providing a yaml input parameters for pyFAT. An example file with all the settings can be printed by runnning 'pyFAT print_examples=True' after installing pyFAT-astro. This command prints both an example yaml file and an example catalogue input file with the names FAT_defaults.yml and FAT_Input_Catalogue.txt, respectively
Below we explain what each section in the yaml file does.

The is divided  in four different sections and has 4 independent keywords. All keywords are optional but pyFAT does require at least a cube name or an input catalogure to work.  The section advanced and the independent keywords are not included in the example file as they require a in depth understanding of pyFAT.
All these options can also be called directly from the command line when calling pyFAT. For example the number of cores can easily be adapted by calling 'pyFAT fitting.ncpu=5'. In the case of a list the option has to be bracketed in apostrophes, i.e. 'pyFAT "fitting.fitting_stages=[Create_FAT_Cube,Run_Sofia]"'.

input keywords
--------
**main_directory**:
  *str, optional, default = os.getcwd()*

  The directory from where the input catalogue paths start. This normally is set to the directory from where pyFAT is called from.

**channel_dependency**:

  *str, optional, default = independent*

  How the channels of the input cubes overlap. Possible options are independent (default), sinusoidal, hanning.

  catalogue: null

  tirific: tirific

  sofia2: sofia2

  sofia_basename: null

  sofia_dir: null
output:
  log_directory: Logs/19-10-2021
  log_file: Log.txt
  catalogue: FAT_results.txt
  new_output: true
  output_quantity: 3
  warp_output: false
  debug: false
  timing: false
fitting:
  catalogue_start_id: '-1'
  catalogue_end_id: '-1'
  fitting_stages:
  - Create_FAT_Cube
  - Run_Sofia
  - Fit_Tirific_OSC
  ring_size: 1.1
  fixed_parameters:
  - Z0
  - XPOS
  - YPOS
  - VSYS
  opt_pixel_beam: 4
  ncpu: 6
  distance: -1.0
advanced:
    start_directory: str =f'{os.getcwd()}'
    max_iterations: int=15 #The maximum number of iterations that FAT tries to calls trific bfeore it call it quits
    loops: int =10 #The number of full loops set for tirific in a  single iteration
    minimum_warp_size: float = 3. # if the number of beams across the major axis/2. is less than this size we will only fit a flat disc,set here.
    minimum_rings: int = 3  # we need at least this amount of rings (Including 0 and 1/5 beam), set here
    too_small_galaxy: float = 1. # if the number of beams across the major axis/2 is less than this we will not fit the galaxy, set here
    unreliable_size: float = 2. #If the final diameter is smaller than this the fit is considered unreliable
    unreliable_inclination: float = 10. #If the final inclination is below this the fit is considered unreliable
    shaker_iterations: int = 20
