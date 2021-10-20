Setting your fit preferences through a yaml file.
=================================

Introduction
--------

In comparision to FAT pyFAT delivers much more control over the fitting procedure. This is made possible by the new modular setup and the use of omegaconf for the input and default settings.
In this page we explain how the fitting in FAT can be adapted by providing a yaml input parameters for pyFAT. An example file with all the settings can be printed by runnning 'pyFAT print_examples=True' after installing pyFAT-astro. This command prints both an example yaml file and an example catalogue input file with the names FAT_defaults.yml and FAT_Input_Catalogue.txt, respectively
Below we explain what each section in the yaml file does.

The is divided  in four different sections and has 4 independent keywords. All keywords are optional but pyFAT does require at least a cube name or an input catalogure to work.  The section advanced and the independent keywords are not included in the example file as they require a in depth understanding of pyFAT.
All these options can also be called directly from the command line when calling pyFAT. For example the number of cores can easily be adapted by calling 'pyFAT fitting.ncpu=5'. In the case of a list the option has to be bracketed in apostrophes, i.e. 'pyFAT "fitting.fitting_stages=[Create_FAT_Cube,Run_Sofia]"'.

Input Keywords
--------
*Specified with input:*

**main_directory**:
  *str, optional, default = os.getcwd()*

  The directory from where the input catalogue paths start. This normally is set to the directory from where pyFAT is called from.

**channel_dependency**:

  *str, optional, default = independent*

  How the channels of the input cubes overlap. Possible options are independent, sinusoidal, hanning.

**catalogue**:

  *str, optional, no default*

  The input catalogue for batch fitting.

**tirific**:

  *str, optional, default = tirific*

  Command used to call tirific from the python subprocess

**sofia2**:

  *str, optional, default = sofia2

  Command to call sofia 2 from the python subprocess

**sofia_basename**:

  *str, optional, no default*

  The basename used in creating the sofia products. This will only be used if the input catalogue is a sofia catalogue or pre-processed sofia products are sepecified.

**sofia_dir**:

  *str, optional, no default*

  Directory of the existing Sofia output. This will only be used if the input catalogue is a sofia catalogue or pre-processed sofia products are sepecified.

Output Keywords
--------
*Specified with output:*

**log_directory**:

  *str, optional, default = Logs/{datetime.now().strftime("%d-%m-%Y")*

  Name of the directory where log files will be stored.

**log_file**:

  *str, optional, default = Log.txt*

  Name of the log file where all print messages are produced.

**catalogue**:

  *str, optional, default = pyFAT_results.txt

  The output catalogue with the results per galaxy indicating a success or failure with possible failure reason.

**new_output**:

  *bool, optional, default = True *

  Create all output anew. Set this to false in case you are looking for a final catalogue and the fitting got interrupted.

**output_quantity**:

  *int, optional, default = 3*

  A numerical indicator, in the range 0 - 5, of how much output you would like to maintain for each galaxy:

    - 0: just organize the output and keep all (This will also happen when a fit is unsuccesful, this can be a lot of files).

    - 1: remove optimized files, log files and input files.

    - 2: remove optimized files, log files, input files, ps files and unsmoothed files.

    - 3: remove optimized files, log files, input files, ps files, unsmoothed files and all model fits files except the final model.

    - 4: keep only the def files and remove all other output.

    - 5: indicates a failed fit clean up.

    - >6 is the same as 0.

    Residuals are created for all cases where the fits files are maintained.

**warp_output**:

  *bool, optional, default = False*

  Switch for whether you want pyFAT top produce a warp radius and tiltograms.

**debug**:

  *bool, optional, default =False*

  Switch for printing debug messages in the log. If you are posting an issue with a log on the github please run once with this turned on.

**timing**:

 *bool, optional, default = False

 Switch for tracking fitting time, CPU usage and RAM usage. This helps a lot to keep track on which stages are taking resources.


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
