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

  *str, optional, default = sofia2*

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

  *str, optional, default = pyFAT_results.txt*

  The output catalogue with the results per galaxy indicating a success or failure with possible failure reason.

**new_output**:

  *bool, optional, default = True *

  Create all output anew. Set this to false in case you are looking for a final catalogue and the fitting got interrupted.

**output_quantity**:

  *int, optional, default = 3*

  A numerical indicator, in the range 0 - 5, of how much output you would like to maintain for each galaxy:

    0: just organize the output and keep all (This will also happen when a fit is unsuccesful, this can be a lot of files).

    1: remove optimized files, log files and input files.

    2: remove optimized files, log files, input files, ps files and unsmoothed files.

    3: remove optimized files, log files, input files, ps files, unsmoothed files and all model fits files except the final model.

    4: keep only the def files and remove all other output.

    5: indicates a failed fit clean up.

    >6 is the same as 0.

  Residuals are created for all cases where the fits files are maintained.

**warp_output**:

  *bool, optional, default = False*

  Switch for whether you want pyFAT top produce a warp radius and tiltograms.

**debug**:

  *bool, optional, default =False*

  Switch for printing debug messages in the log. If you are posting an issue with a log on the github please run once with this turned on.

**timing**:

 *bool, optional, default = False*

 Switch for tracking fitting time, CPU usage and RAM usage. This helps a lot to keep track on which stages are taking resources.


**font_file**

 *str, optional, default =  "/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf"*

  As fonts are a nightmare in matplotlib one can set the location of their preferred font for the plotting.
  On Ubuntu the default can be obtained by installing apt-get install ttf-mscorefonts-installer. This should point to the actual font, if the file is not fond we will fall back to DejaVu Sans.


Fitting Keywords
--------
 *Specified with fitting:*

**catalogue_start_id**:

  *str, optional, default = '-1'*

  Catalogue ID of the first galaxy to be fitted. -1 Means start at the beginning. Note that this is a string so it does not ae to be numerical.

**catalogue_end_id**:

  *str, optional, default = '-1'*

  Catalogue ID the last galaxy to be fitted, if set to -1 the whole catalogue will be fitted.

**fitting_stages**:

  *List, optional, default = ['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC']*

  Possible stages are:

    Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory.

    Catalogue_Sofia: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted.

    Run_Sofia: Run Sofia on the FAT cube and  process the output

    Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet.

    Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth.

    Tirshaker: Bootstrap errors for the final model by running Tirific multiple times with scrambled input. This can take a long time

**ring_size**:

  *float, optional, default = 1.1*

  The size of the rings in number of beams. The minimum is 0.5. A negative value indicates that FAT is not allowed to vary the ring size (although values smaller than 0.5 will still be set to 0.5 and for large galaxies the outer wings will still be doubled)

**fixed_parameters**:

 *List, optional, default = ['Z0','XPOS','YPOS','VSYS']

 A list of the parameters that should stay fixed with radius, i.e. all rings fitted as a single value, in the fitting. The rotation curve (VROT) can not be fixed at the moment. If the surface brightness is fixed it is fitted with a Gaussian after every iterations.
 XPOS, YPOS, and VSYS are always fitted as singular.
**opt_pixel_beam**:

  *int, optional, default=4*

  FAT can regrid the input cubes to have lesser pixels per FWHM of the gaussian clean beam. This can be useful to speed up the fitting process (See Kamphuis et al. 2015). If you want to prevent this set a high number of pixels. FAT will never increase the amount of pixels per FWHM.

**distance**:

  *float, optional, default = -1.*

  Distance to the galaxy. Normally for batch fitting this is taken from the input catalogue. Howeever for individual galaxies this can be set through the yaml.
  -1 means that the distance is derived from the sofia extracted vsys and the hublle flow.



Advanced Keywords
--------
*Specified with advanced:*

**start_directory**:

  *str, optional, default =f'{os.getcwd()}'*

  Vary rarely FAT will change the working directory. This keyword ensures that upon exiting one is back in the directory where started.

**max_iterations**:

  *int, optional, default=15*

  The maximum number of iterations that FAT assumes before it determines a galaxy to be unfittable. call it quits

**loops**:

  *int, optional,  int=10*

  The number of big loops set in tirific in a  single iteration. For small or irregular galaxies a smaller number can be beneficial as the FAT stabelizing routines will be called more often.

**minimum_warp_size**:

  *float, optinal, default = 3.*

  The minimum number of beams in radius that is required to allow the PA and INCL to vary. If the the extend of the model is less than this the PA and INCL will e fitted as a single value, i.e a flat disk.

**minimum_rings**:

  *int, optional, default = 3*

  The minimum amount of rings in the model (including 0 and 1/5 beam). If the galaxy is smaller than this the ring size will be reduced.
  If the ring size is <0.5*FWHM and there are still not enough rings the fit is failed, i.e. models with diameter < 1.4* FWHM wil fail.

**too_small_galaxy**:

  *float, optional, default = 1.*

  If the radius of the becomes less than this in FWHM we will not continue to fit the galaxy.

**unreliable_size**:

  *float, optional, default = 2.*

  If the final diameter of the model is smaller than this number x FWHM the fit will be flagged as unreliable.

**unreliable_inclination**:

  *float, optional, default = 10.*

  If the final inclinaion of the model is smaller than this the fit will be flagged as unreliable.

**shaker_iterations**:

  *int, optional, default = 20*

  If the Tirshaker model is set this keyword controls the amount of iterations.

**multiprocessing**:

  *bool, optional, default = True*

  Use multiprocessing

**per_galaxy_ncpu**

  *int, optional, default = 4*

  When multiprocessing is on pyFAT will attempt to find a good balence between the number of simultaneuous processes, the total allowed number of cpus and the number of cpus available in tirific.
  when multiprocessing is off this number will be set to the global ncpu parameter.
  A high number of cores per galaxy can be beneficial for speeding up the fitting of larger galaxies however, for smaller galaxies less so.

**catalogue_split_character**

  *str, optional, default = '|'*

  The character used to split the columns in the input catalogue. If left unset pyFAT asumes the default and if it can not find all columns then it tries to was a space as a seperation character.
  If it still fails it will thow a bad catalogue error.

**pa_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the PA need to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**incl_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the inclination needs to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**sdis_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the dispersion need to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**z0_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the scale height need to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**vsys_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the systemic need to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**xpos_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the right ascension needs to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**ypos_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the declination needs to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**vrot_input_boundary**

  *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

  The boundaries that the rotation curve needs to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
  Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.

**sbr_input_boundary**

    *list, optional, default = [[0., 0.], [0., 0.], [0., 0.]]*

    The boundaries that the surface brightness profile needs to remain within. Too small boundaries can lead to FAT not finding a a succesfull model.
    Given as min,max for the three areas of the fit: the central part, the approaching side warp, the receding side warp.




Individual Keywords
 --------
*No specifier*

**ncpu**:

  *int, optional, default = number of cores -1*

  Number CPUs used for fitting. In the default mode pyFAT will distribute the input across several calls to the main code.
  to set the number of cores used by tirific you can use the per_galaxy_ncpu parameter in the advanced section.

**print_examples**:

  *bool, optional, default = False*

  Print an example input yaml file and an example catalogue.

**installation_check**:

  *bool, optional , default=False*

  Run the installation check

**cube_name**:

  *str, optional, default = None*

  individual cube name when not batch fitting

**configuration_file**:

  *str, optional, default = None*

  configuration input file
