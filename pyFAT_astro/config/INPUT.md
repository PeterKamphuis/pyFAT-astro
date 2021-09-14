#The Configuration File
=====

Introduction
----

pyFAT uses OmegaConf (https://github.com/omry/omegaconf) to handle the input settings. pyFAT can be ran with default settings or settings from a yml configuration file or with command line input. For batch fitting one needs to set an input catalogue with input for each galaxy.

In a run pyFAT first checks the defaults, then the configuration yaml and finally the command line input. This mean that if a value is set in all three input methods the one from the command line is used.

The yaml file has five main sections globals, input, output, fitting and advanced that contain several parameters that can adjust how pyFAT runs.

The globals that should only be used from the command line.
----

PyFAT uses four globals that are not under any section and can be called directly from the command line. It is not recommended to add these to the yaml file as they alter how the yaml file is dealt with.

These globals and their defaults are:

  print_examples: False
  installation_check: False
  cube_name: None
  configuration_file: None

A configuration file with all default values and an example catalogue file can be printed by

  pyFAT print_examples=True

and the installation check can be ran similiarly

  pyFAT print_installation=True

The cube name can be directly specified from the command line. PyFAT will then simply fit this cube and nothing else. The configuration file for fitting can then still be provided but any catologue information specified in there will be ignored. This can be useful to fit a single galaxy. Or you want to refit a single galaxy from a batch with slightly different parameters.

For most fitting the defaults should suffice however pyFAT is now much more flexible then the previous GDL version and hence many parameters can be set. This is easiest done through a yaml which can be provided to pyFAT by:

  pyFAT configuration_file=FAT_Input.yml

The different sections in the yaml file are explained below.

The Input Section
----
The main input parameters can be set by input.parameter on the command line or in the input: section in the configuration file.

input accepts the following input parameters with their default indicated:

  main_directory: directory from where pyFAT is ran.
This is directory where the cube to be fitted or the directories from the catalogue are residing.

  channel_dependency: independent
This indicates the instrumental velocity resolution which pyFAT attempt to separate from the intrinsic dispersion. The options are independent,sinusoidal, and hanning. Most modern correlators have practically independent channels and thus this is the default. Sinusoidal assumes a a dependency between the channels based on a sinus overlap where the sine wave spans 1.2 x the channelwidth. Hanning assumes that the input is hanning smoothed. Even though smoothing the cube should in principle not affect the TRM fitting there are indications that it can speed up the default fitting procedure of pyFAT.

  catalogue: None
The input catalogue for batch fitting. See the Input Catalogue section for further details.

  tirific: tirific
  sofia2: sofia2/sofia
The names of the commands on the system to run TiRiFiC and SoFiA. If you have installed these programs under a different name you can specify them here. In case of SoFiA pyFAt first looks for sofia2 and if that is not found it tries sofia.

  sofia_basename: None
  sofia_dir: None
In case you want to start the fitting from a pre-made sofia run these parameters indicate the base name that was used in the sofia run. These are output.filename
and output.directory in the sofia2 .par file.

The Output Section
----
The main output parameters can be set by output.parameter on the command line or in the output: section in the configuration file.

output accepts the following input parameters with their default indicated:

  log_directory: Logs/Current Date
This specifies the directory where the log is written. When the debug option is set several other tracking files are written here as well. This directory is also used for the timing output and any Overview files and Log files from previous runs are moved here.  

  log_file: Log.txt
Name of the log file for the run.

  catalogue: 'FAT_results.txt'
Name of the catalogue for batch fitting. See the section Catalogue Input file for how to set up such a file. If you want to batch fit with all defaults you can run pyFAT with only the catalogue set from the command line. It is possible to simple set a SoFiA2 catalogue such that. pyFAT will then automatically organize the SoFiA output into a pyFAT catalogue to fitted where the source finding stage is omitted. This requires output.writeMask = true and output.writeCubelets = true in the SoFiA2 run. To trigger this use the stage Catalogue_Sofia in fitting.fitting_stages (See Below)

  new_output: True
If set to true the output catalogue will be reinitialized and the old catalogue will be moved to a prev.txt file. If False the code will simply append the results to the output catalogue.

  output_quantity: 3
The amount of output you would like to retain for each galaxy. With:
0 = The output will merely be organzied but nothing is deleted. This typically not necessary and can result in a lot files. If a fit fails this will automatically happen such that it easier to trace where the fit has gone wrong.
1 = Remove optimized files, log files and input files.
2 = Remove optimized files, log files, input files, ps files and unsmoothed files.
3 = Remove optimized files, log files, input files, ps files, unsmoothed files and all model fits files except the final model.
4 = keep only the def files and remove all other output.
5 = indicates a failed fit clean up.
>6 is the same as 0.

Residuals are normally created for all cases where the fits files are maintained.   

  warp_output: False 
If you want pyFAT to output a warp radius and tiltograms this should be set to True.

  debug: False
Additonal output is produced for debugging. This typically makes the log file impossible to read except for developers. If you get a mystical message why your fit failed it is useful to run with this option set to true and post the log to the github in order to get a quick turn aoround on a answer.

  timing: False
For optimization pyFAT can track the time it takes to fit the galaxies by setting this option to True. The ouput will be produced into the log file.

Fitting settings
----
The main fitting parameters can be set by fitting.parameter on the command line or in the fitting: section in the configuration file.

fitting accepts the following input parameters with their default indicated:

  catalogue_start_id: -1
The catalogue ID of the first galaxy to be fitted. -1 Means start at the beginning
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



In the fitting section of the configuration file one should set the stages that one want to fit.

  fitting_stages: ['']

The possible stages are

  Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory

  Sofia_Catalogue: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted

  Run_Sofia: Run Sofia on the FAT cube and  process the output

  Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet

  Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth

Note that if pre created sofia output is used several optimization routines are skipped as the sofia output is assumed to be reasonable.
If you want to use these routines it is adviced to first run create_FAT_cube then run Sofia and then feed the output to pyFAT.
