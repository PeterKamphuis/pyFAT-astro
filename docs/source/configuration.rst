The Configuration File
=======================

Introduction
------------

pyFAT uses OmegaConf (https://github.com/omry/omegaconf) to handle the input settings. pyFAT can be ran with default settings or settings from a yml configuration file or with command line input. For batch fitting one needs to set an input catalogue with input for each galaxy.

In a run pyFAT first checks the defaults, then the configuration yaml and finally the command line input. This mean that if a value is set in all three input methods the one from the command line is used.

The yml file has five main sections globals, input, output, fitting and advanced that contain several parameters that can adjust how pyFAT runs. All parameters are described in in the section Advanced Settings <advanced.rst> one by one. Here we give a more general overview of setting up a yml configuration file.

Individual Keywords
-------------------

PyFAT uses four individual keywords  that are not under any section and can be called directly from the command line. It is not recommended to add these to the yml file as they alter how the yml file is dealt with.

These globals and their defaults are:

  print_examples: False

  installation_check: False

  cube_name: None

  configuration_file: None

A configuration file with all default values and an example catalogue file can be printed by

  pyFAT print_examples=True

and the installation check can be ran similarly

  pyFAT installation_check=True

The cube name can be directly specified from the command line. PyFAT will then simply fit this cube and nothing else. The configuration file for fitting can then still be provided but any catalogue information specified in there will be ignored. This can be useful to fit a single galaxy or when you want to refit a single galaxy from a batch with slightly different parameters.

For most fitting the defaults should suffice however pyFAT is now much more flexible then the previous GDL version and hence many parameters can be set. This is easiest done through a yml file which can be provided to pyFAT by:

  pyFAT configuration_file=FAT_Input.yml

The different sections in the yml file are explained in general terms below and all individual keywords can be found in Advanced Settings <advanced.rst>.

The Input Section
-----------------

The main input parameters can be set by input.parameter on the command line or in the input: section in the configuration file. These keywords relate to input file locations such as the input catalogue or specifics about the input cubes.


The Output Section
------------------

The main output parameters can be set by output.parameter on the command line or in the output: section in the configuration file. These keywords relate to names and locations of the output files such as the log name and directory and how much of the pyFAT should be saved after the run ends.
Additionally one can specify whether new output should be appended to old files or start with new files. Finally in this section one can trigger pyFAT to produce debug and timing results for developers.

Fitting settings
----------------

The main fitting parameters can be set by fitting.parameter on the command line or in the fitting: section in the configuration file. These include start and end id's in the catalogue, the number of cpu's to use and other such settings.
It is in this section that the different stages of pyFAT (See The Different Stages <stage.rst>) are set as well throught the fitting_stages keyword. This keyword takes a list of stages to be executed. When specified from the command line the list should be specified between apostrophe's, e.g.:

  pyFAT fitting.fitting_stages='['Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC']'
