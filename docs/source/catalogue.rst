Setting your fit catalogue.
=================================

Introduction
------------
pyFAT can use 3 different input catalogues. If you have already ran SoFiA2 you can feed the output catalogue directly to pyFAT and you can also run with a single cube. The last option is the classical FAT mode where you run with a user provided catalogue. This section of the documentation explains how to run in these different modes.

A Single Cube
-------------
Even though pyFAT is indented for batch fitting, many people find it an easy tool to setup a simple TiRiFiC start/base model which can the be expanded upon by running TiRiFiC and all its capabilities.
If you simply want to fit a single galaxy you can skip the hassle of setting up the catalogue input by simply providing the name HI fits file.
For example if you want to fit your M_83.fits HI cube you start pyFAT from the directory with the cube and simply type:

  pyFAT cube_name=M_83.fits

If you are not happy with the default setting of rings of 1.1 FWHM you can type:

  pyFAT cube_name=M_83.fits fitting.ring_size=0.5

or if you want the inclination to remain flat:

  pyFAT cube_name=M_83.fits 'fitting.fixed_parameters=[INCL]'

or you can provide a yaml file with many non-default settings (see Advanced Settings <advanced.rst> )

  pyFAT cube_name=M_83.fits configuration_file=input.yaml

if cube_name is set pyFAT will ignore all keywords relating to the catalogue.

A SoFiA2 Catalogue
------------------
If you have already ran SoFiA2 on your data cube and extracted the cubelets there is no need to have pyFAT run SoFia2 again, you can simply provide the catalogue.
To run in this mode the fitting.fitting_stages should contain the stage 'Sofia_Catalogue' and the sofia basename (in later sofia2 versions this is specified as output.filename) should be provided, e.g:

  pyFAT input.catalogue="my_sofia_run_cat.txt" fitting.fitting_stages="[Sofia_Catalogue,Fit_Tirific_OSC]" input.sofia_basename='my_sofia_run'

If the sofia_dir is left unset the sofia products are assumed to be in the main directory in a directory called basename_cubelets where basename is replaced with the sofia_basename input. If one sets the sofia_dir this should refer to the cubelets directory.


A User Made Catalogue
---------------------
The classical way to run FAT is with a user constructed catalogue which has the full sample of galaxies that you want to fit. You can either use SoFiA preprocessed data or run SoFiA from within FAT. A catalogue that does not contain SoFiA data should be order in 4 columns separated by | and start with the line,

id|distance|directoryname|cubename

Where id is a string identifier for the galaxy being fit, distance is the distance to the galaxy in Mpc, directoryname is the directory where the galaxy cube with the name specified in cubename is located. The order of the columns does not matter for pyFAT as long as the are properly specified in the first line.
pyFAT always tries to fit the brightest object SoFiA finds in the cube unless this object is very close to the edges. Therefore the best way to ensure that the right object is fitted is to make sure that it is the only bright object in the data cube. Alternatively, one can run sofia manually to make sure it is the only object detected and then provide pyFAT with this input.
In that case the input columns are:

id|distance|directoryname|cubename|basename

Where the last column corresponds to the basename used in the sofia run. pyFAT assumes the sofia products are located in the main fitting directory specified under directoryname. If this is not the case, e.g. because all sofia output is in a single directory, the sofia directory can be specified with the configuration keyeword output.sofia_dir.
pyFAT never touches original input files and thus the sofia input will be copied into Sofia_Output directory in the fitting directory. When providing the SoFiA input the keyword fitting.fitting_stages should include Existing_Sofia instead of Run_Sofia.
