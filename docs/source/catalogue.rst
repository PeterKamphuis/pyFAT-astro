Setting your fit catalogue.
=================================

Introduction
--------

pyFAT can use 3 different input catalogues. If you have already ran SoFiA2 you can feed the output catalogue directly to pyFAT and you can also run with a single cube.
 The last option is the classical FAT mode where you run with a user provided catalogue. This section of the documentation explains how to run in these different modes.

A Single Cube
--------
Even though pyFAT is indented for batch fitting, many people find it an easy tool to setup a simple TiRiFiC start/base model which can the be expanded upon by running TiRiFiC and all its capabilities.
If you simply want to fit a single galaxy you can skip the hassle of setting up the catalogue input by simply providing the name HI fits file.
For example if you want to fit your M_83.fits HI cube you start pyFAT from the directory with the cube and simply type:

pyFAT cube_name=M_83.fits

If you are not happy with the default setting of rings of 1.1 FWHM you can type:

pyFAT cube_name=M_83.fits fitting.ring_size=0.5

or if you want the inclination to remain flat:

pyFAT cube_name=M_83.fits 'fitting.fixed_parameters=[INCL]'

or you can provide a yaml file with many non-default settings (see Advanced Settings<advanced.rst>)

pyFAT cube_name=M_83.fits configuration_file=input.yaml

if cube_name is set pyFAT will ignore all keywords relating to the catalogue.

A SoFiA2 catalogue
--------
If you have already ran SoFiA2 on your data cube and extracted the cubelets there is no need to have pyFAT run SoFia2 again but you can simply provide the catalogue.
To run in this mode the fitting.fitting_stages should contain the stage 'Catalogue_Sofia' and the sofia basename should be provided, e.g:

pyFAT




A user made catalogue
--------
