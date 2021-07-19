# Python Fully Automated TiRiFiC (pyFAT)
=====

Introduction
------------

The Python Fully Automated TiRiFiC is a python (> 3.6) wrapper around the tilted ring fitting code ([TiRiFiC](http://gigjozsa.github.io/tirific/)) that aims to fully automate the process of fitting simple tilted ring models to line emission cubes. This python version is in continuous development and hence errors and bugs can be present. Nevertheless, the code has extensively been tested and the results and a more extensive description of the code are documented in [Kamphuis et al. 2015](http://arxiv.org/abs/1507.00413) and [Kamphuis et al. in prep]

The file [Current_Status.pdf](./Current_Status.pdf)  provides an overview of the default test sets provided by the [HI TRM Database](https://github.com/PeterKamphuis/HI_TRM_Test_Database) which illustrate the performance of the current release version of FAT. This fully tested version is available on the master branch. In release 2.0 the code has permanently switched to python development in order to maintain easier use and a larger development base. The final IDL/GDL version is release 1.5. The python version strives to be an improvement on the IDL/GDL version and has been extensively tested to at the very least match the results of the older code. However, should you find that your fits have significantly degraded please open an issue on the GitHub and let us know. From release 2.0 the intention is to maintain a more regular release policy when the code has undergone major development. For the very latest version of FAT one can always check the available branches but these come without any quality guarantees.

If you are looking for specific functionality or find that FAT is not performing well despite the galaxy having regular rotation or just want to chat about tilted ring modelling pipelines please do not hesitate to contact me.

It is important to remember that FAT is meant for batch fitting. Hence, the aim of the code is to provide tilted ring models that are accurate for a large fraction of galaxies. Ideally, FAT should identify galaxies that are not fitted well however this feature is not optimal yet. When fitting individual galaxies it is recommended to run FAT and then fine tune the model by hand in [TiRiFiC](http://gigjozsa.github.io/tirific/). In most cases such fine tuning will be limited to a few outer rings but in the case of complex galaxies with significant non-cylindrically symmetric motions the models can fail (Or in the case of very bad data but that is not a FAT issue).

FAT is not an automated version of the extended functionality of [TiRiFiC](http://gigjozsa.github.io/tirific/). FAT fits simple rotationally symmetric discs with asymmetric warps and surface brightness distributions. However, [TiRiFiC](http://gigjozsa.github.io/tirific/) itself provides a much more extended functionality and should be used for identifying significant non-cylindrically symmetric motions, thick discs, bars and other such HI features. When modelling such galaxies ideally FAT can provide a base model and setup a .def file with merely a thin disc. These can then be used in [TiRiFiC](http://gigjozsa.github.io/tirific/) in order to explore large scale motions not captured by FAT’s simple model. For examples of such modelling please see [Kamphuis et al. (2011)](http://adsabs.harvard.edu/abs/2011MNRAS.414.3444K), [Zschaechner et al. (2011)](http://adsabs.harvard.edu/abs/2011ApJ...740...35Z), [Kamphuis et al. (2013)](http://adsabs.harvard.edu/abs/2013MNRAS.434.2069K), [Gentile et al. (2013)](http://adsabs.harvard.edu/abs/2013A%26A...554A.125G).

Requirements
------------
The code requires full installation of:

    python v3.6 or higher
    numpy>=1.14, scipy, astropy, omegaconf, matplotlib, future-fstrings, importlib_resources>=3.3.0 (These should be managed through a pip install)
    TiRiFiC v2.2.3 or higher
    SoFiA2  

[python](https://www.python.org/),[TiRiFiC](http://gigjozsa.github.io/tirific/download_and_installation.html), [SoFiA2](https://github.com/SoFiA-Admin/SoFiA-2)

TiRiFiC and SoFiA2 should be accessible for subproccess calls. This normally means that it should be possible to invoke them properly from the command line.
TiRiFiC should be installed as tirific (standard) and can be installed through the kern-suite (https://kernsuite.info)
SoFiA2 should be installed as sofia2 or sofia, pyFAT will first check for the existence of the sofia2 command, if this does not exist it will look for the sofia command and assume it is SoFiA2 if it can run.

Installation
------------

Download the source code from the Github. Once downloaded pyFAT can be installed with:

  	pip install pyfat-astro

This should also install all required python dependencies.
We recommend the use of python virtual environments. If so desired a pyFAT installation would look like:

  	python3 -m venv FAT_venv
  	source FAT_venv/bin/activate
  	pip install pyFAT
  	deactivate
  	source FAT_venv/bin/activate

Once you have installed FAT you can check that it has been installed properly by running FAT as.

  	FAT>  pyFAT installation_check=True

This should take typically 10 min and should finish with the message:

	!!!!--------------------------------------------!!!!!
	!!!! As far as we can tell FAT is installed     !!!!!
	!!!! properly and runs smoothly.                !!!!!
	!!!!--------------------------------------------!!!!!

The check consists of fitting a flat disk on NGC 2903. The data for this galaxy were take as part of the WHISP program.
This survey is decribed in [van der Hulst et al. (2001)](http://adsabs.harvard.edu/abs/2001ASPC..240..451V) and the data can be found at [Westerbork on the Web](http://wow.astron.nl/) or the [WHISP page](https://www.astro.rug.nl/~whisp/).

If you get any other message please do not hesitate to file an issue on the GitHub (https://github.com/PeterKamphuis/pyFAT). Do note however that if you perform this check on a unreleased version/branch it might not perform well. So always check with the master branch.

The Overview.png will contain a comparison with the fit performed by you. These should be the same (the correct fit is classified as Input.)

The plots should look like this:

![Overview plot after running installation check.](Installation_Check/Overview.png)

Sometimes, due to updates in SoFiA2 or TiRiFiC, the check might show differences beyond the tolerance limits. If these are small and you have checked the individual installations of SoFiA2, TiRiFiC and the Installation Check files are older than the latest SoFiA2 or TiRiFiC update, then the installation is probably correct and you can continue. Please do post an issue about the outdated installation check.



Running FAT
-----------
FAT is currently run under python and can be run from the command line

    FAT> pyFAT -h

Will provide an overview of call options. For the most basic usage one can call FAT with a configuration file.

    FAT> pyFAT configuration_file=FAT_Input.yml

or a single Cube

    FAT> pyFAT cube_name=Input_Cube.fits

Where Input_Cube.fits is the observation to be fitted. In this mode configuration_file can still be used to specify fit settings but catalogue and location setting will be ignored. !! If cube_name is set in either the command line or the configuration file this always will always trigger the singular fitting instead of batch fitting.

FAT is intended for batch fitting and as such it is recommended to have all source in separate directories

Configuration File
------

pyFAT uses OmegaConf (https://github.com/omry/omegaconf) to handle the input settings. pyFAT can be ran with default settings or settings from a yml configuration file or with command line input. For batch fitting one needs to set an input catalogue with input for each galaxy.

A configuration file with all default values and an example catalogue file can be printed by

  pyFAT print_examples=True

In a run pyFAT first checks the defaults, then configuration yaml and finally the command line input. This mean that if a value is set in all three input methods the one from the command line is used.

The yaml file has three main sections input, output, fitting that contain several parameters that can adjust how pyFAT runs the specifics of each parameter is explained in read the docs.

In the input section one need to set the catalogue that contains all the galaxies that need to be fitted.

  catalogue = Path_to_catalog_dir/Your_Catalog.txt

How to arrange you catalogue is specified below. If you have already ran sofia 2 the resultsing output catalogue can be fed to pyFAT directly if you have produced cubelets and masks. If so pyFAT will automatically create a directory output_FAT_cubelets in the main directory and create directories for all sources in the sofia catalogue. To use this stage it is crucial that in the fitting section you add 'Sofia_Catalogue' to the fitting stages.

If you are using hanning smoothed cube you should also set

  hanning_smoothing=True

in the input section.

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

Input Catalog
-----------

The input catalog should have at least 4 columns named as

        number|distance|directoryname|cubename

and seperated by |
The number is an easy identifier to keep track of which galaxy is being fitted.
the distance is the distance to the galaxy in Mpc. This is used to make some initial guesses for the structure of the galaxy. If it is unknown it should be set to 1.
The directory name is the name of the directory of the galaxy to be fitted. This directory should be located in the specified maindir in the config file.
cubename is the name of the cube to be fitted. This should be without the fits extension.

An example catalog is included in the distribution. This also gives examples for how to set up a catalog when using pre-made sofia input.
