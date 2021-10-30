The Different Stages.
=================================

Introduction
--------

One major differences between FAT and pyFAT is that when coding up pyFAT an attempt was made to keep major different stages, such as source finding or fitting, as different python modules such that one can easily different combinations or add their own module to replace the existing one.
This are currently 6 different stages which can be applied which together form the pyFAT pipeline in sequence. These stage are explained here. The stages are set through the fitting.fit_stages keyword.
The default order is 'Create_FAT_Cube','Run_Sofia','Fit_Tirific_OSC'.

Create_FAT_Cube
--------
Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory. Running this stage will ensure that the pipeline will not trip up over header keywords deep into the fitting process.


Run Sofia
--------
Run Sofia on the FAT cube and  process the output. The brightest source is selected as the target.

Sofia_Catalogue
--------
The input catalogue is an sofia catalogue. In this case all objects in the catalogue are deemed separate targets and pyFAT will setup the require directory structure itself. When using this stage the keyword sofia_basename is required.
If used in combination with Create_FAT_Cube  all cubelets are checked, however it is more efficient to first run pyFAT with only the Create_FAT_Cube stage on the large cube and then run sofia before running pyFAT with fitting modules.
Existing_Sofia
--------
It is assumed the Sofia output exist and specified in the fitting catalogue, this means a sofia catalogue exists for every cubelet.
Note that if pre created sofia output is used several optimization routines are skipped as the sofia output is assumed to be reasonable.

Fit_Tirific_OSC
--------
Run FAT using the Tirific program in multiple iterations and smooth.

Tirshaker
--------
Bootstrap errors for the final model by running Tirific multiple times with scrambled input. This can take a long time
