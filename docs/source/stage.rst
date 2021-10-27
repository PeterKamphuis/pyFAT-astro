Possible stages are:

  Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory.

  Catalogue_Sofia: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted.

  Run_Sofia: Run Sofia on the FAT cube and  process the output

  Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet.

  Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth.

  Tirshaker: Bootstrap errors for the final model by running Tirific multiple times with scrambled input. This can take a long time


  The possible stages are

    Create_FAT_Cube: Create a FAT compatible cube from the original cube. This will overwrite any previous FAT cube present. If omitted it is assumed the Cube is present in the fitting directory

    Sofia_Catalogue: The catalogue is assumed to be a sofia catalogue and all sources in the catalogue are to be fitted

    Run_Sofia: Run Sofia on the FAT cube and  process the output

    Existing_Sofia: It is assumed the Sofia output exist and specified in the fitting catalogue, this means a catalogue exists for every cubelet

    Fit_Tirific_OSC: Run FAT using the Tirific program in multiple iterations and smooth

  Note that if pre created sofia output is used several optimization routines are skipped as the sofia output is assumed to be reasonable.
  If you want to use these routines it is adviced to first run create_FAT_cube then run Sofia and then feed the output to pyFAT.
