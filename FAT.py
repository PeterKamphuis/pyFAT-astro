#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
# This is the python version of FAT
__version__ = 'pyFAT 1.0.0'
import sys
import os
import copy
import numpy as np
from optparse import OptionParser
import traceback
from datetime import datetime
from astropy.io import fits
from astropy.wcs import WCS


def FAT(argv):
    #Get the directory we are running from
    start_dir = os.getcwd()

    #Constants that are used in the code
    H_0 = 70. #km/s/Mpc #Hubble constant

    #Then check the input options
    parser  = OptionParser()
    parser.add_option('-c','--cf','--configuration_file', action ="store" ,dest = "configfile", default = 'FAT_INPUT.config', help = 'Define the input configuration file.',metavar='CONFIGURATION_FILE')
    parser.add_option('-d','--debug', action ="store_true" ,dest = "debug", default = False, help = 'Print debug messages',metavar = '')
    parser.add_option('-s','--sd','--support_directory', action ="store" ,dest = "supportdir", default=f'{start_dir}/Support', help = 'location where the support files reside.',metavar='SUPPORT_DIR')
    parser.add_option('-i','--ic','--installation_check', action ="store_true" ,dest = "installation_check", default = False, help = 'Run the installation _check.',metavar = '')
    parser.add_option('--LVT','--LVHIS_TEST', action ="store_true" ,dest = "lvhis_test", default = False, help = 'Run the LVHIS Test. Developer Only.')
    parser.add_option('--PT','--PAPER_TEST', action ="store_true" ,dest = "paper_test", default = False, help = 'Run the PAPER Test. Developer Only.')
    parser.add_option('--FD','--FULL_DATABASE', action ="store_true" ,dest = "full_test", default = False, help = 'Run the Full Database Test. Developer Only.')
    parser.add_option('-p','--problems', action ="store_true" ,dest = "problems", default = False, help = 'Run the Problem test set. Developer Only.')
    parser.add_option('-t','--timing', action ="store_true" ,dest = "timing", default = False, help = 'Create a file in the maindir that provides start and stop times for each galaxy.')
    input_parameters,args = parser.parse_args()

    basic_info  = 'BasicInfo'
    if input_parameters.installation_check:
        input_parameters.configfile = 'Installation_Check/FAT_INPUT.config'
    if input_parameters.lvhis_test:
        fat_main_test_dir = os.environ["FAT_TEST_DIR"]
        input_parameters.configfile=fat_main_test_dir+'/LVHIS-26_3/Input.config'
    if input_parameters.paper_test:
        fat_main_test_dir = os.environ["FAT_TEST_DIR"]
        input_parameters.configfile=fat_main_test_dir+'/SmallCat_Warps/Input_1.config'
    if input_parameters.full_test:
        fat_main_test_dir = os.environ["FAT_TEST_DIR"]
        input_parameters.configfile=fat_main_test_dir+'/Full_Database/FAT_INPUT.config'
    if input_parameters.problems:
        fat_main_test_dir = os.environ["FAT_TEST_DIR"]
        input_parameters.configfile=fat_main_test_dir+'/Problems/FAT_INPUT.config'

    #Add the support dir to the system path and read the functions from there
    sys.path.insert(1, input_parameters.supportdir)
    # Functions that read files
    import read_functions as rf
    # functions that are used often for menial tasks
    import support_functions as sf
    # function that keep things orderly and nicely
    import clean_functions as cf
    # Functions that modify or produce fat fits file
    import fits_functions as ff
    # Functions that run external programs such as tirific and sofia
    import run_functions as runf
    # Functions that check on the presence of output or input
    import exist_check as ec
    #Check the existence of the config file and read it
    read_succes = False
    while not read_succes:
        try:
            Original_Configuration = rf.config_file(input_parameters,start_dir)
            read_succes = True
        except:
            print(traceback.print_exc())
            input_parameters.configfile = input('''
                    You have provided a config file but it can't be found.
                    If you want to provide a config file please give the correct name.
                    Else press CTRL-C to abort.
            ''')
    # Add the starting directory to the Configuration
    Original_Configuration['START_DIR'] = start_dir
    # Also add the timing input and some other recurring parameters
    Original_Configuration['TIMING'] = input_parameters.timing
    Original_Configuration['DEBUG'] = input_parameters.debug
    Original_Configuration['FINAL_COMMENT'] = "This fitting stopped with an unregistered exit."
    Original_Configuration['CONVERGENCE_1'] = 'Not converged'
    Original_Configuration['AC1'] = 0
    Original_Configuration['AC2'] = 0
    #Then read the input Catalogue
    Catalogue = rf.catalogue(Original_Configuration['CATALOGUE'])
    # Get the longest directory name to format the output directory properly
    dirname = 'Directory Name'
    maximum_directory_length = len(dirname)
    for directory in Catalogue['DIRECTORYNAME']:
        if directory == './':
            maximum_directory_length = len(Original_Configuration['MAINDIR'].split('/')[-2])
        if len(directory) > maximum_directory_length:
            maximum_directory_length = len(directory)

    # Create a file to write the results to if if required
    if not os.path.exists(Original_Configuration['OUTPUTCATALOGUE']) or Original_Configuration['NEW_OUTPUT']:
        output_catalogue = open(Original_Configuration['OUTPUTCATALOGUE'],'w')
        comment = 'Comments on Fit Result'
        AC1 = 'AC1'
        AC2 = 'AC2'
        output_catalogue.write(f"{dirname:<{maximum_directory_length}s} {AC1:>6s} {AC2:>6s}   {comment}\n")
        output_catalogue.close()

    if Original_Configuration['TIMING']:
        timing_result = open(Original_Configuration['MAINDIR']+'/Timing_Result.txt','w')
        timing_result.write("This file contains the system start and end times for the fitting of each galaxy")
        timing_result.close()
    #  Negative start galaxy is not allowed
    if Original_Configuration['STARTGALAXY'] < 0:
        Original_Configuration['STARTGALAXY'] = 0
    # If the end galaxy is -1 fit the whole catalogue
    if Original_Configuration['ENDGALAXY'] == -1:
        Original_Configuration['ENDGALAXY'] = len(Catalogue['NUMBER'])
    # start the main fitting loop
    for i in range(Original_Configuration['STARTGALAXY'],Original_Configuration['ENDGALAXY']):
        # First check the starttime
        Configuration = copy.deepcopy(Original_Configuration)
        Configuration['START_TIME'] = datetime.now()
        doubled = 0
        ini_mode_factor =25
        #Make a dictionary for the fitsfiles we use
        Fits_Files = {'ORIGINAL_CUBE': f"{Catalogue['CUBENAME'][i]}.fits"}

        # If we have a fitting log we start writing
        log_statement = f'''{sf.linenumber():7s} This file is a log of the fitting process run at {Configuration ['START_TIME']}.
{"":7s} This is version {__version__} of the program.
'''
        # Adapt configuration to hold some specifics to this galaxy
        if Configuration['OUTPUTLOG']:
            if Catalogue['DIRECTORYNAME'][i] == './':
                Configuration['OUTPUTLOG'] = f"{Configuration['MAINDIR']}/{Configuration['OUTPUTLOG']}"
            else:
                Configuration['OUTPUTLOG'] = f"{Configuration['MAINDIR']}/{Catalogue['DIRECTORYNAME'][i]}/{Configuration['OUTPUTLOG']}"
            #If it exists move the previous Log
            if os.path.exists(Configuration['OUTPUTLOG']):
                os.rename(Configuration['OUTPUTLOG'],f"{Configuration['MAINDIR']}/{Catalogue['DIRECTORYNAME'][i]}/Prev_Log.txt")
        sf.print_log(log_statement,Configuration['OUTPUTLOG'])

        # Adapt configuration to hold some specifics to this galaxy
        #Never use the original cube only a fat modified one
        if 'BASENAME' in Catalogue['ENTRIES']:
            Configuration['BASE_NAME'] = Catalogue['BASENAME'][i]
        else:
            Configuration['BASE_NAME'] = Catalogue['CUBENAME'][i]+'_FAT'
        Fits_Files['NOISEMAP'] = f"{Configuration['BASE_NAME']}_noisemap.fits"
        Fits_Files['FITTING_CUBE'] = f"{Catalogue['CUBENAME'][i]}_FAT.fits"
        Fits_Files['MOMENT0'] = f"{Configuration['BASE_NAME']}_mom0.fits"
        Fits_Files['MOMENT1'] = f"{Configuration['BASE_NAME']}_mom1.fits"
        Fits_Files['MOMENT2'] = f"{Configuration['BASE_NAME']}_mom2.fits"
        Fits_Files['MASK'] = f"{Configuration['BASE_NAME']}_mask.fits"


        #Add our fitting directory to the Configuration
        if Catalogue['DIRECTORYNAME'][i] == './':
            Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/"
        else:
            Configuration['FITTING_DIR'] = f"{Configuration['MAINDIR']}/{Catalogue['DIRECTORYNAME'][i]}/"
        if Configuration['FITTING_DIR'][-2:] == '//':
            Configuration['FITTING_DIR'] = Configuration['FITTING_DIR'][:-2]+'/'
        # then we want to read the template
        Tirific_Template = rf.tirific_template(f'{input_parameters.supportdir}/template.def')

        log_statement = f'''We are in loop {i}. This is catalogue number {Catalogue['NUMBER'][i]} and the directory {Catalogue['DIRECTORYNAME'][i]}.'''
        sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)
        #Run cleanup
        cf.cleanup(Configuration)

        # Check if the input cube exists
        if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['ORIGINAL_CUBE']}"):
            log_statement = f'''We cannot find the cube {Fits_Files['ORIGINAL_CUBE']} in the directory {Configuration['FITTING_DIR']}.
         We skip this galaxy.
'''
            sf.print_log(log_statement,Configuration['OUTPUTLOG'], screen =True)
            Configuration['FINAL_COMMENT'] = "This galaxy has no fits cube to work with, it is skipped."
            cf.finish_galaxy(Configuration,maximum_directory_length)
            traceback.print_exc()
            continue

        #Define a dictionary with a set of tracking triggers
        Tracking = {'MAX_RINGS' : 0,
                    'OPTIMIZED' : False,
                    'TOO_SMALL' : False
                    }

        # Let's see if our base cube exists, Note that cleanup removes it if we want to start from the original dir so no need to check start_point
        if not os.path.exists(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}"):
            try:
                ff.create_fat_cube(Configuration, Fits_Files)
            except Exception as e:
                Configuration['FINAL_COMMENT'] = e
                Configuration['MAPS_OUTPUT'] = 5
                cf.finish_galaxy(Configuration,maximum_directory_length)
                traceback.print_exc()
                continue
        # We open the header of the fitting cube and get some parameters and make a header wcs structure
        cube_hdr = fits.getheader(f"{Configuration['FITTING_DIR']}/{Fits_Files['FITTING_CUBE']}")
        cube_wcs = WCS(cube_hdr)
        #If we have Sofia Preprocessed Output request make sure it all exists
        if Configuration['START_POINT'] == 3:
            try:
                ec.sofia_output(Configuration,Fits_Files)
            except Exception as e:
                Configuration['FINAL_COMMENT'] = e + 'Please make sure your Sofia input is in the right location'
                Configuration['MAPS_OUTPUT'] = 5
                cf.finish_galaxy(Configuration,maximum_directory_length)
                traceback.print_exc()
                continue

        if Configuration['START_POINT'] != 3:
            # Run sofia2
            try:
                runf.sofia(Configuration, Fits_Files,cube_hdr,input_parameters.supportdir)
                # check that all is well
                ec.sofia_output(Configuration,Fits_Files)
            except Exception as e:
                Configuration['FINAL_COMMENT'] = e
                Configuration['MAPS_OUTPUT'] = 5
                cf.finish_galaxy(Configuration,maximum_directory_length)
                traceback.print_exc()
                continue
            # We assume sofia is ran and created the proper files
        try:
            # These come out as strings
            name,x,x_min,x_max,y,y_min,y_max,z,z_min,z_max,ra,dec,v_app,f_sum,kin_pa = rf.sofia_catalogue(Configuration,header=cube_hdr)

        except Exception as e:
            Configuration['FINAL_COMMENT'] = e
            Configuration['MAPS_OUTPUT'] = 5
            cf.finish_galaxy(Configuration,maximum_directory_length)
            traceback.print_exc()
            continue
        print([x,x_min,x_max,y,y_min,y_max,z,z_min,z_max])
        x_min,x_max,y_min,y_max,z_min,z_max = sf.convert_type([x_min,x_max,y_min,y_max,z_min,z_max], type = 'int')
        x,y,z,ra,dec,v_app,f_sum,kin_pa = sf.convert_type([x,y,z,ra,dec,v_app,f_sum,kin_pa])
        #How does sofia 2 deal with the fully flagged channels?
        # Need to check for that here if NaNs are included

        # convert the boundaries to real coordinates
        ralow,declow,vellow = cube_wcs.wcs_pix2world(x_min,y_min,z_min,0.)
        rahigh,dechigh,velhigh = cube_wcs.wcs_pix2world(x_max,y_max,z_max,0.)
        DECboun = np.sort([float(declow),float(dechigh)])
        RAboun = np.sort([float(ralow),float(rahigh)])
        VELboun = np.sort([float(vellow),float(velhigh)])
        # If the provided distance  = -1 we assume a Hubble follow
        if Catalogue['DISTANCE'][i] == -1:
            Catalogue['DISTANCE'][i] == v_app/(1000.*H_0)
        if Catalogue['DISTANCE'][i] < 0.5:
            Catalogue['DISTANCE'][i] == 0.5

        print(DECboun)
        print(RAboun)
        print(VELboun)
'''





                                ;let's check the size of this cube to see whether it is reasonable
     IF maxrings EQ 0. then begin
        foundmin=0.
        minx=intarr(sxpar(headermask,'NAXIS1'))
        maxx=intarr(sxpar(headermask,'NAXIS1'))
        countmin=-1
        for pix=0,n_elements(mask[*,0,0])-1 do begin
           tmp=WHERE(mask[pix,*,*] GT 0.)
           IF tmp[0] NE -1 AND foundmin EQ 0. then begin
              countmin++
              minx[countmin]=pix
              foundmin=1
           ENDIF
           IF tmp[0] EQ -1 AND foundmin EQ 1. then begin
              maxx[countmin]=pix
              foundmin=0.
           ENDIF
        endfor
        IF maxx[countmin] EQ 0 then maxx[countmin]=n_elements(mask[*,0,0])-1
        foundmin=0.
        countmin=-1
        miny=intarr(sxpar(headermask,'NAXIS2'))
        maxy=intarr(sxpar(headermask,'NAXIS2'))
        for pix=0,n_elements(mask[0,*,0])-1 do begin
           tmp=WHERE(mask[*,pix,*] GT 0.)
           IF tmp[0] NE -1 AND foundmin EQ 0. then begin
              countmin++
              miny[countmin]=pix
              foundmin=1
           ENDIF
           IF tmp[0] EQ -1 AND foundmin EQ 1. then begin
              maxy[countmin]=pix
              foundmin=0.
           ENDIF
        endfor
        IF maxy[countmin] EQ 0 then maxy[countmin]=n_elements(mask[0,*,0])-1
        xhandle=MAX(maxx-minx)
        yhandle=MAX(maxy-miny)
        ringbuffer=fix(3./ring_spacing)
        if ringbuffer GT 6 then ringbuffer=6
        if ringbuffer LT 2 then ringbuffer=2
        maxrings=round(SQRT((xhandle/2.)^2+(yhandle/2.)^2)/((catmajbeam[i]*ring_spacing)/(ABS(sxpar(header,'cdelt1'))*3600.))+ringbuffer)

        IF NOT smallexists then begin
           newsize=fix((2.*(maxrings+1.)*(catmajbeam[i]*ring_spacing))/(ABS(sxpar(header,'cdelt1'))*3600.))
           IF floor(RApix[0]-newsize/2.) GE 0. AND floor(RApix[0]+newsize/2.) LT  sxpar(header,'NAXIS1') AND $
              floor(DECpix[0]-newsize/2.) GE 0. AND floor(DECpix[0]+newsize/2.) LT sxpar(header,'NAXIS2') AND $
              not smallexists then begin
              newdummy=fltarr(newsize+1,newsize+1,n_elements(dummy[0,0,*]))
              newdummy[*,*,*]=dummy[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.),*]
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"In order to speed up the fitting process significantly we cut the cube to "+string(newsize)+ " pixels."
                 close,66
              ENDIF
              sxaddpar,header,'NAXIS1',newsize+1
              sxaddpar,header,'NAXIS2',newsize+1
              sxaddpar,header,'CRPIX1',sxpar(header,'CRPIX1')-floor(RApix[0]-newsize/2.)
              sxaddpar,header,'CRPIX2',sxpar(header,'CRPIX2')-floor(DECpix[0]-newsize/2.)
              writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+'_small.fits',float(newdummy),header
              dummy=newdummy
              newdummy[*,*,*]=0.
              newdummy[*,*,*]=mask[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.),*]
              writefits,maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'_small.fits',float(newdummy),header
              mask=newdummy
              moment0map=readfits(maindir+'/'+catdirname[i]+'/'+catMom0Name[i]+'.fits',headermap,/SILENT)
              newdummy=fltarr(newsize+1,newsize+1)
              sxaddpar,headermap,'NAXIS1',newsize+1
              sxaddpar,headermap,'NAXIS2',newsize+1
              sxaddpar,headermap,'CRPIX1',sxpar(headermap,'CRPIX1')-floor(RApix[0]-newsize/2.)
              sxaddpar,headermap,'CRPIX2',sxpar(headermap,'CRPIX2')-floor(DECpix[0]-newsize/2.)
              newdummy[*,*]=moment0map[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.)]
              writefits,maindir+'/'+catdirname[i]+'/'+catMom0name[i]+'_small.fits',float(newdummy),headermap
              moment0map=newdummy
              moment1map=readfits(maindir+'/'+catdirname[i]+'/'+catMom1Name[i]+'.fits',headervel,/SILENT)
              newdummy=fltarr(newsize+1,newsize+1)
              sxaddpar,headervel,'NAXIS1',newsize+1
              sxaddpar,headervel,'NAXIS2',newsize+1
              sxaddpar,headervel,'CRPIX1',sxpar(headervel,'CRPIX1')-floor(RApix[0]-newsize/2.)
              sxaddpar,headervel,'CRPIX2',sxpar(headervel,'CRPIX2')-floor(DECpix[0]-newsize/2.)
              newdummy[*,*]=moment1map[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.)]
              writefits,maindir+'/'+catdirname[i]+'/'+catMom1name[i]+'_small.fits',float(newdummy),headervel
              moment1map=newdummy
              moment2map=readfits(maindir+'/'+catdirname[i]+'/'+catMom2Name[i]+'.fits',headerwidth,/SILENT)
              newdummy=fltarr(newsize+1,newsize+1)
              sxaddpar,headerwidth,'NAXIS1',newsize+1
              sxaddpar,headerwidth,'NAXIS2',newsize+1
              sxaddpar,headerwidth,'CRPIX1',sxpar(headerwidth,'CRPIX1')-floor(RApix[0]-newsize/2.)
              sxaddpar,headerwidth,'CRPIX2',sxpar(headerwidth,'CRPIX2')-floor(DECpix[0]-newsize/2.)
              newdummy[*,*]=moment2map[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.)]
              writefits,maindir+'/'+catdirname[i]+'/'+catMom2name[i]+'_small.fits',float(newdummy),headerwidth
              moment2map=newdummy
              noisemap=readfits(maindir+'/'+catdirname[i]+'/'+currentfitcube+'_noisemap.fits',headernoise,/SILENT)
              newdummy=fltarr(newsize+1,newsize+1)
              sxaddpar,headernoise,'NAXIS1',newsize+1
              sxaddpar,headernoise,'NAXIS2',newsize+1
              sxaddpar,headernoise,'CRPIX1',sxpar(headernoise,'CRPIX1')-floor(RApix[0]-newsize/2.)
              sxaddpar,headernoise,'CRPIX2',sxpar(headernoise,'CRPIX2')-floor(DECpix[0]-newsize/2.)
              newdummy[*,*]=noisemap[floor(RApix[0]-newsize/2.):floor(RApix[0]+newsize/2.),floor(DECpix[0]-newsize/2.):floor(DECpix[0]+newsize/2.)]
              writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+'_noisemap_small.fits',float(newdummy),headernoise
              noisemap=newdummy
              noisemapname=currentfitcube+'_noisemap_small'
              cubeext='.fits'
              CD, new_dir, CURRENT=old_dir
              spawn,'rm -f '+currentfitcube+'_cont.png '+currentfitcube+'_scat.png',isthere
              IF allnew LT 2 then begin
                 spawn,'rm -f '+catmaskname[i]+'.fits '+catMom1name[i]+'.fits '+catMom0name[i]+'.fits '+catMom2name[i]+'.fits '+currentfitcube+'_noisemap.fits',isthere
              ENDIF
              CD,old_dir
              currentfitcube=currentfitcube+'_small'
              smallexists=1
              catMom0name[i]=catMom0name[i]+'_small'
              catMom1name[i]=catMom1name[i]+'_small'
              catMom2name[i]=catMom2name[i]+'_small'
              catmaskname[i]=catmaskname[i]+'_small'
              Old=RApix[0]
              RApix[0]=floor(newsize/2.)+(RApix[0]-floor(RApix[0]))
              pixshift=Old-RApix[0]
              RApix[1:2]=RApix[1:2]-pixshift
              Old=DECpix[0]
              DECpix[0]=floor(newsize/2.)+(DECpix[0]-floor(DECpix[0]))
              pixshift=Old-DECpix[0]
              DECpix[1:2]=DECpix[1:2]-pixshift

           ENDIF
        ENDIF
     ENDIF
     centralflux=dummy[fix(RApix[0]),fix(DECpix[0]),fix(VSYSpix[0])]
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We have this value in the central pixel "+string(centralflux)+" at the pixel (x,y,z) "+string(fix(RApix[0]))+","+string(fix(DECpix[0]))+","+string(fix(VSYSpix[0]))
        IF not FINITE(centralflux) then  printf,66,linenumber()+"Central exclude = 1" else   printf,66,linenumber()+"Central exclude = 0"
        close,66
     ENDIF
     IF not FINITE(centralflux) then centralexclude=1 else centralexclude=0.
                                ; Print the source finder results
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The source finder found the following center in pixels."
        printf,66,linenumber()+"DEC center="+string(DECpix[0])+"DEC Boundaries"+string(DECpix[1])+","+string(DECpix[2])
        printf,66,linenumber()+"RA center="+string(RApix[0])+"RA Boundaries"+string(RApix[1])+","+string(RApix[2])
        printf,66,linenumber()+"V_sys center="+string(VSYSpix[0])+"V_sys Boundaries"+string(VSYSpix[1])+","+string(VSYSpix[2])
        printf,66,linenumber()+"Which converts to real coordinates as."
        printf,66,linenumber()+"DEC center="+string(DECdeg)+"DEC Boundaries"+string(DECboundeg[0])+","+string(DECboundeg[1])
        printf,66,linenumber()+"RA center="+string(RAdeg)+"RA Boundaries"+string(RAboundeg[0])+","+string(RAboundeg[1])
        printf,66,linenumber()+"V_sys center="+string(catVSYS[i])+"V_sys Boundaries"+string(ROTboun[0])+","+string(ROTboun[1])
        close,66
     ENDIF

                                ;We check whether the galaxy is bright
                                ;enough
     tmp=WHERE(mask EQ 1.)
     maxbright=MAX(dummy[tmp])
     maxSN=maxbright/catnoise[i]
     IF maxSN LT 2.5 then begin
        comment = 'The maximum Signal to Noise in this cube is '+strtrim(string(MaxSN),2)+' that is not enough for a fit.'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+catDirname[i]+' The maximum Signal to Noise in this cube is '+string(MaxSN)+' that is not enough for a fit.'
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+catDirname[i]+' The maximum Signal to Noise in this cube is '+string(MaxSN)+' that is not enough for a fit.'
        ENDELSE
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF
     pixelarea=ABS(3600.^2*sxpar(header,'CDELT2')*sxpar(header,'CDELT1'))
     pixelsizeRA=ABS(sxpar(header,'CDELT1'))
     pixelsizeDEC=ABS(sxpar(header,'CDELT2'))
     pixelsize=(pixelsizeRA+pixelsizeDEC)/2.
                                ;Let's get the size of the cube
     imagesize=sxpar(header,'CDELT2')*(sxpar(header,'NAXIS2')/2.)*3600.
                                ;IF it is not in km/s convert
     IF strupcase(veltype) EQ 'M/S' then channelwidth=channelwidth/1000.

     beamarea=(!pi*ABS(double(catminbeam[i])*double(catmajbeam[i])))/(4.*ALOG(2.))
     pixperbeam=beamarea/(ABS(pixelsizeRA*3600.)*ABS(pixelsizeDEC*3600.))
                                ;Obtain the cutoff values without
                                ;inclination corrections.
     rad=[0.,((findgen(maxrings+2.))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.)]
     calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i],channelwidth],cutoffor,vsys=catVSYS[i]
                                ;And print the values

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We use an intrinsic cutoff value of: "
        printf,66,linenumber()+"Radius           Minimum SBR"
        for j=0,n_elements(rad)-1 do begin
           printf,66,rad[j],cutoffor[j]
        endfor
        close,66
     ENDIF ELSE begin
        print,linenumber()+"We calculated the cutoff  values. "
     endelse
                                ;create the rings to be fitted with major beam
                                ;Let's calculate how many rings
                                ;we need to investigate the full HI
                                ;image. But let's reduce it by 1 rings
                                ;as to not run into edge problems
                                ;Initial size guess is
     norings=fix(round((imagesize)/(catmajbeam[i]*ring_spacing)))

                                ; If there is lt 1.5 beams in the cube we can stop here

     ;IF norings LT 1.5/ring_spacing then begin
     ;OR maxrings LE 2/ring_spacing then begin
     ;   openu,1,outputcatalogue,/APPEND
     ;   printf,1,format='(A60,2A12,A120)',catDirname[i],0.,0.,"This Cube is too small"
     ;   close,1
     ;   IF size(log,/TYPE) EQ 7 then begin
     ;      openu,66,log,/APPEND
     ;      printf,66,linenumber()+"We are aborting this fit as the located source is too small"
     ;      printf,66,linenumber()+"We have"+strtrim(string(norings[0]/ring_spacing))+" rings."
     ;      close,66
     ;   ENDIF
     ;   print,linenumber()+"We are aborting this fit as the located source is too small"
     ;   bookkeeping=5
     ;   goto,finishthisgalaxy
     ;ENDIF
                                ;Now let's obtain a initial pa
                                ;and inclination for the galaxy
     obtain_pa_incl,moment0map,newPA,newinclination,[RApix[0],DECpix[0]],NOISE=momnoise,BEAM=catmajbeam[i]/(pixelsize*3600.),gdlidl=gdlidl,extent=noringspix

     IF TOTAL([newPA,newinclination]) EQ 0. then begin
        comment = "No initial estimates. Likely the source is too faint."
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We are aborting this fit as we cannot establish proper initial estimates."
           close,66
        ENDIF
        print,linenumber()+"We are aborting this fit as we cannot establish proper initial estimates."
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF
                                ;Let's check wether the PA
                                ;lines up with the kinematical
                                ;or needs a 180  degree shift
                                ;since we have velocities we need to
                                ;check which half our PA should be
                                ;IF negative rotation values are in
                                ;the east the PA should be between 180
                                ;to 360 else between 0-180
                                ;This is so that the rotation curve is always positive
     dummyvel=moment1map
     tmpfinite=WHERE(FINITE(dummyvel) EQ 0)
     if tmpfinite[0] GT -1 then dummyvel[tmpfinite]=0.
     twobeampixel=3.*catmajbeam[i]/(ABS(sxpar(headermap,'cdelt1'))*3600.)
     case (1) of
        (sxpar(headervel,'cdelt1') LT 0):begin
           elowRA=fix(RApix[0]-twobeampixel)
           IF elowRA LT 0 then elowRA=0
           ehighRA=fix(RApix[0])
           IF ehighRA GT n_elements(dummyvel[*,0])-1 then ehighRA=n_elements(dummyvel[*,0])-1
           elowDEC=fix(DECpix[0]-twobeampixel)
           IF elowDEC LT 0 then elowDEC=0
           ehighDEC=fix(DECpix[0]+twobeampixel)
           IF ehighDEC GT n_elements(dummyvel[0,*])-1 then ehighDEC=n_elements(dummyvel[0,*])-1
           IF elowRA GE ehighRA then elowRA=0
           IF elowDEC GE ehighDEC then elowDEC=0
           wlowRA=fix(RApix[0])
           IF wlowRA LT 0 then wlowRA=0
           whighRA=fix(RApix[0]+twobeampixel)
           IF whighRA GT n_elements(dummyvel[*,0])-1 then whighRA=n_elements(dummyvel[*,0])-1
           wlowDEC=fix(DECpix[0]-twobeampixel)
           IF wlowDEC LT 0 then wlowDEC=0
           whighDEC=fix(DECpix[0]+twobeampixel)
           IF whighDEC GT n_elements(dummyvel[0,*])-1 then whighDEC=n_elements(dummyvel[0,*])-1
           IF wlowRA GE whighRA then whighRA=n_elements(dummyvel[*,0])-1
           IF wlowDEC GE whighDEC then whighDEC=n_elements(dummyvel[0,*])-1
        end
        (sxpar(headervel,'cdelt1') GT 0):begin
           elowRA=fix(RApix[0])
           IF elowRA LT 0 then elowRA=0
           ehighRA=fix(RApix[0]+twobeampixel)
           IF ehighRA GT n_elements(dummyvel[*,0])-1 then ehighRA=n_elements(dummyvel[*,0])-1
           elowDEC=fix(DECpix[0]-twobeampixel)
           IF elowDEC LT 0 then elowDEC=0
           ehighDEC=fix(DECpix[0]+twobeampixel)
           IF ehighDEC GT n_elements(dummyvel[0,*])-1 then ehighDEC=n_elements(dummyvel[0,*])-1
           IF elowRA GE ehighRA then ehighRA=n_elements(dummyvel[*,0])-1
           IF elowDEC GE ehighDEC then ehighDEC=n_elements(dummyvel[0,*])-1
           wlowRA=fix(RApix[0]-twobeampixel)
           IF wlowRA LT 0 then wlowRA=0
           whighRA=fix(RApix[0])
           IF whighRA GT n_elements(dummyvel[*,0])-1 then whighRA=n_elements(dummyvel[*,0])-1
           wlowDEC=fix(DECpix[0]-twobeampixel)
           IF wlowDEC LT 0 then wlowDEC=0
           whighDEC=fix(DECpix[0]+twobeampixel)
           IF whighDEC GT n_elements(dummyvel[0,*])-1 then whighDEC=n_elements(dummyvel[0,*])-1
           IF wlowRA GE whighRA then wlowRA=0
           IF wlowDEC GE whighDEC then wlowDEC=0
        end
        else:begin
           elowRA=0
           elowDEC=0
           ehighDEC=n_elements(dummyvel[*,0])-1
           ehighRA=n_elements(dummyvel[0,*])-1
           wlowRA=0
           wlowDEC=0
           whighDEC=n_elements(dummyvel[*,0])-1
           whighRA=n_elements(dummyvel[0,*])-1
        end
     endcase
     IF elowRA gt ehighRA then elowRA=0.
     IF elowDEC gt ehighDEC then elowDEC=0.
     tmpvelmap=dummyvel[elowRA:ehighRA,elowDEC:ehighDEC]
     nozervel=WHERE(tmpvelmap NE 0)
     East=TOTAL(tmpvelmap[nozervel])/n_elements(nozervel)
     tmpvelmap=dummyvel[wlowRA:whighRA,wlowDEC:whighDEC]
     nozervel=WHERE(tmpvelmap NE 0)
     West=TOTAL(tmpvelmap[nozervel])/n_elements(nozervel)
     IF ABS(EAST-WEST) GT channelwidth*(1.+vresolution) then begin
        IF East GT West then begin
           IF newPA[0] GT 180. then begin
              newPA[0]=newPA[0]-180.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We flipped the PA because the East side had higher velocties East="+string(East)+' West='+string(West)
                 close,66
              ENDIF
           ENDIF
           IF newPA[0] LT 170 and newPA[0] GT 10 and ABS(east-west) GT 20. then PAfixboun='e' else PAfixboun='n'
        endif else begin
           IF newPA[0] lT 180. then begin
              newPA[0]=newPA[0]+180.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We flipped the PA because the West side had higher velocties East="+string(East)+' West='+string(West)
                 close,66
              ENDIF

           ENDIF
           IF newPA[0] LT 350 and newPA[0] GT 190 AND ABS(east-west) GT 20. then PAfixboun='w'  else PAfixboun='n'
        endelse
     ENDIF ELSE BEGIN
        newPA[1]=360.
        PAfixboun='n'
     ENDELSE
     avPA=dblarr(2)
     mismatch=0.
                                ;get the kinematical PA
     obtain_velpa,moment1map,velPA,center=[RApix[0],DECpix[0]],intensity=moment0map
     IF ABS(velPA[0]-newPA[0]) LT 25. then begin
        avPA[0]=(velPA[0]/velPA[1]+newPA[0]/newPA[1])/(1./velPA[1]+1./newPA[1])
        avPA[1]=SQRT(velPA[1]^2+newPA[1]^2)/2.
     endif else begin
        IF FINITE(velPA[0]) then begin
           avPA[0]=velPA[0]
           avPA[1]=velPA[1]
;recalculate the inclination based on the velPA in case it is very
;different from the found PA
        ENDIF else avPA=newPA
        mismatch=1
     ENDELSE
                                ;recalculate the inclination



     obtain_w50,dummy,mask,header,W50
                                ;We want to adjust the cutoff values
                                ;with the inclination as edge-on
                                ;galaxies are naturally integrated
                                ;The test galaxy in Jozsa 2007 has an
                                ;inclination of 75 deg, hence above
                                ;that the cutoff value should be
                                ;declined and below that increased
     cutoffcorrection=SIN(75.*!DtoR)/SIN((newinclination[0]+newinclination[1]/2.)*!DtoR)
     IF newinclination[0] GT 80. then cutoffcorrection=0.9
     IF newinclination[0]+newinclination[1]/2. GT 50 then cutoffcorrection=1.
     IF newinclination[0]+newinclination[1]/2. LT 50 AND newinclination[0]+newinclination[1]/2. GT 40  then cutoffcorrection=1.+(50-(newinclination[0]+newinclination[1]/2.))*0.05
     IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
     IF doubled and cutoffcorrection GT 1. then cutoffcorrection=SQRT(cutoffcorrection)

     ringcorrection=cutoffcorrection
     WHILE maxrings-4-floor(ringcorrection/2.) LT 4 AND ringcorrection GT 1.03 do ringcorrection=ringcorrection-0.25

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with"+string(cutoffcorrection)
        printf,66,linenumber()+"And the rings are reduced by "+string(floor(ringcorrection/2.))
        close,66
     ENDIF
     cutoff=cutoffor*cutoffcorrection

                                ;Also our initial ring estimates
                                ;should take this inclination into
                                ;account
     case 1 of
        cutoffcorrection GT 1.03 :norings=maxrings-round((4.+floor(ringcorrection/2.))/ring_spacing)
        newinclination[0] GT 80:norings=maxrings-round(2./ring_spacing)
        else:norings=maxrings-round(4./ring_spacing)
     endcase
     IF norings[0] LT 1  OR norings[0]-round(1./ring_spacing) LT 1. then begin
        norings[0] = maxrings-ringbuffer
        if norings[0] LT 1. then norings[0]=1.
     endif else begin
        IF newinclination[0] LT 60 AND norings[0] GT 4 then norings=norings[0]-round(1./ring_spacing)
     endelse
     ;IF norings LT 3 then norings=3
     noringspix=norings[0]*(catmajbeam[i]*ring_spacing)/(ABS(sxpar(headermap,'cdelt1'))*3600.)

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'We find a W50 of'+string(W50)
        printf,66,linenumber()+'We find an extend of '+strtrim(string(fix(norings)),2)+' at this point from the inclination estimates.'
        close,66
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The original estimates are:"
        printf,66,linenumber()+"maxrings:"+string(fix(maxrings))
        printf,66,linenumber()+"norings:"+strtrim(string(fix(norings)),2)
        printf,66,linenumber()+"Kinematical PA:"+string(velPA[0])+'+/-'+string(velPA[1])
        printf,66,linenumber()+"Morphological PA:"+string(newPA[0])+'+/-'+string(newPA[1])
        printf,66,linenumber()+"noringspix:"+string(noringspix)
        printf,66,linenumber()+"newinclination:"+string(newinclination[0])+'+/-'+string(newinclination[1])
        printf,66,linenumber()+"PA:"+string(avPA[0])+'+/-'+string(avPA[1])
        close,66
     ENDIF
     newPA=avPA
     setfinishafter=0.
                                ; if the galaxy is too small we will
                                ; only fit a flat disk
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"With a ring size of "+strtrim(string(ring_spacing),2)+" we get "+strtrim(string(norings[0]),2)+" number of rings."
        close,66
     ENDIF
     if norings[0] LE 4. then begin
        ring_spacing_new = ring_spacing
        ring_spacing_up=ring_spacing
        ;we always want at least 3 beams in the model
        WHILE ring_spacing_new GT 0.5 AND norings[0] LE 4. do begin

           ring_spacing_new = ring_spacing_new/1.5
           IF ring_spacing_new LT 0.5 then ring_spacing_new=0.5
           if norings[0] LE 2. then begin
              norings[0]=round(norings[0]*ring_spacing_up/ring_spacing_new)
           endif else norings[0]=round(norings[0]*ring_spacing_up/ring_spacing_new-0.66)
           ring_spacing_up=ring_spacing_new
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Therefore we will reduce the ring size."
              printf,66,linenumber()+"With a ring size of "+strtrim(string(ring_spacing_new),2)+" we get "+strtrim(string(norings[0]),2)+"number of rings."
              close,66
           ENDIF
        ENDWHILE

        IF norings[0] LT 3. then begin
           comment = "This Cube is too small"
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We are aborting this fit as the located source is too small"
              close,66
           ENDIF
           print,linenumber()+"We are aborting this fit as the located source is too small"
           bookkeeping=5
           goto,finishthisgalaxy
        ENDIF
        maxrings=maxrings*ring_spacing/ring_spacing_new
        noringspix=norings[0]*(catmajbeam[i]*ring_spacing_new)/(ABS(sxpar(headermap,'cdelt1'))*3600.)
        rad=[0.,((findgen(round(maxrings)))*(catmajbeam[i]*ring_spacing_new)+catmajbeam[i]/5)]
        calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i],channelwidth],cutoffor,vsys=catVSYS[i]
        ring_spacing = ring_spacing_new

        ;if the model is less than 6 beams across we will not fit a warp
        if norings[0]*ring_spacing LT 3. then begin
                                ;we don't want to do individual
                                ;ring fits on this cube
           cutoffor=cutoffor/(norings[0])^0.4
           IF finishafter GT 1. then finishafter=1.1
           maxrings=maxrings-1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"This galaxy is smaller than 6 beams across hence we fit only a flat disc."
              close,66
           ENDIF
                                ;we also want to make sure that we are in the cube

        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Having not enough rings across the major axis we reduced the rings sizes "+catDirname[i]+"in galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           printf,66,linenumber()+"The new cutoff values are: "
           printf,66,linenumber()+"Radius           Minimum SBR"
           for j=0,n_elements(rad)-1 do begin
              printf,66,rad[j],cutoffor[j]
           endfor
           close,66
        ENDIF

        cutoff=cutoffor*cutoffcorrection
     ENDIF ELSE BEGIN
        if norings[0] GT maxrings then begin
           norings[0]=maxrings
           noringspix=norings[0]*catmajbeam[i]*ring_spacing/(ABS(sxpar(headermap,'cdelt1'))*3600.)
        ENDIF


     ENDELSE


     IF maxrings GT 25 AND norings[0] GT 20 then begin
        doubled=1
        tmpring=norings[0]-10.
        IF newinclination[0] LT 40 then tmpring=tmpring+2
        norings=10.+fix(tmpring/2.)
        rings=dblarr(norings[0])
        rings[0:9]=(findgen(10))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
        rings[10:norings[0]-1]=(findgen(fix(tmpring/2.)))*catmajbeam[i]*ring_spacing*2+catmajbeam[i]/5.+11.*catmajbeam[i]*ring_spacing
        norings[0]=norings[0]+1.
        rad=[0.,rings[0:9],(findgen(fix((maxrings-10.)/2.)))*catmajbeam[i]*ring_spacing*2+catmajbeam[i]/5.+11.*catmajbeam[i]]*ring_spacing
        IF finishafter GT 1 then finishafter=2.1
        maxrings=10.+fix((maxrings-10.)/2.)

        calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i],channelwidth],cutoffor,vsys=catVSYS[i]

        cutoff=cutoffor*cutoffcorrection
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Being a large galaxy we have doubled the outer the ring sizes "+catDirname[i]+"in galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           printf,66,linenumber()+"The new cutoff values are: "
           printf,66,linenumber()+"Radius           Minimum SBR"
           for j=0,n_elements(rad)-1 do begin
              printf,66,rad[j],cutoffor[j]
           endfor
           printf,66,linenumber()+"These are corrected with "+strtrim(string(cutoffcorrection),2)
           close,66
        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           printf,66,linenumber()+"Fitting with "+strtrim(string(ring_spacing*2),2)+" times the beam major axis FWHM beyond ring 10."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           print,linenumber()+"Fitting with "+strtrim(string(ring_spacing*2),2)+" times the beam major axis FWHM beyond ring 10."
        Endelse


     ENDIF ELSE BEGIN
        rings=(findgen(norings[0]))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
        norings[0]=norings[0]+1.
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           printf,66,linenumber()+"Fitting with "+strtrim(string(ring_spacing),2)+" times the beam major axis FWHM as ringsize."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The number of rings for the fitting = "+strtrim(string(fix(norings[0])),2)
           print,linenumber()+"Fitting with "+strtrim(string(ring_spacing),2)+"times the beam major axis FWHM as ringsize."
        Endelse
     ENDELSE

                                ;print the inclination and PA we found to the output log
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We are using an inclination of "+string(newinclination[0])+" +/- "+string(newinclination[1])
        printf,66,linenumber()+"We are using an PA of "+string(newPA[0])+" +/- "+string(newPA[1])
        close,66
     ENDIF
                                ; store the original estimates to compare against later on
     IF catPA[i] NE 0 AND catinc[i] NE 0 then begin
        UndatedPA=catPA[i]
        UndatedINCL=catinc[i]
        UndatedINCLdev=catincdev[i]
        Undatedmaxrot=catmaxrot[i]
     ENDIF ELSE BEGIN
        UndatedPA=newPA[0]
        UndatedINCL=newinclination[0]
        UndatedINCLdev=newinclination[1]
        IF catinc[i] LT 40 then begin
           Undatedmaxrot=W50/2./SIN(ABS(newinclination[0]-2.5)*!pi/180.)
        ENDIF ELSE Undatedmaxrot=W50/2./SIN(newinclination[0]*!pi/180.)

     ENDELSE
                                ;Write the new estimates as catalogue value
     catPA[i]=newPA[0]
     catPAdev[i]=newPA[1]
     catinc[i]=newinclination[0]
     catincdev[i]=newinclination[1]
                                ;calculate the vmax based on inclination and W50
     IF catinc[i] LT 40 then begin
        catmaxrot[i]=W50/2./SIN(ABS(newinclination[0]+5.)*!pi/180.)
     ENDIF ELSE catmaxrot[i]=W50/2./SIN(newinclination[0]*!pi/180.)
                                ;Let's optimize the cube for fitting
     IF catminbeam[i]/(cubecdelt*3600.) GT optpixelbeam  then begin
        newcdelt=catminbeam[i]/double(optpixelbeam)
        ratio=double(newcdelt/double((cubecdelt*3600.)))
        new_header=header
        new_Cube=CONGRID(dummy,fix(sxpar(header,'NAXIS1')/ratio),fix(sxpar(header,'NAXIS2')/ratio),sxpar(header,'NAXIS3'))
        ratio=double(sxpar(header,'NAXIS1'))/fix(sxpar(header,'NAXIS1')/ratio)
        sxaddpar,new_header,'CDELT1',sxpar(header,'CDELT1')*(double(sxpar(header,'NAXIS1'))/(fix(sxpar(header,'NAXIS1')/ratio))),format='(E20.12)'
        sxaddpar,new_header,'CDELT2',sxpar(header,'CDELT2')*(double(sxpar(header,'NAXIS2'))/(fix(sxpar(header,'NAXIS2')/ratio))),format='(E20.12)'
        sxaddpar,new_header,'CRPIX1',sxpar(header,'CRPIX1')/(double(sxpar(header,'NAXIS1'))/(fix(sxpar(header,'NAXIS1')/ratio)))+(1.-1./ratio),format='(E20.12)'
        sxaddpar,new_header,'CRPIX2',sxpar(header,'CRPIX2')/(double(sxpar(header,'NAXIS2'))/(fix(sxpar(header,'NAXIS2')/ratio)))+(1.-1./ratio),format='(E20.12)'
        sxaddpar,new_header,'NAXIS1',fix(sxpar(header,'NAXIS1')/ratio)
        sxaddpar,new_header,'NAXIS2',fix(sxpar(header,'NAXIS2')/ratio)
        CD,new_dir,CURRENT=old_dir
        currentfitcube=currentfitcube+'_opt'
        writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+cubeext,float(new_Cube),new_header
        cd,old_dir
        dummy=new_Cube
        header=new_header
        optimized=1
     ENDIF
                                ;Let's write a file with the basic initial parameters
     oldRA=RAdeg
     oldDEC=DECdeg
     oldVSYS=catVSYS[i]
     RAhr=RAdeg
     RAdiff=TOTAL(ABS(RADeg-RAboundeg))/2.*3600./15.
     DEChr=DECdeg
     DECdiff=TOTAL(ABS(DECDeg-DECboundeg))/2.*3600.
     VSYSdiff=TOTAL(ABS(catVSYS[i]-ROTboun))
     totflux[0]=totflux[0]/pixperbeam
     IF FINITE(totflux[0]) NE 1 then begin
        tmp=WHERE(FINITE(moment0map))
        totflux=TOTAL(moment0map[tmp])/pixperbeam/ABS(channelwidth)
     ENDIF
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr
     getDHI,moment0map,headermap,catPA[i],[RADeg,DECdeg,catVSYS[i]],DHI
     catmaxrotdev[i]=(VSYSdiff/4.)/SIN((newinclination[0])*!pi/180.)
     IF catmaxrot[i] LT 20*channelwidth AND catmaxrotdev[i] LT 25.*channelwidth then  catmaxrotdev[i]=25.*channelwidth
     IF catmaxrotdev[i] LT 4*channelwidth then catmaxrotdev[i]=4*channelwidth
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We are using maxrot of "+string(catmaxrot[i])+" +/- "+string(catmaxrotdev[i])
        close,66
     ENDIF
     basicinfofile=maindir+'/'+catdirname[i]+'/'+basicinfo+'_'+catcubename[i]+'.txt'
     openw,1,basicinfofile
     printf,1,"#This file contains the basic parameters of the Galaxy"
     printf,1,"#The total flux is determined by the source finder in the initial estimates and the total emission in the masked cube of the final model in the fits"
     printf,1,"#D_HI is determined as the diameter where the major axis with given PA cuts the 10^20 column of the moment0 map"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)', 'RA','DEC','VSYS','PA','Inclination','Max VRot','V_mask', 'Tot FLux','D_HI', 'Distance','HI Mass','D_HI'
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)', 'hh:mm:ss','dd:mm:ss','km/s','deg','deg','km/s','km/s', 'Jy','arcsec', 'Mpc','M_solar','kpc'
     printf,1,"#The initial input"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catVSYS[i])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catPA[i])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catinc[i])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(catmaxrot[i])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(VSYSdiff)),2)+'+/-'+strtrim(strcompress(string(channelwidth*(1.+vresolution),format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1

     IF Totflux[0] LT 0. then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This galaxy has negative total flux. That will not work. Aborting "
           close,66
        ENDIF
        comment = 'We found an initial negative total flux.'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF

                                ;and we make a pv-diagram based on these parameters
     IF optimized then begin
        noptname=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        nooptcube=readfits(maindir+'/'+catdirname[i]+'/'+noptname[0]+'.fits',noptheader,/NOSCALE,/SILENT)
        extract_pv,nooptcube,noptheader,catpa[i],xv,center=[RAdeg,DECdeg],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_0_xv.fits',float(xv),new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,header,catpa[i],xv,center=[RAdeg,DECdeg],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+'_0_xv.fits',float(xv),new_header
     ENDELSE
                                ;break if we do not want to use tirific
     IF finishafter EQ 0. then begin
        if setfinishafter NE 1 then begin
           comment = 'You have chosen to skip the fitting process after all preparations for the fit'
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
        ENDIF
        if optimized then begin
           tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
           currentfitcube=tmp[0]
        endif
                                ;       bookkeeping=0
        goto,finishthisgalaxy
     ENDIF
                                ;build up moment 0 axis
     ;buildaxii,headermap,xaxmom0,yaxmom0
                                ;some counters for keeping track
     prevmodification=0.
     overwrite=0.
     lastsbr=100
     counter=0.
     countsbr=0.
     restart_counter = 0
     continue_tirific = 'initialized'
     curr_run_id = 'stopped'
     case 1 of
       norings[0] LT 5: inimode = 1
       norings[0] LT 15: inimode = 2
       else: inimode = 3
     endcase
     ringmodifier=0
     plus2beamshift=0
     shiftcentercounter=0.
                                ;CD to the correct directory
     CD, new_dir, CURRENT=old_dir
     newrings=fix(norings[0])
                                ;Calculate the beam in pixels  and get an estimate for the peak sbr
     beaminpixels=fix(catmajbeam[i]/((ABS(pixelsizeRA)+ABS(pixelsizeDEC))/2.*3600.))
     centralarea=moment0map[fix(RApix[0]-beaminpixels):fix(RApix[0]+beaminpixels),fix(DECpix[0]-beaminpixels):fix(DECpix[0]+beaminpixels)]
     cenav=TOTAL(centralarea[WHERE(FINITE(centralarea))])/n_elements(centralarea[WHERE(FINITE(centralarea))])
     peaksbr=cenav*channelwidth/(catmajbeam[i]*catminbeam[i])
                                ; at high inclinations rings overlap
                                ; so the peaksbr needs to be reduced
                                ; according to the inclination and the
                                ; physical size of the galaxy
     IF catinc[i] GT 70. then begin
        IF convertskyanglefunction(DHI,catDistance[i]) GT 1. then $
           peaksbr=peaksbr/(convertskyanglefunction(DHI,catDistance[i])*SIN(catinc[i]*!DtoR)*0.5) else $
              peaksbr=peaksbr/(convertskyanglefunction(3.*catmajbeam[i],catDistance[i])*SIN(catinc[i]*!DtoR)*0.5)
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We estimate the peaksbr to be "+string(peaksbr)
        close,66
     ENDIF

                                ;Calculate an initial surface
                                ;brightness profile based on the
                                ;central brightness and an exponential
     velaverage=catmaxrot[i]
     signvelaverage=catmaxrot[i]
     stringsbr='SBR= '+string(peaksbr)+' '+string(peaksbr)
     scalelength=convertskyanglefunction(7.5,catDistance[i],/PHYSICAL)
     for j=1,norings[0]-2 do begin
        exposbr=peaksbr*exp(-(rings[j]-rings[0])/(scalelength))
        stringsbr=stringsbr+string(exposbr)+' '
     endfor
     tmppos=where('SBR' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringsbr
                                ;setting the values for the 2nd disk
     tmppos=where('SBR_2' EQ tirificfirstvars)
     tmp=str_sep(strtrim(strcompress(stringsbr),2),'=')
     tirificfirst[tmppos]='SBR_2='+tmp[1]
     tmp2=str_sep(strtrim(strcompress(tmp[1]),2),' ')


                                ;but that does mean we'd like
                                ;to change SBR a lot so we set quite
                                ;large starting steps
     string1=''
     string2=''
     string3=''
     string4=''
     string5=''
     string6=''
     string7=''

                                ;The minimum is based on the cutoff
                                ;values and therefore we need to set
                                ;every ring in the tirific file
     IF norings[0]*ring_spacing LT 7 then begin
        for j=norings[0],3,-1 do begin
           string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+' SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
           string2=string2+' 1'
           string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
           string4=string4+' 7.5E-5'
           string5=string5+' 2E-6'
           string6=string6+' 5E-6'
           string7=string7+' 3'
        endfor
     ENDIF else begin
        for j=norings[0],4,-1 do begin
           string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+', SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
           string2=string2+' 1 1'
           string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
           string4=string4+' 7.5E-5 7.5E-5'
           string5=string5+' 2E-6 2E-6'
           string6=string6+' 5E-6 5E-6'
           string7=string7+' 3 3'
        endfor
     ENDELSE
     string1 = STRMID(string1,1,STRLEN(string1)-1)
     SBRinput1=[string1,string2,string3,string4,string5,string6,string7]
     SBRinput2=[' SBR 1 2 3 SBR_2 1 2 3',$
                strtrim(strcompress(string(peaksbr,format='(E12.5)'))),'0','1E-5','1E-6','5E-5','1E-6','3','70','70']
                                ;Now setting some general values
     tmppos=where('INSET' EQ tirificfirstvars)
     tirificfirst[tmppos]='INSET='+strtrim(strcompress(string(currentfitcube+'.fits')))
     tmppos=where('BMIN' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMIN='+strtrim(strcompress(string(catminbeam[i])))
     tmppos=where('BMAJ' EQ tirificfirstvars)
     tirificfirst[tmppos]='BMAJ='+strtrim(strcompress(string(catmajbeam[i])))
     tmppos=where('BPA' EQ tirificfirstvars)
     tirificfirst[tmppos]='BPA='+strtrim(strcompress(string(catbpa[i])))
     tmppos=where('SDIS' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS=  8.'
     tmppos=where('SDIS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='SDIS_2=  8.'
     tmppos=where('CONDISP' EQ tirificfirstvars)
     IF vresolution EQ 0. then begin
        tirificfirst[tmppos]='CONDISP=  '+string(channelwidth*1.2/(2*SQRT(2*ALOG(2.))))
     ENDIF ELSE BEGIN
        tirificfirst[tmppos]='CONDISP=  '+string((channelwidth*(1.+vresolution))/(2*SQRT(2*ALOG(2.))))
     ENDELSE
     tmppos=where('RMS' EQ tirificfirstvars)
     tirificfirst[tmppos]='RMS='+strtrim(strcompress(string(catnoise[i])))
     tmppos=where('DISTANCE' EQ tirificfirstvars)
     tirificfirst[tmppos]='DISTANCE='+strtrim(strcompress(string(catDistance[i])))
     prevrings=norings[0]
                                ; defining some flags that we will use
                                ; for checks on what is going on in
                                ; the fitting
                                ;the First estimate of Z0 is 0.2 kpc
                                ; if 0.2 kpc is less
                                ;than 1/4th of the beam we want to
                                ;make it a quarter of the beam at high inclinations
     If catdistance[i] EQ 1. then begin
        IF catinc[i] GT 80. then inpzval=MAX([(catmajbeam[i]*ring_spacing*norings[0])/150.,catmajbeam[i]/4.]) else  inpzval=(catmajbeam[i]*ring_spacing*norings[0])/150.
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string((inpzval))))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string((inpzval))))
     ENDIF ELSE BEGIN
        IF catinc[i] GT 80. then inpzval=MAX([convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL),catmajbeam[i]/4.]) else inpzval=convertskyanglefunction(0.2,double(catDistance[i]),/PHYSICAL)
        tmppos=where('Z0' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0='+strtrim(strcompress(string(inpzval)))
        tmppos=where('Z0_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='Z0_2='+strtrim(strcompress(string(inpzval)))
     ENDELSE
     constring=0.
     forcedring=0.
     lastcutrings=0.
     lastaddrings=0.
     secondtime=0.
     sofiarings=norings[0]
     stringvelocities='VROT= 0. '+string(ABS(double(catmaxrot[i])))
     tmppos=where('CFLUX' EQ tirificfirstvars)
     IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX= 1e-5'
     tmppos=where('CFLUX_2' EQ tirificfirstvars)
     IF float(string((totflux[0])/7.5E5)) NE 0. then tirificfirst[tmppos]='CFLUX_2= '+string((totflux[0])/7.5E5) else tirificfirst[tmppos]='CFLUX_2= 1e-5'
                                ;If we change the amount of rings we need to come back here
     sbrshift:
     case finishafter of
        2.1:begin
           tmpring=norings[0]-10.
           rings=dblarr(norings[0])
           rings[0:9]=(findgen(10))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
           rings[10:norings[0]-1]=(findgen(fix(tmpring)))*catmajbeam[i]*ring_spacing*2+catmajbeam[i]/5.+11.*catmajbeam[i]*ring_spacing
        end
       else:rings=(findgen(norings[0]))*catmajbeam[i]*ring_spacing+catmajbeam[i]/5.
     endcase


     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
        close,66
     ENDIF ELSE begin
        print,linenumber()+"The number of rings used = "+strtrim(string(fix(norings[0])),2)
     endelse
     fluxadjust=0.
                                ;let's write the number of rings to the tirific file
     tmppos=where('NUR' EQ tirificfirstvars)
     tirificfirst[tmppos]='NUR='+strtrim(strcompress(string(norings[0])),1)
                                ;and cflux


                                ;Now we write the radii of our rings
     stringring='RADI=0.0 '
     for j=0,norings[0]-2 do begin
        stringring=stringring+string(rings[j])+' '
     endfor
     tmppos=where('RADI' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringring
                                ;Using the parameters from sofia or a
                                ;previous fit
     tmppos=where('VROT' EQ tirificfirstvars)
     tirificfirst[tmppos]=stringvelocities
     tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
     tmppos=where('VROT_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VROT_2='+tmp[1]
     tmppos=where('INCL' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
     tmppos=where('PA' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA='+strtrim(strcompress(string(catPA[i])))
     tmppos=where('XPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS='+strtrim(strcompress(string(catvsys[i])))
     tmppos=where('INCL_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
     tmppos=where('PA_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='PA_2='+strtrim(strcompress(string(catPA[i])))
     tmppos=where('XPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='XPOS_2='+strtrim(strcompress(string(RAdeg)))
     tmppos=where('YPOS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='YPOS_2='+strtrim(strcompress(string(DECdeg)))
     tmppos=where('VSYS_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='VSYS_2='+strtrim(strcompress(string(catvsys[i])))

                                ; Set the tirific fitting parameters for the inclination
                                ;If the inclination is above 75 the
                                ;initial estimates must be good so we
                                ;want less change
     case 1 of
        catinc[i] LT 40: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.2),string(0.1),string(5.0),string(0.1),'5','70','70']
           ;We put this inside might not work then take it out again
           IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
           IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
        end
        catinc[i] LT 75: begin
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(5.),string(0.1),string(1.0),string(0.1),'3','70','70']

           IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
           IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
        end
        else:INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                         string(catinc[i]+catincdev[i]),string(catinc[i]-catincdev[i]),string(1.0),string(0.01),string(1.0),string(0.001),'3','70','70']
     endcase
                                ; have to ensure that the parameters are within the limits
     IF INCLinput1[1] GT 90. then INCLinput1[1]='90.'
     IF INCLinput1[1] LT 30.then INCLinput1[1]='30.'

                                ;Set the input for the PA

     PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string((catPA[i])+catPAdev[i]+10.),string((catPA[i])-catPAdev[i]-10),string(1.),string(0.1),string(1.),string(0.1),'3','70','70']

                                ;If the estimate PA isn't close
                                ;to 180 or 360 we want to block
                                ;changes that can flip the rotation curve
     IF PAfixboun NE 'n' then begin
        IF PAfixboun EQ 'w' then begin
           IF catPA[i]+catPAdev[i]+10. GT 360 then PAinput1[1]='360'
           IF catPA[i]-catPAdev[i]-10. LT 180 then PAinput1[2]='180'
        ENDIF
        IF PAfixboun EQ 'e' then begin
           IF catPA[i]+catPAdev[i]+10. GT 180. then PAinput1[1]='180'
           IF catPA[i]-catPAdev[i]-10. LT 0. then PAinput1[2]='0'
        ENDIF
     ENDIF
                                ;Also the minimum and maximum of Z0 which have to be based on physical
                                ;values
     Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.1,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(1.0,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.01,double(catDistance[i]),/PHYSICAL)),'3','70','70']
                                ;if the inclination is high we want to
                                ;make sure that the scale height can
                                ;be at least half a beam or 1 kpc
                                ;whichever is larger because otherwise
                                ;it might be pushing the inclination
                                ;down.

     IF catinc[i] GT 80 then begin
        Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
     ENDIF
     IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([catmajbeam[i]*norings[0]*ring_spacing/66.,catmajbeam[i]/2.]))
        ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*ring_spacing*norings[0]/66.)
        Z0input1[2]='0.'
        Z0input1[3]=string(-1*catmajbeam[i]*ring_spacing*norings[0]/1.5E5)
        Z0input1[4]=string(catmajbeam[i]*ring_spacing*norings[0]/1.5E6)
        Z0input1[5]=string(catmajbeam[i]*ring_spacing*norings[0]/1500)
        Z0input1[6]=string(catmajbeam[i]*ring_spacing*norings[0]/1.5E6)
     ENDIF





                                ;Now The rotation

     VROTmax=catmaxrot[i]+catmaxrotdev[i]
     VROTmin=channelwidth
                                ;ensure reasonable vrotmax dependent
                                ;on inclination as more unsure at low inclination
     IF VROTmax LT 80. then VROTmax=80.
     IF VROTmax GT 600. then VROTmax=600.

     string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
     string2=string(VROTmax)
     string3=string(VROTmin)
     string4=string(channelwidth)
                                ;At low inclinations we want to be
                                ;very accepting of the rotational
                                ;values as they are very unsure.
     case 1 of
        catinc[i] LT 30:begin
           string4=string(2*channelwidth)
           string5=string(0.02*channelwidth)
           string6=string(0.2*channelwidth)
           string7=string(0.02*channelwidth)
        end
        else:begin
           string4=string(channelwidth)
           string5=string(0.01*channelwidth)
           string6=string(0.1*channelwidth)
           string7=string(0.01*channelwidth)
        end
     ENDCASE
     string8='3'
     string9='70'
                                ;If the model is large enough we fit
                                ;only a slope to the outer quarter of
                                ;the model at low inclinations we want
                                ;this slope to be at least half of the
                                ;rings and be flat
     IF centralexclude then begin
        IF catinc[i] LT 40 then begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDIF ELSE BEGIN
        IF catinc[i] LT 40 then begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        ENDIF else begin
           string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
        ENDELSE
     ENDELSE
     IF norings[0] GT 4 then begin
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
     ENDIF ELSE BEGIN
        VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
     ENDELSE
     vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                  strtrim(strcompress(string(catvsys[i]+100.)),1),strtrim(strcompress(string(catvsys[i]-100.)),1),'10','0.01','0.5','0.01','3','70','70']

                                ;IF we have a small cube we will only accept very small changes to the
                                ;central position
     IF norings[0] LT 4 then begin
        maxxpos=strcompress(string(RAdeg+catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-catmajbeam[i]/3600.))
        ;We'll be satified with a change less than a tenth of a pixel
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), +$
                      maxxpos,minxpos,strtrim(string(pixelsizeRA)),strtrim(string(pixelsizeRA/20.)),+$
                      strtrim(string(pixelsizeRA/10.)),strtrim(string(pixelsizeRA/20.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                     maxypos,minypos,strtrim(string(pixelsizeDEC)),strtrim(string(pixelsizeDEC/20.)),+$
                      strtrim(string(pixelsizeDEC/10.)),strtrim(string(pixelsizeDEC/20.)),'3','70','70']
     ENDIF ELSE BEGIN
        maxxpos=strcompress(string(RAdeg+3.*catmajbeam[i]/3600.))
        minxpos=strcompress(string(RAdeg-3.*catmajbeam[i]/3600.))
        maxypos=strcompress(string(DECdeg+3.*catmajbeam[i]/3600.))
        minypos=strcompress(string(DECdeg-3.*catmajbeam[i]/3600.))
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1), +$
                      maxxpos,minxpos,strtrim(string(pixelsizeRA*2.)),strtrim(string(pixelsizeRA/10.)),+$
                      strtrim(string(pixelsizeRA/2.)),strtrim(string(pixelsizeRA/10.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),+$
                     maxypos,minypos,strtrim(string(pixelsizeDEC*2.)),strtrim(string(pixelsizeDEC/10.)),+$
                      strtrim(string(pixelsizeDEC/2.)),strtrim(string(pixelsizeDEC/10.)),'3','70','70']
     ENDELSE
     tryone=0.
                                ;   IF we have a low inclination we first want to fit the PA by itself
     PAest=0.
     INCLest=0.
     INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                 ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                 string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(2),string(0.1),string(0.1),'3']
     IF (catinc[i]  LT 50 OR mismatch EQ 1 OR ceil(norings[0]*COS(catinc[i]*!DtoR))*ring_spacing LE 5) AND counter EQ 0. then begin
                                ;If we have a small number of beams
                                ;acros the minor axis the inclination
                                ;is very unsure and we first want to
                                ;fit the inclination. If the KinPA and VelPA do not match we want to fit the PA
        PAinincl=catPA[i]
        INCLinincl=catinc[i]
        maxrotinincl=catmaxrot[i]
        fixstring='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
        IF ceil(norings[0]*COS(catinc[i]*!DtoR))*ring_spacing LE 5 OR catinc[i] LT 30. OR mismatch EQ 1 then begin

                                ; Avoid runaway rotation
           IF mismatch EQ 1 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Because of a difference between the kinematic PA and the morphological PA we will scramble the inclination."
                 printf,66,linenumber()+"Original inclination ="+string(catinc[i])+" New = "+string(catinc[i]-5.)
                 close,66
              ENDIF
              catinc[i]=catinc[i]-5.
              catincdev[i]=catincdev[i]+10.
           ENDIF ELSE BEGIN
              catinc[i]=catinc[i]-3.
              catincdev[i]=catincdev[i]+5.
           ENDELSE
           IF catinc[i] LT 5 then catinc[i]=5.

           tmppos=where('INCL' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
           tmppos=where('INCL_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
           TMP_VEL=W50/2./SIN(ABS(INCLinincl)*!pi/180.)
           IF TMP_VEL GT 300. then TMP_VEL= 300.
           tmppos=where('VROT' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT= 0. '+strtrim(strcompress(string(tmp_vel)))
           tmppos=where('VROT_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT_2= 0. '+strtrim(strcompress(string(tmp_vel)))
           INCLest=1
           ;INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
           ;             ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
           ;             '90.',string(catinc[i]-catincdev[i]-5),string(1.),string(0.1),string(0.5),string(0.05),'5','70','70']
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                        ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                        '90.',string(catinc[i]-catincdev[i]-5),string(5.),string(0.5),string(0.1),'5']
           ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
           ;           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
           ;           string((catPA[i])+catPAdev[i]+40.),string((catPA[i])-catPAdev[i]-40),string(6),string(0.1),string(1.0),string(0.01),'3','70','70']
           IF mismatch EQ 1 THEN BEGIN

             PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                      string((catPA[i])+catPAdev[i]+40.),string((catPA[i])-catPAdev[i]-40),string(10),string(1.),string(0.5),'3']
           ENDIF ELSE BEGIN
             PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                      string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(1.),string(0.25),'3']
           ENDELSE

           IF INCLinput1[2] LT 2. then INCLinput1[2]='2'

           IF PAfixboun NE 'n' then begin
              IF PAfixboun EQ 'w' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
                 IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
              ENDIF
              IF PAfixboun EQ 'e' then begin
                 IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
                 IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
              ENDIF
           ENDIF

           ;VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
          ;                ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
          ;                ,string(VROTmax),string(0.),$
          ;                string(channelwidth),string(0.1*channelwidth),string(0.5*channelwidth),string(0.01*channelwidth),'3','70','70',fixstring]
           VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
                          ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
                          ,string(VROTmax),string(0.), string(4*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth),'5',fixstring]
           Writefittingvariables,tirificfirst,inclinput1,VROTinputINCL,painput1,sbrinput1,sbrinput2,INIMODE=inimode

           againINCLestimate:
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because there are only a few beams across the minor axis or we have low inclination we first adjust the inclination "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           restart_counter++
           tmppos=where('RESTARTID' EQ tirificfirstvars)
           tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))

           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificfirst)-1 do begin
              printf,1,tirificfirst[index]
           endfor
           close,1

           IF testing GE 1 then goto,testing1INCL
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           print,linenumber()+"Starting tirific the INCL estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()


           run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name='1stfit.', $
                  store_name='1stfitold.', log=log, loops=loops,nopoints=nopoints, $
                  AC=AC1,toymodels=toymodels,run_unit=run_unit
           IF bookkeeping EQ 5 THEN goto,finishthisgalaxy

           ;If not accepted we try again
           ;IF AC1 EQ 0. and INCLest LT 1 then begin
            ;  IF size(log,/TYPE) EQ 7 then begin
            ;     openu,66,log,/APPEND
            ;     printf,66,linenumber()+"The INCL estimate was not accepted try again. Tries: "+string(INCLest)
            ;     close,66
            ;  ENDIF
            ;  VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2','SBR','SBR_2']
            ;  firstfitvalues=0.
            ;  writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
            ;  catPA[i]=firstfitvalues[0,0]
            ;  catinc[i]=firstfitvalues[0,2]
            ;  catmaxrot[i]=firstfitvalues[1,4]
            ;  IF catinc[i] LT 5 then catinc[i]=5.
            ;  PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(1.),string(PAinput1[5]*(1+INCLest*0.25)),'5']
            ;  ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;  ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;  ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
            ;  IF PAfixboun NE 'n' then begin
            ;     IF PAfixboun EQ 'w' then begin
            ;        IF catPA[i]+catPAdev[i]+40. GT 360 then PAinput1[1]='360'
            ;        IF catPA[i]-catPAdev[i]-40. LT 180 then PAinput1[2]='180'
            ;     ENDIF
            ;     IF PAfixboun EQ 'e' then begin
            ;        IF catPA[i]+catPAdev[i]+40. GT 180. then PAinput1[1]='180'
            ;        IF catPA[i]-catPAdev[i]-40. LT 0. then PAinput1[2]='0'
            ;     ENDIF
            ;  ENDIF
            ;  INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
            ;            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
            ;            '90.',string(catinc[i]-catincdev[i]-5),string(5.),string(0.5),string(0.1*(1.+INCLest*0.5)),'5']
            ;  VROTinputINCL=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
            ;              ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
            ;              ,string(VROTmax),string(0.), string(4*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth*(1.+INCLest*0.5)),'5',fixstring]
              ;INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
              ;            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
              ;            '90.',string(catinc[i]-catincdev[i]-20),string(1.),string(0.5),string(0.1),string(0.01),'3','70','70']
            ;  Writefittingvariables,tirificfirst,painput1,inclinput1,vrotinputINCL,INIMODE=inimode


            ;  INCLest=INCLest+1.
            ;  goto,againINCLestimate
           ;ENDIF
           ;ELSE BEGIN

              ;When accepted update our values
            ;  VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2','SBR','SBR_2']
            ;  firstfitvalues=0.
          ;    writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted

          ;    IF size(log,/TYPE) EQ 7 then begin
          ;       openu,66,log,/APPEND
          ;       printf,66,linenumber()+"Because of low inclination (<30) or small number of beams across minor axis (< 5):"
          ;       printf,66,linenumber()+"We have adjusted the PA from "+string(PAinincl)+" to "+string(firstfitvalues[0,0])+'+/-'+string(newPA[1])
          ;       printf,66,linenumber()+"The inclination from "+string(INCLinincl)+" to "+string(firstfitvalues[0,2])+'+/-'+string(newinclination[1])
          ;       printf,66,linenumber()+"The max rotation from "+string(maxrotinincl)+" to "+string(firstfitvalues[n_elements(firstfitvalues[*,0])-2,4])+'+/-'+string(catmaxrotdev[1])
          ;       printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT."
          ;       close,66
          ;    ENDIF
          ;    catPA[i]=firstfitvalues[0,0]
          ;    catinc[i]=firstfitvalues[0,2]
          ;    IF catinc[i] LT 5 then catinc[i]=5.
          ;    catmaxrot[i]=firstfitvalues[1,4]


           ;ENDELSE
        ENDIF
        testing1INCL:


        ;We only want to check the PA if we have not checked the INCL
        IF INCLest EQ 0 THEN BEGIN
          ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
          ;          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
          ;          string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(2.5),string(0.1),string(0.05),string(0.01),'3','70','70']
          PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.5),string(0.1),'3']
          IF PAfixboun NE 'n' then begin
             IF PAfixboun EQ 'w' then begin
                IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
             ENDIF
             IF PAfixboun EQ 'e' then begin
                IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
             ENDIF
          ENDIF
          ;VROTinputPA=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
          ;             ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
          ;             ,string(VROTmax),string(0.),$
          ;             string(channelwidth),string(0.01*channelwidth),string(channelwidth),string(0.01*channelwidth),'3','70','70',fixstring]
          VROTinputPA=['!VROT '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2 '+$
                       ' VROT_2 '+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+':2'$
                       ,string(VROTmax),string(0.), string(3*channelwidth),string(0.1*channelwidth),string(0.1*channelwidth),'5',fixstring]
          Writefittingvariables,tirificfirst,painput1,sbrinput1,sbrinput2,VROTinputPA,INIMODE=inimode


          againPAestimate:
          restart_counter++
          tmppos=where('RESTARTID' EQ tirificfirstvars)
          tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))
          openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
          for index=0,n_elements(tirificfirst)-1 do begin
             printf,1,tirificfirst[index]
          endfor
          close,1

          IF testing GE 1 then goto,testing1PA
          IF size(log,/TYPE) EQ 7 then begin
             openu,66,log,/APPEND
             printf,66,linenumber()+"Starting tirific the PA estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
             close,66
          ENDIF
          print,linenumber()+"Starting tirific the PA estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()

          run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name='1stfit.', $
                    store_name='1stfitold.', log=log,loops=loops,nopoints=nopoints,AC=AC1, $
                    toymodels=toymodels,run_unit=run_unit
          IF bookkeeping EQ 5 THEN goto,finishthisgalaxy

          ;IF AC1 EQ 0. and PAest LT 1 then begin
          ;   IF size(log,/TYPE) EQ 7 then begin
          ;      openu,66,log,/APPEND
          ;      printf,66,linenumber()+"The Pa estimate was not accepted try again. Tries: "+string(PAest)
          ;      close,66
        ;     ENDIF
        ;     VariablesWanted=['PA','PA_2','VROT','VROT_2','SBR','SBR_2']
        ;     firstfitvalues=0.
        ;     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
        ;     catPA[i]=firstfitvalues[0,0]
        ;     IF catPA[i] GT 360. then begin
        ;        WHILE catPA[i] GT 360. do catPA[i]=catPA[i]-360.
        ;     ENDIF
        ;     PAinput1[5]=string(PAinput1[5]*(1+PAest*0.5))
             ;PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
              ;         ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
              ;         string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
        ;     IF PAfixboun NE 'n' then begin
        ;        IF PAfixboun EQ 'w' then begin
        ;           IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
        ;           IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
        ;        ENDIF
        ;        IF PAfixboun EQ 'e' then begin
        ;           IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
        ;           IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
        ;        ENDIF
        ;     ENDIF
        ;     IF size(log,/TYPE) EQ 7 then begin
        ;         openu,66,log,/APPEND
        ;         printf,66,linenumber()+"Trying the PA estimate again with PA="+strtrim(string(catPA[i]))+" and the boundaries ["+STRJOIN([strtrim(string(PAinput1[2]),2),strtrim(string(PAinput1[1]),2)],',')+"]"
        ;         close,66
        ;      ENDIF
        ;     Writefittingvariables,tirificfirst,painput1,vrotinputPA,INIMODE=inimode


        ;     PAest=PAest+1.
        ;     goto,againPAestimate
        ;  ENDIF
        ENDIF
        testing1PA:
                                ;get the values from these adapted
                                ;models and write them into the
                                ;tirific array
        VariablesWanted=['PA','PA_2','INCL','INCL_2','VROT','VROT_2','SBR','SBR_2']
        firstfitvalues=0.
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
                                ;if the change is large we want to
                                ;adjust our inc and vrot estimates as
                                ;well
        newPA=[firstfitvalues[0,0],ABS(firstfitvalues[0,0]-newPA[0])+newPA[1]]
        newinclination=[firstfitvalues[0,2],catincdev[i]]
        newrot=firstfitvalues[1,4]

                                ;when we get a new inclination we need
                                ;to reset our cutoff values.
        cutoffcorrection=SIN(75.*!DtoR)/SIN((newinclination[0]+newinclination[1]/2.)*!DtoR)
        IF newinclination[0] GT 80. then cutoffcorrection=0.9
        IF newinclination[0]+newinclination[1]/2. GT 50 then cutoffcorrection=1.
        IF newinclination[0]+newinclination[1]/2. LT 50 AND newinclination[0]+newinclination[1]/2. GT 40  then cutoffcorrection=1.+(50-(newinclination[0]+newinclination[1]/2.))*0.05

        IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
        IF doubled and cutoffcorrection GT 1. then cutoffcorrection=SQRT(cutoffcorrection)

        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with"+string(cutoffcorrection)
           close,66
        ENDIF
        cutoff=cutoffor*cutoffcorrection

        IF INCLest GT 0. then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because of low inclination (<30) or small number of beams across minor axis (< 5):"
              printf,66,linenumber()+"We have adjusted the PA from "+string(PAinincl)+' to '+string(newPA[0])+'+/-'+string(newPA[1])
              printf,66,linenumber()+"The inclination from "+string(INCLinincl)+" to "+string(newinclination[0])+'+/-'+string(newinclination[1])
              printf,66,linenumber()+"The max rotation from "+string(maxrotinincl)+" to "+string(newrot[0])+'+/-'+string(catmaxrotdev[i])
              printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT."
              close,66
           ENDIF
        ENDIF ELSE BEGIN
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Because of low inclination (<50) or a mismatch between the KinPA and VelPA estimate :"
              printf,66,linenumber()+"We have adjusted the PA from "+string(PAinincl)+' to '+string(newPA[0])+'+/-'+string(newPA[1])
              printf,66,linenumber()+"This should not change The inclination from "+string(INCLinincl)+" to "+string(newinclination[0])+'+/-'+string(newinclination[1])
              printf,66,linenumber()+"The max rotation from "+string(maxrotinincl)+" to "+string(newrot[0])+'+/-'+string(catmaxrotdev[i])
              printf,66,linenumber()+"We have adjusted the fit setting of the inclination and VROT."
              close,66
           ENDIF

        ENDELSE
        catPAdev[i]=15.
        catinc[i]=newinclination[0]
        catincdev[i]=newinclination[1]
        catPA[i]=newPA[0]
        catmaxrot[i]=newrot[0]
                                ;adapt the inc values
        ;tmppos=where('INCL' EQ tirificfirstvars)
        ;tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
        ;tmppos=where('INCL_2' EQ tirificfirstvars)
        ;tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
        ;INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
        ;            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
        ;            string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']
        INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string(catinc[i]+catincdev[i]+10),string(catinc[i]-catincdev[i]-10),string(0.5),string(0.1),string(0.1),'3']

        IF INCLinput1[1] GT 90. then INCLinput1[1]='90'
        IF INCLinput1[1] LT 30 then INCLinput1[1]='30'

        IF INCLinput1[2] LT 0.1 then INCLinput1[2]='5'
        IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
                                ;and the vrot values

        VROTmax=catmaxrot[i]+catmaxrotdev[i]
        VROTmin=channelwidth
                                ;WE want to ensure a reasonable vrot max
        IF VROTmax LT 80 then VROTmax=80.
        IF VROTmax GT 600. then VROTmax=600.
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
        string2=string(VROTmax)
        string3=string(VROTmin)
        string4=string(2*channelwidth)
        string5=string(0.05*channelwidth)
        ;string6=string(0.2*channelwidth)
        string7=string(0.05*channelwidth)
        string8='5'
        IF centralexclude then begin
           IF catinc[i] LT 40 then begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDIF ELSE BEGIN
           IF catinc[i] LT 40 then begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDELSE

        IF norings[0] GT 4 then begin
           VROTinput1=[string1,string2,string3,string4,string5,string7,string8,string10]
        ENDIF ELSE BEGIN
           VROTinput1=[string1,string2,string3,string4,string5,string7,string8,' ']
        ENDELSE
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+string(catmaxrot[i])
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+string(catmaxrot[i])
     ENDIF
                                ;Since our first SBR estimate is very
                                ;crude we want to adapt that by itself first
     IF counter EQ 0 then begin
        Writefittingvariables,tirificfirst,sbrinput1,sbrinput2,painput1,vrotinput1,INIMODE=inimode
        restart_counter++
        tmppos=where('RESTARTID' EQ tirificfirstvars)
        tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))
        openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
        for index=0,n_elements(tirificfirst)-1 do begin
           printf,1,tirificfirst[index]
        endfor
        close,1

        IF testing GE 1 then goto,testing1SBR
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Starting tirific SBR estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
           close,66
        ENDIF
        print,linenumber()+"Starting tirific SBR estimate in  "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()

        run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name='1stfit.', $
                  store_name='1stfitold.', log=log,loops=loops,nopoints=nopoints,AC=AC1,$
                  toymodels=toymodels,run_unit=run_unit
        IF bookkeeping EQ 5 THEN goto,finishthisgalaxy
        testing1SBR:


        VariablesWanted=['RADI','SBR','SBR_2','PA','PA_2','VROT','VROT_2']
        firstfitvalues=0.
        check_cflux,nopoints,norings[0],tirificfirst,tirificfirstvars,cfluxadjusted,log=log
        IF cfluxadjusted then  VariablesWanted=['RADI','SBR','SBR_2','VROT','VROT_2']
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
        IF not cfluxadjusted then begin
           newPA=[firstfitvalues[0,4],ABS(firstfitvalues[0,4]-newPA[0])+newPA[1]]
        ENDIF ELSE BEGIN
           tmp=firstfitvalues
           firstfitvalues=dblarr(n_elements(tmp[*,0]),n_elements(tmp[0,*])+2)
           firstfitvalues[*,0:2]=tmp[*,0:2]
           firstfitvalues[*,3:4]=newPA[0]
           firstfitvalues[*,5:6]=tmp[*,3:4]
           VariablesWanted=['RADI','SBR','SBR_2','PA','PA_2','VROT','VROT_2']
        ENDELSE



        tmppos=where('SBR' EQ VariablesWanted)
        SBRarr=firstfitvalues[*,tmppos]

        tmppos=where('SBR_2' EQ VariablesWanted)
        SBRarr2=firstfitvalues[*,tmppos]
                                ;We always want to smooth the surface
                                ;brightnes. Added 16-06-2017
        tmppos=where('RADI' EQ VariablesWanted)
        SBRarr=fat_savgol(SBRarr,firstfitvalues[*,tmppos],/HALF)
        SBRarr2=fat_savgol(SBRarr2,firstfitvalues[*,tmppos],/HALF)
                                ;checking the surface brightness
        sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff
        IF newinclination[0] LT 40 then newrot=W50/2./SIN(ABS(newinclination[0]+5)*!pi/180.) else newrot=W50/2./SIN(newinclination[0]*!pi/180.)
        catmaxrotdev[i]=(VSYSdiff/4.)/SIN(newinclination[0]*!pi/180.)
        IF catmaxrotdev[i] LT 3*channelwidth then catmaxrotdev[i]=3*channelwidth
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Because of changing the SBR,"
           printf,66,linenumber()+"We have adjusted the PA from "+string(catPA[i])+" to "+string(firstfitvalues[0,4])
           printf,66,linenumber()+"The inclination from "+string(catInc[i])+" to "+string(newinclination[0])+'+/-'+string(newinclination[1])
           printf,66,linenumber()+"Maxrot from "+string(catmaxrot[i])+" to "+string(newrot[0])+'+/-'+string(catmaxrotdev[i])
           close,66
        ENDIF
                                ;write new values to overall arrays
        catPAdev[i]=15.
        catinc[i]=newinclination[0]
        catincdev[i]=newinclination[1]
        catmaxrot[i]=newrot
        catPA[i]=firstfitvalues[0,4]
                                ;adapt the inc values
        tmppos=where('INCL' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL='+strtrim(strcompress(string(catinc[i])))
        tmppos=where('INCL_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL_2='+strtrim(strcompress(string(catinc[i])))
                                ;Inclination is tricky and the settings
                                ;depend greatly on whether high or low inclination
        case 1 of
           catinc[i] GE 75:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]),string(catinc[i]-catincdev[i]),string(0.1),string(0.01),string(0.5),string(0.01),'3','70','70']
           END
           catinc[i] GE 40:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']
           END
           catinc[i] LT 30:begin
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.05),string(3),string(0.05),'5','70','70']
           end
           catinc[i] LT 40:begin
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.05),string(3),string(0.05),'5','70','70']
           end
           ELSE:BEGIN
              INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                          string(catinc[i]+catincdev[i]+20),string(catinc[i]-catincdev[i]-20),string(0.2),string(0.1),string(1.),string(0.1),'3','70','70']
           END
        endcase
        IF INCLinput1[1] GT 90. then INCLinput1[1]='90'
        IF INCLinput1[1] LT 30 then INCLinput1[1]='30'
        IF catinc[i] LT 75 then begin
           IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
           IF INCLinput1[2] GT 60 then INCLinput1[2]='60'
        endif
                                ;and the vrot values
        VROTmax=catmaxrot[i]+catmaxrotdev[i]
        VROTmin=channelwidth
                                ;We want to ensure that we have a
                                ;reasonable Vrot max
        IF VROTmax LT 80. then VROTmax=80.
        IF VROTmax GT 600. then VROTmax=600.
        string1='!VROT '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2 VROT_2 '+strtrim(strcompress(string(norings[0],format='(I3)')),1)+':2'
        string2=string(VROTmax)
        string3=string(VROTmin)

        IF catinc[i] LT 40 then begin
           string4=string(channelwidth)
           string5=string(0.1*channelwidth)
           string6=string(channelwidth)
           string7=string(0.1*channelwidth)
        ENDIF ELSE BEGIN
           string4=string(2*channelwidth)
           string5=string(0.01*channelwidth)
           string6=string(0.2*channelwidth)
           string7=string(0.01*channelwidth)
        ENDELSE
        string8='3'
        string9='70'
        IF centralexclude then begin
           IF catinc[i] LT 40 then begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDIF ELSE BEGIN
           IF catinc[i] LT 40 then begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*0.5+2),format='(I3)')),1)
           ENDIF else begin
              string10='VROT '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)+' VROT_2 '+strtrim(strcompress(string(norings[0]-1,format='(I3)')),1)+':'+strtrim(strcompress(string(fix(norings[0]*3./4.),format='(I3)')),1)
           ENDELSE
        ENDELSE

        IF norings[0] GT 4 then begin
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,string10]
        ENDIF ELSE BEGIN
           VROTinput1=[string1,string2,string3,string4,string5,string6,string7,string8,string9,string9,' ']
        ENDELSE
        tmppos=where('VROT' EQ VariablesWanted)
        VROTarr=firstfitvalues[*,tmppos]
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
     ENDIF
                                ;Setting some SBR values for the fitting
     minstep=strtrim(strcompress(string(cutoff[norings[0]-1]/3.,format='(E8.1)')),1)
     IF countsbr EQ 0. then begin
        startstep='1E-5'
        minstep='1E-6'
     ENDIF else begin
        startstep=strtrim(strcompress(string(5.*cutoff[norings[0]-1],format='(E8.1)')),1)
        minstep=strtrim(strcompress(string(cutoff[norings[0]-1]/2.,format='(E8.1)')),1)
     ENDELSE

     satlevel=strtrim(strcompress(string(5.*double(startstep)/10.,format='(E8.1)')),1)

     string1=''
     string2=''
     string3=''
     string4=''
     string5=''
     string6=''
     string7=''

                                ;The minimum is based on the cutoff
                                ;values and therefore we need to set
                                ;every ring in the tirific file

     for j=norings[0],3,-1 do begin
        string1=string1+','+'SBR '+strtrim(strcompress(string(j,format='(I3)')),1)+', SBR_2 '+strtrim(strcompress(string(j,format='(I3)')),1)
        string2=string2+' 1 1'
        IF doubled then $
           string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/4.,format='(E12.5)')),1)+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/4.,format='(E12.5)')),1) $
           else string3=string3+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)+' '+strtrim(strcompress(string(cutoff[fix(j-1)]/2.,format='(E12.5)')),1)
        string4=string4+' '+startstep+' '+startstep
        string5=string5+' '+minstep+' '+minstep
        string6=string6+' '+strtrim(strcompress(string(satlevel,format='(E8.1)')),1)+' '+strtrim(strcompress(string(satlevel,format='(E8.1)')),1)
        string7=string7+' 3 3'
     endfor
     ;ENDELSE
     string1 = STRMID(string1,1,STRLEN(string1)-1)

     SBRinput1=[string1,string2,string3,string4,string5,string6,string7]
                                ;The inner rings should be fitted as one
     SBRinput2=[' SBR 1 2  SBR_2 1 2 ',$
                strtrim(strcompress(string(ABS(SBRarr[2])*2.,format='(E12.5)'))) ,'0','1E-5','1E-6','5E-6','3']
     smoothrotation=0.
     keepcenter=0.
     fixedcenter=0.
     shiftcenter:
     notacceptedone:
     IF smoothrotation then vrotinput1=vrotinputoriginal
     smoothrotation=0.
     smoothingone:
                                ;write the paramter to the array based
                                ;on inclination
     if not keepcenter and not fixedcenter and not cfluxadjusted then begin
        case 1 of
           catinc[i] LT 30: Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,INIMODE=inimode
           catinc[i] LT 50:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,inclinput1,sbrinput1,sbrinput2,INIMODE=inimode
           catinc[i] GT 75:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,z0input1,INIMODE=inimode
           else:Writefittingvariables,tirificfirst,xposinput1,yposinput1,vsysinput1,painput1,inclinput1,vrotinput1,sbrinput1,sbrinput2,INIMODE=inimode
        endcase
    ENDIF ELSE BEGIN
        case 1 of
           catinc[i] LT 30:Writefittingvariables,tirificfirst,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,INIMODE=inimode
           catinc[i] LT 50:Writefittingvariables,tirificfirst,painput1,vrotinput1,inclinput1,sbrinput1,sbrinput2,INIMODE=inimode
           catinc[i] GT 75:Writefittingvariables,tirificfirst,painput1,vrotinput1,sbrinput1,sbrinput2,inclinput1,z0input1,INIMODE=inimode
           else:Writefittingvariables,tirificfirst,painput1,inclinput1,vrotinput1,sbrinput1,sbrinput2,INIMODE=inimode
        endcase
     ENDELSE
                                ;If we tried fifty times and have not
                                ;found a reasonable solution then we
                                ;give up
     IF counter GT 30 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have rerun "+string(counter)+" times and still haven't found an acceptable solution."
           printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
           printf,66,linenumber()+"We are aborting this galaxy."
           close,66
        ENDIF
        comment = 'We could not find a proper center.'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        bookkeeping=5
        goto,finishthisgalaxy
     ENDIF
     ;Write the current array to the def file
     restart_counter++
     tmppos=where('RESTARTID' EQ tirificfirstvars)
     tirificfirst[tmppos]='RESTARTID= '+string(fix(restart_counter))
     print,linenumber()+'We used '+tirificfirst[tmppos]+' to restart'
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificfirst)-1 do begin
        printf,1,tirificfirst[index]
     endfor
     close,1
                                ;And starting the first fit while
                                ;moving any previous fits to old.
                                ;if this is the second fit we do want
                                ;to keep the first one not being old
     IF testing GE 1 then goto,testing1
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific first estimate in  "+catDirname[i]+"  which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
        printf,66,linenumber()+"We have rerun "+string(counter)+" times."
        printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
        close,66
     ENDIF
     print,linenumber()+"Starting tirific first estimate in  "+catDirname[i]+"  which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
     IF tryone EQ 1 then begin
        output_name = '1stfit.'
        store_name = '1stfitall.'
     ENDIF ELSE BEGIN
        output_name = '1stfit.'
        store_name = '1stfitold.'
     ENDELSE
     counter++
     run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name=output_name, $
                  store_name=store_name, log=log, loops=loops,nopoints=nopoints,AC=AC1, $
                  toymodels=toymodels,run_unit=run_unit

     testing1:

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC1 EQ 1 then begin
           printf,66,linenumber()+"The first estimate was accepted in this run."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The First estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF
     IF AC1 EQ 1 then begin
        print,linenumber()+"The first estimate is accepted."
     ENDIF ELSE begin
        print,linenumber()+"The first estimate is not accepted."
     ENDELSE
                                ;Check the number of points in the models
     check_cflux,nopoints,norings[0],tirificfirst,tirificfirstvars,cfluxadjusted,log=log
     IF  cfluxadjusted then goto,notacceptedone

     IF testing GE 1 then begin
        maxchangeRA=ABS(0.15*catmajbeam[i])/3600.
        maxchangeDEC=ABS(0.15*catmajbeam[i])/3600.
        IF maxchangeRA LT ABS(0.5*pixelsizeRA) then maxchangeRA=ABS(pixelsizeRA)
        IF maxchangeDEC LT ABS(0.5*pixelsizeDEC) then maxchangeDEC=ABS(pixelsizeDEC)
        maxchangevel=ABS(0.5*channelwidth)
        IF maxchangeRA LT  1./3600. then maxchangeRA=1./3600.
        IF maxchangeDEC LT  1./3600. then maxchangeDEC=1./3600.
        IF maxchangevel LT 2.5 then  maxchangevel=2.5
        goto,testing1skip
     ENDIF
                               ;Let's extract some fitting parameters from the file we just fitted
     VariablesWanted=['INCL','PA','SBR','XPOS','YPOS','VSYS','CFLUX','VROT','SBR_2','RADI','Z0']
     firstfitvalues=0.
     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted,/EXTRACT
     tmppos=where('SBR' EQ VariablesWanted)
     SBRarr=firstfitvalues[*,tmppos]
     tmppos=where('SBR_2' EQ VariablesWanted)
     SBRarr2=firstfitvalues[*,tmppos]
                                ;We always want to smooth the surface
                                ;brightnes. Added 16-06-2017
     tmppos=where('RADI' EQ VariablesWanted)
     ;To judge whether we want extend we want to use the non-smoothed profiles
     SBRarrunmod=SBRarr
     SBRarr2unmod=SBRarr2
     SBRarr=fat_savgol(SBRarr,firstfitvalues[*,tmppos],/half)
     SBRarr2=fat_savgol(SBRarr2,firstfitvalues[*,tmppos],/half)

     tmppos=where('VROT' EQ VariablesWanted)
     VROTarr=firstfitvalues[*,tmppos]
     stringvelocities='VROT= 0. '+STRJOIN(VROTarr[1:n_elements(VROTarr)-1],' ')
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"SBRs before they are checked."
        printf,66,SBRarr,SBRarr2
        close,66
     endif
     sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"SBRs after they are checked."
        printf,66,SBRarr,SBRarr2
        close,66
     endif
                                ; Let's see if the values are too close to the boundaries or not
                                ;First the inclination
     diffmax=ABS(double(inclinput1[1])-firstfitvalues[0,0])
     diffmin=ABS(firstfitvalues[0,0]-double(inclinput1[2]))
     IF inclinput1[2] EQ 5. then diffmin = 0.2*catincdev[i]
     IF inclinput1[1] EQ 90 then diffmax = 0.2*catincdev[i]
     oldinc=catinc[i]
     oldpa=catpa[i]
     paraised=0.
     inclraised=0
     IF (diffmax GT 0.1*catincdev[i] and diffmin GT 0.1*catincdev[i]) then begin
                                ; if the  value is  more than 10% of
                                ; the error removed from the boundary
                                ; values it is ok and we update the
                                ; inclination or we have high
                                ; inclination
        catinc[i]=ABS(firstfitvalues[0,0])
        inclraised=1
        IF catinc[i] GT 90. then catinc[i]=90.-(catinc[i]-90.)
        tmppos=where('INCL' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL= '+strtrim(string(catinc[i]))
        tmppos=where('INCL_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='INCL_2= '+strtrim(string(catinc[i]))
        catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.)
        catmaxrotdev[i]=(VSYSdiff/4.)/SIN(catinc[i]*!pi/180.)
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]=stringvelocities
        tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2='+tmp[1]
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The incl was adjusted to"+string(catinc[i])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The incl was adjusted to"+string(catinc[i])
        ENDELSE
     ENDIF ELSE Begin
        ; Otherwise we will adjust the boundaries if they are insecure
        IF ceil(norings[0]*ring_spacing*COS(firstfitvalues[0,0]*!DtoR)) GT 4 OR (norings[0]*ring_spacing GT 6. AND catinc[i] GT 80.) then begin
                                ;If we are on the boundary we reject
                                ;the new value unless we are in very
                                ;uncertain conditions
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The incl was on the boundary"+string(firstfitvalues[0,0])+'+'+string(inclinput1[1])+'-'+string(inclinput1[2])
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The incl was on the boundary"+string(firstfitvalues[0,0])+'+'+string(inclinput1[1])+'-'+string(inclinput1[2])
           ENDELSE
                                ;However if this is due to high
                                ;inclination we still want to reset vrot
           IF  (norings[0]*ring_spacing GT 6. AND catinc[i] GT 80.) then begin
              catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.)
              catmaxrotdev[i]=(VSYSdiff/4.)/SIN(catinc[i]*!pi/180.)
              tmppos=where('VROT' EQ tirificfirstvars)
              tirificfirst[tmppos]=stringvelocities
              tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
              tmppos=where('VROT_2' EQ tirificfirstvars)
              tirificfirst[tmppos]='VROT_2='+tmp[1]
          ENDIF ELSE AC1=0
        ENDIF ELSE BEGIN
                                ; Otherwise we will adjust the boundaries if they are insecure
           inclraised=1
           catinc[i]=ABS(firstfitvalues[0,0])
           AC1=0
           IF catinc[i] GT 90. then catinc[i]=90.-(catinc[i]-90.)
           tmppos=where('INCL' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL= '+strtrim(string(catinc[i]))
           tmppos=where('INCL_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='INCL_2= '+strtrim(string(catinc[i]))
           catmaxrot[i]=W50/2./SIN(catinc[i]*!pi/180.)
           catmaxrotdev[i]=(VSYSdiff/4.)/SIN(catinc[i]*!pi/180.)
           tmppos=where('VROT' EQ tirificfirstvars)
           tirificfirst[tmppos]=stringvelocities
           tmp=str_sep(strtrim(strcompress(stringvelocities),2),'=')
           tmppos=where('VROT_2' EQ tirificfirstvars)
           tirificfirst[tmppos]='VROT_2='+tmp[1]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Although the incl hit its boundary we adjusted to"+string(catinc[i])+"due to the small ring size. "
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"Although the incl hit its boundary we adjusted to"+string(catinc[i])+"due to the small ring size. "
           ENDELSE
           VROTmax=catmaxrot[i]+catmaxrotdev[i]
           INCLmin=catinc[i]-catincdev[i]-20
                                ;We base the minimum vrot max on
                                ;inclination however the maximum of
                                ;600 is a physical limit.
           IF VROTmax LT 80 then VROTmax=80.
           IF VROTmax GT 600. then VROTmax=600.
           VROTinput1[1]=strtrim(strcompress(string(VROTmax)),2)
           INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(90.),string(INCLmin),string(1.0),string(0.1),string(1.),string(0.1),'5','70','70']
        ENDELSE
     ENDELSE
                                ;when we get a new inclination we need
                                ;to reset our cutoff values.
     cutoffcorrection=SIN(75.*!DtoR)/SIN((catinc[i]+newinclination[1]/2.)*!DtoR)
     IF newinclination[0] GT 80. then cutoffcorrection=0.9
     IF catinc[i]+newinclination[1]/2. GT 50 then cutoffcorrection=1.
     IF catinc[i]+newinclination[1]/2. LT 50 AND catinc[i]+newinclination[1]/2. GT 40  then cutoffcorrection=1.+(50-(catinc[i]+newinclination[1]/2.))*0.05

     IF cutoffcorrection GT 2.5 then cutoffcorrection=2.5
     IF doubled and cutoffcorrection GT 1. then cutoffcorrection=SQRT(cutoffcorrection)

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Cutoff values will be adjusted for inclination by multiplying them with "+string(cutoffcorrection)
        close,66
     ENDIF
     cutoff=cutoffor*cutoffcorrection
                                ;Then we look at the PA
     diffmax=ABS(firstfitvalues[0,1]-double(painput1[1]))
     diffmin=ABS(firstfitvalues[0,1]-double(painput1[2]))
     IF diffmax GT 0.1*catpadev[i] and diffmin GT 0.1*catpadev[i] then begin
        catpa[i]=firstfitvalues[0,1]
        paraised=1
        ;if not between -360 and 360 adapt
        WHILE ABS(catpa[i])/360. GT 1. DO BEGIN
           IF catpa[i] GT 360. then catpa[i]=catpa[i]-360.
           IF catpa[i] LT -360. then catpa[i]=catpa[i]+360.
        ENDWHILE
        IF catpa[i] LT 0. then catpa[i]=catpa[i]+360.
        tmppos=where('PA' EQ tirificfirstvars)
        tirificfirst[tmppos]='PA= '+strtrim(string(catpa[i]))
        tmppos=where('PA_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='PA_2= '+strtrim(string(catpa[i]))
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The pa was adjusted to "+string(catpa[i])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The pa was adjusted to "+string(catpa[i])
        ENDELSE
     ENDIF ELSE begin
        IF shiftcentercounter EQ 1. then begin
           IF UndatedPA LT 180 then begin
              firstfitvalues[*,1]=firstfitvalues[*,1]+180.
              catpa[i]=catpa[i]+180.
              UndatedPA=UndatedPA+180.
           endif else begin
              firstfitvalues[*,1]=firstfitvalues[*,1]-180.
              catpa[i]=catpa[i]-180.
              UndatedPA=UndatedPA-180.
           endelse

        ENDIF ELSE BEGIN
           IF oldinc LT 50 AND UndatedINCL LT 65 AND catinc[i] LT 55 then begin

              paraised=1
              catpa[i]=firstfitvalues[0,1]
              AC1=0
              WHILE ABS(catpa[i])/360. GT 1. DO BEGIN
                 IF catpa[i] GT 360. then catpa[i]=catpa[i]-360.
                 IF catpa[i] LT -360. then catpa[i]=catpa[i]+360.
              ENDWHILE
              IF catpa[i] LT 0. then catpa[i]=catpa[i]+360.
              tmppos=where('PA' EQ tirificfirstvars)
              tirificfirst[tmppos]='PA= '+strtrim(string(catpa[i]))
              tmppos=where('PA_2' EQ tirificfirstvars)
              tirificfirst[tmppos]='PA_2= '+strtrim(string(catpa[i]))
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Although the pa hit its boundary we adjusted to"+string(catpa[i])+"because of low inclination."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"Although the pa hit its boundary we adjusted to"+string(catpa[i])+"because of low inclination."
              ENDELSE
              IF PAest LT 5 then begin
                 PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                           string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(5),string(0.1),string(0.05),string(0.01),'3','70','70']
                 IF PAfixboun NE 'n' then begin
                    IF PAfixboun EQ 'w' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                       IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
                    ENDIF
                    IF PAfixboun EQ 'e' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                       IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
                    ENDIF
                 ENDIF
                 Writefittingvariables,tirificfirst,painput1,INIMODE=inimode
                 PAest++
                 goto,againPAestimate
              ENDIF else begin
                 PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                           string((catPA[i])+catPAdev[i]+30.),string((catPA[i])-catPAdev[i]-30),string(2.5),string(0.1),string(0.05),string(0.01),'3','70','70']
                 IF PAfixboun NE 'n' then begin
                    IF PAfixboun EQ 'w' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 360 then PAinput1[1]='360'
                       IF catPA[i]-catPAdev[i]-30. LT 180 then PAinput1[2]='180'
                    ENDIF
                    IF PAfixboun EQ 'e' then begin
                       IF catPA[i]+catPAdev[i]+30. GT 180. then PAinput1[1]='180'
                       IF catPA[i]-catPAdev[i]-30. LT 0. then PAinput1[2]='0'
                    ENDIF
                 ENDIF
              ENDELSE
           ENDIF
        ENDELSE
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The pa was on the boundary"+string(firstfitvalues[0,1])+'+'+string(painput1[1])+'-'+string(painput1[2])
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The pa was on the boundary"+string(firstfitvalues[0,1])+'+'+string(painput1[1])+'-'+string(painput1[2])
        ENDELSE
        AC1=0
     ENDELSE

     tmppos=where('Z0' EQ tirificfirstvars)
     tirificfirst[tmppos]='Z0= '+strtrim(string(firstfitvalues[0,10]))
     tmppos=where('Z0_2' EQ tirificfirstvars)
     tirificfirst[tmppos]='Z0_2= '+strtrim(string(firstfitvalues[0,10]))
                                ; let's check wether the center
                                ; shifted a lot since this affects everything
                                ;if there is a major change in central
                                ;position as it determines many parameters

     newxpos=firstfitvalues[0,3]
     newypos=firstfitvalues[0,4]
     newvsys=firstfitvalues[0,5]
     IF doubled then beamfrac=0.3 else beamfrac=0.15
     maxchangeRA=ABS(beamfrac*catmajbeam[i])/3600.
     maxchangeDEC=ABS(beamfrac*catmajbeam[i])/3600.
     IF maxchangeRA LT ABS(0.5*pixelsizeRA) then maxchangeRA=ABS(pixelsizeRA)
     IF maxchangeDEC LT ABS(0.5*pixelsizeDEC) then maxchangeDEC=ABS(pixelsizeDEC)
     maxchangevel=ABS(0.5*channelwidth)
     IF maxchangeRA LT  1./3600. then maxchangeRA=1./3600.
     IF maxchangeDEC LT  1./3600. then maxchangeDEC=1./3600.
     IF maxchangevel LT 2.5 then  maxchangevel=2.5
                                ;if the change between this fit and
                                ;the previous one is bigger that the
                                ;limits, refit
     IF ABS(RADeg-newxpos) GT maxchangeRA OR ABS(DECDeg-newypos) GT maxchangeDEC OR ABS(catvsys[i]-newvsys) GT maxchangevel OR fixedcenter EQ 1. then begin
        ; If the center was fixed refit
        IF fixedcenter EQ 1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center was not fitted."
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center was not fitted."
           ENDELSE
           fixedcenter=0
           goto,shiftcenter
        ENDIF ELSE BEGIN
                                ;if the shift is more than two beams
                                ;from the initial guess something went
                                ;wrong
           IF norings*ring_spacing LE 25 then resetlimit=catmajbeam[i]/1800. else begin
              resetlimit=catmajbeam[i]/3600.*norings[0]*ring_spacing*0.08
           ENDELSE
           IF    ABS(RADeg-newxpos) GT resetlimit OR ABS(DECDeg-newypos) GT resetlimit then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center shifted more than "+string(resetlimit*3600/catmajbeam[i])+" major beams. Not applying this shift."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"The center shifted more than "+string(resetlimit*3600/catmajbeam[i])+" major beams. Not applying this shift."
              ENDELSE
              IF paraised then begin
                 tmppos=where('PA' EQ tirificfirstvars)
                 tirificfirst[tmppos]='PA= '+strtrim(string(oldpa))
                 tmppos=where('PA_2' EQ tirificfirstvars)
                 tirificfirst[tmppos]='PA_2= '+strtrim(string(oldpa))
              ENDIF
              IF inclraised then begin
                 tmppos=where('INCL' EQ tirificfirstvars)
                 tirificfirst[tmppos]='INCL= '+strtrim(string(oldinc))
                 tmppos=where('INCL_2' EQ tirificfirstvars)
                 tirificfirst[tmppos]='INCL_2= '+strtrim(string(oldinc))
              ENDIF
              fixedcenter=1
              tryone=0.
              plus2beamshift++
                                ; if the center has tried to shift by 2 beams
                                ; more than 10 times we are just going
                                ; to keep this center
              IF plus2beamshift GT 10 then begin
                 keepcenter=1.
                 fixedcenter=0.
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"The center keeps trying to change unreasonably. Accepting the current center and fixing it."
                    close,66
                 ENDIF ELSE BEGIN
                    print,linenumber()+"The center keeps trying to change unreasonably. Accepting the current center and fixing it."
                 ENDELSE
              ENDIF
              goto,shiftcenter
           ENDIF
                                ;If the shift in center is jumping
                                ;back and forth between two poistions
                                ;accept the first/last position
           IF ABS(oldRA-newxpos) LT maxchangeRA AND ABS(oldDEC-newypos) LT maxchangeDEC AND ABS(oldVSYS-newvsys) LT maxchangevel then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center shifted back to the previous fit. Accepting this center and fixing it."
                 close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"The center shifted back to the previous fit. Accepting this center and fixing it."
              ENDELSE
              keepcenter=1.
              goto,shiftback
           ENDIF
                                ;If the shift was to large try again
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted too much trying again with new center."
              printf,66,linenumber()+"The RA has shifted from "+string(RAdeg)+" to "+string(newxpos)+" which is a difference of"$
                     +string(ABS(RADeg-newxpos))+" needed ="+string(maxchangeRA)
              printf,66,linenumber()+"The DEC has shifted from "+string(DECdeg)+" to "+string(newypos)+" which is a difference of"$
                     +string(ABS(DECDeg-newypos))+" needed ="+string(maxchangeDEC)
              printf,66,linenumber()+"The systemic has shifted from "+string(catvsys[i])+" to "+string(newvsys)+" which is a difference of"$
                     +string(ABS(catvsys[i]-newvsys))+" needed ="+string(maxchangevel)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center shifted too much trying again with new center."
           ENDELSE
        ENDELSE
        oldRA=RADeg
        oldDEC=DECdeg
        oldVSYS=catvsys[i]
        RADeg=newxpos
        DECdeg=newypos
        catvsys[i]=newvsys
        VariablesWanted=['XPOS','XPOS_2','YPOS','YPOS_2','VSYS','VSYS_2']
        firstfitvalues=0.
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=VariablesWanted
        tryone=0.
        goto,shiftcenter
     ENDIF  ELSE BEGIN
        IF newxpos GT RAboundeg[0] AND newxpos LT RAboundeg[1] AND newypos GT DECboundeg[0] AND newypos LT DECboundeg[1] AND $
           newvsys GT ROTboun[0] AND newvsys LT ROTboun[1] AND testing LT 1 then begin
           shiftback:
           RADeg=newxpos
           DECdeg=newypos
           catvsys[i]=newvsys
           if fixedcenter EQ 1. then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The center was accepted because it was fixed."
                 close,66
              ENDIF
              goto,shiftcenter
           endif


                                ;let's see if the model has the right size
           IF secondtime OR doubled THEN norings=newrings ELSE get_newringsv9,SBRarrunmod,SBRarr2unmod,cutoff,newrings
           ;Let' see what happens when we fix the rings for the first step, That is a bad idea as the last rings throw things off

                                ;cannot have newsize smaller than 4
           if newrings LT 3 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to 3."
                 close,66
              ENDIF
              newrings=3
           ENDIF
                                ;see if the newsize is not too big
           IF newrings GT maxrings then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to maxrings."
                 close,66
              ENDIF
              newrings=maxrings
           ENDIF
                                ;If we have a small size we do not
                                ;vary by more that a single ring
           IF norings[0] LE 8 OR newrings LE 8 then begin
              IF newrings LT norings[0]-1 then newrings=norings[0]-1
              IF newrings GT norings[0]+1 then newrings=norings[0]+1
           ENDIF



           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We found this as rings="+strtrim(string(fix(norings[0])),2)+"  new="+strtrim(string(fix(newrings)),2)
              close,66
           ENDIF
                                ;  Let's see whether the
                                ;  improved version can allow
                                ;  large estimates of the initial
                                ;  version changed -2 to -3 (15-05-2015)
           IF newrings GT sofiarings+(2/ring_spacing) OR newrings LT sofiarings-(3./ring_spacing) then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The new amount of rings ("+strtrim(string(fix(newrings)),2)+") deviates too much from the sofia estimate ("+string(sofiarings)+")."
                 close,66
              ENDIF
              IF newrings LT norings[0] then begin
                 IF norings[0] GT sofiarings-round(2/ring_spacing) then newrings=norings[0]-1 else newrings=sofiarings-round(2/ring_spacing)
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"We set the number of rings to "+strtrim(string(fix(newrings)),2)
                    close,66
                 ENDIF
              endif Else begin
                 IF norings[0] LT sofiarings+round(2/ring_spacing) then newrings=norings[0]+1 else newrings=sofiarings+round(2/ring_spacing)
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"We set the number of rings to "+strtrim(string(fix(newrings)),2)
                    close,66
                 ENDIF
              ENDELSE

           ENDIF
           ;ENDIF
           IF newrings LT norings[0] AND constring NE 1 AND smoothrotation EQ 0. then begin
              IF newrings LE forcedring then begin
                 IF norings[0] GT forcedring then begin
                    newrings=forcedring
                    constring=1
                    prevmodification=0.
                 ENDIF ELSE BEGIN
                    IF size(log,/TYPE) EQ 7 then begin
                       openu,66,log,/APPEND
                       printf,66,linenumber()+"As the cut would modify the rings to less than the maximum forced addition of a ring we do not apply the cut."
                       printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)+"forcedring="+string(forcedring)
                       Close,66
                    ENDIF
                    goto,nocutfirst
                 ENDELSE
              ENDIF




                                ;check that we do not go back and forth
              IF prevmodification EQ 1 then begin
                 IF newrings GT norings[0]-1 then begin
                    prevmodification=-1
                 ENDIF else BEGIN
                    IF newrings EQ norings[0]-1 then begin
                       IF size(log,/TYPE) EQ 7 then begin
                          openu,66,log,/APPEND
                          printf,66,linenumber()+"As the last modification was the addition of a ring we will fix this ring number."
                          printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)
                          Close,66
                       ENDIF
                       constring=1
                       secondtime=1
                       newrings=newrings+1
                    ENDIF
                    prevmodification=-2
                ENDELSE
              ENDIF
              oldrings=norings[0]
                                ;If no change then go on with the fitting process
              IF newrings EQ norings[0] THEN goto,nocutfirst ELSE norings[0]=newrings

                                ;check if we fitted this size before
              tmp=WHERE(prevrings EQ newrings)
              IF tmp[0] NE -1 AND newrings NE oldrings then begin
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                    close,66
                 ENDIF
                 secondtime=1
              ENDIF
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
                 Close,66
              ENDIF ELSE begin
                 print,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
              ENDELSE
              lastcutrings=norings[0]
              changeradii,tirificfirst,norings[0],continue_tirific
              fluxadjust=0.
              sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff
              countsbr++
              overwrite=0.
              prevmodification=-1
              ringmodifier=ringmodifier+1
              prevrings=[prevrings,norings[0]]
              ;refit
              goto,sbrshift
           ENDIF
           nocutfirst:
           overwrite=0
           ;if the last rings are really bright force add an addition
           IF (SBRarr[newrings-2] GT 7.5*cutoff[newrings-2] OR SBRarr2[newrings-2] GT 7.5*cutoff[newrings-2]) AND newrings GT norings[0] AND prevmodification EQ -1 then begin
              prevmodification=0.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We force added a ring to the profile."
                 printf,66,linenumber()+"Newrings"+strtrim(string(fix(newrings)),2)+", Rings"+strtrim(string(fix(norings[0])),2)
                 Close,66
              ENDIF
              IF newrings GT forcedring then forcedring=newrings
              overwrite=1
           ENDIF
                                ;If we previously subtracted we will not add
           IF prevmodification EQ -1 AND overwrite NE 1 AND newrings GT norings[0] then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We wanted to add a ring to the model but the previous step was a subtraction."
                 Close,66
              ENDIF
           ENDIF
                                ;If we want to add a ring we need to
                                ;set some parameters and make estimate
           IF newrings GT norings[0] AND (prevmodification NE -1 OR overwrite EQ 1) AND secondtime NE 1 AND smoothrotation EQ 0 then begin
              tmp=WHERE(prevrings EQ newrings)
              IF tmp[0] NE -1 AND newrings NE norings[0] then begin
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                    close,66
                 ENDIF
                 secondtime=1
              ENDIF
              norings[0]=newrings
              lastaddrings=norings[0]
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We added a ring to the model. "+string(secondtime)
                 Close,66
              ENDIF ELSE BEGIN
                 print,linenumber()+"We added a ring to the model."
              ENDELSE
              changeradii,tirificfirst,norings[0],continue_tirific
              ;if norings[0]*ring_spacing LT 6. then catinc[i]=catinc[i]+2.
              prevmodification=1.
              ringmodifier=ringmodifier+1
              countsbr++
              prevrings=[prevrings,norings[0]]
              goto,sbrshift
           ENDIF
        ENDIF ELSE BEGIN
           ; If we really go crazy then just end the fitting process
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The fit diverged out of the boundaries set by the parameterization. That can't be good."
              printf,66,linenumber()+"Boundaries are; DEC current="+string(newypos)+"DEC Boundaries"+string(DECboundeg[0])+","+string(DECboundeg[1])
              printf,66,linenumber()+"RA current="+string(newxpos)+"RA Boundaries"+string(RAboundeg[0])+","+string(RAboundeg[1])
              printf,66,linenumber()+"V_sys current="+string(newvsys)+"V_sys Boundaries"+string(ROTboun[0])+","+string(ROTboun[1])
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The fit diverged out of the boundaries set by the parameterization. That can't be good."
           ENDELSE
           comment = 'This galaxy diverged out of the set boundaries.'
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(0),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
           bookkeeping=5
           goto,finishthisgalaxy
        ENDELSE
     ENDELSE
     IF AC1 NE 1 And tryone LT 0. then begin
                                ;If not accepted let's try again


        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The fit is not accepted trying once more."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"The fit is not accepted trying once more."
        ENDELSE
                                ;but from the newvalues
        writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def'
        tryone=tryone+1.
        goto,notacceptedone
     endif
                                ;let's do a final check with a
                                ;smoothed rotation curve when all is well
     IF NOT smoothrotation then begin
        smoothrotation=1
        tmppos=where('VROT' EQ VariablesWanted)
        VROTarr=firstfitvalues[*,tmppos[0]]
        sbrpos=where('SBR' EQ VariablesWanted)
        sbrpos2=where('SBR_2' EQ VariablesWanted)
        SBRarrcom=(firstfitvalues[*,sbrpos[0]]+firstfitvalues[*,sbrpos2])/2.
          ;We always want to smooth the surface brightnes. Added 16-06-2017
        SBRarrcom=fat_savgol(SBRarrcom,firstfitvalues[*,9])
        VROTarr[0]=0.
        vmaxdev=MAX([30,7.5*channelwidth*(1.+vresolution)])
        verror=MAX([5.,channelwidth/2.*(1.+vresolution)/SQRT(sin(catinc[i]*!DtoR))])
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We shall now smooth the rotation curve."
           close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"We shall now smooth the rotation curve."
        ENDELSE
        sigmarot=verror
        IF VROTarr[1] LT 120. AND VROTarr[2] LT 120. then begin
           revised_regularisation_rot,VROTarr,SBRarrcom, firstfitvalues[*,9],/REVERSE,fixedrings=norings[0]-fix(norings[0]*3./4.),difference=verror,cutoff=cutoff,arctan=0.,order=polorder,max_par=VROTmax,min_par=channelwidth,log=log ,error=sigmarot,ring_spacing=ring_spacing
        ENDIF ELSE Begin
           revised_regularisation_rot,VROTarr,SBRarrcom, firstfitvalues[*,9],/REVERSE,fixedrings=norings[0]-fix(norings[0]*3./4.),difference=verror,cutoff=cutoff,arctan=0.,order=polorder,max_par=VROTmax,min_par=channelwidth,/NOCENTRAL,log=log ,error=sigmarot,ring_spacing=ring_spacing
        ENDELSE
        tmp0check=WHERE(VROTarr LT 0)
        IF tmp0check[0] NE -1 then VROTarr[tmp0check]=3.*channelwidth
        tmppos=where('VROT' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT= 0. '+STRJOIN(strtrim(string(VROTarr[1:n_elements(VROTarr[*])-1]),2),' ')
        tmppos=where('VROT_2' EQ tirificfirstvars)
        tirificfirst[tmppos]='VROT_2= 0. '+STRJOIN(strtrim(string(VROTarr[1:n_elements(VROTarr[*])-1]),2),' ')
                                ;We also want a full smoothing on the
                                ;profile
        tmppos=where('RADI' EQ VariablesWanted)
        sbrarr=fat_savgol(firstfitvalues[*,sbrpos],firstfitvalues[*,tmppos])
        sbrarr2=fat_savgol(firstfitvalues[*,sbrpos2],firstfitvalues[*,tmppos])

        sbr_check,tirificfirst, tirificfirstvars,sbrarr,sbrarr2,cutoff


        vrotinputoriginal=vrotinput1
        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    string(VROTmax),string(VROTmin),string(channelwidth),string(0.1*channelwidth),string(channelwidth),string(0.1*channelwidth),'3','70','70', ' ']
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have smoothed the rotation curve in the first fit."
           close,66
           print,linenumber()+"We have smoothed the rotation curve in the first fit."
        ENDIF ELSE BEGIN
           print,linenumber()+"We have smoothed the rotation curve in the first fit."
        ENDELSE

        goto,smoothingone
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We finished the first fit of "+catDirname[i]+" at "+systime()
        printf,66,linenumber()+"We have rerun "+string(counter)+" times."
        printf,66,linenumber()+"We have modified the rings "+string(countsbr)+" times."
        IF AC1 EQ 0 then  printf,66,linenumber()+"The fit was not accepted."
        IF AC1 EQ 1 then  printf,66,linenumber()+"The fit was accepted."
        close,66
     ENDIF
     testing1skip:
                                ;Reading out the fit values
     Basicinfovars=['XPOS','YPOS','VSYS','PA','INCL','VROT']
     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,Variables=tirificfirstvars,/EXTRACT
     if testing EQ 1 then begin
        newxpos=Basicinfovalues[0,0]
        newypos=Basicinfovalues[0,1]
     ENDIF
     RAhr=Basicinfovalues[0,0]
     RAdiff=maxchangeRA*3600./15.
     DEChr=Basicinfovalues[0,1]
     DECdiff=maxchangeDEC*3600.

                               ;If the cube is optimized then make a model at full resolution else make sure to kill tirific
     IF optimized then begin
        tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        IF n_elements(tmp) EQ 2 then begin
           currentfitcube=tmp[0]
           tmppos=where('INSET' EQ tirificfirstvars)
           tirificfirst[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
           writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def'
           tmppos=where('LOOPS' EQ tirificfirstvars)
           tirificfirst[tmppos]='LOOPS=  0'
           tmppos=where('INIMODE' EQ tirificfirstvars)
           tirificfirst[tmppos]='INIMODE=  0'
           tmppos=where('RESTARTNAME' EQ tirificfirstvars)
           tirificfirst[tmppos]='RESTARTNAME= '
           openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
           for index=0,n_elements(tirificfirst)-1 do begin
              printf,1,tirificfirst[index]
           endfor
           close,1
           IF testing GE 1 then goto,skiporver
           output_name = '1stfit.'
           store_name = '1stfit_opt.'
           continue_tirific ='optimized'
           IF size(log,/TYPE) EQ 7 then begin
             openu,66,log,/APPEND
             printf,66,linenumber()+"As we are using a optimized resolution cube, we are making a full resolution model"
             close,66
           ENDIF
           run_tirific, continue_tirific, curr_run_id, bookkeeping, output_name=output_name, $
                  store_name=store_name, log=log,loops=loops,nopoints=nopoints,AC=AC1, $
                  toymodels=toymodels,run_unit=run_unit
           IF bookkeeping EQ 5 THEN goto,finishthisgalaxy

           skiporver:
           currentfitcube=currentfitcube+'_opt'
        ENDIF
     ENDIF ELSE BEGIN
       kill_tirific,curr_run_id,run_unit
     ENDELSE
     continue_tirific = 're-initialized for 2ndfit'


     tmpcube=readfits(maindir+'/'+catdirname[i]+'/1stfit.fits',hedtmp1stcube,/SILENT)
                                ;and we make a pv-diagram based on these parameters
     IF optimized then begin
        extract_pv,nooptcube,noptheader,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_1_xv.fits',float(xv),new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,header,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+'_1_xv.fits',float(xv),new_header
     ENDELSE
     extract_pv,tmpcube,hedtmp1stcube,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
     writefits,maindir+'/'+catdirname[i]+'/1stfit_xv.fits',float(xv),new_header
                                ;We check that the model actually has enough flux to be reasonable
     hedtmp1st=hedtmp1stcube
     tmpix=WHERE(tmpcube GT catnoise[i])
     IF tmpix[0] EQ -1 then begin
        comment = 'The first fit does not have flux above the noise level, this means a mis fit.'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"No flux in final fit, aborting."
           close,66
        ENDIF
        bookkeeping=5
        goto,finishthisgalaxy
     endif
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+'we use this mask for the moment maps '+catmaskname[i]
        close,66
     ENDIF
     mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',headermask,/NOSCALE,/SILENT)
                                ;mask the data cube
     tmpmask=fltarr(n_elements(tmpcube[*,0,0]),n_elements(tmpcube[0,*,0]),n_elements(tmpcube[0,0,*]))
     tmpmask[WHERE(mask GT 0.)]=tmpcube[WHERE(mask GT 0.)]
     momentsv2,tmpmask,tmpmap,hedtmp1st,0.
     writefits,maindir+'/'+catdirname[i]+'/1stfit_mom0.fits',float(tmpmap),hedtmp1st
     hedtmp1stv=hedtmp1stcube
     momentsv2,tmpmask,tmpmapv,hedtmp1stv,1.,gdlidl=gdlidl
     writefits,maindir+'/'+catdirname[i]+'/1stfit_mom1.fits',float(tmpmapv),hedtmp1stv
     hedtmp1stw=hedtmp1stcube
     momentsv2,tmpmask,tmpmapw,hedtmp1stw,2.,gdlidl=gdlidl
     writefits,maindir+'/'+catdirname[i]+'/1stfit_mom2.fits',float(tmpmapw),hedtmp1stw
     getDHI,tmpmap,hedtmp1st,Basicinfovalues[0,3],[RAhr,DEChr,Basicinfovalues[0,4]],DHI
     totflux=[TOTAL(tmpcube[tmpix])/pixperbeam,(TOTAL(2.*cutoff[0:norings[0]-1]))/(n_elements(tmpix)/pixperbeam)]
     VSYSdiff=maxchangevel
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr

                                ;write our general overview parameters
     openu,1,basicinfofile,/APPEND
     printf,1,"#After the first fit"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,2])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,3])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,4])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[n_elements(Basicinfovalues[*,5])-1,5])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(W50)),2)+'+/-'+strtrim(strcompress(string(channelwidth,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1

                                ;Let's add errors to the file
     firstfitvaluesnames=0.
     writenewtotemplate,tirificfirst,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=firstfitvaluesnames,Variables=tirificfirstvars
     errorsadd=['VROT','VROT_2','PA','PA_2','INCL','INCL_2','SDIS','SDIS_2']
     tmpfile=strarr(n_elements(tirificfirst)+8)
     added=0
                                ;Add the errors to the final tirific
                                ;file
     newrings=norings[0]
     for j=0,n_elements(tirificfirst)-1 do begin
        tmp=str_sep(strtrim(strcompress(tirificfirst[j]),2),'=')
        tmpex=WHERE(tmp[0] EQ errorsadd)
        if tmpex[0] NE -1 then begin
           tmpfile[j+added]=tirificfirst[j]
           added++
           case tmpex[0] of
              0:begin
                 IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings)
                 errs=dblarr(n_elements(sigmarot))
                 errs=sigmarot
                                ;/SIN(incarr1*!DtoR)
                 errs[0]=double(channelwidth)
                 tmpposv=WHERE(firstfitvaluesnames EQ 'VROT')
                 avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])
                 tmpavrot=WHERE(errs GT avrot/2.)

                 IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
              end
              1: begin
                 IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings)
                 errs=sigmarot
                 errs[0]=double(channelwidth)
                 tmpposv=WHERE(firstfitvaluesnames EQ 'VROT_2')
                 avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])
                 tmpavrot=WHERE(errs GT avrot/2.)

                 IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
              end
              2: begin
                 sigmapa1=replicate(1.,newrings)
                 IF n_elements(sigmapa1) LT newrings then sigmapa1=replicate(2,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa1[0:newrings-1]))),' ')
              end
              3: begin
                 sigmapa2=replicate(1.,newrings)
                 IF n_elements(sigmapa2) LT newrings then  sigmapa2=replicate(2,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa2[0:newrings-1]))),' ')
              end
              4: begin
                 sigmaincl1=replicate(4.*exp(-catinc[i]^2.5/10^3.5)+1.5,newrings)
                 IF n_elements(sigmaincl1) LT newrings then  sigmaincl1=replicate(4,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl1[0:newrings-1]))),' ')
              end
              5: begin
                 sigmaincl2=replicate(4.*exp(-catinc[i]^2.5/10^3.5)+1.5,newrings)
                 IF n_elements(sigmaincl2) LT newrings then   sigmaincl2=replicate(4,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl2[0:newrings-1]))),' ')
              end
              6: begin
                 sigmasdis = replicate(channelwidth/2.,newrings)
                 IF n_elements(sigmasdis) LT newrings then   sigmasdis=replicate(channelwidth/2.,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmasdis[0:newrings-1]))),' ')
              end
              7: begin
                 sigmasdis = replicate(channelwidth/2.,newrings)
                 IF n_elements(sigmasdis) LT newrings then   sigmasdis=replicate(channelwidth/2.,newrings)
                 tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmasdis[0:newrings-1]))),' ')
              end
              else:print,'odd'
           endcase
        endif else tmpfile[j+added]=tirificfirst[j]
     endfor
     tirificfirst=tmpfile
                                ;write the final file and produce the final model
     openw,1,maindir+'/'+catdirname[i]+'/1stfit.def'
     for index=0,n_elements(tirificfirst)-1 do begin
        printf,1,tirificfirst[index]
     endfor
     close,1






     IF finishafter EQ 1 then begin
        comment = 'You have chosen to skip the fitting process after the first fit'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        bookkeeping=bookkeeping+0.5
        goto,finishthisgalaxy
     ENDIF
     prevmodification=0
;******************************************This is the end of the first fit******************************
     endtime_firststep = systime()
     sigmapa1=0.
     sigmapa2=0.
     sigmaincl1=0.
     sigmaincl2=0.
     sigmasdis=0.
     sigmarot=0.
     velfixrings=1
     tmppos=where('INSET' EQ tirificsecondvars)
     tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
                                ;then open the previous fit
     firstfitvaluesnames=0.
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/1stfit.def',Arrays=firstfitvalues,VariableChange=firstfitvaluesnames,Variables=tirificsecondvars
     tmppos=where('XPOS' EQ firstfitvaluesnames)
     Final1stXPOS=firstfitvalues[0,tmppos]
     tmppos=where('YPOS' EQ firstfitvaluesnames)
     Final1stYPOS=firstfitvalues[0,tmppos]

     tmppos=where('RADI' EQ firstfitvaluesnames)
     RADarr=firstfitvalues[*,tmppos]
     tmppos=where('VROT' EQ firstfitvaluesnames)
     VROTarr=dblarr(n_elements(firstfitvalues[*,tmppos]))
     Vfrom1=dblarr(n_elements(firstfitvalues[*,tmppos]))
     Vfrom1=firstfitvalues[*,tmppos]
     half=fix(n_elements(firstfitvalues[*,tmppos])/2.)
     VROTarr[*]=firstfitvalues[*,tmppos]
     tmppos=where('SDIS' EQ firstfitvaluesnames)
     SDISarr=dblarr(n_elements(firstfitvalues[*,tmppos]))
     tmppos=where('SBR' EQ firstfitvaluesnames)
     SBRarr=firstfitvalues[*,tmppos]
     tmppos=where('INCL' EQ firstfitvaluesnames)
     INCLang=firstfitvalues[*,tmppos]
     tmppos=where('PA' EQ firstfitvaluesnames)
     PAang=firstfitvalues[*,tmppos]
     VROTarr2=VROTarr
     tmppos=where('SBR_2' EQ firstfitvaluesnames)
     SBRarr2=firstfitvalues[*,tmppos]
                                ;We always want to smooth the surface brightnes. Added 16-06-2017
     SBRarr=fat_savgol(SBRarr,RADarr)
     SBRarr2=fat_savgol(SBRarr2,RADarr)
     SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.
     SBRarr2[0:1]=SBRarr[0:1]
                                ;To make sure that we fit a warp we
                                ;want to increase the brightness of
                                ;the last two rings after the first
                                ;fit
     SBRarr[n_elements(SBRarr)-3:n_elements(SBRarr)-1]=SBRarr[n_elements(SBRarr)-3:n_elements(SBRarr)-1]*1.2
     SBRarr2[n_elements(SBRarr2)-3:n_elements(SBRarr2)-1]=SBRarr[n_elements(SBRarr2)-3:n_elements(SBRarr2)-1]*1.2

     tmppos=where('SBR' EQ tirificsecondvars)
     tirificsecond[tmppos]='SBR= '+STRJOIN(SBRarr[0:n_elements(SBRarr)-1],' ')
     tmppos=where('SBR_2' EQ tirificsecondvars)
     tirificsecond[tmppos]='SBR_2= '+STRJOIN(SBRarr2[0:n_elements(SBRarr2)-1],' ')
     tmppos=where('INCL_2' EQ firstfitvaluesnames)
     INCLang2=firstfitvalues[*,tmppos]
     IF catinc[i] GT 80. then begin
        INCLang=catinc[i]
        INCLang2=catinc[i]
     ENDIF
     tmppos=where('PA_2' EQ firstfitvaluesnames)
     PAang2=firstfitvalues[*,tmppos]
     tmppos=where('NUR' EQ firstfitvaluesnames)
     norings=firstfitvalues[0,tmppos]
     IF norings*ring_spacing LE 4 then begin
        IF norings*ring_spacing GE 2 and finishafter NE 1.1 and ring_spacing GT 0.75 then begin
           ring_spacing_new = ring_spacing/1.5
           norings[0]=round(norings[0]*ring_spacing/ring_spacing_new)
        ;we always want at least 3 beams in the model
           WHILE ring_spacing_new GT 0.75 AND norings[0] LT 3 do begin
              ring_spacing_new = ring_spacing_new/2.
              norings[0]=round(norings[0]*ring_spacing/ring_spacing_new)
           ENDWHILE

           IF norings[0] LT 3. then begin
              comment = 'The first fit model is too small to fit variations.'
              commentlen='A'+strtrim(string(strlen(comment)),2)
              openu,1,outputcatalogue,/APPEND
              printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
              close,1
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+'The first fit model is too small to fit variations.'
                 printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
                 close,66
              ENDIF
              IF optimized then begin
                 currentfitcube = noptname[0]
              ENDIF
              bookkeeping=bookkeeping+0.5
              goto,finishthisgalaxy
           ENDIF
           maxrings=maxrings*ring_spacing/ring_spacing_new
           noringspix=norings[0]*(catmajbeam[i]*ring_spacing_new)/(ABS(sxpar(headermap,'cdelt1'))*3600.)
           rad=[0.,((findgen(round(maxrings)))*(catmajbeam[i]*ring_spacing_new)+catmajbeam[i]/5)]
           calc_edge,catnoise[i],rad,[catmajbeam[i],catminbeam[i],channelwidth],cutoffor,vsys=catVSYS[i]


           ring_spacing = ring_spacing_new
           finishafter=1.1
           cutoff=cutoffor*cutoffcorrection
           tmp=rad[0:norings[0]-1]
           rad=tmp
           tmp=1
           interpolate,VROTarr,RADarr,newradii=rad,output=tmp
           VROTarr=tmp
           VROTarr2=tmp
           tmp=1
           interpolate,SBRarr,RADarr,newradii=rad,output=tmp
           SBRarr=tmp
           tmp=1
           interpolate,SBRarr2,RADarr,newradii=rad,output=tmp
           SBRarr2=tmp
           SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.
           SBRarr2[0:1]=SBRarr[0:1]
           RADarr=rad
           tmppos=where('RADI' EQ tirificsecondvars)
           tirificsecond[tmppos]='RADI= '+STRJOIN(strtrim(strcompress(string(RADarr))),' ')
           tmppos=where('SBR' EQ tirificsecondvars)
           tirificsecond[tmppos]='SBR= '+STRJOIN(strtrim(strcompress(string(SBRarr))),' ')
           tmppos=where('SBR_2' EQ tirificsecondvars)
           tirificsecond[tmppos]='SBR_2= '+STRJOIN(strtrim(strcompress(string(SBRarr2))),' ')
           tmppos=where('NUR' EQ tirificsecondvars)
           tirificsecond[tmppos]='NUR= '+STRJOIN(strtrim(strcompress(string(norings[0]))),' ')


        endif else begin
           comment = 'The first fit model is too small to fit variations.'
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'The first fit model is too small to fit variations.'
              printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           IF optimized then begin
              currentfitcube = noptname[0]
           ENDIF
           IF finishafter GT 1.75 then finishafter=finishafter-1
           bookkeeping=bookkeeping+0.5
           goto,finishthisgalaxy
        ENDELSE
     ENDIF
;Let's see how many of the inner rings we want to fix
     tmppos=where('VSYS' EQ firstfitvaluesnames)
     vsys=firstfitvalues[0,tmppos]
     levels=(sbrarr+sbrarr2)/2.*1000.
     Columndensity,levels,vsys,[catmajbeam[i],catminbeam[i]],/ARCSQUARE
     tmp=WHERE(levels GT 2E20)

     tmp=WHERE(levels GT 2E20)
     IF tmp[0] NE -1 then innerfix=floor(tmp[n_elements(tmp)-1]/1.5)-1. else innerfix=4
     IF innerfix LT 4 OR innerfix GE norings[0] OR finishafter EQ 1.1 then innerfix=4
     IF centralexclude then begin
        cen=0
        WHILE levels[cen] LT 1E20 AND cen LT n_elements(levels)-1 DO cen++
        IF cen GT innerfix then innerfix=cen else innerfix++
     ENDIF




     lowring=norings[0]-(3/ring_spacing)
     highring=norings[0]+(2/ring_spacing)

     tmppos=where('VROT' EQ tirificsecondvars)
     tirificsecond[tmppos]='VROT= 0.'+STRJOIN(strtrim(strcompress(string(VROTarr[1:n_elements(VROTarr)-1]))),' ')
     tmppos=where('VROT_2' EQ tirificsecondvars)
     tirificsecond[tmppos]='VROT_2= 0. '+STRJOIN(strtrim(strcompress(string(VROTarr[1:n_elements(VROTarr)-1]))),' ')
     tmppos=where('DISTANCE' EQ tirificsecondvars)
     tirificsecond[tmppos]='DISTANCE='+strtrim(strcompress(string(catDistance[i])))
                                ;check the sbr
     sbr_check,tirificsecond, tirificsecondvars,sbrarr,sbrarr2,cutoff

     IF norings[0] LT 4 then begin
        norings[0]=4.
        changeradii,tirificsecond,norings[0],continue_tirific
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We added a ring to the model cause there were too few the new number of rings = "+strtrim(string(fix(norings[0])),2)
           Close,66
        ENDIF ELSE BEGIN
           print,linenumber()+"We added a ring to the model cause there were too few the new number of rings = "+strtrim(string(fix(norings[0])),2)
        ENDELSE
     ENDIF

     tmppos=where('VSYS' EQ firstfitvaluesnames)
     catvsys[i]=firstfitvalues[0,tmppos]
     newrings=norings[0]
                                ;Let's set some limits for the
                                ;first make the best case scenario arrays
                                ;INCL
     INCLtir=(TOTAL(INCLang)+TOTAL(INCLang2))/(n_elements(INCLang)+n_elements(INCLang2))
     INCLmax=INCLtir+30
     IF INCLmax GT 90. then INCLmax=90
     INCLmin=INCLtir-30
     IF INCLmin LT 0.1 then INCLmin=0.1
     INCLinput1=['INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1),string(INCLmax),string(INCLmin),string(1.),string(0.1),string(0.5),string(0.1),'3','70','70']
                                ;PA
     PAtir=(TOTAL(PAang)+TOTAL(PAang2))/(n_elements(PAang)+n_elements(PAang2))
     PAinput1=['PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1),string(PAtir+40),string(PAtir-40),string(5),string(0.1),string(0.1),string(0.1),'5','70','70']
                                ;VROT
     VROTmax=MAX(Vfrom1[1:n_elements(Vfrom1)-1],Min=VROTmin)
     VROTmax=VROTmax*1.5
     VROTmin=VROTmin/4.
                                ;the minimum can never be less than
                                ;the channelwidth
     IF VROTmin LT channelwidth then VROTmin=channelwidth
                                ; The maximum should not be less than
                                ; 80 /sin(inclination) as it becomes
                                ; more uncartain at lower inclination
     IF VROTmax LT 60. then VROTmax=60.
     IF vrotmax GT 600. then VROTmax=600.
                                ;SDIS
     SDISmax=MAX(SDISarr,Min=SDISmin)
     SDISmax=SDISmax*1.5
     SDISmin=SDISmin/4.
                                ;the minimum can never be less than
                                ;the channelwidth
     IF SDISmin LT channelwidth then SDISmin=channelwidth
                                ; The maximum should not be less than
                                ; 80 /sin(inclination) as it becomes
                                ; more uncartain at lower inclination
     IF SDISmax LT 25. then SDISmax=25.
     IF sdismax GT 50. then SDISmax=50.
                                ; See how much of the rotation curve we want to fit as a slope
     get_newringsv9,SBRarr,SBRarr2,2.5*cutoff,velconstused
     velconstused--
     IF norings[0]*ring_spacing GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1
     set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
     set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,avinner=avinner,finish_after=finishafter,fix=fix_sdis
                                ;set the surface brightness values
     set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,/initial,doubled=doubled
                                ;SDIS 14-09-2018 let's make an
                                ;attempt at fitting the dispersion in
                                ;a similar fashion to the pa
                                ;Z0
     Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
               ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.005,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
               string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),'3','70','70']
     IF catinc[i] GT 80 then begin
        Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
     ENDIF
     IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([catmajbeam[i]*ring_spacing*norings[0]/60.,catmajbeam[i]/2.]))
        ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*ring_spacing*norings[0]/60.)
        Z0input1[2]='0.'
        Z0input1[3]=string(-1*catmajbeam[i]*ring_spacing*norings[0]/6E4)
        Z0input1[4]=string(catmajbeam[i]*ring_spacing*norings[0]/6E5)
        Z0input1[5]=string(catmajbeam[i]*ring_spacing*norings[0]/60.)
        Z0input1[6]=string(catmajbeam[i]*ring_spacing*norings[0]/6E5)
     ENDIF
                                ;And then make the other string input variables with the same
                                ;parameters and an additional string
                                ;where we can put slope fitting for
                                ;low SBR rings
     INCLinput2=[INCLinput1,' ']
     INCLinput3=[INCLinput1,' ']
                                ;If we have a high inclination we put stringent limits on the inner
                                ;inclination
     IF catinc[i] GT 80 then begin
        INCLinput1[1]=catinc[i]+catincdev[i]
        IF INCLinput1[1] GT 90. then INCLinput1[1]='90'
        INCLinput1[2]=catinc[i]-catincdev[i]
        IF INCLinput1[2] LT 5 then INCLinput1[2]='5'
     ENDIF
                                ;PA
     PAinput2=[PAinput1,' ']
     PAinput3=[PAinput1,' ']
                                ;IF we have a decent amount of rings
                                ;and we are not fitting a flat disk
                                ;with half beams than we let the rings
                                ;beyond 4 free
     IF norings[0] GT 4 AND finishafter NE 1.1 then begin
                                ;And some different fitting parameters
                                ;for the free rings
        INCLinput2[5]='0.5'
        INCLinput2[3]='1.0'
        INCLinput3[5]='0.5'
        INCLinput3[3]='1.0'
        PAinput2[5]='2.0'
        PAinput2[3]='0.5'
        PAinput3[5]='2.0'
        PAinput3[3]='0.5'
                                                          ;update INCL and PA
        set_warp_slopev3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix
        if fix_incl EQ 0. then INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                             ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        if fix_pa EQ 0. then PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                         ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
     endif else begin
                                ;Otherwise we just fit  a single flat disk
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        IF doubled then begin

        ENDIF else SBRinput1[2]=strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/4.,format='(E12.5)')),1)
     endelse
     ;Set the inimode with new newradii
     case 1 of
       norings[0] LT 5: inimode = 2
       else: inimode = 3
     endcase

                                ;If the first fit is accepted we only change the central position minimally
     IF AC1 EQ 1 then begin
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),'360','0', +$
                      strtrim(string(pixelsizeRA)),strtrim(string(pixelsizeRA/20.)),+$
                      strtrim(string(pixelsizeRA/10.)),strtrim(string(pixelsizeRA/20.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),'90','-90',+$
                     strtrim(string(pixelsizeDEC)),strtrim(string(pixelsizeDEC/20.)),+$
                     strtrim(string(pixelsizeDEC/10.)),strtrim(string(pixelsizeDEC/20.)),'3','70','70']

        vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),strtrim(strcompress(string(catvsys[i]+100.)),1),$
                     strtrim(strcompress(string(catvsys[i]-100.)),1),'0.5','0.01','2','0.5','3','70','70']
        case 1 of
           norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
              Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           fix_incl EQ 0. and fix_pa EQ 1.: begin
              Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           fix_incl EQ 1. and fix_pa EQ 0.: begin
              Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,$
                                 vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           else: begin
              Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
        endcase

     endif else begin
                                ;Else we allow more variation in the
                                ;center and fit them first
        xposinput1=[ ' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),'360','0', +$
                      strtrim(string(pixelsizeRA*3.)),strtrim(string(pixelsizeRA/10.)),+$
                      strtrim(string(pixelsizeRA)),strtrim(string(pixelsizeRA/10.)),'3','70','70']
        yposinput1=[ ' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),'90','-90',+$
                     strtrim(string(pixelsizeDEC*3.)),strtrim(string(pixelsizeDEC/10.)),+$
                     strtrim(string(pixelsizeDEC)),strtrim(string(pixelsizeDEC/10.)),'3','70','70']
        vsysinput1=[ ' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                     ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),strtrim(strcompress(string(catvsys[i]+100.)),1),strtrim(strcompress(string(catvsys[i]-100.)),1),'2','0.5','2','0.5','3','70','70']
        case 1 of
           norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
              Writefittingvariables,tirificsecond, xposinput1,yposinput1,vsysinput1,inclinput1,painput1,vrotinput1,sdisinput1,$
              sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
           end
           fix_incl EQ 0. and fix_pa EQ 1.: begin
              Writefittingvariables,tirificsecond, xposinput1,yposinput1,vsysinput1,inclinput1,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
           end
           fix_incl EQ 1. and fix_pa EQ 0.: begin
              Writefittingvariables,tirificsecond, xposinput1,yposinput1,vsysinput1,inclinput1,inclinput2,inclinput3,painput1,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
           end
           else: begin
              Writefittingvariables,tirificsecond, xposinput1,yposinput1,vsysinput1,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
           end
        endcase
     endelse
                                ;Setting a bunch of fit tracking variables
     forcedring=0.
     letvrotvary=0.
     counter=0.
     trytwo=0.
     constring=0.
     sbrmodify=0.
     slopedrings=[0,0]
     prevslopedrings=[0,0]
     lastcutrings=0.
     noacceptbeforesmooth=0.
     lastaddrings=0.
     polorder1=[!values.f_nan,!values.f_nan]
     polorder2=[!values.f_nan,!values.f_nan]
     centrejump=0.
     prevrings=norings[0]
     secondtime=0.
     prevXPOS=Final1stXPOS
     prevYPOS=Final1stYPOS
                                ;If the second fit is not accepted we
                                ;come back to this point to do the refitting
     notacceptedtwo:

     prevXPOS=newXPOS
     prevYPOS=newYPOS
     prevslopedrings=slopedrings
                                ;If we have optimized the cube we want
                                ;to set the _opt extension to the
                                ;input cubes name. Otherwise the input
                                ;cube should have properly copied from
                                ;the 1st fit
     IF optimized then begin
        tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        if n_elements(tmp) LT 2 then begin
           currentfitcube=tmp[0]+'_opt'
           tmppos=where('INSET' EQ tirificsecondvars)
           tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
        ENDIF
     ENDIF
     ;case 1 of
    ;    norings[0] LT 5: memode=0
    ;    norings[0] LT 15: memode=1
    ;    norings[0] LT 30: memode=2
    ;    else: memode=3
     ;endcase
     ;tmppos=where('INIMODE' EQ tirificsecondvars)
     ;tirificsecond[tmppos]='INIMODE=  '+strtrim(strcompress(string(memode)))
     tmppos=where('LOOPS' EQ tirificsecondvars)
     case 1 of
        trytwo LT 5: begin
           tirificsecond[tmppos]='LOOPS=  12'
        end
        trytwo LT 10: begin
           tirificsecond[tmppos]='LOOPS=  8'
        end
        else: begin
           tirificsecond[tmppos]='LOOPS=  5'
        end
     endcase



;If we are testing than skip the fitting
     IF testing GE 2 then goto,testing2
                                ;Write the tirific file and update log
     restart_counter++
     tmppos=where('RESTARTID' EQ tirificsecondvars)
     tirificsecond[tmppos]='RESTARTID= '+string(fix(restart_counter))
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificsecond)-1 do begin
        printf,1,tirificsecond[index]
     endfor
     close,1
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific Second estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        printf,66,linenumber()+"We have changed the ring numbers "+string(sbrmodify)+" times."
        printf,66,linenumber()+"We have changed the fitting parameters "+string(trytwo)+" times."
        close,66
     ENDIF
     print,linenumber()+"Starting tirific Second estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
                                ;rename the files
     counter++
     IF MEAN(VROTarr) GT 1000. then begin
        print,'something has gone horribly wrong'
        stop
     ENDIF
     nopoints=0
     run_tirific, continue_tirific, curr_run_id,bookkeeping,output_name='2ndfit.', $
                  store_name='2ndfitold.',log=log,loops=loops,nopoints=nopoints,AC=AC2, $
                  toymodels=toymodels,run_unit=run_unit



     testing2:
                                ;Then we check wether the fit is accepted if not we see if the PA and
                                ;INCL can be accepted IF the positions were not accepted before we try
                                ;then as well

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC2 EQ 1 then begin
           printf,66,linenumber()+"The Second estimate was accepted in this run."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models"
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The Second estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models"
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF
                                ;let's moderate the cflux to ensure between 1.5 and 3 million points
     check_cflux,nopoints,norings[0],tirificsecond,tirificsecondvars,cfluxadjusted,log=log
     IF AC2 EQ 1 then print,linenumber()+"The second estimate was accepted in this run." $
     else print,linenumber()+"The second estimate was not accepted."
     secondfitvaluesnames=0.
     secondfitvalues=0.
     finalsmooth=0.
                                ; If the fit is fully accepted and we
                                ; only want to do a smoothing we come
                                ; back to this point
     lastadjust:

                                ;We get the previous fit and read it into the fitting template
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=secondfitvalues,VariableChange=secondfitvaluesnames,Variables=tirificsecondvars
     tmppos=where('XPOS' EQ secondfitvaluesnames)
     newXPOS=secondfitvalues[0,tmppos]
     tmppos=where('YPOS' EQ secondfitvaluesnames)
     newYPOS=secondfitvalues[0,tmppos]
     tmppos=where('RADI' EQ secondfitvaluesnames)
     RADarr=secondfitvalues[*,tmppos]
     maxrad=MAX(RADarr)
     tmppos=where('VROT' EQ secondfitvaluesnames)
     VROTarr=secondfitvalues[*,tmppos]
                                ;We check that it is not out off bound
     tmp=WHERE(VROTarr GT VROTmax)
     IF tmp[0] ne -1 then VROTarr[tmp]=VROTmax*0.8
     tmp=WHERE(VROTarr LT VROTmin)
     IF tmp[0] ne -1 then begin
        IF tmp[0] EQ 0 then begin
           IF n_elements(tmp) GT 1 then VROTarr[tmp[1:n_elements(tmp)-1]]=VROTmin*2.
        ENDIF ELSE VROTarr[tmp]=VROTmin*2.
     ENDIF
     tmppos=where('SDIS' EQ secondfitvaluesnames)
     SDISarr=secondfitvalues[*,tmppos]
                                ;We check that it is not out off bound
     tmp=WHERE(SDISarr GT SDISmax)
     IF tmp[0] ne -1 then SDISarr[tmp]=SDISmax*0.8
     tmp=WHERE(SDISarr LT SDISmin)
     IF tmp[0] ne -1 then SDISarr[tmp]=SDISmin*1.5

     tmppos=where('SBR' EQ secondfitvaluesnames)
     SBRarr=secondfitvalues[*,tmppos]
     SBRarrunmod=SBRarr
                                ;We always want to smooth the surface brightnes. Added 16-06-2017


                                ;and take the central point as the
                                ;extension of the previous two. Added
                                ;18-10-2018. This is done in
                                ;fat_savgol already
     ;inSBR=SBRarr[1:2]
     ;inRad=RADarr[1:2]
     ;newRAD=RADarr[0:2]
     ;newSBR=1
     ;interpolate,inSBR,inRad,newradii=newRAD,output=newSBR
     ;SBRarr[0]=newSBR[0]
    ;

     tmppos=where('INCL' EQ secondfitvaluesnames)
     INCLang=secondfitvalues[*,tmppos]
     tmppos=where('PA' EQ secondfitvaluesnames)
     PAang=secondfitvalues[*,tmppos]
     tmppos=where('VROT_2' EQ secondfitvaluesnames)
     VROTarr2=secondfitvalues[*,tmppos]
                                   ;We check that it is not out off bound
     tmp=WHERE(VROTarr2 GT VROTmax)
     IF tmp[0] ne -1 then VROTarr2[tmp]=VROTmax*0.8
     tmp=WHERE(VROTarr2 LT VROTmin)
     IF tmp[0] ne -1 then begin
        IF tmp[0] EQ 0 then begin
           IF n_elements(tmp) GT 1 then VROTarr2[tmp[1:n_elements(tmp)-1]]=VROTmin*2.
        ENDIF ELSE VROTarr2[tmp]=VROTmin*2.
     ENDIF

     tmppos=where('SDIS_2' EQ secondfitvaluesnames)
     SDISarr2=secondfitvalues[*,tmppos]
                                ;We check that it is not out off bound
     tmp=WHERE(SDISarr2 GT SDISmax)
     IF tmp[0] ne -1 then SDISarr2[tmp]=SDISmax*0.8
     tmp=WHERE(SDISarr2 LT SDISmin)
     IF tmp[0] ne -1 then SDISarr2[tmp]=SDISmin*1.5


     tmppos=where('SBR_2' EQ secondfitvaluesnames)
     SBRarr2=secondfitvalues[*,tmppos]
     SBRarr2unmod=SBRarr2
                                ;We always want to smooth the surface brightnes. Added 16-06-2017
    ; IF finalsmooth GE 1 OR norings[0]*ring_spacing GT 4 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Smoothing the full SBR profile."
           Close,66
        ENDIF
        SBRarr=fat_savgol(SBRarr,RADarr)
        SBRarr2=fat_savgol(SBRarr2,RADarr)
     ;endif else begin
     ;   IF size(log,/TYPE) EQ 7 then begin
     ;      openu,66,log,/APPEND
     ;      printf,66,linenumber()+"Smoothing half the SBR profile."
     ;      Close,66
     ;   ENDIF
     ;   SBRarr=fat_savgol(SBRarr,RADarr,/Half)
     ;   SBRarr2=fat_savgol(SBRarr2,RADarr,/Half)
     ;endif

     SBRarr[0]=(SBRarr[0]+SBRarr2[0])/2.
     SBRarr2[0]=SBRarr[0]
     tmppos=where('SBR' EQ tirificsecondvars)
     tirificsecond[tmppos]='SBR= '+STRJOIN(SBRarr[0:n_elements(SBRarr)-1],' ')
     tmppos=where('SBR_2' EQ tirificsecondvars)
     tirificsecond[tmppos]='SBR_2= '+STRJOIN(SBRarr2[0:n_elements(SBRarr2)-1],' ')
     SBRarror=SBRarr
     SBRarr2or=SBRarr2

     tmppos=where('INCL_2' EQ secondfitvaluesnames)
     INCLang2=secondfitvalues[*,tmppos]
     tmppos=where('PA_2' EQ secondfitvaluesnames)
     PAang2=secondfitvalues[*,tmppos]
     IF PAang2[0] NE PAang[0] then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is where the inner PA starts to deviate."
           Close,66
        ENDIF
     ENDIF
                                ;First we check that the center has
                                ;not jumped. If it has we reset the
                                ;central position. And make sure that
                                ;SBR is at least 2 the cutoff and
                                ;refit

                                ; We do not want the center to deviate
                                ; more than 2 beams from the previous
                                ; fit, 4 beams of the 1st fit or more
                                ; than a 20th of the diameter of the galaxy

     IF ABS(newXPOS-prevXPOS) GT catmajbeam[i]/1800. OR  ABS(newYPOS-prevYPOS) GT catmajbeam[i]/1800. OR $
        ABS(newXPOS-Final1stXPOS) GT catmajbeam[i]/900. OR ABS(newYPOS-Final1stYPOS) GT catmajbeam[i]/900. OR $
        ABS(newXPOS-prevXPOS) GT maxrad/36000. OR ABS(newYPOS-prevYPOS) GT  maxrad/36000. then begin
        IF ABS(newXPOS-prevXPOS) GT catmajbeam[i]/1800. OR  ABS(newYPOS-prevYPOS) GT catmajbeam[i]/1800. then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted more than 2 major beams. Not applying this shift and refitting."
              printf,66,linenumber()+"The center is bad at "+strtrim(string(newXPOS),2)+', '+strtrim(string(newYPOS),2)
              printf,66,linenumber()+"It was at "+strtrim(string(prevXPOS),2)+', '+strtrim(string(prevYPOS),2)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center shifted than 2 major beams. Not applying this shift and refitting."
           ENDELSE
        ENDIF
        IF ABS(newXPOS-Final1stXPOS) GT catmajbeam[i]/900. OR ABS(newYPOS-Final1stYPOS) GT catmajbeam[i]/900.  then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted more than 4 major beams from the original fit. Not applying this shift and refitting."
              printf,66,linenumber()+"The center is bad at "+strtrim(string(newXPOS),2)+', '+strtrim(string(newYPOS),2)
              printf,66,linenumber()+"It was at "+strtrim(string(prevXPOS),2)+', '+strtrim(string(prevYPOS),2)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center shifted than 4 major beams from the original fit. Not applying this shift and refitting."
           ENDELSE
        ENDIF
        IF ABS(newXPOS-prevXPOS) GT maxrad/36000. OR ABS(newYPOS-prevYPOS) GT  maxrad/36000. then begin
        IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted more 10% of the galaxy maximum radius from the original fit. Not applying this shift and refitting."
              printf,66,linenumber()+"The center is bad at "+strtrim(string(newXPOS),2)+', '+strtrim(string(newYPOS),2)
              printf,66,linenumber()+"It was at "+strtrim(string(prevXPOS),2)+', '+strtrim(string(prevYPOS),2)
              close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"The center  shifted more 10% of the galaxy maximum radius from the original fit. Not applying this shift and refitting."
           ENDELSE
        ENDIF
;If we have jumped from the centre too often then abort
        IF centrejump GT 10 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The center shifted too often out of the set boundaries."
              printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy #  "+strtrim(string(fix(i)),2)+" at "+systime()
              close,66
           ENDIF
           comment = 'The first fit centre kept jumping away.'
           commentlen='A'+strtrim(string(strlen(comment)),2)
           openu,1,outputcatalogue,/APPEND
           printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(0),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
           close,1
           bookkeeping=5
           goto,finishthisgalaxy
        ENDIF

        tmppos=where('XPOS' EQ tirificsecondvars)
        tirificsecond[tmppos]=' XPOS= '+string(prevXPOS)
        tmppos=where('XPOS_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=' XPOS_2= '+string(prevXPOS)
        tmppos=where('YPOS' EQ tirificsecondvars)
        tirificsecond[tmppos]=' YPOS= '+string(prevYPOS)
        tmppos=where('YPOS_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=' YPOS_2= '+string(prevYPOS)
        for j=2,n_elements(SBRarr)-1 do begin
           IF SBRarr[j] LT cutoff[j] then SBRarr[j]=cutoff[j]*3.
           IF SBRarr2[j] LT cutoff[j] then SBRarr2[j]=cutoff[j]*3.
        endfor
        AC2=0
        centrejump++
     ENDIF ELSE BEGIN
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The center is swell it is at "+strtrim(string(newXPOS),2)+', '+strtrim(string(newYPOS),2)
           printf,66,linenumber()+"It was at "+strtrim(string(prevXPOS),2)+', '+strtrim(string(prevYPOS),2)
           close,66
        ENDIF
     ENDELSE
                                ;Sometimes TiRiFiC can place values
                                ;outside the boundaries leading to
                                ;peculiar effects. Hence we first check
                                ;that the PA and INCL are within our
                                ;set boundaries. This seems to make
                                ;matters worse. because of messing up
                                ;the PA. Let's turn off PA for now (15-05-2017)
     check=0.
     WHILE check EQ 0 DO BEGIN
        tmp=WHERE(INCLang GT 90.)
        IF tmp[0] NE -1 then INCLang[tmp]=180-INCLang[tmp] else check=1
     ENDWHILE
     check=0.
     WHILE check EQ 0 DO BEGIN
        tmp=WHERE(INCLang2 GT 90.)
        IF tmp[0] NE -1 then INCLang2[tmp]=180-INCLang2[tmp] else check=1
     ENDWHILE
     check=0.
     WHILE check EQ 0 DO BEGIN
        tmp=WHERE(INCLang LT 0.)
        IF tmp[0] NE -1 then INCLang[tmp]=ABS(INCLang[tmp]) else check=1
     ENDWHILE
     check=0.
     WHILE check EQ 0 DO BEGIN
        tmp=WHERE(INCLang2 LT 0.)
        IF tmp[0] NE -1 then INCLang2[tmp]=ABS(INCLang2[tmp]) else check=1
     ENDWHILE

;Let's see how many of the inner rings we want to fix
     tmppos=where('VSYS' EQ firstfitvaluesnames)
     vsys=firstfitvalues[0,tmppos]
     levels=(sbrarr+sbrarr2)/2.*1000.
     Columndensity,levels,vsys,[catmajbeam[i],catminbeam[i]],/ARCSQUARE
     tmp=WHERE(levels GT 2E20)
     IF tmp[0] NE -1 then innerfix=floor(tmp[n_elements(tmp)-1]/1.5)-1.
     IF innerfix LT 4 OR innerfix GE norings[0] OR finishafter EQ 1.1 then innerfix=4
     IF centralexclude then begin
        cen=0
        WHILE levels[cen] LT 1E20 AND cen LT n_elements(levels)-1 DO cen++
        IF cen GT innerfix then innerfix=cen else innerfix++
     ENDIF

     tmppos=where('NUR' EQ secondfitvaluesnames)
     norings[0]=secondfitvalues[0,tmppos]
     sbr_check,tirificsecond, tirificsecondvars,sbrarr,sbrarr2,cutoff
                                ;We get the rings for which we only
                                ;want to fit a slope in the rotation curve
     get_newringsv9,SBRarr,SBRarr2,2.5*cutoff,velconstused
     velconstused--
     IF norings[0]*ring_spacing GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1
     xind=0
     ;IF (TOTAL(VROTarr[1:2])/2. GT 120 AND VROTarr[1] GT VROTarr[2] AND VROTarr[1] GT VROTarr[3]) OR $
     ;   (MEAN(VROTarr[2:n_elements(vrotarr)-1]) GT 200.) then begin
     IF TOTAL(VROTarr[1:2])/2. GT 180. OR $
        (MEAN(VROTarr[2:n_elements(vrotarr)-1]) GT 120.) then begin
        x=n_elements(VROTarr)-1
        WHILE VROTarr[x] GT  VROTarr[x-1] AND x GT fix(n_elements(VROTarr)/2) DO x--


        IF x LT n_elements(VROTarr)-1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+'This is a massive galaxy hence we flatten the outer part if '+strtrim(string(x),2)+' is less '+$
                     strtrim(string(n_elements(VROTarr)-1),2)
              printf,66,linenumber()+'From ring '+strtrim(string(x),2)+' on we have the value '+strtrim(string(VROTarr[x]),2)
              close,66
           ENDIF
           VROTarr[x:n_elements(VROTarr)-1]=VROTarr[x]
           stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
           tmppos=where('VROT' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
           stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
           tmppos=where('VROT_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
           slope=1
           xind=x
           IF norings[0]-x GT velconstused then velconstused=norings[0]-x
        ENDIF ELSE BEGIN
           IF MEAN(VROTarr[1:n_elements(vrotarr)-1]) GT 200. then begin
              min=MAX(VROTarr)
              xind=n_elements(VROTarr)-1
              for x=fix(n_elements(VROTarr)/2),n_elements(VROTarr)-1 do begin
                 IF VROTarr[x] LT min then begin
                    min=VROTarr[x]
                    xind=x
                 ENDIF
              ENDFOR
              IF xind LT n_elements(VROTarr)-1 then begin
                 VROTarr[xind:n_elements(VROTarr)-1]=VROTarr[xind]
                 stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
                 tmppos=where('VROT' EQ tirificsecondvars)
                 tirificsecond[tmppos]=stringVROT
                 stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
                 tmppos=where('VROT_2' EQ tirificsecondvars)
                 tirificsecond[tmppos]=stringVROT
                 slope=1
                                ;            IF norings[0]-xind GT velconstused then velconstused=norings[0]-xind
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+'This is a massive galaxy hence we flatten the outer part if xind is less '+$
                           strtrim(string(n_elements(VROTarr)-1),2)
                    printf,66,linenumber()+'From ring '+strtrim(string(xind),2)+' on we have the value '+strtrim(string(VROTarr[xind]),2)
                    close,66
                 ENDIF
              ENDIF
           ENDIF
        ENDELSE
     ENDIF


     locmax=WHERE(MAX(VROTarr) EQ VROTarr)
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Vmax occurs at ring no "+string(locmax[n_elements(locmax)-1]+1)
        close,66
     ENDIF
     slope=0
     if locmax[n_elements(locmax)-1] LT n_elements(VROTarr)/2. AND MEAN(VROTarr[fix(n_elements(VROTarr)/2.):n_elements(VROTarr)-1]) GT 120. then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This curve is declining so we flatten the outerpart instead of fitting a slope"
           close,66
        ENDIF
        slope=1
        velfixrings=norings[0]-velconstused
        IF velfixrings LT 1 then velfixrings=1
        IF velfixrings EQ 1 AND norings[0] GT 5 then velfixrings=2
        IF velfixrings GT 1 then VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1]=TOTAL(VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1])/n_elements(VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1])
     ENDIF
     set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,finish_after=finishafter,slope=slope
     set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,finish_after=finishafter,slope=slope,fix=fix_sdis
                                ;Then set the surface brighness profile parameters
     set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,doubled=doubled
                                ;Update the rings to be fitted in the other parameters

     IF norings[0] LE 4 or finishafter EQ 1.1 then begin
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
     ENDIF ELSE BEGIN
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;update INCL and PA
        set_warp_slopev3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix,sloped=slopedrings
        if fix_incl EQ 0. then   INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                               ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        if fix_pa EQ 0. then   PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                           ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
     ENDELSE

                                ;additionallly we should updat the min
                                ;and max in the fitting if
                                ;two rings reach it
                                ;This doesn't happen when the
                                ;fit is accepted
                                ;when this happens all rings after are
                                ;unreliable and should be set to the
                                ;last ring that isn't affected
     boundaryadjustment=0.
     lastreliablerings=n_elements(INCLang)

     PAinput1old=PAinput1
     INCLinput1old = INCLinput1
     PAinput2old=PAinput2
     INCLinput2old = INCLinput2
     PAinput3old=PAinput3
     INCLinput3old = INCLinput3


     IF double(INCLinput1[1]) LT 90. OR double(INCLinput2[1]) LT 90. OR  double(INCLinput3[1]) LT 90.  AND norings[0] GT 4 then begin
        ;tmp=WHERE(INCLang GE (double(INCLinput2[1])-double(INCLinput2[3])))
        ;tmp2=WHERE(INCLang2 GE (double(INCLinput3[1])-double(INCLinput3[3])))
        IF finishafter EQ 1.1 or fix_incl EQ 0. then begin
           tmp=WHERE(INCLang GE (double(INCLinput1[1])-double(INCLinput1[4])))
           tmp2=WHERE(INCLang2 GE (double(INCLinput1[1])-double(INCLinput1[4])))
        ENDIF ELSE BEGIN
           tmp=WHERE(INCLang GE (double(INCLinput2[1])-double(INCLinput2[4])))
           tmp2=WHERE(INCLang2 GE (double(INCLinput3[1])-double(INCLinput3[4])))
        ENDELSE
        IF tmp[0] EQ n_elements(INCLang)-2 OR tmp[0] EQ n_elements(INCLang)-1 then tmp=[-1]
        IF tmp2[0] EQ n_elements(INCLang2)-2  OR tmp2[0] EQ n_elements(INCLang2)-1 then tmp2=[-1]
        case 1 of
           (n_elements(tmp) GE 2 OR n_elements(tmp) GE 2) AND  (finishafter EQ 1.1 or fix_incl EQ 0.): begin
              IF INCLinput1[1] LT 90. then begin
                 INCLinput1[1]=strtrim(strcompress(string(double(INCLinput1[1]+5.))),2)
                 INCLinput2[1]=INCLinput1[1]
                 INCLinput3[1]=INCLinput1[1]
                 boundaryadjustment=1
              ENDIF
           end
           tmp[0] NE -1 or  tmp2[0] NE -1:begin
              IF INCLinput2[1] LT 90 then begin
                 INCLinput2[1]=strtrim(strcompress(string(double(INCLinput2[1]+5.))),2)
                 boundaryadjustment=1
              ENDIF
              IF INCLinput3[1] LT 90 then begin
                 INCLinput3[1]=strtrim(strcompress(string(double(INCLinput3[1]+5.))),2)
                 boundaryadjustment=1
              ENDIF

           END


           else:
        endcase

        IF boundaryadjustment EQ 1 then begin
           IF tmp[0] EQ -1 then tmp[0]=n_elements(INCLang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(INCLang)
           tmp=MAX([tmp,tmp2],min=lastreliablerings)
        ENDIF
     ENDIF
     IF double(INCLinput1[2]) GT 5 OR double(INCLinput2[2]) GT 5 OR double(INCLinput3[2]) GT 5 AND norings[0] GT 4 then begin
        IF finishafter EQ 1.1 or fix_incl EQ 0. then begin
           tmp=WHERE(INCLang LE (double(INCLinput1[2])+double(INCLinput1[3])))
           tmp2=WHERE(INCLang2 LE (double(INCLinput1[2])+double(INCLinput1[3])))
        ENDIF ELSE BEGIN
           tmp=WHERE(INCLang LE (double(INCLinput2[2])+double(INCLinput2[3])))
           tmp2=WHERE(INCLang2 LE (double(INCLinput3[2])+double(INCLinput3[3])))
        ENDELSE

        IF tmp[0] EQ n_elements(INCLang)-2 OR  tmp[0] EQ n_elements(INCLang)-1 then tmp=[-1]
        IF tmp2[0] EQ n_elements(INCLang2)-2 OR tmp2[0] EQ n_elements(INCLang2)-1 then tmp2=[-1]
        case 1 of
           (n_elements(tmp) GE 2 OR n_elements(tmp) GE 2) AND  (finishafter EQ 1.1 or fix_incl EQ 0.): begin
              IF INCLinput1[2] GT 5. then begin
                 INCLinput1[2]=strtrim(strcompress(string(double(INCLinput1[2]-5.))),2)
                 INCLinput2[2]=INCLinput1[2]
                 INCLinput3[2]=INCLinput1[2]
                 boundaryadjustment=1
              ENDIF
           end
           tmp[0] NE -1 or tmp2[0] NE -1:begin
              IF INCLinput2[2] GT 5 then begin
                 INCLinput2[2]=strtrim(strcompress(string(double(INCLinput2[2]-5.))),2)
                 boundaryadjustment=1
              ENDIF
              IF INCLinput3[2] GT 5 then begin
                 INCLinput3[2]=strtrim(strcompress(string(double(INCLinput3[2]-5.))),2)
                 boundaryadjustment=1
              ENDIF
           END

           else:
        endcase


        IF boundaryadjustment EQ 1 then begin
           IF tmp[0] EQ -1 then tmp[0]=n_elements(INCLang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(INCLang)
           tmp=MAX([tmp,tmp2,lastreliablerings],min=lastreliableringsINCL)
           if lastreliableringsINCL LT lastreliablerings then lastreliablerings=lastreliableringsINCL
        ENDIF
     ENDIF
     IF ABS(double(PAinput1[1])-double(PAinput1[2])) LT 400 OR ABS(double(PAinput2[1])-double(PAinput2[2])) LT 400 OR  ABS(double(PAinput3[1])-double(PAinput3[2])) LT 400 AND norings[0] GT 4 then begin
        IF finishafter EQ 1.1 or fix_pa EQ 0. then begin
           tmp=WHERE(PAang GE (double(PAinput1[1])-double(PAinput1[3])))
           tmp2=WHERE(PAang2 GE (double(PAinput1[1])-double(PAinput1[3])))
        ENDIF ELSE BEGIN
           tmp=WHERE(PAang GE (double(PAinput2[1])-double(PAinput2[3])))
           tmp2=WHERE(PAang2 GE (double(PAinput3[1])-double(PAinput3[3])))
        ENDELSE

        IF tmp[0] EQ n_elements(PAang)-2 OR tmp[0] EQ n_elements(PAang)-1 then tmp=[-1]
        IF tmp2[0] EQ n_elements(PAang2)-2 OR tmp2[0] EQ n_elements(PAang2)-1 then tmp2=[-1]
        case 1 of
           (n_elements(tmp) GE 2 OR n_elements(tmp) GE 2) AND  (finishafter EQ 1.1 or fix_pa EQ 0.): begin
              IF PAinput1[1] LT 400 then begin
                 PAinput1[1]=strtrim(strcompress(string(double(PAinput1[1]+25.))),2)
                 PAinput2[1]=PAinput1[1]
                 PAinput3[1]=PAinput1[1]
                 boundaryadjustment=1
              ENDIF
           end
           tmp[0] NE -1 or tmp2[0] NE -1: begin
              IF PAinput2[1] LT 400 then begin
                 PAinput2[1]=strtrim(strcompress(string(double(PAinput2[1]+25.))),2)
                 boundaryadjustment=1
              ENDIF
              IF PAinput3[1] LT 400 then begin
                 PAinput3[1]=strtrim(strcompress(string(double(PAinput3[1]+25.))),2)
                 boundaryadjustment=1
              ENDIF
           end
           else:
        endcase

        IF finishafter EQ 1.1 or fix_pa EQ 0. then begin
           tmp=WHERE(PAang LE (double(PAinput1[2])+double(PAinput1[3])))
           tmp2=WHERE(PAang2 LE (double(PAinput1[2])+double(PAinput1[3])))
        ENDIF ELSE BEGIN
           tmp=WHERE(PAang LE (double(PAinput2[2])+double(PAinput2[3])))
           tmp2=WHERE(PAang2 LE (double(PAinput3[2])+double(PAinput3[3])))
        ENDELSE
        IF tmp[0] EQ n_elements(PAang)-2 OR  tmp[0] EQ n_elements(PAang)-1 then tmp=[-1]
        IF tmp2[0] EQ n_elements(PAang2)-2 OR tmp2[0] EQ n_elements(PAang2)-1 then tmp2=[-1]
        case 1 of
           (n_elements(tmp) GE 2 OR n_elements(tmp) GE 2) AND  (finishafter EQ 1.1 or fix_pa EQ 0.): begin
              IF PAinput1[2] GT -40 then begin
                 PAinput1[2]=strtrim(strcompress(string(double(PAinput1[2]-25.))),2)
                 PAinput2[2]=PAinput1[2]
                 PAinput3[2]=PAinput1[2]
                 boundaryadjustment=1
              ENDIF
           end
           tmp[0] NE -1 OR tmp2[0] NE -1: begin ;We want to change both limits as warps are predominatly mostly symmetric
              IF PAinput2[2] GT -40 then begin
                 PAinput2[2]=strtrim(strcompress(string(double(PAinput2[2]-25.))),2)
                 boundaryadjustment=1
              ENDIF
              IF PAinput3[2] GT -40 then begin
                 PAinput3[2]=strtrim(strcompress(string(double(PAinput3[2]-25.))),2)
                 boundaryadjustment=1
              ENDIF
           end
           else:
        endcase
        if boundaryadjustment EQ 1 then begin
           IF tmp[0] EQ -1 then tmp[0]=n_elements(PAang)
           IF tmp2[0] EQ -1 then tmp2[0]=n_elements(PAang2)
           tmp=MAX([tmp,tmp2],min=lastreliableringsPA)
           if lastreliableringsPA LT lastreliablerings then lastreliablerings=lastreliableringsPA
        endif
     ENDIF
     INCLangor=INCLang
     INCLang2or=INCLang2
     PAangor=PAang
     PAang2or=PAang2
     IF float(INCLinput1old[1]) NE float(INCLinput1[1]) OR $
        float(INCLinput2old[1]) NE float(INCLinput2[1]) OR $
        float(INCLinput3old[1]) NE float(INCLinput3[1]) OR $
        float(INCLinput1old[2]) NE float(INCLinput1[2]) OR $
        float(INCLinput2old[2]) NE float(INCLinput2[2]) OR $
        float(INCLinput3old[2]) NE float(INCLinput3[2]) OR $
        float(PAinput1old[1]) NE float(PAinput1[1]) OR $
        float(PAinput2old[1]) NE float(PAinput2[1]) OR $
        float(PAinput3old[1]) NE float(PAinput3[1]) OR $
        float(PAinput1old[2]) NE float(PAinput1[2]) OR $
        float(PAinput2old[2]) NE float(PAinput2[2]) OR $
        float(PAinput3old[2]) NE float(PAinput3[2]) THEN boundaryadjustment = 1 ELSE boundaryadjustment = 0


                                ;If we determined that the boundaries should be updated than we need
                                ;to check SBR and reset the values, as
                                ;well as smooth the INCL, PAand VROT
     IF boundaryadjustment then begin
        IF lastreliablerings EQ 0 then lastreliablerings=1
        for j=n_elements(SBRarr)-1,2,-1 do begin
           IF SBRarr[j] LE cutoff[j] then begin
              SBRarr[j]=cutoff[j]*2.
           ENDIF else break
        endfor

        IF lastreliablerings LT n_elements(SBRarr)-1 then begin
           ;It is also important to give the outer rings a lot of flux
           SBRarr[lastreliablerings:n_elements(SBRarr)-1]=cutoff[lastreliablerings:n_elements(SBRarr)-1]*6.
           INCLang[*]=MEAN(INCLang)
           PAang[lastreliablerings:n_elements(SBRarr)-1]=PAang[lastreliablerings]
           IF lastreliablerings EQ 1 then begin
              VROTarr[lastreliablerings:n_elements(SBRarr)-1]=MEAN(VROTarr[1:n_elements(VROTarr)-1])
              SDISarr[lastreliablerings:n_elements(SBRarr)-1]=MEAN(SDISarr[1:n_elements(SDISarr)-1])
           ENDIF else BEGIN
              VROTarr[lastreliablerings:n_elements(SBRarr)-1]=VROTarr[lastreliablerings]
              SDISarr[lastreliablerings:n_elements(SDISarr)-1]=SDISarr[lastreliablerings]
           ENDELSE
           comin=[[PAang],[INCLang]]
           errors=[[0.],[0.]]

           revised_regularisation_com,comin,SBRarror,RADarr,fixedrings=3,difference=[1.,4.*exp(-catinc[i]^2.5/10^3.5)+1.5],cutoff=cutoff,arctan=1,order=polorder,max_par=[PAinput2[1],INCLinput2[1]],min_par=[PAinput2[2],INCLinput2[2]],accuracy=[1./4.,1.],error=errors ,gdlidl=gdlidl,log=log ,ring_spacing=ring_spacing
           PAang=comin[*,0]
           INCLang=comin[*,1]
           sigmapa1=errors[*,0]
           sigmaincl1=errors[*,1]
           stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
           tmppos=where('INCL' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
           stringINCL='PA= '+STRJOIN(string(PAang),' ')
           tmppos=where('PA' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
        ENDIF
        for j=n_elements(SBRarr2)-1,2,-1 do begin
           IF SBRarr2[j] LE cutoff[j] then begin
              SBRarr2[j]=cutoff[j]*2.
           ENDIF else break
        endfor

        IF lastreliablerings LT n_elements(SBRarr)-1 then begin
             ;It is also important to give the outer rings a lot of flux
           SBRarr2[lastreliablerings:n_elements(SBRarr2)-1]=cutoff[lastreliablerings:n_elements(SBRarr2)-1]*6.
                                ;INCLang2[lastreliablerings:n_elements(SBRarr2)-1]=INCLang2[lastreliablerings]
           INCLang2[*]=MEAN(INCLang2)
           PAang2[lastreliablerings:n_elements(SBRarr2)-1]=PAang2[lastreliablerings]
           IF lastreliablerings EQ 1 then BEGIN
              VROTarr2[lastreliablerings:n_elements(SBRarr)-1]=MEAN(VROTarr2[1:n_elements(VROTarr2)-1])
              SDISarr2[lastreliablerings:n_elements(SBRarr)-1]=MEAN(SDISarr2[1:n_elements(SDISarr2)-1])
           ENDIF else BEGIN
              VROTarr2[lastreliablerings:n_elements(SBRarr2)-1]=VROTarr2[lastreliablerings]
              SDISarr2[lastreliablerings:n_elements(SBRarr2)-1]=SDISarr2[lastreliablerings]
           ENDELSE
            comin=[[PAang2],[INCLang2]]
           errors=[[0.],[0.]]
           revised_regularisation_com,comin,SBRarr2or,RADarr,fixedrings=3,difference=[1.,4.*exp(-catinc[i]^2.5/10^3.5)+1.5],cutoff=cutoff,arctan=1,order=polorder,max_par=[PAinput3[1],INCLinput3[1]],min_par=[PAinput3[2],INCLinput3[2]],accuracy=[1./4.,1.],error=errors ,gdlidl=gdlidl,log=log,ring_spacing=ring_spacing
           PAang2=comin[*,0]
           INCLang2=comin[*,1]
           sigmapa2=errors[*,0]
           sigmainc2=errors[*,1]
           stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
           tmppos=where('INCL_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
           stringINCL='PA_2= '+STRJOIN(string(PAang2),' ')
           tmppos=where('PA_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
        ENDIF
        VROTarr=(VROTarr+VROTarr2)/2.
        tmpSBR=(SBRarror+SBRarr2or)/2.
        verror=MAX([5.,channelwidth/2.*(1.+vresolution)/SQRT(sin(catinc[i]*!DtoR))])
        IF double(VROTarr[1]) LT 120. AND  double(VROTarr[2]) LT 120. then begin

           revised_regularisation_rot,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=1,order=polorder,error=sigmarot,log=log,max_par=VROTmax,min_par=channelwidth,ring_spacing=ring_spacing
        ENDIF ELSE BEGIN
           revised_regularisation_rot,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=1,/NOCENTRAL,error=sigmarot,log=log,max_par=VROTmax,min_par=channelwidth,ring_spacing=ring_spacing
        ENDELSE

        stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
        tmppos=where('VROT' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringVROT
        stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
        tmppos=where('VROT_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringVROT
                                ;If we do this we also need to check the vrot min
        if VROTarr[lastreliablerings-1] LT avinner then vrotinput1[2]=string(VROTmin)+' '+string(VROTarr[lastreliablerings-1]*0.6)
                                ;SDIS
        SDISarr=(SDISarr+SDISarr2)/2.
        if norings[0] LT 4 OR finishafter EQ 1.1 then begin
           SDISarr[*]=MEAN(SDISarr)
           sigmasdis=replicate(channelwidth/4.,n_elements(SDISarr))
        endif else begin
           regularisation_sdis,SDISarr,tmpSBR,RADarr,log=log,max_par=SDISmax,min_par=SDISmin,arctan=1,fixedrings=velfixrings,difference=channelwidth/4.,cutoff=cutoff,gdlidl=gdlidl ,error= sigmasdis
        endelse
        stringSDIS='SDIS= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
        tmppos=where('SDIS' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSDIS
        stringSDIS='SDIS_2= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
        tmppos=where('SDIS_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSDIS
        SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.
        SBRarr2[0:1]=SBRarr[0:1]
        stringSBR='SBR= '+STRJOIN(string(SBRarr),' ')
        tmppos=where('SBR' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        stringSBR='SBR_2= '+STRJOIN(string(SBRarr2),' ')
        tmppos=where('SBR_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
                                ;Write the fitting variables to files
        case 1 of
           norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
              Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           fix_incl EQ 0. and fix_pa EQ 1.: begin
              Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           fix_incl EQ 1. and fix_pa EQ 0.: begin
              Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,$
                                 vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
           else: begin
              Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                    vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                    xposinput1,yposinput1,vsysinput1,INIMODE=inimode
           end
        endcase



        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We changed the boundaries. These are the new boundaries: "
           printf,66,linenumber()+"PA upper limits="+PAinput1[1]+" "+PAinput2[1]+" "+PAinput3[1]+" "
           printf,66,linenumber()+"PA lower limits="+PAinput1[2]+" "+PAinput2[2]+" "+PAinput3[2]+" "
           printf,66,linenumber()+"INCL upper limits="+INCLinput1[1]+" "+INCLinput2[1]+" "+INCLinput3[1]+" "
           printf,66,linenumber()+"INCL lower limits="+INCLinput1[2]+" "+INCLinput2[2]+" "+INCLinput3[2]+" "
           printf,66,linenumber()+"Wit the the following PAs "
           printf,66,linenumber()+STRJOIN(strtrim(string(PAangor[*]),2),' ')
           printf,66,linenumber()+STRJOIN(strtrim(string(PAang2or[*]),2),' ')
           printf,66,linenumber()+"And the following Inclinations "
           printf,66,linenumber()+STRJOIN(strtrim(string(INCLangor[*]),2),' ')
           printf,66,linenumber()+STRJOIN(strtrim(string(INCLang2or[*]),2),' ')

           Close,66
        ENDIF
                                ;goto the point where the the second
                                ;fit is not accepted and try again
        goto,notacceptedtwo
     ENDIF
                                ;If we are smoothing with polynomials
                                ;we want to make sure that the minimum
                                ;sbr is actually of some value for the
                                ;weighing in the smoothing.
     IF finalsmooth EQ 1 then begin
        for j=n_elements(SBRarr)-1,2,-1 do begin
           IF SBRarr[j] LE cutoff[j] then begin
              SBRarr[j]=cutoff[j]*0.5
           ENDIF else break
        endfor
        tmp=WHERE(SBRarr LT 0.)
        IF tmp[0] NE -1 then SBRarr[tmp]=cutoff[tmp]*0.95
        for j=n_elements(SBRarr2)-1,2,-1 do begin
           IF SBRarr2[j] LE cutoff[j] then begin
              SBRarr2[j]=cutoff[j]*0.5
           ENDIF else break
        endfor
        tmp=WHERE(SBRarr2 LT 0.)
        IF tmp[0] NE -1 then SBRarr2[tmp]=cutoff[tmp]*0.95
     ENDIF
                                ;We will correct the output by fitting a polynomial
                                ;first we set the amount of rings that
                                ;are clearly one value

     IF finalsmooth EQ 1 then begin
        get_fixedringsv9,[[PAang],[PAang2],[INCLang],[INCLang2]],fixedrings,/smooth
     ENDIF ELSE get_fixedringsv9,[[PAang],[PAang2],[INCLang],[INCLang2]],fixedrings
     IF fixedrings GE norings[0] then begin
        IF norings[0] GT 5 then fixedrings=norings[0]-2 else fixedrings=norings[0]-1
     ENDIF
                                ;if we are not smoothing with a
                                ;polynomial we still want to simply
                                ;smooth the profile by averaging such
                                ;that the jigsaw pattern
                                ;doesn't amplify every iteration
     IF finalsmooth EQ 1 then prefunc=0 else prefunc=1
     centralincl=(INCLang[0]+INCLang2[0])/2.
     IF norings[0]*ring_spacing GT 15 then accuracy=0.5+COS(centralincl*!DtoR)*0.5*15./(norings[0]*ring_spacing) else accuracy=1

                                ;Smoothing PA_1
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin
        comin=[[PAang],[INCLang]]
        errors=[[0.],[0.]]
        padiv=1.5-(ATAN((norings[0]*ring_spacing-5.)*2.))*5.25/!pi
        IF padiv LT 1.0 then padiv=1.0
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We are Fitting the first PA and INCL"
           printf,66,linenumber()+"The PA profile is forced to remain between "+PAinput2[1]+" and "+PAinput2[2]
           printf,66,linenumber()+"The INC profile is forced to remain between "+INCLinput2[1]+" and "+INCLinput2[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
           Close,66
        ENDIF
        IF keyword_set(debug) then print,linenumber()+'making sure it is here before'
        if fix_pa EQ 1 or fix_incl EQ 1 then begin
           revised_regularisation_com,comin,SBRarror,RADarr,fixedrings=fixedrings,difference=[padiv,4.*exp(-catinc[i]^2.5/10^3.5)+1.5],cutoff=cutoff,arctan=prefunc,order=polorder1,max_par=[PAinput2[1],INCLinput2[1]],min_par=[PAinput2[2],INCLinput2[2]],accuracy=[accuracy/4.,accuracy],error=errors ,gdlidl=gdlidl,log=log,sloped=prevslopedrings[0],ring_spacing=ring_spacing
        ENDIF
        IF keyword_set(debug) then print,linenumber()+'making sure it is here after'
        IF fix_pa EQ 1 then PAang=comin[*,0]
        IF fix_incl EQ 1 then INCLang=comin[*,1]
        sigmapa1=errors[*,0]
        sigmaincl1=errors[*,1]

     ENDIF
     padiff=PAang[*]-PAang[0]
                                ;Smoothing PA_2
     if finalsmooth LE 1 AND norings[0] GT 4 and finishafter NE 1.1 then begin
        comin=[[PAang2],[INCLang2]]
        errors=[[0.],[0.]]
        padiv=1.5-(ATAN((norings[0]*ring_spacing-5.)*2.))*5.25/!pi
        IF padiv LT 1.0 then padiv=1.0
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We are Fitting the second PA and INCL"
           printf,66,linenumber()+"The PA profile is forced to remain between "+PAinput3[1]+" and "+PAinput3[2]
           printf,66,linenumber()+"The INCL profile is forced to remain between "+INCLinput3[1]+" and "+INCLinput3[2]
           printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
           Close,66
        ENDIF
        if fix_pa EQ 1 or fix_incl EQ 1 then begin
           revised_regularisation_com,comin,SBRarr2or,RADarr,fixedrings=fixedrings,difference=[padiv,4.*exp(-catinc[i]^2.5/10^3.5)+1.5],cutoff=cutoff,arctan=prefunc,order=polorder2,max_par=[PAinput3[1],INCLinput3[1]],min_par=[PAinput3[2],INCLinput3[2]],accuracy=[accuracy/4.,accuracy],error=errors ,gdlidl=gdlidl,log=log,sloped=prevslopedrings[1],ring_spacing=ring_spacing
        endif
        if fix_pa EQ 1 then PAang2=comin[*,0]
        if fix_incl EQ 1 then INCLang2=comin[*,1]
        sigmapa2=errors[*,0]
        sigmaincl2=errors[*,1]


     ENDIF
                                ;making sure that the rings that are
                                ;fixed are the average of both sides
                                ;and the same
     IF polorder1[0] EQ 0 then begin
        PAang[*]=(double(PAang[0])+double(PAang2[0]))/2.
     ENDIF ELSE PAang[0:fixedrings-1]=(double(PAang[0])+double(PAang2[0]))/2.
     IF polorder2[0] EQ 0 then begin
        PAang2[*]=PAang[0]
     ENDIF ELSE PAang2[0:fixedrings-1]=PAang[0]
                                ;And writing them to the file template
     stringPA='PA= '+STRJOIN(string(PAang),' ')
     tmppos=where('PA' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringPA
     stringPA='PA_2= '+STRJOIN(string(PAang2),' ')
     tmppos=where('PA_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringPA

                                ;making sure that the rings that are
                                ;fixed are the average of both sides
                                ;and the same

     IF polorder1[1] EQ 0 then begin
        INCLang[*]=(double(INCLang[0])+double(INCLang2[0]))/2.
     ENDIF ELSE INCLang[0:fixedrings-1]=(double(INCLang[0])+double(INCLang2[0]))/2.
     IF polorder2[1] EQ 0 then begin
        INCLang2[*]=INCLang[0]
     ENDIF ELSE INCLang2[0:fixedrings-1]=INCLang[0]


     stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
     tmppos=where('INCL' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringINCL
     stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
     tmppos=where('INCL_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringINCL
                                ;rotation curve
                                ;THis is a bit more tricky since we
                                ;want the flat parts on the outside we
                                ;need to inverse them
                                ;First let's see if the size of the
                                ;galaxy is reasonable
                                ;As we fix the rotation curve to be
                                ;the same for both sides we want to
                                ;use an average surface brightness
                                ;profile
                                ;Not dividing the profile by two seems
                                ;like a mistake but as this
                                ;combination gives good results
                                ;let's leave it like this.
     SBRav=(SBRarr+SBRarr2)/2.
     get_newringsv9,SBRav,SBRav,2.5*cutoff,velconstused
     velconstused--
     IF double(norings[0])*ring_spacing GT 15. then begin
        IF double(norings[0])*ring_spacing LT 25. then fact = 10 else fact= 20-(norings[0]*ring_spacing/2.5)
        if fact LT 6 then fact=6

        IF velconstused GT norings[0]-ceil(norings[0]/fact) then velconstused=norings[0]-ceil(norings[0]/fact)

     ENDIF
     prefunc=0.
     IF norings[0]-velconstused LT 2 then velconstused=norings[0]-1
     IF norings[0]*ring_spacing GT 8 AND not finishafter EQ 2.1 then velconstused=velconstused-1


     if finalsmooth LE 1 then begin

        IF finalsmooth EQ 1 AND norings[0] GT 4 then prefunc=0 else prefunc=1
                                ; IF finalsmooth EQ 1 then we first
                                ; want to refit the rotation curve
                                ; with the smoothed inclination as it
                                ; often drives eachother
        IF finalsmooth EQ 1 then begin
           tmppos=where('VARINDX' EQ tirificsecondvars)
           tmp=strsplit(tirificsecond[tmppos],': ',/extract)
           IF n_elements(tmp) GT 3 then begin
                                ;IF we have a warp slope the location
                                ;of VROT is different hence we need to
                                ;look for the VROT keyword.
              locv=WHERE('VROT' EQ strtrim(STRUPCASE(tmp),2))
              IF n_elements(tmp)-1 GT locv[0]+2 then begin
                 IF isnumeric(tmp[locv[0]+2]) then velconstused=tmp[locv[0]+2]-1 else velconstused=norings[0]
              ENDIF else velconstused=norings[0]
           ENDIF else velconstused=norings[0]
        ENDIF

        ;Check that or fixedrings are larger then the flattened outer part
        IF norings[0]-xind GT velconstused then velconstused=norings[0]-xind


        velfixrings=norings[0]-velconstused
        IF velfixrings LT 1 then velfixrings=1
        IF velfixrings EQ 1 AND norings[0] GT 5 then velfixrings=2


        if locmax[n_elements(locmax)-1] LT n_elements(VROTarr)/2.  AND MEAN(VROTarr[fix(n_elements(VROTarr)/2.):n_elements(VROTarr)-1]) GT 120. then begin
            IF size(log,/TYPE) EQ 7 then begin
               openu,66,log,/APPEND
               printf,66,linenumber()+"This curve is declining so we fit a flat outer part"
               close,66
            ENDIF
           slope=1
           IF velfixrings GT 1 then VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1]=TOTAL(VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1])/n_elements(VROTarr[n_elements(VROTarr)-velfixrings:n_elements(VROTarr)-1])
           stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
           tmppos=where('VROT' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
           stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
           tmppos=where('VROT_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
        ENDIF
        tmpSBR=(SBRarror+SBRarr2or)/2.
        SDISarr=(SDISarr+SDISarr2)/2.
        if norings[0] LT 4*ring_spacing OR finishafter EQ 1.1 then begin
           SDISarr[*]=MEAN(SDISarr)
           sigmasdis=replicate(channelwidth/4.,n_elements(SDISarr))
        endif else begin
           regularisation_sdis,SDISarr,tmpSBR,RADarr,log=log,max_par=SDISmax,min_par=SDISmin,arctan=prefunc,fixedrings=velfixrings,difference=channelwidth/4.,cutoff=cutoff,gdlidl=gdlidl ,error= sigmasdis
        endelse
        stringSDIS='SDIS= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
        tmppos=where('SDIS' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSDIS
        stringSDIS='SDIS_2= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
        tmppos=where('SDIS_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSDIS

                                ;it is important to use the original
                                ;non centrallly modified profiles here
                                ;for the error weighing. For PA and
                                ;INCL it doesn't matter as
                                ;those errors are set to 1

        tmpSBR=(SBRarror+SBRarr2or)/2.
        tmphigh=str_sep(strtrim(strcompress(VROTinput1[1]),2),' ')
                                ;This way things would be set to
                                ;avinner which would heavily penalize
                                ;slowly rising curves.
                                ;tmplow=str_sep(strtrim(strcompress(VROTinput1[2]),2),' ')

        tmplow=vrotmin
        vmaxdev=MAX([30.,7.5*channelwidth*(1.+vresolution)])
        verror=MAX([4.,channelwidth/2.*(1.+vresolution)/SQRT(sin(INCLang[0]*!DtoR))])
        IF ~(FINITE(verror)) then verror=channelwidth/2.*(1.+vresolution)/SQRT(sin(5.*!DtoR))
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is VROT."
           printf,66,VROTarr
           printf,66,linenumber()+"This is SBR."
           printf,66,tmpSBR
           printf,66,linenumber()+"The number of rings= "+strtrim(string(fix(norings[0])),2)
           printf,66,linenumber()+"The maximum/minimum= "+strtrim(string(tmphigh[0]),2)+"/"+strtrim(string(tmplow[0]),2)
           printf,66,linenumber()+"The number of fixed rings= "+strtrim(string(fix(velfixrings)),2)
           printf,66,linenumber()+"The velocity error = "+strtrim(string(verror),2)
           Close,66
        ENDIF
                                ; if we do not have more freering than
                                ; max order we only want to smooth
                                ;if n_elements(VROTarr)-velfixrings LT
                                ;8 AND SUM([polorder1,polorder2]) EQ
                                ;0. then velprefunc=0 else
                                ;velprefunc=prefunc
        ;if n_elements(VROTarr)-velfixrings LT 8 then velprefunc=1 else velprefunc=prefunc
                                ;If the central parameters are LT 120
                                ;we want to take the central ring into
                                ;account otherwise not

        velprefunc=prefunc
        IF (double(VROTarr[1]) LT 120. AND  double(VROTarr[2]) LT 120.) then begin
           revised_regularisation_rot,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=velprefunc,order=polorder,max_par=tmphigh[0],min_par=tmplow[0],accuracy=accuracy,error=sigmarot,log=log,ring_spacing=ring_spacing
        ENDIF ELSE BEGIN
           revised_regularisation_rot,VROTarr,tmpSBR, RADarr,/REVERSE,fixedrings=velfixrings,difference=verror,cutoff=cutoff,arctan=velprefunc,/NOCENTRAL,order=polorder,max_par=tmphigh[0],min_par=tmplow[0],accuracy=accuracy,error=sigmarot,log=log,ring_spacing=ring_spacing
        ENDELSE
        tmp0check=WHERE(VROTarr LT 0)
        IF tmp0check[0] NE -1 then VROTarr[tmp0check]=tmplow[0]
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"This is VROT after smoothing."
           printf,66,VROTarr
           Close,66
        ENDIF
        VROTarr[0]=0.
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We used a polynomial of the "+string(polorder)+" order to smooth the VROT profile."
           printf,66,linenumber()+"And fixed "+string(velfixrings)+" rings."
           Close,66
        ENDIF
         ;SDIS
       ; SDISarr=(SDISarr+SDISarr2)/2.
       ; if norings[0] LT 4 OR finishafter EQ 1.1 then begin
       ;    SDISarr[*]=MEAN(SDISarr)
       ;    sigmasdis=replicate(channelwidth/4.,n_elements(SDISarr))
       ; endif else begin
       ;    regularisation_sdis,SDISarr,tmpSBR,RADarr,log=log,max_par=SDISmax,min_par=SDISmin,arctan=prefunc,fixedrings=velfixrings,difference=channelwidth/4.,cutoff=cutoff,gdlidl=gdlidl ,error= sigmasdis
       ; endelse


     ENDIF
                                ;As the SBR used is the average and
                                ;the rotation curve is the same for
                                ;both side we only need to fit one
                                ;side and can copy the
                                ;other. Velconstused also determines
                                ;the amount of rings only fitted in a slope
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We fit a slope to the rotation curve from ring "+string(velconstused)+" to ring"+strtrim(string(fix(norings[0])),2)
        for x=0,n_elements(SBRarr)-1 do begin
           printf,66,SBRarr[x],SBRarr2[x],2.*cutoff[x]
        endfor
        Close,66
     ENDIF ELSE BEGIN
        print,linenumber()+"We fit a slope to the rotation curve from ring "+string(velconstused)+" to ring"+strtrim(string(fix(norings[0])),2)
     ENDELSE
     VROTarr[0]=0.
     VROTarr2=VROTarr
     IF velconstused LT norings[0]-1 AND finalsmooth LT 1 then begin
        VROTarr[velconstused-1:norings[0]-1]=VROTarr[velconstused-2]
        VROTarr2=VROTarr
     ENDIF
     newindrings=[n_elements(SBRarr),n_elements(SBRarr2)]

                                ;IF we are on the last smoothing we
                                ;want to set the rings in the SBR
                                ;profile that are below the cutoff to
                                ;a very small number so they do not
                                ;affect the fit
     IF finalsmooth GE 1 then begin
        get_newringsv9,SBRarr,SBRarr2,cutoff,newindrings,/individual
        IF newindrings[0] EQ n_elements(SBRarr)-2 and SBRarr[n_elements(SBRarr)-1] GT cutoff[n_elements(SBRarr)-1] then newindrings[0]=n_elements(SBRarr)
        IF newindrings[0] LT n_elements(SBRarr) then begin
           SBRarr[newindrings[0]-1:n_elements(SBRarr)-1]=1E-16
        ENDIF
        IF newindrings[1] EQ n_elements(SBRarr2)-2 and SBRarr2[n_elements(SBRarr2)-1] GT cutoff[n_elements(SBRarr2)-1] then newindrings[1]=n_elements(SBRarr2)
        IF newindrings[1] LT n_elements(SBRarr2) then begin
           SBRarr2[newindrings[1]-1:n_elements(SBRarr2)-1]=1E-16
        ENDIF
        SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.
        SBRarr2[0:1]=SBRarr[0:1]
        stringSBR='SBR= '+STRJOIN(string(SBRarr),' ')
        tmppos=where('SBR' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        stringSBR='SBR_2= '+STRJOIN(string(SBRarr2),' ')
        tmppos=where('SBR_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"The last SBR 1 is "+string(STRJOIN(string(SBRarr),' '))
           printf,66,linenumber()+"The last SBR 2 is "+string(STRJOIN(string(SBRarr2),' '))
           close,66
        ENDIF
     ENDIF
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"The last SBR ring 1 is "+string(newindrings[0])
        printf,66,linenumber()+"The last SBR ring 2 is "+string(newindrings[1])
        printf,66,linenumber()+"The amount of rings is "+strtrim(string(fix(norings[0])),2)
        close,66
     ENDIF
     MAXindrings=MAX(newindrings)
     IF MAXindrings LT n_elements(VROTarr) then begin
        VROTarr[MAXindrings:n_elements(VROTarr)-1]=VROTarr[MAXindrings-1]
     ENDIF
    ;Write all to the template
     stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
     tmppos=where('VROT' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringVROT
     stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')
     tmppos=where('VROT_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringVROT

     stringSDIS='SDIS= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
     tmppos=where('SDIS' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringSDIS
     stringSDIS='SDIS_2= '+STRJOIN(string(SDISarr[0:n_elements(SDISarr)-1]),' ')
     tmppos=where('SDIS_2' EQ tirificsecondvars)
     tirificsecond[tmppos]=stringSDIS


                                ;IF we have not done too many
                                ;iterations and we are not on the
                                ;final smoothing than we want to check
                                ;the extend of the models
     IF sbrmodify LT 10. AND finalsmooth LT 1 AND secondtime NE 1  then begin
                                ;check the SBR profiles to see if we
                                ;need to cut/add a ring.
        get_newringsv9,SBRarrunmod,SBRarr2unmod,cutoff,newrings
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Newrings should always be bigger or equal to velconstused "+strtrim(string(fix(newrings)),2)+" =>"+string(velconstused)
           close,66
        ENDIF
                                ;Checking the new amount of rings against various criteria
                                ;Not too small
        if newrings LT 4 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to 4."
              close,66
           ENDIF
           newrings=4
        ENDIF
                                ;not too big
        IF newrings GT maxrings then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" we set it to maxrings."
              close,66
           ENDIF
           newrings=maxrings
        ENDIF
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We found this as rings = "+strtrim(string(fix(norings[0])),2)+"  new="+strtrim(string(fix(newrings)),2)
           close,66
        ENDIF
      ;  IF NOT sofiafail then begin
        IF newrings GT sofiarings+(3/ring_spacing) OR newrings LT sofiarings-(3/ring_spacing) then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"The new amount of rings ("+strtrim(string(fix(newrings)),2)+") deviates too much from the sofia estimate ("+string(sofiarings)+")."
              close,66
           ENDIF
           IF newrings LT norings[0] then begin
              IF norings[0] GT sofiarings-round(3/ring_spacing) then newrings=norings[0]-1 else newrings=sofiarings-round(3/ring_spacing)
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We set the number of rings to "+strtrim(string(fix(newrings)),2)
                 close,66
              ENDIF
           endif Else begin
              IF norings[0]*ring_spacing LT sofiarings+round(3/ring_spacing) then newrings=norings[0]+1 else newrings=sofiarings+round(3/ring_spacing)
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We set the number of rings to "+strtrim(string(fix(newrings)),2)
                 close,66
              ENDIF
           ENDELSE

        ENDIF
        ;ENDIF

                                ;Not previously fitted. If it has been
                                ;than fix to this number of rings.
        tmp=WHERE(prevrings EQ newrings)
        IF tmp[0] NE -1 AND newrings NE norings[0] then begin
           IF finalsmooth EQ 2 then begin
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before and we are smoothing we are going to fix the rings at "+strtrim(string(fix(norings[0])),2)
                 close,66
              ENDIF
              newrings=norings[0]
           ENDIF ELSE BEGIN
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Cause the newrings was "+strtrim(string(fix(newrings)),2)+" which we have fitted before we gonna fix the rings at this value."
                 close,66
              ENDIF
           ENDELSE
           secondtime=1
        ENDIF
                                ;IF we want to cut a ring we check
                                ;what we did previously and whether it
                                ;is acceptable.
        IF newrings LT norings[0] AND constring NE 1 AND newrings GE lowring then begin
           IF newrings LE forcedring then begin
              IF norings[0] GT forcedring then begin
                 newrings=forcedring
                 constring=1
                 prevmodification=0.
              ENDIF ELSE BEGIN
                 IF size(log,/TYPE) EQ 7 then begin
                    openu,66,log,/APPEND
                    printf,66,linenumber()+"As the cut would modify the rings to less than the maximum forced addition of a ring we do not apply the cut."
                    printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)+" forcedring="+string(forcedring)
                    Close,66
                 ENDIF
                 goto,nocut
              ENDELSE
           ENDIF
           IF prevmodification EQ 1 then begin
              IF newrings GE norings[0]-1 then begin
                 prevmodification=-1
              ENDIF else BEGIN
                 prevmodification=-2
              ENDELSE
           ENDIF
           IF finalsmooth EQ 2 then begin
              newrings=norings[0]-1
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"As we are smoothing we will only subtract one ring."
                 printf,66,linenumber()+"new rings "+strtrim(string(fix(newrings)),2)+" old rings "+strtrim(string(fix(norings[0])),2)
                 Close,66
              ENDIF
           ENDIF
           oldrings=norings[0]
           IF newrings EQ norings[0] then goto,nocut else norings[0]=newrings

           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
              Close,66
           ENDIF ELSE begin
              print,linenumber()+"We cut a ring! lastcut"+string(lastcutrings)+" lastadd "+string(lastaddrings)+" new "+strtrim(string(fix(newrings)),2)+" old "+string(oldrings)
           ENDELSE
                                ;change the template to have a ring less
           changeradii,tirificsecond,norings[0],continue_tirific
           lastcutrings=norings[0]
                                ;Set the fitting parameters
           set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
           set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,finish_after=finishafter,slope=slope,fix=fix_sdis

                                ;Adapt the other fitting parameters
                                ;first the SBRarr
           set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log,doubled=doubled
           IF norings[0] LE 4 OR finishafter EQ 1.1 then begin
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)

           ENDIF ELSE BEGIN
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;update INCL and PA
              set_warp_slopev3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix,sloped=slopedrings
              if fix_pa EQ 0. then   PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                                 ' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              if fix_incl EQ 0. then INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                                                   ' INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
                                ;And the central parameters
              z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                          ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                            ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)

           ENDELSE
           case 1 of
              norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
                 Writefittingvariables,tirificsecond,inclinput1,painput1,$
                                       vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                       xposinput1,yposinput1,vsysinput1,INIMODE=inimode

              end
              fix_incl EQ 0. and fix_pa EQ 1.: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE
              end
              fix_incl EQ 1. and fix_pa EQ 0.: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE

              end
              else: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE

              end
           endcase



                                ;Update counters and go to a new
                                ;iteration of the model.
           sbrmodify++
           finalsmooth=0.
           overwrite=0.
           prevrings=[prevrings,norings[0]]
           goto,notacceptedtwo
        ENDIF
        nocut:
        overwrite=0
                                ;Under certain conditions we want to
                                ;add a previously cut ring back on
        IF newrings GT norings[0] then begin
           IF (SBRarr[newrings-2] GT 7.5*cutoff[newrings-2] OR SBRarr2[newrings-2] GT 7.5*cutoff[newrings-2] OR $
               (SBRarr[newrings-3] GT 5*cutoff[newrings-3] AND SBRarr[newrings-2] GT 3.*SBRarr[newrings-3] ) OR $
               (SBRarr2[newrings-3] GT 5*cutoff[newrings-3] AND SBRarr2[newrings-2] GT 3.*SBRarr2[newrings-3] )) $
               AND prevmodification EQ -1 then begin
              prevmodification=0.
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"We force added a ring to the profile."
                 printf,66,linenumber()+"Newrings "+strtrim(string(fix(newrings)),2)+", Rings "+strtrim(string(fix(norings[0])),2)
                 Close,66
              ENDIF
              IF newrings GT forcedring then forcedring=newrings
              overwrite=1
           ENDIF
        ENDIF

                                ;If we do not have those conditions
                                ;and the previous iteration was a cut
                                ;we do not add a ring
        IF prevmodification EQ -1 AND overwrite NE 1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We wanted to add a ring to the model but the previous step was a subtraction."
              Close,66
           ENDIF
        ENDIF
                                ;If all conditions for adding a ring
                                ;are satified then update the tirific file
        IF newrings GT norings[0] AND (prevmodification NE -1 OR overwrite EQ 1) AND newrings LT highring then begin
           norings[0]=newrings
           addring:
           changeradii,tirificsecond,norings[0],continue_tirific
                                ;we need to estimate some values for the added rings
           prefunc=0.
           tmppos=where('RADI' EQ tirificsecondvars)
           tmp=str_sep(strtrim(strcompress(tirificsecond[tmppos]),2),'=')
           RADarr=double(str_sep(strtrim(strcompress(tmp[1]),2),' '))

           tmpsbr=[SBRarror,cutoff[n_elements(SBRarr)]]
                                ;PA_1
           comin=[[PAang,PAang[n_elements(PAang)-1]],[INCLang,INCLang[n_elements(INCLang)-1]]]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We are extending the first PA and INCL"
              printf,66,linenumber()+"The profile is forced to remain between "+PAinput2[1]+" and "+PAinput2[2]
              printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
              Close,66
           ENDIF
           if fix_pa eq 1 or fix_incl EQ 1 then begin
              revised_regularisation_com,comin,tmpsbr,RADarr,fixedrings=fixedrings,difference=[2.,6.*exp(-catinc[i]^2.5/10^3.5)+4.],cutoff=cutoff,arctan=prefunc,order=polorder,max_par=[PAinput2[1],INCLinput2[1]],min_par=[PAinput2[2],INCLinput2[2]],gdlidl=gdlidl,log=log,/extending,ring_spacing=ring_spacing
           endif
           tmp2=PAang
           PAang=dblarr(n_elements(RADarr))
           PAang[0:n_elements(PAang)-2]=tmp2
           IF fix_pa EQ 1 then PAang[n_elements(PAang)-1]=comin[n_elements(comin[*,0])-1,0] else  PAang[n_elements(PAang)-1]=PAang[0]

           tmp2=INCLang
           INCLang=dblarr(n_elements(RADarr))
           INCLang[0:n_elements(INCLang)-2]=tmp2
           IF fix_incl EQ 1 then INCLang[n_elements(INCLang)-1]=comin[n_elements(comin[*,1])-1,1] else  INCLang[n_elements(INCLang)-1]=INCLang[0]
                                ;PA_2

           tmpsbr2=[SBRarr2or,cutoff[n_elements(SBRarr2)]]
           comin=[[PAang2,PAang2[n_elements(PAang2)-1]],[INCLang2,INCLang2[n_elements(INCLang2)-1]]]
            IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We are extending the second PA and INCL"
              printf,66,linenumber()+"The profile is forced to remain between "+PAinput3[1]+" and "+PAinput3[2]
              printf,66,linenumber()+"The inner" +string(fixedrings)+" rings are fixed."
              Close,66
           ENDIF
           if fix_pa eq 1 or fix_incl EQ 1 then begin
              revised_regularisation_com,comin,tmpsbr2,RADarr,fixedrings=fixedrings,difference=[2.,6.*exp(-catinc[i]^2.5/10^3.5)+4.],cutoff=cutoff,arctan=prefunc,order=polorder2,max_par=[PAinput3[1],INCLinput3[1]],min_par=[PAinput3[2],INCLinput3[2]],gdlidl=gdlidl,log=log,/extending,ring_spacing=ring_spacing
           endif
           tmp2=PAang2
           PAang2=dblarr(n_elements(RADarr))
           PAang2[0:n_elements(PAang2)-2]=tmp2
           IF fix_pa EQ 1 then PAang2[n_elements(PAang2)-1]=comin[n_elements(comin[*,0])-1,0] else  PAang2[n_elements(PAang)-1]=PAang2[0]
           tmp2=INCLang2
           INCLang2=dblarr(n_elements(RADarr))
           INCLang2[0:n_elements(INCLang2)-2]=tmp2
           IF fix_incl EQ 1 then INCLang2[n_elements(INCLang2)-1]=comin[n_elements(comin[*,1])-1,1] else INCLang2[n_elements(INCLang)-1]=INCLang2[0]

           PAang[0:fixedrings]=(PAang[0]+PAang2[0])/2.
           PAang2[0:fixedrings]=(PAang[0]+PAang2[0])/2.
           stringPA='PA= '+STRJOIN(string(PAang),' ')
           tmppos=where('PA' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringPA
           stringPA='PA_2= '+STRJOIN(string(PAang2),' ')
           tmppos=where('PA_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringPA
                                ;inclination


           INCLang[0:fixedrings]=(INCLang[0]+INCLang2[0])/2.
           INCLang2[0:fixedrings]=(INCLang[0]+INCLang2[0])/2.
           tmp=WHERE(Inclang GT 90.)
           IF tmp[0] NE -1 then INCLang[tmp]=90.
           tmp=WHERE(Inclang LT 0.)
           IF tmp[0] NE -1 then INCLang[tmp]=0.
           tmp=WHERE(Inclang2 GT 90.)
           IF tmp[0] NE -1 then INCLang2[tmp]=90.
           tmp=WHERE(Inclang2 LT 0.)
           IF tmp[0] NE -1 then INCLang2[tmp]=0.
           stringINCL='INCL= '+STRJOIN(string(INCLang),' ')
           tmppos=where('INCL' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
           stringINCL='INCL_2= '+STRJOIN(string(INCLang2),' ')
           tmppos=where('INCL_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringINCL
                                ;VROT
           IF n_elements(VROTarr) GE 4 then begin
              vrotext=string(strtrim(strcompress(TOTAL([VROTarr[n_elements(VROTarr)-4:n_elements(VROTarr)-1],VROTarr2[n_elements(VROTarr2)-4:n_elements(VROTarr2)-1]])/8.),2))
           ENDIF ELSE BEGIN
              vrotext=string(strtrim(strcompress(TOTAL([VROTarr,VROTarr2])/(n_elements(VROTarr)+n_elements(VROTarr))),2))
              vrotext2=string(strtrim(strcompress(TOTAL(VROTarr2)/n_elements(VROTarr2)),2))
           ENDELSE
           stringVROT='VROT= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')+' '+vrotext
           tmppos=where('VROT' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
           stringVROT='VROT_2= 0. '+STRJOIN(string(VROTarr[1:n_elements(VROTarr)-1]),' ')+' '+vrotext
           tmppos=where('VROT_2' EQ tirificsecondvars)
           tirificsecond[tmppos]=stringVROT
                                ;Update fitting setting and Log file
           lastaddrings=norings[0]
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"We added a ring to the model."
              Close,66
           ENDIF ELSE BEGIN
              print,linenumber()+"We added a ring to the model."
           ENDELSE
           prevmodification=1.
           INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
           PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
           IF norings[0] LE 4 OR finishafter EQ 1.1 then begin
              INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           ENDIF ELSE BEGIN
                                ;update INCL and PA
              set_warp_slopev3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix,sloped=slopedrings
              if fix_pa eq 0 then PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
              if fix_incl EQ 0. then INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           ENDELSE
                                ;Updating surface brightness settings
           set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log,doubled=doubled
                                ;Central parameter setting and Z0
           z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                       ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                         ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
                               ;Rotation settings
           set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
           set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,finish_after=finishafter,slope=slope,fix=fix_sdis

                                ;Write the parameters to the tirific
                                ;array
           case 1 of
              norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
                 Writefittingvariables,tirificsecond,inclinput1,painput1,vrotinput1,sdisinput1, $
                                       sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,xposinput1,yposinput1,vsysinput1,INIMODE=inimode
              end
              fix_incl EQ 0. and fix_pa EQ 1.: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE
              end
              fix_incl EQ 1. and fix_pa EQ 0.: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE

              end
              else: begin
                 IF norings[0] LE 6 then begin
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDIF ELSE BEGIN
                    Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                          vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                 ENDELSE

              end
           endcase

                                ;Update counters and go back to the fitting proces
           sbrmodify++
           finalsmooth=0.
           prevrings=[prevrings,norings[0]]
           goto,notacceptedtwo
        ENDIF
     ENDIF
     noadd:
                                ;If we have tried the second fit three
                                ;times without changing the rings or
                                ;boundaries than accept that the
                                ;second fit failed.
     IF AC2 EQ 1 and trytwo EQ 3. then begin
        AC2=0
        AC1=1
     endif
                                ;IF the second fit is not accepted and
                                ;we have not tried three times without
                                ;changing anything already then update
                                ;parameters and try again
     IF AC2 NE 1  And trytwo LT 3.  then begin
        INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' INCL_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
        PAinput1[0]='PA 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(innerfix,format='(F7.4)')),1)
                                ;Need to update the arrays otherwise might use
                                ;them later with wrong ring numbers
        IF norings[0] LE 4 OR finishafter EQ 1.1  then begin
           INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
           PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0]-1,format='(F7.4)')),1)
        ENDIF ELSE BEGIN
                                ;update INCL and PA
           set_warp_slopev3,SBRarr,SBRarr2,cutoff,INCLinput2,PAinput2,INCLinput3,PAinput3,norings,log=log,innerfix=innerfix,sloped=slopedrings
           if fix_pa eq 0 then PAinput1[0]='PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' PA_2 1:'+strtrim(strcompress(string(norings[0]-1,format='(F7.4)')),1)
           if fix_incl eq 0. then INCLinput1[0]='INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' INCL_2 1:' +strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        ENDELSE
                                ;set the SBR fitting parameters
        set_sbr,SBRinput1,SBRinput2,SBRinput3,SBRinput4,SBRinput5,SBRinput6,SBRarr,cutoff,norings,finishafter,log=log,doubled=doubled
                                ;and some other parameters
        z0input1[0]='Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                    ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        xposinput1[0]=' XPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' XPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        yposinput1[0]=' YPOS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' YPOS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        vsysinput1[0]=' VSYS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                      ' VSYS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
       ; sdisinput1[0]=' SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
        ;              ' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)
        set_vrotv6,vrotinput1,VROTarr,velconstused,vrotmax,vrotmin,norings,channelwidth,avinner=avinner,centralexclude=centralexclude,finish_after=finishafter
        set_sdis,sdisinput1,SDISarr,velconstused,sdismax,sdismin,norings,channelwidth,finish_after=finishafter,slope=slope,fix=fix_sdis

                                ;But we need to reset all the fitting parameters
        IF trytwo LT 2. then begin
           IF AC1 NE 1 then begin
              ;If the first estimate is not accepted either

              case 1 of
                 norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
                    Writefittingvariables,tirificsecond,inclinput1,painput1,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                          xposinput1,yposinput1,vsysinput1,INIMODE=inimode

                 end
                 fix_incl EQ 0. and fix_pa EQ 1.: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDELSE
                 end
                 fix_incl EQ 1. and fix_pa EQ 0.: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDELSE

                 end
                 else: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,$
                                             xposinput1,yposinput1,vsysinput1,INIMODE=inimode
                    ENDELSE

                 end
              endcase
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Both the first and the second estimate weren't accepted, retrying."
                 close,66
              ENDIF
           endif else begin
                                ;IF only the second model is not
                                ;accepted try without fitting the center

              case 1 of
                 norings[0] LE 4 OR finishafter EQ 1.1 OR (fix_incl EQ 0. and fix_pa EQ 0.): begin
                    Writefittingvariables,tirificsecond,inclinput1,painput1,$
                                          vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                 end
                 fix_incl EQ 0. and fix_pa EQ 1.: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDELSE
                 end
                 fix_incl EQ 1. and fix_pa EQ 0.: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDELSE

                 end
                 else: begin
                    IF norings[0] LE 6 then begin
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput1,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDIF ELSE BEGIN
                       Writefittingvariables,tirificsecond,inclinput1,inclinput2,inclinput3,painput1,painput2,painput3,$
                                             vrotinput1,sdisinput1,sbrinput5,sbrinput1,sbrinput6,sbrinput2,sbrinput4,sbrinput3,z0input1,INIMODE=inimode
                    ENDELSE

                 end
              endcase
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"The second estimate wasn't accepted. Running with fixed center."
                 close,66
              ENDIF
           endelse
        endif else begin
                                ;if we are on the third try
           IF AC1 NE 1 then begin
                                ;if the center is not yet accepted try
                                ;just changing the center
              xposinput1[3]='1E3'
              yposinput1[3]='1E3'
              vsysinput1[3]='4'
              Writefittingvariables,tirificsecond,xposinput1,yposinput1,vsysinput1,INIMODE=inimode
              IF size(log,/TYPE) EQ 7 then begin
                 openu,66,log,/APPEND
                 printf,66,linenumber()+"Last try on the second fit."
                 close,66
              ENDIF
           endif else begin
                                ;Else wrap up a failed fit
              goto,noretrytwo
           endelse
        endelse
                                ; update counter and try again
        trytwo++
        goto,notacceptedtwo
     endif
     noretrytwo:
                                ;if we get here without polynomial
                                ;smoothing then set smoothing and run parameters
     IF finalsmooth LT 1 then begin
        finalsmooth=1
        IF trytwo EQ 3 and AC2 NE 1 then noacceptbeforesmooth=1
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"We have a succesful model and we will smooth with polynomials now."
           close,66
        END
        goto,lastadjust
     ENDIF
                                ;IF the fitting rings are very small set them to the minimum of 4
     IF newindrings[0] LT 4 then newindrings[0]=4
     IF newindrings[1] LT 4 then newindrings[1]=4
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We'll fit from "+strtrim(strcompress(STRJOIN(string(newindrings,format='(F7.1)'),' ')),1)
        close,66
     ENDIF
                                ;set the parameters for single fitting
     SBRinput1=['SBR '+strtrim(strcompress(string(newindrings[0],format='(F7.4)')),1)+':1','3',strtrim(strcompress(string(cutoff[n_elements(cutoff)-1]/2.,format='(E12.5)')),1),'1E-6','1E-7','1E-5','1E-7','3','70','70']
     SBRinput2=SBRinput1
     SBRinput2[0]='SBR_2 '+strtrim(strcompress(string(newindrings[1],format='(F7.4)')),1)+':3'
     IF norings GT 4 and finishafter NE 1.1 then begin
                                ;PA
        PAinput1=['PA 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' '+$
                  'PA_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                  '360','0',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70']
                                ;INCL
        INCLinput1=['INCL 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' '+$
                    'INCL_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                    '90','5',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70']

        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                    ,'500',string(channelwidth),string(channelwidth*0.5),string(0.1*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70',' ']
                                ;SDIS
        If INCLang[0] GT 60 then begin
           SDISinput1=['SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                       string(MAX(SDISarr)),'3.',string(channelwidth*0.5),string(0.1*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70']
        ENDIF ELSE BEGIN
          SDISinput1=['SDIS 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' SDIS_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                      '40','3.',string(channelwidth*0.5),string(0.1*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70']
       ENDELSE
                                ;Z0
        Z0input1=['Z0 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+$
                  ' Z0_2 1:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1),$
                  string(convertskyanglefunction(0.4,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.075,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.005,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.5,double(catDistance[i]),/PHYSICAL)),$
                  string(convertskyanglefunction(0.0005,double(catDistance[i]),/PHYSICAL)),'3','70','70']

        IF catinc[i] GT 80 then begin
           Z0input1[1]=string(MAX([convertskyanglefunction(1.0,double(catDistance[i])),catmajbeam[i]/2.]))
        ENDIF
        IF catDistance[i] EQ 1. then begin
                                ;ok if we do not know the distance
                                ;then let us assume each disk is about
                                ;30 kpc which means that
                                ;0.5=norings[0]/60. and so on
           IF catinc[i] GT 80 then begin
              Z0input1[1]=string(MAX([catmajbeam[i]*ring_spacing*norings[0]/60.,catmajbeam[i]/2.]))
           ENDIF ELSE  Z0input1[1]=string(catmajbeam[i]*ring_spacing*norings[0]/60.)
           Z0input1[2]='0.'
           Z0input1[3]=string(-1*catmajbeam[i]*ring_spacing*norings[0]/6E4)
           Z0input1[4]=string(catmajbeam[i]*ring_spacing*norings[0]/6E5)
           Z0input1[5]=string(catmajbeam[i]*ring_spacing*norings[0]/60.)
           Z0input1[6]=string(catmajbeam[i]*ring_spacing*norings[0]/6E5)
        ENDIF
     ENDIF ELSE BEGIN
        VROTinput1=['VROT 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)+' VROT_2 2:'+strtrim(strcompress(string(norings[0],format='(F7.4)')),1)$
                    ,'500',string(channelwidth),string(channelwidth*0.1),string(0.05*channelwidth),string(0.5*channelwidth),string(0.1*channelwidth),'3','70','70',' ']
     ENDELSE
                                ;write the parameters to the tirific array
     Writefittingvariables,tirificsecond,PAinput1,INCLinput1,VROTinput1,SDISinput1,SBRinput1,SBRinput2,Z0input1,INIMODE=inimode
     IF testing GE 2 OR finalsmooth EQ 2. then goto,testing2check
                                ;Write to file
     restart_counter++
     tmppos=where('RESTARTID' EQ tirificsecondvars)
     tirificsecond[tmppos]='RESTARTID= '+string(fix(restart_counter))
     openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
     for index=0,n_elements(tirificsecond)-1 do begin
        printf,1,tirificsecond[index]
     endfor
     close,1
                                ;Perform tirific check with only changes as a whole to the parameters
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Starting tirific check of the smoothed second estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        printf,66,linenumber()+"We have changed the ring numbers "+string(sbrmodify)+" times."
        printf,66,linenumber()+"We have changed the fitting parameters "+string(trytwo)+" times."
        close,66
     ENDIF
     IF MEAN(VROTarr) GT 1000. then begin
        print,'something has gone horribly wrong'
        stop
     ENDIF
     print,linenumber()+"Starting tirific check of second estimate in "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
     run_tirific, continue_tirific, curr_run_id,bookkeeping,output_name='2ndfit.', $
                  store_name='2ndfitunsmooth.',log=log,loops=loops,nopoints=nopoints,AC=AC2, $
                  toymodels=toymodels,run_unit=run_unit


     spawn,'cp tirific.def the_last_input.def'
     finalsmooth=2.
     overwrite=0.

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        IF AC2 EQ 1 then begin
           printf,66,linenumber()+"The smoothed second estimate was accepted in this run."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDIF ELSE BEGIN
           printf,66,linenumber()+"The smoothed second estimate was not accepted."
           printf,66,linenumber()+"Tirific ran through "+strtrim(strcompress(string(loops)))+" loops and produced "+strtrim(strcompress(string(toymodels)))+" models."
           printf,66,linenumber()+"Disk 1 had "+strtrim(strcompress(string(nopoints[0])))+" point sources and Disk 2 "+strtrim(strcompress(string(nopoints[1])))
        ENDELSE
        close,66
     ENDIF
                                ;Go and see wether all parameters are satisfied
     goto,lastadjust

     testing2check:





                           ;Clean up the model
     IF finalsmooth LT 3 then begin
        firstfitvaluesnames=0.
        writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=firstfitvalues,VariableChange=firstfitvaluesnames,Variables=tirificsecondvars
        tmp=WHERE(firstfitvaluesnames EQ 'SBR')
        SBRarr=firstfitvalues[*,tmp]
        tmp=WHERE(firstfitvaluesnames EQ 'SBR_2')
        SBRarr2=firstfitvalues[*,tmp]
        get_newringsv9,SBRarr,SBRarr2,cutoff,finalrings,/individual
        IF finalrings[0]+1 LT n_elements(SBRarr) then SBRarr[finalrings[0]+1:n_elements(SBRarr)-1]=1e-16
        IF finalrings[1]+1 LT n_elements(SBRarr2) then SBRarr2[finalrings[1]+1:n_elements(SBRarr2)-1]=1e-16

        SBRarr[0:1]=(SBRarr[0:1]+SBRarr2[0:1])/2.
        SBRarr2[0:1]=SBRarr[0:1]

        stringSBR='SBR= '+STRJOIN(string(SBRarr),' ')
        tmppos=where('SBR' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR
        stringSBR='SBR_2= '+STRJOIN(string(SBRarr2),' ')
        tmppos=where('SBR_2' EQ tirificsecondvars)
        tirificsecond[tmppos]=stringSBR

                                ;In this case we want to add them back
                                ;in
        tmpSBR=(SBRarr+SBRarr2)/2.
        IF finishafter EQ 1.1 then begin
           get_newringsv9,tmpSBR,tmpSBR,cutoff,newend
           IF newend LT 4 then newend=4
           rename,'2ndfit.','2ndfituncor.'
           tmp=WHERE(firstfitvaluesnames EQ 'NUR')
           rings=firstfitvalues[0,tmp]
           IF newend LT rings then begin
              changeradii,tirificsecond,newend,continue_tirific
              INCLinput1=['INCL 1','90','89',string(0.5),string(0.1),string(0.5),string(0.1),'3','70','70',' ']
              writefittingvariables,tirificsecond,INCLinput1,INIMODE=inimode
              newrings=newend
           ENDIF ELSE newrings=rings
        ENDIF
        IF optimized then begin
           tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
           IF n_elements(tmp) EQ 2 then begin
              currentfitcube=tmp[0]
              tmppos=where('INSET' EQ tirificsecondvars)
              tirificsecond[tmppos]='INSET=  '+strtrim(strcompress(string(currentfitcube+'.fits')))
           ENDIF

        ENDIF
        tmp1=WHERE(SBRarr LT 1.1E-16)
        tmp2=WHERE(SBRarr2 LT 1.1E-16)
        IF n_elements(tmp1) GT 1 AND n_elements(tmp2) GT 1 then begin
           IF tmp1[0] GT tmp2[0] then newrings=tmp1[0]+1 else newrings=tmp2[0]+1
           IF newrings LT 4 then newrings=4
           changeradii,tirificsecond,newrings,continue_tirific
        ENDIF ELSE BEGIN
           IF finishafter NE 1.1 then newrings=n_elements(SBRarr)
        ENDELSE
        tmppos=where('LOOPS' EQ tirificsecondvars)
        tirificsecond[tmppos]='LOOPS=  0'
        tmppos=where('RESTARTNAME' EQ tirificsecondvars)
        tirificsecond[tmppos]='RESTARTNAME= '
        print,'!!!!!!!!! Rhis is being weird'
        print,  tirificsecond[tmppos]
        print,tirificsecondvars[tmppos]
        continue_tirific = 'optimized'
        tmppos=where('INIMODE' EQ tirificsecondvars)
        tirificsecond[tmppos]='INIMODE=  0'
        errorsadd=['VROT','VROT_2','PA','PA_2','INCL','INCL_2','SDIS','SDIS_2']
        tmpfile=strarr(n_elements(tirificsecond)+8)
        added=0
                                ;Add the errors to the final tirific
                                ;file

        for j=0,n_elements(tirificsecond)-1 do begin
           tmp=str_sep(strtrim(strcompress(tirificsecond[j]),2),'=')
           tmpex=WHERE(tmp[0] EQ errorsadd)
           if tmpex[0] NE -1 then begin
              tmpfile[j+added]=tirificsecond[j]
              added++
              case tmpex[0] of
                 0:begin
                    IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings)
                    errs=dblarr(n_elements(sigmarot))
                    errs=sigmarot
                    ;/SIN(incarr1*!DtoR)
                    errs[0]=double(channelwidth)
                    tmpposv=WHERE(firstfitvaluesnames EQ 'VROT')
                    avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])
                    tmpavrot=WHERE(errs GT avrot/2.)

                    IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
                 end
                 1: begin
                    IF n_elements(sigmarot) LT newrings then sigmarot=replicate(channelwidth,newrings)
                    errs=sigmarot
                    errs[0]=double(channelwidth)
                    tmpposv=WHERE(firstfitvaluesnames EQ 'VROT_2')
                    avrot=TOTAL(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])/n_elements(firstfitvalues[1:n_elements(firstfitvalues[*,tmpposv])-1,tmpposv])
                    tmpavrot=WHERE(errs GT avrot/2.)

                    IF tmpavrot[0] NE -1 AND avrot NE 0 then errs[WHERE(errs GT avrot/2.)]=avrot/2.
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(errs[0:newrings-1]))),' ')
                 end
                 2: begin
                    IF n_elements(sigmapa1) LT newrings then sigmapa1=replicate(2,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa1[0:newrings-1]))),' ')
                 end
                 3: begin
                     IF n_elements(sigmapa2) LT newrings then  sigmapa2=replicate(2,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmapa2[0:newrings-1]))),' ')
                 end
                 4: begin
                     IF n_elements(sigmaincl1) LT newrings then  sigmaincl1=replicate(4,newrings)
                     tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl1[0:newrings-1]))),' ')
                  end
                 5: begin
                     IF n_elements(sigmaincl2) LT newrings then   sigmaincl2=replicate(4,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmaincl2[0:newrings-1]))),' ')
                 end
                 6: begin
                    IF n_elements(sigmasdis) LT newrings then   sigmasdis=replicate(channelwidth/2.,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmasdis[0:newrings-1]))),' ')
                 end
                 7: begin
                     IF n_elements(sigmasdis) LT newrings then   sigmasdis=replicate(channelwidth/2.,newrings)
                    tmpfile[j+added]='# '+tmp[0]+'_ERR='+STRJOIN(strtrim(strcompress(string(sigmasdis[0:newrings-1]))),' ')
                 end
                 else:print,'odd'
              endcase
           endif else tmpfile[j+added]=tirificsecond[j]
        endfor
        tirificsecond=tmpfile
                                ;write the final file and produce the final model
        restart_counter++
        tmppos=where('RESTARTID' EQ tirificsecondvars)
        tirificsecond[tmppos]='RESTARTID= '+string(fix(restart_counter))
        openw,1,maindir+'/'+catdirname[i]+'/tirific.def'
        for index=0,n_elements(tirificsecond)-1 do begin
           printf,1,tirificsecond[index]
        endfor
        close,1
        if optimized then begin
           output_name= '2ndfit.'
           store_name = '2ndfit_opt.'
        endif else BEGIN
           output_name= '2ndfit.'
           store_name = '2ndfitold.'
        ENDELSE
        run_tirific, continue_tirific, curr_run_id,bookkeeping,output_name=output_name, $
                     store_name=store_name,log=log,run_unit=run_unit



        spawn,'cp tirific.def 2ndfit.def',isthere
        finalsmooth=3
        secondtime=1
        IF size(log,/TYPE) EQ 7 then begin
          openu,66,log,/APPEND
          printf,66,linenumber()+"This should be the very last thing and not happening again."
          close,66
        ENDIF
     ENDIF
                                ;Write the final info and pv-diagrams
     Basicinfovars=['XPOS','YPOS','VSYS','PA','INCL','VROT','SDIS']
     writenewtotemplate,tirificsecond,maindir+'/'+catdirname[i]+'/2ndfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,Variables=tirificsecondvars,/EXTRACT
     RAhr=Basicinfovalues[0,0]
     RAdiff=maxchangeRA*3600./15.
     DEChr=Basicinfovalues[0,1]
     DECdiff=maxchangeDEC*3600.
     tmpcube=readfits(maindir+'/'+catdirname[i]+'/2ndfit.fits',hedtmp1stcube,/SILENT)
     IF optimized then begin
        extract_pv,nooptcube,noptheader,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+noptname[0]+'_2_xv.fits',float(xv),new_header
     ENDIF ELSE BEGIN
        extract_pv,dummy,header,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
        writefits,maindir+'/'+catdirname[i]+'/'+currentfitcube+'_2_xv.fits',float(xv),new_header
     ENDELSE
     extract_pv,tmpcube,hedtmp1stcube,Basicinfovalues[0,3],xv,center=[RAhr,DEChr],xvheader=new_header
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_xv.fits',float(xv),new_header
     hedtmp1st=hedtmp1stcube
     tmpix=WHERE(tmpcube GT catnoise[i])
     totflux=[TOTAL(tmpcube[tmpix])/pixperbeam,(2.*cutoff*norings[0])/(n_elements(tmpix)/pixperbeam)]
     mask=readfits(maindir+'/'+catdirname[i]+'/'+catmaskname[i]+'.fits',headermask,/NOSCALE,/SILENT)
                                ;mask the data cube
     tmpmask=fltarr(n_elements(tmpcube[*,0,0]),n_elements(tmpcube[0,*,0]),n_elements(tmpcube[0,0,*]))
     tmpmask[WHERE(mask GT 0.)]=tmpcube[WHERE(mask GT 0.)]
     momentsv2,tmpmask,tmpmap,hedtmp1st,0.
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_mom0.fits',float(tmpmap),hedtmp1st
     hedtmp1stv=hedtmp1stcube
     momentsv2,tmpmask,tmpmapv,hedtmp1stv,1.,gdlidl=gdlidl
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_mom1.fits',float(tmpmapv),hedtmp1stv
     hedtmp1stw=hedtmp1stcube
     momentsv2,tmpmask,tmpmapw,hedtmp1stw,2.,gdlidl=gdlidl
     writefits,maindir+'/'+catdirname[i]+'/2ndfit_mom2.fits',float(tmpmapw),hedtmp1stw
     getDHI,tmpmap,hedtmp1st,Basicinfovalues[0,3],[RAhr,DEChr,Basicinfovalues[0,4]],DHI
     VSYSdiff=maxchangevel
     HIMASS=2.36E5*catDistance[i]^2*totflux*ABS(channelwidth)
     convertradec,RAhr,DEChr
     openu,1,basicinfofile,/APPEND
     printf,1,"#After the second fit"
     printf,1,format='(A25,A25,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20,A20)',string(RAhr+'+/-'+strtrim(strcompress(string(RAdiff,format='(F6.1)')),2)),$
            string(DEChr+'+/-'+strtrim(strcompress(string(DECdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,2])),2)+'+/-'+strtrim(strcompress(string(VSYSdiff,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,3])),2)+'+/-'+strtrim(strcompress(string(catPAdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[0,4])),2)+'+/-'+strtrim(strcompress(string(catincdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(Basicinfovalues[n_elements(Basicinfovalues[*,5])-1,5])),2)+'+/-'+strtrim(strcompress(string(catmaxrotdev[i],format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(W50)),2)+'+/-'+strtrim(strcompress(string(channelwidth,format='(F6.1)')),2)),$
            string(strtrim(strcompress(string(totflux[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(DHI,format='(F8.1)')),2)),$
            string(strtrim(strcompress(string(catDistance[i])),2)),$
            string(strtrim(strcompress(string(HIMass[0],format='(E10.3)')),2)),$
            string(strtrim(strcompress(string(convertskyanglefunction(DHI,catDistance[i]),format='(F8.1)')),2))
     close,1

     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"We finished the second fit of "+catDirname[i]+"  at "+systime()
        case prevmodification of
           0:lastringstr=" no change"
           -1:lastringstr=" a subtraction"
           1:lastringstr=" adding a ring"
           else: lastringstr=" something weird happened"
        endcase
        printf,66,linenumber()+"We have rerun "+string(counter)+" times and the last modification was "+lastringstr
        close,66
     ENDIF
     IF noacceptbeforesmooth EQ 1 then begin
        IF AC2 EQ 1 then AC2=2
     ENDIF
     IF finishafter EQ 2 OR finishafter EQ 2.1 then begin
                                ;Below 8 beams in diameter FAT is not
                                ;reliable so if in that range or when
                                ;the inclination is below 20 the fit
                                ;is failed
        case 1 of
           norings[0]*ring_spacing LE 5:begin
              comment = 'The finalmodel is less than 8 beams in diameter. FAT is not necessarily reliable in this range.'
              commentlen='A'+strtrim(string(strlen(comment)),2)
              openu,1,outputcatalogue,/APPEND
              printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(2),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
              close,1
           end
           Basicinfovalues[0,4] LT 20:begin
              comment = 'The final inclination is below 20 degrees. FAT is not necessarily reliable in this range.'
              commentlen='A'+strtrim(string(strlen(comment)),2)
              openu,1,outputcatalogue,/APPEND
              printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(2),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
              close,1
           end
           else:begin
              comment = ' You have chosen to skip the fitting process after the second fit'
              commentlen='A'+strtrim(string(strlen(comment)),2)
              openu,1,outputcatalogue,/APPEND
              printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(fix(AC2)),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
              close,1
           end
        endcase
        goto,finishthisgalaxy
     ENDIF ELSE begin
        IF finishafter EQ 1.1 then begin
                                ;Below 8 beams in diameter FAT is not
                                ;reliable so if in that range or when
                                ;the inclination is below 20 the fit
                                ;is failed
           case 1 of
              norings[0]*ring_spacing LE 5:begin
                 comment = 'The finalmodel is less than 8 beams in diameter. FAT is not necessarily reliable in this range.'
                 commentlen='A'+strtrim(string(strlen(comment)),2)
                 openu,1,outputcatalogue,/APPEND
                 printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(2),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
                 close,1
              end
              Basicinfovalues[0,4] LT 20:begin
                 comment = 'The final inclination is below 20 degrees. FAT is not necessarily reliable in this range.'
                 commentlen='A'+strtrim(string(strlen(comment)),2)
                 openu,1,outputcatalogue,/APPEND
                 printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(2),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
                 close,1
              end
              else:begin
                 comment = 'As the galaxy radius is already halved we stop here'
                 commentlen='A'+strtrim(string(strlen(comment)),2)
                 openu,1,outputcatalogue,/APPEND
                 printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(fix(AC2)),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
                 close,1
              end
           endcase
           goto,finishthisgalaxy
        ENDIF ELSE BEGIN
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              IF AC2 EQ 1 then printf,66,linenumber()+"The Second estimate was accepted." $
              else printf,66,linenumber()+"The Second estimate was not accepted."
              close,66
           ENDIF ELSE BEGIN
              IF AC2 EQ 1 then print,linenumber()+"The Second estimate was accepted." $
              else print,linenumber()+"The Second estimate was not accepted."
           ENDELSE
        ENDELSE
     ENDELSE
     IF AC1 EQ 0 AND AC2 EQ 0 then begin
        comment = 'Both AC1 and AC2 could not get accepted. Aborting the fit.'
        commentlen='A'+strtrim(string(strlen(comment)),2)
        openu,1,outputcatalogue,/APPEND
        printf,1,catDirname[i],strtrim(string(fix(AC1)),2),strtrim(string(fix(AC2)),2),comment,format='("",('+dirformat+')," ",2(A6),"  ",('+commentlen+'))'
        close,1
        goto,finishthisgalaxy
     ENDIF
     IF warpoutput then begin
        Basicinfovars=['RADI','SBR','SBR_2','PA','PA_2','PA_ERR','PA_2_ERR','INCL','INCL_2','INCL_ERR','INCL_2_ERR','VROT','VROT_ERR','RMS','BMAJ','BMIN']
        tirificfirst=1
        writenewtotemplate,tirificfirst,new_dir+'2ndfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,/EXTRACT
        get_fixedringsv9,[[Basicinfovalues[*,3]],[Basicinfovalues[*,4]],[Basicinfovalues[*,7]],[Basicinfovalues[*,8]]],fixedrings,/warp_output,workingdir=new_dir,radii=[Basicinfovalues[*,0]],SBR=[[Basicinfovalues[*,1]],[Basicinfovalues[*,2]]]
     ENDIF

;!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!The current Code ends here!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                                ;This is where the third estimate
                                ;starts in previous versions. Here we change to 0.5 beam
                                ;size radii and we try to determine
                                ;the start of the warp fitting those
                                ;parameters as 1 for greater accuracy
                                ;first the standard changes


     finishthisgalaxy:
     ; Make sure all tirifics are dead
     kill_tirific,curr_run_id,run_unit
     IF optimized then begin
        catcubename[i]= noptname[0]
        tmp=str_sep(strtrim(strcompress(currentfitcube),2),'_opt')
        IF n_elements(tmp) EQ 2 then begin
           currentfitcube=tmp[0]
        ENDIF
     ENDIF


     cd, new_dir
                                ;If we want information about the warp
                                ;and we have a warp then make a
                                ;directory for the warp output and
                                ;write the tiltograms and info. Else
                                ;log why not
     IF warpoutput AND finishafter GE 2 AND bookkeeping LT 5 then begin
        IF size(log,/TYPE) EQ 7 then begin
           openu,66,log,/APPEND
           printf,66,linenumber()+"Creating the Warp_Info Directory for advanced information about the fitted warp"
           close,66
        ENDIF
        Basicinfovars=['RADI','SBR','SBR_2','PA','PA_2','PA_ERR','PA_2_ERR','INCL','INCL_2','INCL_ERR','INCL_2_ERR','VROT','VROT_ERR','RMS','BMAJ','BMIN']
        tirificfirst=1
        writenewtotemplate,tirificfirst,'2ndfit.def',Arrays=Basicinfovalues,VariableChange=Basicinfovars,/EXTRACT
        get_fixedringsv9,[[Basicinfovalues[*,3]],[Basicinfovalues[*,4]],[Basicinfovalues[*,7]],[Basicinfovalues[*,8]]],fixedrings,/warp_output,radii=[Basicinfovalues[*,0]],SBR=[[Basicinfovalues[*,1]],[Basicinfovalues[*,2]]]
     ENDIF ELSE BEGIN



        IF  warpoutput and finishafter EQ 1.1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"As the galaxy is small and we have not allowed for a warp, the directory Warp_Info will not be created."
              close,66
           ENDIF
        ENDIF
        IF warpoutput and finishafter EQ 1 then begin
           IF size(log,/TYPE) EQ 7 then begin
              openu,66,log,/APPEND
              printf,66,linenumber()+"As you have chosen not to fit a warp, the directories Warp_Info and No_Warp will not be created."
              close,66
           ENDIF
        ENDIF
     ENDELSE
     ; We should add some standard errors to the 1st fit as well






                                ;Clearing up the direcory and
                                ;organizing the output in proper names
                                ;and such
     if finishafter LT 1.75 then begin
        fix_pa = 0
        fix_incl = 0.
        fix_sdis = 0
     ENDIF
     names=[currentfitcube,catMom0name[i],catMom1name[i],catmaskname[i],noisemapname,catCatalogname[i],basicinfo,catMom2name[i]]
     book_keeping,names,bookkeeping,catdistance[i],gdlidl,log=log,noise=catnoise[i],finishafter=finishafter,fixedpars=[fix_pa,fix_incl,fix_sdis]
     IF size(log,/TYPE) EQ 7 then begin
        openu,66,log,/APPEND
        printf,66,linenumber()+"Finished "+catDirname[i]+" which is galaxy # "+strtrim(string(fix(i)),2)+" at "+systime()
        close,66
     ENDIF
     IF bookkeeping EQ 5 then print,linenumber()+'FAT did not run the full fitting routines. Please check your output log or screen messages carefully.'
     CD,old_dir
     if KEYWORD_SET(timing) then BEGIN
        openu,22,maindir+'/Timing_Result.txt',/APPEND
        printf,22,"The galaxy in directory  "+catDirname[i]+" started at"+starttime
        printf,22,"Finished the first step at "+endtime_firststep
        printf,22,"And finished the whole process at "+systime()
        close,22
     ENDIF
  endfor
  CD,originaldir
  close,3
  if keyword_set(installation_check) then begin
     check=install_check(gdlidl)
     check_error:
     case 1 of
        check EQ 0:begin
           print, ' '
           print,'!!!!--------------------------------------------!!!!!'
           print,'!!!! As far as we can tell FAT is installed     !!!!!'
           print,'!!!! properly and runs smoothly.                !!!!!'
           print,'!!!!--------------------------------------------!!!!!'
        end
        check EQ 1:begin
           print, ' '
           print,'!!!!---------------------------------------------!!!!!'
           print,'!!!! FAT has crashed somewhere during the        !!!!!'
           print,'!!!! fitting. Please check the error message     !!!!!'
           print,'!!!! and all dependencies. If you are unable to  !!!!!'
           print,'!!!! resolve the issue please file a bug report  !!!!!'
           print,'!!!! at:                                         !!!!!'
           print,'!!!!                                             !!!!!'
           print,'!!!! https://github.com/PeterKamphuis/FAT/issues !!!!!'
           print,'!!!!---------------------------------------------!!!!!'
        end
        check eq 2: begin
           print, ' '
           print,'!!!!---------------------------------------------!!!!!'
           print,'!!!! FAT ran through the fitting process but the !!!!!'
           print,'!!!! fitted values differ too much from their    !!!!!'
           print,'!!!! expectations. Please update SoFiA and other !!!!!'
           print,'!!!! dependencies. If you are unable to resolve  !!!!!'
           print,'!!!! resolve the issue please file a bug report  !!!!!'
           print,'!!!! at:                                         !!!!!'
           print,'!!!!                                             !!!!!'
           print,'!!!! https://github.com/PeterKamphuis/FAT/issues !!!!!'
           print,'!!!!                                             !!!!!'
           print,'!!!! Please add the Log.txt file in the directory!!!!!'
           print,'!!!! Installation_Check and the Finalmodel.def   !!!!!'
           print,'!!!! in the Installation_Check/Finalmodel/       !!!!!'
           print,'!!!! directory to your report.                   !!!!!'
           print,'!!!!---------------------------------------------!!!!!'
        end
        else: begin
           print, ' '
           print,'!!!!---------------------------------------------!!!!!'
           print,'!!!! This should not happen please file a bug    !!!!!'
           print,'!!!! report at:                                  !!!!!'
           print,'!!!!                                             !!!!!'
           print,'!!!! https://github.com/PeterKamphuis/FAT/issues !!!!!'
           print,'!!!!---------------------------------------------!!!!!'
        end
     endcase
                                ;We always will want to clean up the directory
  endif
end

'''

FAT.__doc__ = '''
;+
; NAME:
;      FAT
; PURPOSE:
;      Fit Tilted Ring Models with Tirific in a fully automated manner
; CATEGORY:
;      Main for fitting galaxies. Tirific still requires interactive fitting this code attempts
;      to remedy that
;
; CALLING SEQUENCE:
;      FAT,support='supportdir',configuration_file='configfile'
;
; INPUTS:
;      -
; OPTIONAL INPUTS:
;      SUPPORT  = path to the directory where FAT's support
;      routines are located. The default location is ./Support/
;      CONFIGURATION_FILE = A configuration file for FAT. This file
;      should contain the locations of the galaxies to be fitted. See
;      readme for more detailed info.
;
; OPTIONAL INPUT KEYWORDS
;     /INSTALLATION_CHECK = Flag to run the Installation check.
; ---------------------------------------------------------------------------------
;     The following input keywords are only meant to be used by
;     developers. Except for the /debug flag they will not work for
;     the common user. If you want to know about these please contact Peter
;     Kamphuis.
; ---------------------------------------------------------------------------------
;     /DEBUG = Flag to print debugging information in several routines
;     /LVHIS_TEST = Flag to run the LVHIS Test.
;     /PAPER_TEST = Flag to run the paper artificial galaxies.
;     /RESOLUTION_TEST = Flag to run the additional resolution tests
; OPTIONAL KEYWORD OUTPUT:
;      -
;
; OUTPUTS:
;     See Readme or just run the code
;
; EXAMPLE:
;     python3 FAT.py --cf /home/your_computer/FAT_dir/FAT_INPUT.config'
'''


if __name__ == '__main__':
    error_message = '''
                    Your code has crashed for some reason. If this message completely baffles you then please submit the trace back as a bug report to: \n
                    https://github.com/PeterKamphuis/FAT/issues \n
                    If the error occured while fitting a galaxy, please attach your fitting log as well
                    '''
    try:
        FAT(sys.argv[1:])
    except Exception as exc:
        print('------------When filing a bug report please copy all output  below this line------------')
        traceback.print_exc()
        print(error_message)
        sys.exit(1)
