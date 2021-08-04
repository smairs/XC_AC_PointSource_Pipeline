def TCTrigger(input_data,protocat,diskcat,region, aperture_diam = 0.00083333, trigger_thresh = 5, brightness_thresh = 0.0, sd_thresh = 2,wave='850',mjypbmfactor=537000.0,mjyparcsecfactor=2340.0,fidnoiseterm=14,fidcalterm=0.02,WEIGHTED=False,GOODBOX=False):
    '''
    This program loads all of the functions defined in TCTriggerFunctions.py
    and exectues them in their correct order such that we are able to
    find potential variable sources as new Transient data is collected.

    #########################

    IN:

    input_data        = list or array          # path names to the CROPPED AND SMOOTHED data
    protocat          = string: 'protocat.txt' # the protostar catalogue 
    diskcat           = string: 'diskcat.txt'  # the disk catalogue
    region            = string                 # The region name you are analysing


    aperture_diam     = 0.00083333             # Keyword. 3 arcseconds in degrees - 1 pixel, the area over 
                                                          which to measure the flux for each source
    trigger_thresh    = 5                      # Keyword. The threshold for identifying new variables:
                                                          abs(flux_measure - flux_mean)/SD >= trigger_thresh
                                                          If the above expression is true for any peak flux
							  measurement, the source is marked as a potential variable

    brightness_thresh = keyword. The minimum brightness we will consider sources when measuring their potentially discrepant standard deviation
                                     in peak flux. In Jy/beam

    sd_thresh         = keyword. The source SD is computed and compared to the fiducial model presented in Johnstone et al. 2018. If sd_measured / sd_fiducial > sd_thresh,
                                 the source will be marked as a potential variable. Only sources with brightnesses greater than brightness_thresh will be considered.

    wave              = keyword. Which wavelength? '450' or '850'

    mjypbmfactor      = keyword. The nominal calibration factor for conversion to mJy/beam for this wavelength provided by the JCMT (850 microns = 537000, 450 microns = 491000)

    mjyparcsecfactor  = keyword. The nominal calibration factor for conversion to mJy/arcsec^2 for this wavelength provided by the JCMT (850 microns = 2340, 450 microns = 4710)

    fidnoiseterm      = keyword. The noise term in the fiducial standard deviation model (noiseterm**2+(calterm*flux)**2)**0.5 (Johnstone et al. 2018)

    fidcalterm        = keyword. The calibration term in the fiducial standard deviation model (noiseterm**2+(calterm*flux)**2)**0.5 (Johnstone et al. 2018)

    WEIGHTED          = boolean. Has this input data been processed by the weighted calibration algorithm

    GOODBOX           = boolean. Has this input data been identified as a "Good Map" at 450 microns?

    ###########################

    OUT (placed in a directory called "output_dir"):

    1. coadded image of all input epochs

    2. An astropy table containing the following columns:
       UTdates, julian_dates, scans, elevs, tau225s, noises, noise_units

       UNDER DEVELOPMENT: Adding a column for Relative Flux Calibration Factors

    3. FellWalker Catalogue in the same style as a JSA_Catalogue "Peak" .FIT file


    4. Astropy table with the following columns:

       source index, ra, dec, distance to nearest protostar ,protostar class (Megeath et al. 2012 style classes for Protostars and Candidates), Distance to the nearest Class II source, Disk Class (Megeath et al 2012. style: D)


    5. "REGION_variables_DATE.txt" with today's date, containing information about identified potential variables.
    
    6. Astropy table with the following columns:

       Source ID, ra, dec, distance to nearest protostar, distance to nearest Class II object, Mean Peak Flux, Relative Standard Deviation in Peak Fluxes, Slope of linear fit of Peak fluxes over time, Error in the slope, Intercept of the linear fit of peak fluxes over time, Error in the Intercept of the linear fit of peak fluxes over time                                                                                                                                 

     7. The source IDs of Identified Potential Variables

     8. If EXTRASOURCES = True - the above tables and files will also be created for a set of additional sources found in Johnstone et al. 2018.

     ################################

     EXAMPLE:

     >>>from TCTrigger import *
     >>>TCTrigger(['testdata/SERPM_20160202_00054_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160223_00050_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160317_00051_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160415_00046_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160521_00039_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160722_00023_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160827_00012_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20160929_00012_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170222_00070_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170320_00056_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170403_00053_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170417_00044_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170505_00035_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170519_00030_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170602_00041_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170616_00025_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170705_00027_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170721_00017_850_EA3_cal_crop_smooth_jypbm.sdf','testdata/SERPM_20170810_00024_850_EA3_cal_crop_smooth_jypbm.sdf'],'config/protocat.txt','config/diskcat.txt','SerpensMain',aperture_diam = 0.000833333,trigger_thresh = 4, brightness_thresh = 0.0, sd_thresh = 2,wave='850')

     ################################
     ################################
    '''

    #Measure the time it takes to complete this pipeline
    import time
    start_time = time.time()    
    import datetime
    from astropy.table import Table
    from TCTriggerPS.TCTriggerFunctionsPS import TCCoadd,TCMetadata,TCIdentifySources,TCYSOcompare,TCTrackSources,TCCheck4Variables
    import os
    import numpy as np
    import glob
    
    # Make a directory to store all of the output files, so we can organise them at the end of this pipeline

    output_dir = 'pointsource_results/'

    os.system('mkdir '+output_dir)
    os.system('mkdir '+output_dir+'/'+region)

    #############################
    ####### Make a Co-add #######
    #############################

    if WEIGHTED:
        if GOODBOX:
            prev_coadd_file = '*_Wcal_GoodMaps_coadd.sdf'
        else:
            prev_coadd_file ='*_Wcal_coadd.sdf'
    else:
        prev_coadd_file = '*mJybmsm_coadd.sdf'


    if len(sorted(list(glob.glob(output_dir+'/'+region+'/*'+wave+prev_coadd_file))))>0:
        inputfilelist = []
        inputfilelist.append(sorted(list(glob.glob(output_dir+'/'+region+'/*'+wave+prev_coadd_file)))[-1])
        if not GOODBOX:
            for eachnewfile in sorted(input_data)[1:]:
                inputfilelist.append(eachnewfile)
        else:
            for eachnewfile in sorted(input_data):
                inputfilelist.append(eachnewfile)
        coadd_name = TCCoadd(inputfilelist,wave,GOODBOX=GOODBOX)
    else:
        coadd_name = TCCoadd(input_data,wave,GOODBOX=GOODBOX)

    os.system('mv -f '+coadd_name+' '+output_dir+'/'+region)

    ##############################
    ##############################

#
#

    ##################################
    ###### Build Metadata Table ######
    ##################################

    metadata = TCMetadata(input_data,region,output_dir+'/'+region,wave=wave,mjypbmfactor=mjypbmfactor,mjyparcsecfactor=mjyparcsecfactor,WEIGHTED=WEIGHTED,GOODBOX=GOODBOX)
   
    ###################################
    ###################################

#
#

    ##########################################
    ######### Load Source Catalogue ##########
    ##########################################

    # 20171229 - No longer identifying a peak catalogue each time we do a coadd. This could cause sources
    # to move around or have different peak IDs if one becomes brighter than another!
    # So use the catalogues Doug created for his midproject paper at 850 microns and new cats for 450 microns 
    # and high mass regions

    #peakcat_name = TCIdentifySources(output_dir+coadd_name.split('/')[-1])

    if region not in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']:
        peakcat_name = 'config/'+region+'_'+wave+'_sourcecat_20201201.fits'
    elif wave == '450':
        peakcat_name = 'config/'+region+'_'+wave+'_sourcecat_20200911.fits'
    elif wave == '850':
        peakcat_name = 'config/'+region+'_'+wave+'_sourcecat_20170616.fits'

    if os.path.exists(peakcat_name):

        ##########################################
        ##########################################

#
#

        ###############################################
        # If we haven't already done so in the past,  #
        # associate each source with its nearest YSOs #
        ###############################################


        if not os.path.exists(output_dir+'/'+region+'/'+region+'_YSOcompare_'+wave+'.txt'):

            YSOtable     = TCYSOcompare(peakcat_name,protocat,diskcat,region,wave=wave)

            os.system('mv '+region+'_YSOcompare_'+wave+'.txt '+output_dir+'/'+region)
    
        else:

            YSOtable = Table.read(output_dir+'/'+region+'/'+region+'_YSOcompare_'+wave+'.txt',format='ascii') 

        ##################################################
        ##################################################

#
#

        #########################################
        ### Find the Peak Flux of Each Source ###
        #########################################

        source_dict = TCTrackSources(input_data,peakcat_name,region,output_dir+'/'+region,aperture_diam = aperture_diam,wave=wave,WEIGHTED=WEIGHTED,GOODBOX=GOODBOX)

        # Source_dict is a dictionary with each key representing one source
        # Source ID key contains date and peak flux information

        ##########################################
        ##########################################

#
#

        #########################################
        ### Perform the Tests for Variability ###
        #########################################

        os.system('mkdir '+output_dir+'/light_curves/')
        os.system('mkdir '+output_dir+'/light_curves/'+region)
        
        source_slope_table,triggered_sources = TCCheck4Variables(source_dict,YSOtable,region,trigger_thresh = trigger_thresh,brightness_thresh = brightness_thresh,sd_thresh = sd_thresh,wave=wave,fidnoiseterm=fidnoiseterm,fidcalterm=fidcalterm,WEIGHTED=WEIGHTED,GOODBOX=GOODBOX)


        ##################
        #### Clean Up ####
        ##################

        os.system('mv '+region+'*txt '+output_dir+'/'+region)
        #2021-05-14 -- Changed to simply remove these light curves since we have better ones now further down the main pipeline
        #os.system('mv *_lightcurve.pdf '+output_dir+'/light_curves/'+region+'/')
        os.system('rm -f *_lightcurve.pdf')

        ##############################################################################
        #### Add pound sign to each txt file header if it doesn't have it already ####
        ##############################################################################

        # Sourceinfo file:
        # Get most recent file's date for name
        mostrecentdate = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[1]
        mostrecentscan = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[2]

        if WEIGHTED:
            if GOODBOX:
                sourceinfoname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_Wcal_GoodMaps_sourceinfo.txt'
            else:
                sourceinfoname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_Wcal_sourceinfo.txt'
        else:
            sourceinfoname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_sourceinfo.txt'
        needtooverwrite= 0
        f    = open(sourceinfoname)
        newf = open('newfile.txt','w')
        lines = f.readlines() # read old content
        if lines[0][0]=='#':
            newf.close()
            os.system('rm -f newfile.txt')
        else:
            needtooverwrite = 1
            newf.write('#'+lines[0]) # write new content at the beginning
            for line in lines[1:]: # write old content after new
                newf.write(line)
            newf.close()
        f.close()

        if needtooverwrite==1:
            os.system('mv -f newfile.txt '+sourceinfoname)
        
        # YSOcompare File
        needtooverwrite = 0
        f    = open(output_dir+'/'+region+'/'+region+'_YSOcompare_'+wave+'.txt')
        newf = open('newfile.txt','w')
        lines = f.readlines() # read old content
        if lines[0][0]=='#':
            newf.close()
            os.system('rm -f newfile.txt')
        else:
            needtooverwrite = 1
            newf.write('#'+lines[0]) # write new content at the beginning
            for line in lines[1:]: # write old content after new
                newf.write(line)
            newf.close()
        f.close()
    
        if needtooverwrite==1:
            os.system('mv -f newfile.txt '+output_dir+'/'+region+'/'+region+'_YSOcompare_'+wave+'.txt')

        # Metadata file
        if WEIGHTED:
            if GOODBOX:
                metadataname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_Wcal_GoodMaps_metadata.txt'
            else:
                metadataname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_Wcal_metadata.txt'
        else:
            metadataname = output_dir+'/'+region+'/'+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_metadata.txt'
        needtooverwrite = 0
        f    = open(metadataname)
        newf = open('newfile.txt','w')
        lines = f.readlines() # read old content
        if lines[0][0]=='#':
            newf.close()
            os.system('rm -f newfile.txt')
        else:
            needtooverwrite = 1
            newf.write('#'+lines[0]) # write new content at the beginning
            for line in lines[1:]: # write old content after new
                newf.write(line)
            newf.close()
        f.close()

        if needtooverwrite == 1:
            os.system('mv -f newfile.txt '+metadataname)
        #########################################
        #########################################


    else:
        print(peakcat_name+' DOES NOT EXIST!')

    print("\tTCTriggerPS took ",round((time.time() - start_time)/60.0,3), " minutes to run")
