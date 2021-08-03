#########################
#########################
#########################

def TCCoadd(input_data,wave,GOODBOX=False): 
    '''                                                                                                                                                                      
    This Function takes a list of                                                                                                                                            
    cropped, smoothed sdf files (with their full path)                                                                                                                                         
    and creates a coadd named after 
    the most recent file in the list.

    The reference image used is the
    first chronological file in the list.

    ############

    IN:

    input_data = list or array

    A list containing the full path names of
    the input sdf files to be coadded.

    #############

    OUT:

    coadded image placed in the
    same directory as the input images

    #############
    #############
    '''

    from starlink import kappa
    import subprocess
    import os
    import numpy as np

    # List the input sdf files in a text file to be
    # used by KAPPA's WCSMOSAIC.

    wcsmosaic_in = open('input_data.txt','w')
    for eachfile in input_data:
        wcsmosaic_in.write(eachfile+'\n')
    wcsmosaic_in.close()

    # Define the reference file (first file in input_data, chronologically)
    # and the output file name (most recent file with the suffix "coadd")

    reference_file = sorted(input_data)[0]

    datescans_for_sorting = []
    for i in input_data:
        datescans_for_sorting.append(i.split('_'+wave)[0].split('_')[-2]+'_'+i.split('_'+wave)[0].split('_')[-1])
    #print(i.split('_'+wave)[0].split('_')[-2]+'_'+i.split('_'+wave)[0].split('_')[-1])
    datescans_for_sorting = np.array(datescans_for_sorting)
    input_data = np.array(input_data)

    if GOODBOX:
        output_name     = input_data[np.argsort(datescans_for_sorting)][-1].split('.sdf')[0]+'_GoodMaps_coadd'
    else:
        output_name    =  input_data[np.argsort(datescans_for_sorting)][-1].split('.sdf')[0]+'_coadd'

    #Correct our little "NGC2071 is actually NGC2068" faux pas
    #if output_name.split('/')[-1].split('_')[0]=='NGC2071':
    #    output_name = output_name.split('NGC2071')[0]+'NGC2068'+output_name.split('NGC2071')[1]


    # Run wcsmosaic
    
    kappa.wcsmosaic('^input_data.txt',output_name,ref=reference_file,lbnd='!',ubnd='!')

    #Do not crop the file anymore, Instead, we are giving this function files which are already cropped and smoothed

    ##Make the cropping parameter file.
    #crop_parms = open('crop.ini', 'w')
    #crop_parms.write('[CROP_SCUBA2_IMAGES]\n')
    #crop_parms.write('CROP_METHOD = CIRCLE\n')
    #crop_parms.write('MAP_RADIUS = 1200\n')
    #crop_parms.close()

    ##Perform the cropping.
    #crop_command = '${ORAC_DIR}/etc/picard_start.sh CROP_SCUBA2_IMAGES '
    #crop_command += '-log f -recpars crop.ini '+output_name+'.sdf'+' ${1+"$@"};'
    #subprocess.call(crop_command, shell=True)

    # Clean up
    os.system('rm -f input_data.txt')

    return(output_name+'.sdf')

##############################                                                                                                                                                                
##############################                                                                                                                                                                
##############################                                                                                                                                                                

def TCMetadata(input_data,region,output_dir,wave='850',WEIGHTED=False,GOODBOX=False):
    '''                    
    Read in the essential meta-data for each epoch:
    date, scan, zenith angle, tau225, calibration factor, noise, decimal Julian date
                                                                                    
    This can mostly be accomplished by accessing the fits header.

    ##################

    IN:

    input_data = list or array

    A list containing the full path names of
    the input sdf files.

    region     = Region name, example: 'IC348'

    output_dir = The final location of the metadata table

    wave       = keyword. '850' is default. This is for the name of the output file.

    ###################

    OUT:
     
    An astropy table containing the following columns:
    UTdates, julian_dates, scans, elevs, tau225s, noises, noise_units

    UNDER DEVELOPMENT: Adding a column for Relative Flux Calibration Factors

    ###################
    ###################
    '''                                                                             
    # Import the necessary astropy modules                                             
    from starlink import picard                                                     
    from astropy.io import fits                                                     
    from astropy.time import Time                                                   
    from astropy.table import Table
    from astropy.io import ascii as apascii
    import numpy as np                                                              
    import os                                                                       
    import datetime
    import glob

    # First, check to see if a metadata table already exists. This is one of the longest steps in
    # the TCTrigger pipeline - so if a table already exists, don't reinvent the wheel

    #print('\n\n\n',input_data,'\n\n\n')

    if wave == '450':
        mjypbmfactor = 491000
        mjyparcsecfactor = 4710
    elif wave == '850':
        mjypbmfactor = 537000
        mjyparcsecfactor = 2340

    if WEIGHTED:
        if GOODBOX:
            fileending = wave+"_Wcal_GoodMaps_metadata.txt"
            ind1 = -7
        else:
            fileending = wave+"_Wcal_metadata.txt"
            ind1 =  -6
    else:
        fileending = wave+"_metadata.txt"
        ind1 = -5
    existing_metadata_file = sorted(list(glob.glob(output_dir+"/*"+fileending)))

    if len(existing_metadata_file)>0:
        table_already_exists = True 

        t = Table.read(existing_metadata_file[-1],format='ascii') # The list of metadata files has been sorted, above

        new_input_data = []
        for eachinput in sorted(list(input_data)):
            datescanthisinput = eachinput.split('_')[ind1]+'_'+eachinput.split('_')[ind1+1]
            if region+'_'+datescanthisinput+'_'+wave not in list(t['Name']):
                #print('NEW DATA!',eachinput)
                new_input_data.append(eachinput)

        input_data = new_input_data
        if len(new_input_data)>0:
            anythingnew = 1
        else:
            anythingnew = 0

        ## If the table already exists, check to see if we are adding any new information, if not, skip this whole function!
        #anythingnew = 1
        #t = Table.read(existing_metadata_file[-1],format='ascii') # The list of metadata files has been sorted, above
        #print('\n\n\n',t['Name'][-1],'\n\n\n')
        #print(np.array(input_data)[np.argsort(np.array(input_data))][-1].split('_'))
        #if t['Name'][-1] == region+'_'+np.array(input_data)[np.argsort(np.array(input_data))][-1].split('_')[-6]+'_'+np.array(input_data)[np.argsort(np.array(input_data))][-1].split('_')[-5]+'_'+wave:
        #    anythingnew = 0
        #    input_data = []
	## If the table already exists and we are adding something new, just get information from the files not yet included and append them to the existing table
        #if anythingnew == 1:
        #    input_data = [np.array(input_data)[np.argsort(np.array(input_data))][-1]]
        #existing_metadata_file = np.array(existing_metadata_file)[np.argsort(np.array(existing_metadata_file))]
    else:
        table_already_exists = False

    # Generate empty lists for the metadata:

    IDs          = []
    names        = []
    UTdates      = []
    julian_dates = []
    scans        = []
    elevs        = []
    tau225s      = []
    noises       = []
    noise_units  = []

    dummyID = 0
    # For each file, run SCUBA2_MAPSTATS to get all the info we need, 
    # less the relative FCF                                           

    for eachfile in input_data:

        #print('\n\n\n\n',eachfile,'\n\n\n\n')
        dummyID=dummyID+1
        picard_out = picard.scuba2_mapstats(list([eachfile]))
        #print('\n\n\n\n',picard_out,'\n\n\n\n')
        picard_logfiles = picard_out.logfiles
        for eachlog in picard_logfiles:
            if 'log.mapstats' in eachlog:
                mapstats_logfile = eachlog
        metadata   = np.loadtxt(mapstats_logfile,dtype={'names':('UT','HST','Obs','Source','Mode','Filter',
                                                                     'El','Airmass','Trans','Tau225','Tau','t_elapsed',
                                                                     't_exp','rms','rms_units','nefd','nefd_units','RA','DEC',
                                                                     'mapsize','pixscale','project','recipe','filename'),     
                                                            'formats':('f8','<U23','f8','<U10','<U10','f8','f8','f8','f8','f8',
                                                                       'f8','f8','f8','f8','<U11','f8','<U8','<U11','<U11','f8','f8','<U8','<U11','<U34')})

        # Put the dates in the correct format to be converted to JD
        # or to be easily read in by Python's astropy.time.Time function
        source     = metadata['Source']                                  
        if len(str(int(metadata['Obs']))) == 1:
            scan   = '0000'+str(int(metadata['Obs']))
        elif len(str(int(metadata['Obs']))) == 2:
	            scan   = '000'+str(int(metadata['Obs']))
        elif len(str(int(metadata['Obs']))) == 3:
	            scan   = '00'+str(int(metadata['Obs']))
        elif len(str(int(metadata['Obs']))) == 4:
	            scan   = '0'+str(int(metadata['Obs']))
        elif len(str(int(metadata['Obs']))) == 5:
	            scan   = str(int(metadata['Obs']))
        eachUTdate = str(metadata['UT'])                                                                                                                  
        year       = eachUTdate[0:4]                                                                                                                      
        month      = eachUTdate[4:6]                                                                                                                      
        day        = eachUTdate[6:8]                                                                                                                      
        hours      = float('0.'+eachUTdate.split('.')[-1])*24                                                                                             


        hour       = str(int(np.floor(hours)))                                                                                                            
        if int(hour)>0:
            minutes    = (hours % float(hour))*60                                                                                                             
        else:
            minutes    = hours*60


        minute     = str(int(np.floor(minutes)))
        if int(minute)>0:
            second     = str(round((minutes % float(minute))*60,4))                                                                                           
        else:
            second     = str(round(minutes*60,4))
        if float(hour)<10: 
            isotUThour = '0'+hour
        else:
            isotUThour = hour


        if float(minute)<10:
            isotUTminute = '0'+minute
        else:
            isotUTminute = minute

        if float(second)<10:
            isotUTsecond = '0'+second
        else:
            isotUTsecond = second

        # Construct the astropy.table columns
        IDs.append(dummyID)
        names.append(region+'_'+year+month+day+'_'+scan+'_'+wave)
        UTdates.append(year+'-'+month+'-'+day+'T'+isotUThour+':'+isotUTminute+':'+isotUTsecond)
        #print('\n\n\n',year,month,day,isotUThour,isotUTminute,isotUTsecond,'\n\n\n')
        julian_dates.append(str(Time(year+'-'+month+'-'+day+'T'+isotUThour+':'+isotUTminute+':'+isotUTsecond,format='isot',scale='utc').jd))
        scans.append(str(int(metadata['Obs'])))
        elevs.append(str(int(metadata['El'])))
        tau225s.append(str(metadata['Tau225']))
        #noises.append(str(round(mjypbmfactor*metadata['rms'],2))) # 2021-05-26 -- I guess the variance map doesn't get updated when I cmult?!
        noises.append(str(round(float(metadata['rms']),2))) # 2021-07-21 -- NOISES OFF BY THE mjypbmfactor!!
        noise_units.append(str(metadata['rms_units']))

    # If the table already exists, read it in and append the row
    if table_already_exists:
        if anythingnew==1:
            for eachname,eachUTdate,eachJD,eachscan,eachEL,eachTau225,eachnoise,eachnoiseu in zip(names, UTdates, julian_dates, scans, elevs, tau225s, noises, noise_units):
                t.add_row([t['ID'][-1]+1, [eachname], [eachUTdate], [eachJD], [eachscan], [eachEL], [eachTau225], [eachnoise], [eachnoiseu]])
            #t.add_row([t['ID'][-1]+1, names, UTdates, julian_dates, scans, elevs, tau225s, noises, noise_units])
    # Build the Table if it doesn't already exist:
    else:
        t = Table([IDs, names, UTdates, julian_dates, scans, elevs, tau225s, noises, noise_units], names=('ID', 'Name', 'UT', 'JD', 'Obs', 'Elev', 'Tau225', 'RMS', 'RMS_unit'), meta={'name': 'Meta Data Table'})

    if len(input_data)>0:
        # Get today's date into a string
        now   = datetime.datetime.now()
        today = '{:02d}'.format(now.year)+'{:02d}'.format(now.month)+'{:02d}'.format(now.day)
        #apascii.write(t,region+'_metadata_'+today+'_'+wave+'.txt')

        # Get most recent file's date for name
        mostrecentdate = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[1]
        mostrecentscan = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[2]
        apascii.write(t,region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+fileending)

        os.system('mv -f '+region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+fileending+' '+output_dir)

    # Clean up the working directory the Starlink pywrapper creates
    os.system('rm -rf PICARDworkin*')

    return(t)

########################
########################
########################

def TCIdentifySources(coadd):                                                                       
    '''
    This program runs FellWalker in the same way
    JSA_Catalogue runs FellWalker for the peak sources.

    We will run FellWalker on the coadd in order to determine
    peak flux positions.

    #################

    IN:

    coadd = .sdf coadd file created using TCCoadd

    #################

    OUT:

    FellWalker Catalogue in the same style as a JSA_Catalogue "Peak"
    .FIT file

    #################
    #################

    '''                                                                                             
    from starlink import cupid
    import subprocess
    import numpy as np
    import astropy.io.fits as apfits
    import os

    # derive the noise of the coadd while
    # simultaneously making a FITS version

    #Convert the .sdf to .fits.
    fits_name = coadd.split('.sdf')[0]+'.fits'
    sdf_conv_command = '$CONVERT_DIR/ndf2fits '+coadd+' !'+fits_name
    subprocess.call(sdf_conv_command, shell=True)

    #Read in the variance extension of the file.
    vardata = apfits.getdata(fits_name, 1)

    #Find the noise in the central portion of the map.
    vardata_rad = vardata.shape[1]/2
    #Pick a region in the middle of the map to find the median noise. The X and
    # Y ranges will be the same.
    vardata_lo = int((vardata_rad-(vardata_rad/2)))
    vardata_hi = int((vardata_rad+(vardata_rad/2)))
    vardata_mid = vardata[0, vardata_lo:vardata_hi, vardata_lo:vardata_hi]
    noise = np.sqrt(np.median(vardata_mid))

    outcat_name = coadd.split('_EA')[0]+'_coadd_cat.FIT'

    # Run FellWalker
    cupid.findclumps(coadd,coadd.split('.sdf')[0]+'_out',method='FellWalker',config='^config/FellWalkerIn.txt',rms=noise,outcat=outcat_name,wcspar=True)

    os.system('rm -f '+coadd.split('.sdf')[0]+'_out.sdf')

    return(outcat_name)

###########################
###########################
###########################

def TCYSOcompare(peakcat,protocatalogue,diskcatalogue,region,wave='850'):                    
    '''
    This program associates each peak in the coadded image
    with the nearest Protostar/Protostar Candidate and
    Class II object.

    ####################

    IN:

    peakcat        = The Fellwalker peak catalogue
                     produced by TCIdentifySources

    protocatalogue = A text file describing the positions
                     of protostars and protostar candidates
                     (ID, RA, DEC, CLASS)

    diskcatalogue  = A text file describing the positions
                     of Class II Disc objects (ID, RA, DEC, CLASS)

    wave           = keyword. The wavelength: '850' is default. This is for the name of the output file 

    ######################

    OUT:

    Astropy table with the following columns:

    source index, ra, dec, distance to nearest protostar ,protostar class (Megeath et al. 2012 style classes for Protostars and Candidates), Distance to the nearest Class II source, Disk Class (Megeath et al 2012. style: D)

    ######################
    ######################
    '''                                                                                                                     
    from astropy.io import fits                                                                                             
    import numpy as np                                                                                                      
    from astropy.table import Table                                                                                         
    from astropy.io import ascii as apascii
    import datetime

    # Read in the Fellwalker catalogue and the YSO catalogues
    peak_cat      = fits.getdata(peakcat) 
    proto_cat     = np.loadtxt(protocatalogue,dtype=[('index','i'),('ra','f8'),('dec','f8'),('class','<U2')])
    disk_cat      = np.loadtxt(diskcatalogue,dtype=[('index','i'),('ra','f8'),('dec','f8'),('class','<U2')]) 
    classII_ind   = np.where(disk_cat['class'] == 'D')

    # Generate empty lists for the information we would like to collect
    peak_int_ind   = []
    peak_ind       = []
    ra             = []
    dec            = []
    proto_distance = []
    proto_class    = []
    disk_distance  = []
    disk_class     = []
    
    # For each source in the FellWalker catalogue,
    for eachsource in range(len(peak_cat)):

        peak_int_ind.append(eachsource)

        # Identify the index of the closest protostar in the protostar catalogue
        closest_YSO_ind = (np.sqrt(np.abs((proto_cat['ra'] - peak_cat['RA'][eachsource])*np.cos(peak_cat['DEC'][eachsource]*np.pi/180.0))**2.0+np.abs(proto_cat['dec'] - peak_cat['DEC'][eachsource])**2.0)).argmin()

        # Calculate the distance to each protostar in the protostar catalogue
        closest_YSO_dists   = 3600.0*np.sqrt(np.abs((proto_cat['ra'] - peak_cat['RA'][eachsource])*np.cos(peak_cat['DEC'][eachsource]*np.pi/180.0))**2.0+np.abs(proto_cat['dec'] - peak_cat['DEC'][eachsource])**2.0) # In arcseconds!

        # Save the closest protostar distance and the Megeath et al 2012 style Class (P, F, R)
        proto_distance.append(closest_YSO_dists[closest_YSO_ind])
        proto_class.append(proto_cat['class'][closest_YSO_ind])

        # Save information about the source
        peak_ind.append(peak_cat['ID'][eachsource])
        ra.append(peak_cat['RA'][eachsource])
        dec.append(peak_cat['DEC'][eachsource])


        # Identify the index of the closest disc in the disc catalogue
        closest_YSO_ind = (np.sqrt(np.abs((disk_cat['ra'][classII_ind] - peak_cat['RA'][eachsource])*np.cos(peak_cat['DEC'][eachsource]*np.pi/180.0))**2.0+np.abs(disk_cat['dec'][classII_ind] - peak_cat['DEC'][eachsource])**2.0)).argmin()

        # Calculate the distance to each disc in the disc catalogue
        closest_YSO_dists   = 3600.0*np.sqrt(np.abs((disk_cat['ra'][classII_ind] - peak_cat['RA'][eachsource])*np.cos(peak_cat['DEC'][eachsource]*np.pi/180.0))**2.0+np.abs(disk_cat['dec'][classII_ind] - peak_cat['DEC'][eachsource])**2.0) # In arcseconds!

        # Save the closest disc distance and the Megeath et al 2012 style Class (D)
        disk_distance.append(closest_YSO_dists[closest_YSO_ind])
        disk_class.append(disk_cat['class'][classII_ind][closest_YSO_ind])

    # Construct the astropy table and return it
    t = Table([peak_int_ind,peak_ind,ra,dec,proto_distance,proto_class,disk_distance,disk_class], names=('Index','ID','RA','DEC', 'Proto_Dist', 'Proto_Class', 'Disk_Dist', 'Disk_Class'), meta={'name': 'YSO Table'})

    # Get today's date into a string
    now   = datetime.datetime.now()
    today = '{:02d}'.format(now.year)+'{:02d}'.format(now.month)+'{:02d}'.format(now.day)
    apascii.write(t,region+'_YSOcompare_'+wave+'.txt')

    return(t)

####################
####################
####################

def TCTrackSources(input_data,peakcat,region,output_dir,aperture_diam = 0.00083333,wave='850',WEIGHTED=False,GOODBOX=False):  
    '''                                                                                                                                                                                                                            
    This program uses the results of TCIdentifySources
    to obtain the peak brightness for each source in each individual epoch. 

    ###################

    IN:

    input_data    = list or array containing the full path names of
                    the input sdf files to be coadded.

    peakcat       = The Fellwalker peak catalogue
                    produced by TCIdentifySources

    region        = The Transient field being analysed

    output_dir    = Where all the final files are kept

    aperture_diam = keyword. The diameter of the aperture, in degrees,
                             that is placed at the peak pixel location
	                     for each epoch. This should be set to 1 pixel 
		             (0.0008333 degrees = 3" pixels),
		             as it only returns the total flux. So, we obtain
		             the total flux of the peak pixel.

    ###################

    OUT:

    results_dict = A dictionary where each key corresponds
                   to a source in the FellWalker Peak catalogue.
                   Each key contains its own dictionary, containing
		   the peak fluxes measured on each julian date.

    ###################
    ###################
    '''
    import astropy.io.fits as apfits
    from astropy.coordinates import SkyCoord
    from astropy.time import Time  
    from astropy.table import Table
    from starlink import kappa              
    from starlink import convert
    #from starlink.ndfpack import Ndf        
    import os
    import pickle
    import glob

    # Check to see if a previous sourceinfo dict has been created by TCCheck4Variables:
    previous_results=False
    if WEIGHTED:
        if GOODBOX:
            fileending = wave+'_Wcal_GoodMaps_sourceinfo.txt'
        else:
            fileending = wave+'_Wcal_sourceinfo.txt'
    else:
        fileending = wave+'_sourceinfo.txt'
    previous_results_list = sorted(glob.glob(output_dir+'/*'+fileending))
    previous_results_thisregion = []
    if len(previous_results_list)>0:
        #print('\n\n\n\n\nPREVIOUS SOURCEINFO FILE FOUND\n\n\n\n\n')
        previous_results=True
        previous_results_thisregion.append(previous_results_list[-1])

    datesinfile = []
    alldatescans = []
    if previous_results:
       sourceinfotable = Table.read(sorted(previous_results_thisregion)[-1],format='ascii') 
       for eachcolname in sourceinfotable.colnames:
           if eachcolname[0:2]=='f_':
               datesinfile.append(eachcolname.split('_')[1])
               alldatescans.append(eachcolname.split('_')[1]+'_'+eachcolname.split('_')[2])
       #print('DATES IN PREVIOUS FILE:',datesinfile,'\n\n\n')

    peak_cat    = apfits.getdata(peakcat)

    # Initialise a dictionary to keep track of the data sensically.
    # The dictionary will be organised by source ID. Each Source ID
    # Will contain the peak flux value and the associated dates    

    results_dict = {}

    for eachsource in range(len(peak_cat)):
        results_dict[peak_cat['ID'][eachsource]]               = {}
        results_dict[peak_cat['ID'][eachsource]]['peakfluxes'] = []
        results_dict[peak_cat['ID'][eachsource]]['dates']      = []
        results_dict[peak_cat['ID'][eachsource]]['dates_reg']  = []
        results_dict[peak_cat['ID'][eachsource]]['scan']       = []
        results_dict[peak_cat['ID'][eachsource]]['index']      = eachsource

    if previous_results:
        # Now loop over all the date_scans in order and add previously known results to the results_dict
        for eachdatescan in alldatescans:
            date_reg  = eachdatescan.split('_')[0]
            scan      = eachdatescan.split('_')[1]
            jd        = Time(date_reg[0:4]+'-'+date_reg[4:6]+'-'+date_reg[6:8]+'T00:00:00.00',format='isot').jd
            # Loop over all sources:
            for eachsource in range(len(peak_cat)):
                peak_flux = sourceinfotable['f_'+date_reg+'_'+scan][eachsource]

                # Save the peak fluxes and their associated dates for this source, do this for every source
                results_dict[peak_cat['ID'][eachsource]]['peakfluxes'].append(peak_flux)
                results_dict[peak_cat['ID'][eachsource]]['dates'].append(jd)
                results_dict[peak_cat['ID'][eachsource]]['dates_reg'].append(date_reg)
                results_dict[peak_cat['ID'][eachsource]]['scan'].append(scan)
 
        



    # Now loop over each file --- THIS WILL APPEND COLUMNS AT THE END - SO IF WE ARE RE-RUNNING OLD DATES THIS WILL PUT THE LIST OUT OF ORDER!
    for eachfile in input_data:
#        myndf    = Ndf(eachfile)
#
#        # Make a line break separated string of the header
#        myndf_string_list = []
#        for i in myndf.head['FITS']:
#            myndf_string_list.append(i.decode('UTF-8'))
#
#        # myndf.head['FITS'] is a fits header that is given as a list of strings
#        # so split up the strings into values and comments and load it into a
#        # standard astropy fits header
#
#        hdr = apfits.Header.fromstring('\n'.join(myndf_string_list), sep='\n')

        # "Ndf" not working (bug with unknown fix time), so use convert and get the header from a fits file

        convert.ndf2fits(eachfile,out=eachfile.split('.sdf')[0]+'.fits') 

        hdr = apfits.getheader(eachfile.split('.sdf')[0]+'.fits')

        os.system('rm -f '+eachfile.split('.sdf')[0]+'.fits')

        # Obtain the UT date and the Julian Date
        date = hdr['DATE-OBS']
        scani= hdr['OBSNUM']
        if len(str(scani))==1:
            scan = "0000"+str(scani)
        elif len(str(scani)) == 2:
            scan = "000"+str(scani)
        elif len(str(scani)) == 3:
            scan = "00"+str(scani)
        elif len(str(scani)) == 4:
            scan = "0"+str(scani)
        elif len(str(scani)) == 5:
            scan = str(scani)
        date_reg=date[0:4]+date[5:7]+date[8:10]
        jd  = Time(date,format='isot',scale='utc').jd

        if date_reg not in datesinfile:
    
            # For each soure in the FellWalker catalogue,
            for eachsource in range(len(peak_cat)):
    
                # Get the PeakX and peakY coordinates in a form Starlink's KAPPA APERADD function will understand (sexigesmal)
                peak_coord_degrees = SkyCoord(peak_cat['RA'][eachsource],peak_cat['DEC'][eachsource],frame='icrs',unit='deg')
                peakX              = str(int(peak_coord_degrees.ra.hms.h))+':'+str(int(peak_coord_degrees.ra.hms.m))+':'+str(peak_coord_degrees.ra.hms.s)
                # peak_coord_degrees.dec.dms.d might come out to be: -0.0. But int(-0.0) = 0 and we lose the parity information. Fix that.
                if int(peak_coord_degrees.dec.dms.d) == 0:
                    isneg     = 0
                    firstchar = str(peak_coord_degrees.dec.dms.d)[0]
                    if firstchar == '-':
                        isneg = 1
                    
                    if isneg == 0:
                        peakY              = str(int(peak_coord_degrees.dec.dms.d))+':'+str(abs(int(peak_coord_degrees.dec.dms.m)))+':'+str(abs(peak_coord_degrees.dec.dms.s))
                    else:
                        peakY              = '-'+str(int(peak_coord_degrees.dec.dms.d))+':'+str(abs(int(peak_coord_degrees.dec.dms.m)))+':'+str(abs(peak_coord_degrees.dec.dms.s))
                else:
                    peakY              = str(int(peak_coord_degrees.dec.dms.d))+':'+str(abs(int(peak_coord_degrees.dec.dms.m)))+':'+str(abs(peak_coord_degrees.dec.dms.s))
    
                # Obtain the peak flux by placing an aperture on that location - THIS IS THE BIGGEST TIME SINK OF THE TRIGGER CODE
                peak_flux = kappa.aperadd(eachfile,centre='"'+peakX+','+peakY+'"',diam=aperture_diam).total
                results_dict[peak_cat['ID'][eachsource]]['peakfluxes'].append(peak_flux)
                results_dict[peak_cat['ID'][eachsource]]['dates'].append(jd)
                results_dict[peak_cat['ID'][eachsource]]['dates_reg'].append(date_reg)
                results_dict[peak_cat['ID'][eachsource]]['scan'].append(scan)
        #else:
        #    for eachsource in range(len(peak_cat)):
        #        peak_flux = sourceinfotable['f_'+date_reg+'_'+scan][eachsource]
        #
        #        # Save the peak fluxes and their associated dates for this source, do this for every source
        #        results_dict[peak_cat['ID'][eachsource]]['peakfluxes'].append(peak_flux)
        #        results_dict[peak_cat['ID'][eachsource]]['dates'].append(jd)
        #        results_dict[peak_cat['ID'][eachsource]]['dates_reg'].append(date_reg)
        #        results_dict[peak_cat['ID'][eachsource]]['scan'].append(scan)

    # Return the peak fluxes and dates organised by each source
    if WEIGHTED:
        if GOODBOX:
            sourcecat_name = output_dir+'/'+region+'_'+wave+'_Wcal_GoodMaps_sourcecat.bin'
        else:
            sourcecat_name = output_dir+'/'+region+'_'+wave+'_Wcal_sourcecat.bin'
    else:
        sourcecat_name = output_dir+'/'+region+'_'+wave+'_sourcecat.bin'
    pickle.dump(results_dict,open(sourcecat_name,'wb'))
    return(results_dict)

############################
############################
############################

def TCCheck4Variables(source_dict,YSOtable,region,trigger_thresh = 5,brightness_thresh = 0.0,sd_thresh = 4, wave ='850',WEIGHTED=False,GOODBOX=False): 
    '''
    This code uses data from TCTrackSources and TCYSOcompare
    to hunt for variability in the newly obtained images.

    ####################

    IN:

    source_dict       = The dictionary produced by TCTrackSources
                        containing peak fluxes and dates for every
		        source identified by FellWalker in TCIdentifySources.

    YSOtable          = The associated closest YSOs and their respective distances
                        from each source identified by FellWalker. This is the table
   		        produced by TCYSOcompare

    region            = The region we are analysing. This is to name the variables file.

    trigger_thresh    = Keyword. The threshold defined to identify a source as variable:
    		                 abs(peak flux measurement - average peak flux measurement) / SD in peak flux measurements > trigger_thresh = Variable

    brightness_thresh = keyword. The minimum brightness we will consider sources when measuring their potentially discrepant standard deviation
                                 in peak flux. In Jy/beam

    sd_thresh         = keyword. The source SD is computed and compared to the fiducial model presented in Johnstone et al. 2018. If sd_measured / sd_fiducial > sd_thresh,
                                 the source will be marked as a potential variable. Only sources with brightnesses greater than brightness_thresh will be considered.

    wave              = keyword. Wavelength. Default is '850'. This if for the names of the output file.

    ####################

    OUT:

    This program prints information about each source identified as a variable to a text file called "REGION_variables_DATE_WAVELENGTH.txt" with today's date. It also compiles a table with the following columns:

    Source ID, ra, dec, distance to nearest protostar, distance to nearest Class II object, Mean Peak Flux, Relative Standard Deviation in Peak Fluxes, Slope of linear fit of Peak fluxes over time, Error in the slope, Intercept of the linear fit of peak fluxes over time, Error in the Intercept of the linear fit of peak fluxes over time

    and returns the source IDs of Identified Variables.

    ####################
    ####################
    '''                                                                                                                                              
    from astropy.table import Table                                                                                                                  
    from astropy.io import ascii as apascii
    import numpy as np                                                                                                                               
    import datetime
    import matplotlib.pyplot as plt

    # First, correlate Source IDs with known special names
    # so when we send out the announcement - we'll know, for instance,
    # if the source is EC53

    special_names={}
    special_names['JCMTPP_J182951.2+011638']='EC53'
    special_names['JCMTPP_J053522.4-050111']='HOPS88'
    special_names['JCMTPP_J053525.2-051914']='V1017Ori'
    special_names['JCMTPP_J053527.4-050929']='HOPS370'
    special_names['JCMTPP_J182952.0+011550']='SMM10'
    special_names['JCMTPP_J182949.8+011520']='SMM1'
    special_names['JCMTPP_J162624.4-241616']='OPH162624'
    special_names['JCMTPP_J182937.8-015103']='IRAS18270-0153'
    special_names['JCMTPP_J054631.0-000232']='HOPS373'

    special_names_keys = []
    for i in special_names.keys():
        special_names_keys.append(i)

    # Get the name of one of the sources to help construct the astropy table

    sourcename_fortable = list(source_dict.keys())[0]  
    number_of_epochs    = len(source_dict[sourcename_fortable]['dates'])
    number_of_sources   = len(list(source_dict.keys()))
    #print('\n\n\n\nSOURCE_DICT DATES IN IMPORTANT CODE',source_dict[sourcename_fortable]['dates'],'\n\n\n\n')


    # Get today's date into a string
    now   = datetime.datetime.now()
    today = '{:02d}'.format(now.year)+'{:02d}'.format(now.month)+'{:02d}'.format(now.day)

    # Generate a tuple of the names we need for the table (different every time because we are adding flux values when new epochs are observed)

    name_list=['Index','ID','RA','DEC','proto_dist','disk_dist']
    dtype_list=['i4','<U23','f8','f8','f8','f8',]

    for eachname in ['mean_peak_flux','sd_peak_flux','sd_fiducial','sd/sd_fiducial','slope(%/year)','delta_slope','abs(slope/delta_slope)','intercept','delta_intercept']:
        name_list.append(eachname)
        dtype_list.append('f8')

    for eachflux in range(len(np.array(source_dict[sourcename_fortable]['peakfluxes'])[np.argsort(source_dict[sourcename_fortable]['dates'])])):
        name_list.append('f_'+sorted(source_dict[sourcename_fortable]['dates_reg'])[eachflux]+'_'+np.array(source_dict[sourcename_fortable]['scan'])[np.argsort(source_dict[sourcename_fortable]['dates'])][eachflux])
        dtype_list.append('f8')

    for eachflux in range(len(np.array(source_dict[sourcename_fortable]['peakfluxes'])[np.argsort(source_dict[sourcename_fortable]['dates'])])):
        name_list.append('abs(f_'+sorted(source_dict[sourcename_fortable]['dates_reg'])[eachflux]+'_'+np.array(source_dict[sourcename_fortable]['scan'])[np.argsort(source_dict[sourcename_fortable]['dates'])][eachflux]+'-mean_peak_flux)/sd_peak_flux')
        dtype_list.append('f8')

    name_tuple=tuple(name_list)
    dtype_tuple=tuple(dtype_list)

    # Generate an astropy table:
    t = Table(names=name_tuple,dtype=dtype_tuple)

    trigger_messages_stochastic       = []
    trigger_messages_stochastic_newind= []
    trigger_messages_stochastic_epochs= []

    trigger_messages_fiducial         = []

    triggered_sources = []

    # Initialise the file that will contain information about each identified variable
    # Get most recent file's date for name
    mostrecentdate = sorted(source_dict[sourcename_fortable]['dates_reg'])[-1]#np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[1]
    mostrecentscan = np.array(source_dict[sourcename_fortable]['scan'])[np.argsort(np.array(source_dict[sourcename_fortable]['dates_reg']))][-1]#np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[2]
    variables_file    = open(region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+wave+'_variables.txt','w')

    variables_file.write('Hello Everyone,\n\nAs of '+today+', the '+region+' region has '+str(number_of_epochs)+' Transient Survey epochs.\n\nWe are tracking '+str(number_of_sources)+' sources in this region.\n\nHere are the latest results from the automatic variability detection pipeline:')

    # For each source, calculate the mean flux and standard deviation
    for eachsource in source_dict.keys():                            
        source_dict[eachsource]['trigger']     = []
        source_dict[eachsource]['mean']        = np.average(source_dict[eachsource]['peakfluxes'])
        source_dict[eachsource]['sd']          = np.sqrt(sum((source_dict[eachsource]['peakfluxes']-np.average(source_dict[eachsource]['peakfluxes']))**2.0)/(len(source_dict[eachsource]['peakfluxes'])-1))
        source_dict[eachsource]['sd_rel']      = np.sqrt(sum((source_dict[eachsource]['peakfluxes']-np.average(source_dict[eachsource]['peakfluxes']))**2.0)/(len(source_dict[eachsource]['peakfluxes'])-1))/np.average(source_dict[eachsource]['peakfluxes'])
        if wave == '850':
            source_dict[eachsource]['sd_fiducial'] = np.sqrt(14**2.0+(0.02*np.average(source_dict[eachsource]['peakfluxes']))**2.0) # mJy -- should be 14 rather than 0.014
        if wave == '450':
            source_dict[eachsource]['sd_fiducial'] = np.sqrt((0.05*np.average(source_dict[eachsource]['peakfluxes']))**2.0) # mJy -- 5\% cutoff indicates the faint sources are noisy and we don't know what they are necessarily doing on the sd vs mean_flux (in linear,xlog space) graph

        #The diagonal elements of cov are the variances of the coefficients in z, i.e. np.sqrt(np.diag(cov)) gives you the standard deviations of the coefficients. You can use the standard deviations to estimate the probability that the absolute error exceeds a certain value, e.g. by inserting the standard deviations in the uncertainty propagation calculation. If you use e.g. 3*standard deviations in the uncertainty propagation, you calculate the error which will not be exceeded in 99.7% of the cases.

        p,cov           = np.polyfit(np.array(sorted(source_dict[eachsource]['dates']))-sorted(source_dict[eachsource]['dates'])[0],np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/np.average(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]),1,cov=True,full=False)

        slope                = p[0]
        delta_slope          = np.sqrt(np.diag(cov))[0]
        intercept            = p[1]
        delta_intercept      = np.sqrt(np.diag(cov))[1]

        # Append the slope, intercept, and associated uncertainties to the dictionary
        source_dict[eachsource]['slope']           = 100*slope*365.24 # Now in %/year
        source_dict[eachsource]['delta_slope']     = 100*delta_slope*365.24
        source_dict[eachsource]['intercept']       = intercept
        source_dict[eachsource]['delta_intercept'] = delta_intercept

        YSOtableind = np.where(YSOtable['ID']==eachsource)[0][0]
        ra         = YSOtable['RA'][YSOtableind]
        dec        = YSOtable['DEC'][YSOtableind]
        proto_dist = YSOtable['Proto_Dist'][YSOtableind]
        disk_dist  = YSOtable['Disk_Dist'][YSOtableind]

#########
#########

        # The actual trigger! abs(flux - flux_m)/SD > trigger_thresh
        
        for eachmeasurement in range(len(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))])):

            peakfluxes_not_including_current_measurement = []
            for eachpeakflux in range(len(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))])):
                if eachpeakflux != eachmeasurement:
                    peakfluxes_not_including_current_measurement.append(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachpeakflux])

            peakfluxes_not_including_current_measurement = np.array(peakfluxes_not_including_current_measurement)
            mean_not_including_current_measurement       = np.average(peakfluxes_not_including_current_measurement)
            SD_not_including_current_measurement         = np.sqrt(sum((peakfluxes_not_including_current_measurement-mean_not_including_current_measurement)**2.0)/(len(peakfluxes_not_including_current_measurement)-1))

            trigger = abs(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement] - mean_not_including_current_measurement)/SD_not_including_current_measurement
 
            source_dict[eachsource]['trigger'].append(trigger)

            if wave == '850':
                fiducial_model_SD = np.sqrt(14**2.0+(0.02*source_dict[eachsource]['mean'])**2.0) # mJy - should be 14 rather than 0.014
            if wave == '450':
                fiducial_model_SD = np.sqrt((0.05*source_dict[eachsource]['mean'])**2.0) # mJy - 5\% contraint right now for sd vs mean peak flux graph (linear/log) -- just saying we don't know what is happening for the faint sources at 450 yet, too noisy.

            # Check to see if any of the peak flux measurements on any date have a large variance

            if trigger>=trigger_thresh:
	        # If triggered, save the source ID, print the information to the screen, and add to the variables_REGION_DATE.txt text file

                triggered_sources.append(eachsource)

                if eachsource in special_names_keys:

                    trigger_message = '\n\n####################\n####################\n\nSource '+special_names[eachsource]+' (index = '+str(source_dict[eachsource]['index'])+') has abs(flux - flux_m)/SD = '+str(round(trigger,2))+' on JD: '+str(np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement])+' = '+str(np.array(source_dict[eachsource]['dates_reg'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement])+' (Epoch '+str(eachmeasurement+1)+'/'+str(number_of_epochs)+')\nThis is greater than the current abs(flux - flux_m)/SD threshold: '+str(trigger_thresh)+'.\nMean Source Brightness: '+str(round(source_dict[eachsource]['mean'],4))+' Jy/beam.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\nMean Peak Brightness = '+str(source_dict[eachsource]['mean'])+'\nSD                   = '+str(source_dict[eachsource]['sd'])+'\nSD_fid               = '+str(fiducial_model_SD)+'\nSD/SD_fid            = '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,5))+'\n\n####################\n####################\n\n'
                
                else: 

                    trigger_message = '\n\n####################\n####################\n\nSource '+str(eachsource)+' (index = '+str(source_dict[eachsource]['index'])+') has abs(flux - flux_m)/SD = '+str(round(trigger,2))+' on JD: '+str(np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement])+' = '+str(np.array(source_dict[eachsource]['dates_reg'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement])+' (Epoch '+str(eachmeasurement+1)+'/'+str(number_of_epochs)+')\nThis is greater than the current abs(flux - flux_m)/SD threshold: '+str(trigger_thresh)+'.\nMean Source Brightness: '+str(round(source_dict[eachsource]['mean'],4))+' Jy/beam.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\nMean Peak Brightness = '+str(source_dict[eachsource]['mean'])+'\nSD                   = '+str(source_dict[eachsource]['sd'])+'\nSD_fid               = '+str(fiducial_model_SD)+'\nSD/SD_fid            = '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,5))+'\n\n####################\n####################\n\n'

                trigger_messages_stochastic.append(trigger_message)
                trigger_messages_stochastic_epochs.append(eachmeasurement+1)
                trigger_messages_stochastic_date = np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement]
                if trigger_messages_stochastic_date == np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][-1]:
                    trigger_messages_stochastic_newind.append(1)
                else:
                    trigger_messages_stochastic_newind.append(0)

                plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label=str(source_dict[eachsource]['index'])+':'+eachsource.replace('_',''))
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])+source_dict[eachsource]['sd'],color='k',linestyle='dashed')
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])+source_dict[eachsource]['sd_fiducial'],color='b',linestyle='dotted')

                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes']),color='k',linestyle='solid')
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])-source_dict[eachsource]['sd'],color='k',linestyle='dashed')
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])-source_dict[eachsource]['sd_fiducial'],color='b',linestyle='dotted')
                plt.legend(loc='lower left')
                plt.savefig('index_'+str(source_dict[eachsource]['index'])+'_'+eachsource+'_'+wave+'_ACcal_lightcurve.pdf',format='pdf')
                plt.clf()
##########
##########

        # Now check for the other indicator of variability - a large relative standard deviation relative to the fiducial model in Doug's Midproject paper (like, for instance, EC53)
        # Also, a birghtness_threshold is coded in just in case we get too many spurious non-detections from faint sources and we want to get rid of those

        if wave == '850':
            fiducial_model_SD = np.sqrt(14**2.0+(0.02*source_dict[eachsource]['mean'])**2.0) # mJy -- should be 14 rather than 0.014
        if wave == '450':
            fiducial_model_SD = np.sqrt((0.05*source_dict[eachsource]['mean'])**2.0) # mJy - 5\% contraint right now for sd vs mean peak flux graph (linear/log) -- just saying we don't know what is happening for the faint sources at 450 yet, too noisy.
        
        source_dict[eachsource]['sd_fiducial'] = fiducial_model_SD

        source_dict[eachsource]['sd_fiducial_trigger'] = source_dict[eachsource]['sd']/fiducial_model_SD

        if np.logical_and(source_dict[eachsource]['mean']>brightness_thresh,source_dict[eachsource]['sd']/fiducial_model_SD>sd_thresh):

                triggered_sources.append(eachsource)

                if eachsource in special_names_keys:

                    trigger_message = '\n\n####################\n####################\n\nSource '+special_names[eachsource]+' (index = '+str(source_dict[eachsource]['index'])+') has an SD = '+str(round(source_dict[eachsource]['sd'],4))+' over all peak flux measurements.\n\nThis is greater than the fiducial SD model would predict for a mean brightness of '+str(round(source_dict[eachsource]['mean'],3))+' Jy/beam by a factor of '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,2))+'.\nThe current SD/SD_fiducial threshold is set to '+str(sd_thresh)+'.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\nMean Peak Brightness = '+str(source_dict[eachsource]['mean'])+'\nSD                   = '+str(source_dict[eachsource]['sd'])+'\nSD_fid               = '+str(fiducial_model_SD)+'\nSD/SD_fid            = '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,5))+'\n\n####################\n####################\n\n'

                else:

                    trigger_message = '\n\n####################\n####################\n\nSource '+str(eachsource)+' (index = '+str(source_dict[eachsource]['index'])+') has an SD = '+str(round(source_dict[eachsource]['sd'],4))+' over all peak flux measurements.\n\nThis is greater than the fiducial SD model would predict for a mean brightness of '+str(round(source_dict[eachsource]['mean'],3))+' Jy/beam by a factor of '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,2))+'.\nThe current SD/SD_fiducial threshold is set to '+str(sd_thresh)+'.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\nMean Peak Brightness = '+str(source_dict[eachsource]['mean'])+'\nSD                   = '+str(source_dict[eachsource]['sd'])+'\nSD_fid               = '+str(fiducial_model_SD)+'\nSD/SD_fid            = '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,5))+'\n\n####################\n####################\n\n'

                trigger_messages_fiducial.append(trigger_message)

                plt.scatter(np.array(range(len(np.array(source_dict[eachsource]['dates_reg'])[np.argsort(np.array(source_dict[eachsource]['dates']))])))+1,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label=eachsource.replace('_',''))
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])+source_dict[eachsource]['sd'],color='k',linestyle='dashed')
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes']),color='k',linestyle='solid')
                plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])-source_dict[eachsource]['sd'],color='k',linestyle='dashed')
                plt.legend(loc='upper right')
                plt.savefig('index_'+str(source_dict[eachsource]['index'])+'_'+eachsource+'_'+wave+'_lightcurve.pdf',format='pdf')
                plt.clf()
##########
##########

        # Close the text file, construct the table, and return the table along with the list of variable sources

        peakfluxes_fortable = ' '.join(list(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(source_dict[eachsource]['dates'])].astype('str')))

        list_to_add_to_table = [source_dict[eachsource]['index'],eachsource,ra,dec,proto_dist,disk_dist]

        for eachprop in [source_dict[eachsource]['mean'],source_dict[eachsource]['sd'],source_dict[eachsource]['sd_fiducial'],source_dict[eachsource]['sd_fiducial_trigger'],100*slope*365.24,100*delta_slope*365.24,abs(slope/delta_slope),intercept,delta_intercept]:
            list_to_add_to_table.append(eachprop)

        for peakfluxind in range(len(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(source_dict[eachsource]['dates'])])):
            list_to_add_to_table.append(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(source_dict[eachsource]['dates'])][peakfluxind])

        for peakfluxind in range(len(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(source_dict[eachsource]['dates'])])):
            list_to_add_to_table.append(source_dict[eachsource]['trigger'][peakfluxind]) # ALREADY IN ORDER!! ADDING ARGSORT WOULD PUT THIS OUT OF ORDER!

        t.add_row(list_to_add_to_table)

    # Order the table

    t.sort('Index')

    # Order the Stochastic Variable Message Lists:

    trigger_messages_stochastic = np.array(trigger_messages_stochastic)[np.argsort(trigger_messages_stochastic_epochs)]
    trigger_messages_stochastic_newind = np.array(trigger_messages_stochastic_newind)[np.argsort(trigger_messages_stochastic_epochs)]

    if len(triggered_sources)<1:
        variables_file.write('\nNo potential variables found so far... :(\n\n')
        variables_file.close()
    else:
        variables_file.write('\n\n\nNEW INDIVIDUAL EPOCH OUTLIER DETECTIONS:\n\n')
        if len(np.array(trigger_messages_stochastic)[np.where(np.array(trigger_messages_stochastic_newind)>0)])>0:
            for i in np.array(trigger_messages_stochastic)[np.where(np.array(trigger_messages_stochastic_newind)>0)]:
                variables_file.write(i)
        variables_file.write('\n\n\nTIME SERIES STOCHASTICITY DETECTIONS:\n\n')
        if len(trigger_messages_fiducial)>0:
            for i in trigger_messages_fiducial:
                variables_file.write(i)
        variables_file.write('\n\n\nOLD INDIVIDUAL EPOCH OUTLIER DETECTIONS:\n\n')
        if len(np.array(trigger_messages_stochastic)[np.where(np.array(trigger_messages_stochastic_newind)==0)])>0:
            for i in np.array(trigger_messages_stochastic)[np.where(np.array(trigger_messages_stochastic_newind)==0)]:
                variables_file.write(i)

        variables_file.write('\n\nIf you see any issues with this message, please contact Steve Mairs at s.mairs@eaobservatory.org.\n\nHave a great day!\n\nSteve (via the automated variability detection pipeline)')
        variables_file.close()
        

    # Get today's date into a string
    now   = datetime.datetime.now()
    today = '{:02d}'.format(now.year)+'{:02d}'.format(now.month)+'{:02d}'.format(now.day)
    #apascii.write(t,region+'_sourceinfo_'+today+'_'+wave+'.txt')

    # Get most recent file's date for name - already found, above
    #mostrecentdate = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[1]
    #mostrecentscan = np.array(input_data)[np.argsort(np.array(input_data))][-1].split('/')[-1].split('_')[2]
    if WEIGHTED:
        if GOODBOX:
            fileending = wave+'_Wcal_GoodMaps_sourceinfo.txt'
        else:
            fileending = wave+'_Wcal_sourceinfo.txt'
    else:
        fileending = wave+'_sourceinfo.txt'
    apascii.write(t,region+'_'+mostrecentdate+'_'+mostrecentscan+'_'+fileending)

    return(t,triggered_sources)
