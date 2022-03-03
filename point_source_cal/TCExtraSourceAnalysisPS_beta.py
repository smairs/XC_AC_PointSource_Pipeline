def TCYSOcompareES(peakcat,protocatalogue,diskcatalogue,region): #The ES stands for EXTRA SOURCES - this had to be built because the original YSOcomapre only accepted FellWalker Style Catalogues
    from astropy.io import fits
    from astropy.io import ascii as apascii
    import numpy as np
    from astropy.table import Table
    import datetime

    # Read in the Fellwalker catalogue and the YSO catalogues
    peak_cat      = np.loadtxt(peakcat,dtype=[('name','<U15'),('region','<U13'),('ra','f8'),('dec','f8')]) # Make sure ra and dec are in degrees
    proto_cat     = np.loadtxt(protocatalogue,dtype=[('index','i'),('ra','f8'),('dec','f8'),('class','<U2')])
    disk_cat      = np.loadtxt(diskcatalogue,dtype=[('index','i'),('ra','f8'),('dec','f8'),('class','<U2')])

    # Generate empty lists for the information we would like to collect
    peak_ind       = []
    ra             = []
    dec            = []
    proto_distance = []
    proto_class    = []
    disk_distance  = []
    disk_class     = []

    # For each source in the FellWalker catalogue,
    for eachsource in range(len(peak_cat['name'])):

        # Identify the index of the closest protostar in the protostar catalogue
        closest_YSO_ind = (np.sqrt(np.abs((proto_cat['ra'] - peak_cat['ra'][eachsource])*np.cos(peak_cat['dec'][eachsource]*np.pi/180.0))**2.0+np.abs(proto_cat['dec'] - peak_cat['dec'][eachsource])**2.0)).argmin()

        # Calculate the distance to each protostar in the protostar catalogue
        closest_YSO_dists   = 3600.0*np.sqrt(np.abs((proto_cat['ra'] - peak_cat['ra'][eachsource])*np.cos(peak_cat['dec'][eachsource]*np.pi/180.0))**2.0+np.abs(proto_cat['dec'] - peak_cat['dec'][eachsource])**2.0) # In arcseconds!

        # Save the closest protostar distance and the Megeath et al 2012 style Class (P, F, R)
        proto_distance.append(closest_YSO_dists[closest_YSO_ind])
        proto_class.append(proto_cat['class'][closest_YSO_ind])

        # Save information about the source
        peak_ind.append(peak_cat['name'][eachsource])
        ra.append(peak_cat['ra'][eachsource])
        dec.append(peak_cat['dec'][eachsource])


        # Identify the index of the closest disc in the disc catalogue
        closest_YSO_ind = (np.sqrt(np.abs((disk_cat['ra'] - peak_cat['ra'][eachsource])*np.cos(peak_cat['dec'][eachsource]*np.pi/180.0))**2.0+np.abs(disk_cat['dec'] - peak_cat['dec'][eachsource])**2.0)).argmin()

        # Calculate the distance to each disc in the disc catalogue
        closest_YSO_dists   = 3600.0*np.sqrt(np.abs((disk_cat['ra'] - peak_cat['ra'][eachsource])*np.cos(peak_cat['dec'][eachsource]*np.pi/180.0))**2.0+np.abs(disk_cat['dec'] - peak_cat['dec'][eachsource])**2.0) # In arcseconds!

        # Save the closest disc distance and the Megeath et al 2012 style Class (D)
        disk_distance.append(closest_YSO_dists[closest_YSO_ind])
        disk_class.append(disk_cat['class'][closest_YSO_ind])

    # Construct the astropy table and return it
    t = Table([peak_ind,ra,dec,proto_distance,proto_class,disk_distance,disk_class], names=('ID','RA','DEC', 'Proto_Dist', 'Proto_Class', 'Disk_Dist', 'Disk_Class'), meta={'name': 'Extra Sources YSO Table'})

    # Get today's date into a string
    now   = datetime.datetime.now()
    today = str(now.year)+str(now.month)+str(now.day)
    apascii.write(t,region+'_YSOcompare_ES_'+today+'.txt')

    return(t)


def TCTrackSourcesES(input_data,peakcat,region,aperture_diam = 0.00083333):
    import astropy.io.fits as apfits
    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    from starlink import kappa
    from starlink.ndfpack import Ndf
    import numpy as np

    peak_cat    = np.loadtxt(peakcat,dtype=[('name','<U15'),('region','<U13'),('ra','f8'),('dec','f8')])

    # Initialise a dictionary to keep track of the data sensically.
    # The dictionary will be organised by source ID. Each Source ID
    # Will contain the peak flux value and the associated dates

    results_dict = {}

    #The extra_sources peak catalogue contains sources from many different regions - so we just want to focus on the current one
    for eachsource in range(len(peak_cat['name'])):
        if peak_cat['region'][eachsource] == region:
            results_dict[peak_cat['name'][eachsource]]               = {}
            results_dict[peak_cat['name'][eachsource]]['peakfluxes'] = []
            results_dict[peak_cat['name'][eachsource]]['dates']      = []


    # Now loop over each file
    for eachfile in input_data:
        myndf    = Ndf(eachfile)

        # Make a line break separated string of the header
        myndf_string_list = []
        for i in myndf.head['FITS']:
            myndf_string_list.append(i.decode('UTF-8'))

        # myndf.head['FITS'] is a fits header that is given as a list of strings
        # so split up the strings into values and comments and load it into a
        # standard astropy fits header

        hdr = apfits.Header.fromstring('\n'.join(myndf_string_list), sep='\n')

        # Obtain the UT date and the Julian Date
        date = hdr['DATE-OBS']
        jd  = Time(date,format='isot',scale='utc').jd
        # For each soure in the FellWalker catalogue,
        for eachsource in range(len(peak_cat['name'])):
            if peak_cat['region'][eachsource] == region:
                # Get the PeakX and peakY coordinates in a form Starlink's KAPPA APERADD function will understand (sexigesmal)
                peak_coord_degrees = SkyCoord(peak_cat['ra'][eachsource],peak_cat['dec'][eachsource],frame='icrs',unit='deg')
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

                # Save the peak fluxes and their associated dates for this source, do this for every source
                results_dict[peak_cat['name'][eachsource]]['peakfluxes'].append(peak_flux)
                results_dict[peak_cat['name'][eachsource]]['dates'].append(jd)

    # Return the peak fluxes and dates organised by each source
    return(results_dict)


def TCCheck4VariablesES(source_dict,YSOtable,region,trigger_thresh = 5, brightness_thresh = 0.0,sd_thresh=4):
    from astropy.table import Table
    from astropy.io import ascii as apascii
    import numpy as np
    import datetime

    # Get today's date into a string
    now   = datetime.datetime.now()
    today = str(now.year)+str(now.month)+str(now.day)

    # Generate an astropy table:
    t = Table(names=('Name','RA','DEC','proto_dist','disk_dist','mean_peak_flux','rel_sd_peak_flux','slope','delta_slope','intercept','delta_intercept','abs(f_last - f_mean)/SD'),dtype=('<U15','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8','f8'))

    triggered_sources = []

    # Initialise the file that will contain information about each identified variable
    variables_file    = open(region+'_variables_ES_'+today+'.txt','w')

    # For each source, calculate the mean flux and standard deviation
    for eachsource in source_dict.keys():
        source_dict[eachsource]['mean']   = np.average(source_dict[eachsource]['peakfluxes'])
        source_dict[eachsource]['sd']     = np.sqrt(sum((source_dict[eachsource]['peakfluxes']-np.average(source_dict[eachsource]['peakfluxes']))**2.0)/(len(source_dict[eachsource]['peakfluxes'])-1))
        source_dict[eachsource]['sd_rel'] = np.sqrt(sum((source_dict[eachsource]['peakfluxes']-np.average(source_dict[eachsource]['peakfluxes']))**2.0)/(len(source_dict[eachsource]['peakfluxes'])-1))/np.average(source_dict[eachsource]['peakfluxes'])                                                                                   

        #The diagonal elements of cov are the variances of the coefficients in z, i.e. np.sqrt(np.diag(cov)) gives you the standard deviations of the coefficients. You can use the standard deviations to estimate the probability that the absolute error exceeds a certain value, e.g. by inserting the standard deviations in the uncertainty propagation calculation. If you use e.g. 3*standard deviations in the uncertainty propagation, you calculate the error which will not be exceeded in 99.7% of the cases.

        p,cov           = np.polyfit(sorted(source_dict[eachsource]['dates']),np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],1,cov=True,full=False)
        slope           = p[0]
        delta_slope     = np.sqrt(np.diag(cov))[0]
        intercept       = p[1]
        delta_intercept = np.sqrt(np.diag(cov))[1]

        # Append the slope, intercept, and associated uncertainties to the dictionary
        source_dict[eachsource]['slope']           = slope
        source_dict[eachsource]['delta_slope']     = delta_slope
        source_dict[eachsource]['intercept']       = intercept
        source_dict[eachsource]['delta_intercept'] = delta_intercept

        YSOtableind = np.where(YSOtable['ID']==eachsource)[0][0]
        ra         = YSOtable['RA'][YSOtableind]
        dec        = YSOtable['DEC'][YSOtableind]
        proto_dist = YSOtable['Proto_Dist'][YSOtableind]
        disk_dist  = YSOtable['Disk_Dist'][YSOtableind]

#########
#########

        # The actual trigger! abs(flux_last - flux_m)/SD > trigger_thresh

        for eachmeasurement in range(len(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))])):
            trigger = abs(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement] - source_dict[eachsource]['mean'])/source_dict[eachsource]['sd']

            # Check to see if any of the peak flux measurements on any date have a large variance

            if trigger>=trigger_thresh:
                # If triggered, save the source ID, print the information to the screen, and add to the variables_REGION_DATE.txt text file

                triggered_sources.append(eachsource)

                trigger_message = '\n\n####################\n####################\n\nSource '+str(eachsource)+' has abs(flux - flux_m)/SD = '+str(round(trigger,2))+' on JD: '+str(np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][eachmeasurement])+'\n This is greater than the current threshold: '+str(trigger_thresh)+'.\nMean Source Brightness: '+str(round(source_dict[eachsource]['mean'],4))+' Jy/beam.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\n####################\n####################\n\n'

                print(trigger_message)

                variables_file.write(trigger_message)


##########
##########

        # Now check for the other indicator of variability - a large relative standard deviation relative to the fiducial model in Doug's Midproject paper (like, for instance, EC53)
        # Also, a birghtness_threshold is coded in just in case we get too many spurious non-detections from faint sources and we want to get rid of those

        # From Johnstone et al. 2018:
        fiducial_model_SD = np.sqrt(0.014**2.0+(0.02*source_dict[eachsource]['mean'])**2.0)

        if np.logical_and(source_dict[eachsource]['mean']>brightness_thresh,source_dict[eachsource]['sd']/fiducial_model_SD>sd_thresh):

                triggered_sources.append(eachsource)

                trigger_message = '\n\n####################\n####################\n\nSource '+str(eachsource)+' has an SD = '+str(round(source_dict[eachsource]['sd'],4))+' over all peak flux measurements.\n\nThis is greater than the fiducial SD model would predict for a mean brightness of '+str(round(source_dict[eachsource]['mean'],3))+' Jy/beam by a factor of '+str(round(source_dict[eachsource]['sd']/fiducial_model_SD,2))+'.\nThe current threshold is set to '+str(sd_thresh)+'.\nThis source is located at (RA, dec) = ('+str(ra)+','+str(dec)+')\nThe nearest protostar is '+str(round(proto_dist,2))+'" away and the nearest disc is '+str(round(disk_dist,2))+'" away.\n\n####################\n####################\n\n'

                print(trigger_message)

                variables_file.write(trigger_message)


##########
##########


        # Close the text file, construct the table, and return the table along with the list of variable sources

        t.add_row([eachsource,ra,dec,proto_dist,disk_dist,source_dict[eachsource]['mean'],source_dict[eachsource]['sd_rel'],slope,delta_slope,intercept,delta_intercept,trigger]) # SAVE THE LAST TRIGGER VALUE

    if len(triggered_sources)<1:
        variables_file.write('\nNo potential variables found... :(\n\n')
        variables_file.close()
    else:
        variables_file.close()

    # Get today's date into a string
    now   = datetime.datetime.now()
    today = str(now.year)+str(now.month)+str(now.day)
    apascii.write(t,region+'_sourceinfo_ES_'+today+'.txt')

    return(t,triggered_sources)

