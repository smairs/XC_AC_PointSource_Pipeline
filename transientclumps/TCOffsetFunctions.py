import numpy as np
import scipy.spatial.distance as ssd
import astropy.io.fits as apfits


### cull_matches ###

# Using the offsets of all sources found within 10 arcsecond separation from one another
# Loop over the x and y offsets
#     CODA
#     Zero sources - warning - error - do not use, move on
#     One source - warning - set STD to 0, since that shows the error is infinite (if we divide it later)
#     Two sources - warning - set STD to difference between them ( delta_xi - delta_xj )
#     Three or more:
#      	n objects = n(n-1)/2 pairs
#      	calculate difference between each x offset and y offset value for each source
# 	If a source has one object that is less than 4 arcseconds away - keep it, if not, throw the source away= CULL
#        then go back to CODA
#
#   Then find averages with what's left: avg_delta_x, avg_delta_y and peak/radius ratios
#
#
# Keywords:
#
# source_xoffs - the x offsets of all the sources within 10 arcseconds from one another - make sure this has been corrected by the factor of cos(declination)
# source_yoffs - the y offsets of all the sources within 10 arcseconds from one another
# cutoff       - the maximum distance (in arcseconds) the culled sources are allowed to be apart from one another
#
# Returns:
#
# good_ind     - the index of all the sources which survived the cull - this can be applied to the offsets and peak ratios and radius ratios 
# stdx         - the error in the x offsets
# stdy         - the error in the y offsets
#


def cull_matches(source_xoffs,source_yoffs,cutoff):

    source_xoffs=np.array(source_xoffs)
    source_yoffs=np.array(source_yoffs)

    # In the case of two different sized x offset and y offset arrays:
    
    if len(source_xoffs)!=len(source_yoffs):
        print('\nERROR: THE X OFFSET AND Y OFFSETS DO NOT HAVE THE SAME NUMBER OF ELEMENTS\n')
        good_ind=np.nan
        stdx=0
        stdy=0

    # In the case of 0 sources found, produce an empty set

    if len(source_xoffs)==0:     
        print('\nWARNING: THERE WERE NO SOURCES FOUND\n')
        good_ind=[]
        stdx=0
        stdy=0
    
    # In the case of 1 source found:
    
    if len(source_xoffs)==1:
        print('\nWARNING: THERE WAS ONLY 1 SOURCE FOUND, STD = 0\n')
        good_ind=[0]
        stdx=0
        stdy=0

    # In the case of 2 sources found:
    
    if len(source_xoffs)==2:
        print('\nWARNING: THERE WERE ONLY 2 SOURCES FOUND\n')
        good_ind=[0,1]
        stdx=source_xoffs[1]-source_xoffs[0]
        stdy=source_yoffs[1]-source_yoffs[0]

    # In the case of 3 or more sources found:

    if len(source_xoffs)>2:

        #set up the index of good source
        good_ind=[]

        #loop over every source 
        for i in range(len(source_xoffs)):
            #define a list of the distances between sources on the XOFFSET YOFFSET plane 
            offset_differences=[]
            #for every source
            for j in range(len(source_xoffs)):
                #leave out the current source you are working on - because that will give a difference of zero 
                if j==i: continue
                #calculate the distance between this current source i and all the others on the XOFFSET YOFFSET plane
                offset_differences.append(np.sqrt((source_xoffs[j]-source_xoffs[i])**2.0+(source_yoffs[j]-source_yoffs[i])**2.0))
            #Set a dummy variable that we can use to indicate when we have a source that survives the cull
            dummy='dummy'
            #for each of the distances between the sources on the XOFFSET YOFFSET plane
            for k in offset_differences:
                #Is this distance less than the cutoff?
                if k<cutoff:
                    #If it is, set the dummy variable to indicate the source survives!
                    dummy='SURVIVES'
            #If the source survives, it means it has a neighbour on the offset plane that is within the cutoff distance, so add id to good_ind
            if dummy=='SURVIVES':
                #If the source is already in good_ind, then skip this
                if i in good_ind: continue
                good_ind.append(i)

        #Now we cant to check how many sources survived the cull and repeat the first part of the function for the cases of 0 sources, 1 source, 2 sources, and 3 or more sources

        if len(good_ind)==0:
            print('\nWARNING: THERE WERE NO SOURCES THAT SURVIVED THE CULL\n')
            stdx=0
            stdy=0
  
        if len(good_ind)==1:
            print('\nWARNING: THERE WAS ONLY 1 SOURCE FOUND, THAT SURVIVED THE CULL STD = 0\n')
            stdx=0
            stdy=0        
            
        if len(good_ind)==2:
            print('\nWARNING: THERE WERE ONLY 2 SOURCES THAT SURVIVED THE CULL\n')
            stdx=source_xoffs[1]-source_xoffs[0]
            stdy=source_yoffs[1]-source_yoffs[0]
       
        if len(good_ind)>2: 
          stdx=np.std(source_xoffs[good_ind])
          stdy=np.std(source_yoffs[good_ind])         


    return good_ind,stdx,stdy

### source_match ###
# Takes a gaussclumps catalog for a target image and compares the catalog
# to a reference catalog to find matches and report offsets in position and
# flux.
#
# Keywords:
# cat_name - name of the target gaussclumps catalog in .fits format (string).
#          -- required
# ref_name - name of the reference gaussclumps catalog in .fits format (string).
#          -- required
# minpeak - minimum brightness of a source to consider matching it (Jy/Beam).
# maxrad - maximum effective radius of a source to consider matching it (").
# maxsep - maximum allowed separation between two sources to match them (").
#
# Output:
# off - a 4-element array that gives the x and y offset of the target image
#           with respect to the reference image (ie: target - reference) as well
#           as peak and radius ratios (ie: target / reference). All values are
#           averages.
# err - a 4-element array that gives the standard deviation in each of the above
#           measurements, respectively.

def source_match(cat_name, ref_name, minpeak=0.2, maxrad=15, maxsep=10, cutoff=4, pix_scale=3.0):

    #Read in the reference catalog and prepare the data.
    ref_data = apfits.getdata(ref_name, 0)
    ref_nsources = len(ref_data['Cen1'])
    ref_posxvals = ref_data['Cen1']
    ref_posyvals = ref_data['Cen2']
    ref_peakvals = ref_data['Peak']
    ref_fwhm1vals = ref_data['GCFWHM1']
    ref_fwhm2vals = ref_data['GCFWHM2']

    #Calculate the effective radius, the square root of the two FWHM multiplied
    # together. Multiply by 1.5 to account for FWHM is diameter but want radius (divide by 2),
    # and units are pixels but want arcseconds (multiply by 3 - SHOULD MAKE THIS PIX_SCALE BECAUSE 450 DATA HAS 2 ARCSECOND PIXELS).
    ref_rvals = np.sqrt( np.multiply( ref_fwhm1vals, ref_fwhm2vals ) ) * pix_scale/2.0

    #Keep an array to track index of successful matches in the target catalog.
    matched_sources = np.zeros(ref_nsources)
    matched_sources[:] = -1

    #Read in the target catalog and prepare the data.
    targ_data = apfits.getdata(cat_name, 0)
    targ_nsources = len(targ_data['Cen1'])
    targ_posxvals = targ_data['Cen1']
    targ_posyvals = targ_data['Cen2']
    targ_peakvals = targ_data['Peak']
    targ_fwhm1vals = targ_data['GCFWHM1']
    targ_fwhm2vals = targ_data['GCFWHM2']

    #Calculate effective radius, same as above CHANGE THIS TOO, LIKE ABOVE.
    targ_rvals = np.sqrt( np.multiply( targ_fwhm1vals, targ_fwhm2vals ) ) * pix_scale/2.0

    #Keep an array to track whether or not target catalog sources have found
    # a match, because they can't match more than once. Brightest sources
    # get priority.
    targ_matchcheck = np.zeros(targ_nsources)

    #Loop over all sources in the reference catalog.
    for i in range(ref_nsources):

        #Test a radius and peak condition here:
        if ref_rvals[i] > maxrad: continue
        if ref_fwhm1vals[i]*pix_scale/2.0>maxrad:continue
        if ref_fwhm2vals[i]*pix_scale/2.0>maxrad:continue
        if ref_peakvals[i] < minpeak: continue

        #Declare the reference coordinates.
        ref_coords = np.array([[ref_posxvals[i], ref_posyvals[i]],])

        #Set a minimum distance tracking variable and the matched source
        # indexing variable.
        mindist = -1
        match_ind = -1

        #Loop over the target catalog sources.
        for j in range(targ_nsources):

            #Don't match a source that already has a match.
            if targ_matchcheck[j] == 1: continue

            #Declare the target coordinates
            targ_coords = np.array([[targ_posxvals[j], targ_posyvals[j]],])

            #Calculate distance.
            #dist = ssd.cdist(ref_coords, targ_coords, metric='euclidean')
            dist = 3600.0*np.sqrt(((ref_coords[0][0]-targ_coords[0][0])*np.cos(ref_coords[0][1]*np.pi/180.0))**2.0+(ref_coords[0][1]-targ_coords[0][1])**2.0)
            #Check if this source meets the criterion for matching distance.
            if dist < maxsep:

                #If this is the first iteration then mindist will be -1
                if mindist == -1:
                    #Auto assign mindist and match_ind
                    match_ind = j
                    mindist = dist
                ###fi

                #Now check if a new closest source.
                if dist < mindist:
                    #Remember the new source number and assign dist to mindist.
                    match_ind = j
                    mindist = dist
                    ###fi
            ###fi
        ###j

        #If we are at the end of the reference catalog then assign
        # this source if it found a match.
        if match_ind >= 0:
            targ_matchcheck[match_ind] = 1
            matched_sources[i] = match_ind
        ##fi
    ###i

    #Find out the number of matched sources.
    n_matched = len(np.where(matched_sources != -1)[0])
    print('\n'+str(n_matched)+' sources were matched!\n')

    if n_matched > 0:

        #Where matched_sources is greater than -1, the reference catalogue has a match. 
        #The value of matched_sources gives the index of the target catalogue source to which it is matched
        ref_sources_with_match=np.where(matched_sources != -1)
    

        #Get rid of all the sources that do not match in the reference catalogue
        ref_peakvals=np.array(ref_peakvals)[ref_sources_with_match]
        ref_posxvals=np.array(ref_posxvals)[ref_sources_with_match]
        ref_posyvals=np.array(ref_posyvals)[ref_sources_with_match]
        ref_rvals=np.array(ref_rvals)[ref_sources_with_match]

        #matched_sources is the array containing the index of the target catalogue's matching sources
        corresponding_targ_sources=np.array(matched_sources)[np.where(matched_sources != -1)]
        corresponding_targ_sources=corresponding_targ_sources.astype(int)

        #Get rid of all the sources that do not match in the target catalogue
        targ_peakvals=np.array(targ_peakvals)[corresponding_targ_sources]
        targ_posxvals=np.array(targ_posxvals)[corresponding_targ_sources]
        targ_posyvals=np.array(targ_posyvals)[corresponding_targ_sources]
        targ_rvals=np.array(targ_rvals)[corresponding_targ_sources]

        #Now find the offset between the matched reference and target sources as well as the peak ratios and radius ratios:
        matched_xoff=[]
        matched_yoff=[]
        matched_peakr=[]
        matched_rr=[]
        for i in range(len(np.where(matched_sources != -1)[0])):
            #Convert offsets into arcseconds^2
            matched_xoff.append(((targ_posxvals[i])-(ref_posxvals[i]))*np.cos(ref_posyvals[i]*np.pi/180.0)*3600.0)
            matched_yoff.append((targ_posyvals[i]-ref_posyvals[i])*3600.0)
            matched_peakr.append(targ_peakvals[i]/ref_peakvals[i])
            matched_rr.append(targ_rvals[i]/ref_rvals[i])

        number_of_sources_matched=len(np.where(matched_sources != -1)[0])

        #Now, find the index and the standard deviation for all the sources which survived the cull
        culled_ind,stdx,stdy=cull_matches(matched_xoff,matched_yoff,cutoff)

        print('\n\n'+str(len(culled_ind))+' sources survived the cull!\n\n')


        #Calculate averages.
        avg_xoff = np.average(np.array(matched_xoff)[culled_ind])
        avg_yoff = np.average(np.array(matched_yoff)[culled_ind])
        avg_peakr = np.average(np.array(matched_peakr)[culled_ind])
        avg_rr = np.average(np.array(matched_rr)[culled_ind])

        #Calculate error.
        std_xoff = stdx
        std_yoff = stdy
        
        if len(culled_ind)==0:       
            std_peakr = 0
            std_rr = 0

        if len(culled_ind)==1:
            std_peakr = 0
            std_rr = 0

        if len(culled_ind)==2:
            std_peakr = matched_peakr[1]-matched_peakr[0]
            std_rr = matched_rr[1]-matched_rr[0]

        if len(culled_ind)>2:
            std_peakr = np.std(matched_peakr)
            std_rr = np.std(matched_rr)

        return [avg_xoff,avg_yoff,avg_peakr,avg_rr], [std_xoff,std_yoff,std_peakr,std_rr], [len(culled_ind)], [number_of_sources_matched], np.where(matched_sources != -1), corresponding_targ_sources,np.where(matched_sources != -1)[0][culled_ind],corresponding_targ_sources[culled_ind]
    ##fi

    #If no sources found return None.
    if n_matched == 0:
        print('\nNo sources found.\n')
        return [None,None,None,None], [None,None,None,None],[0], [0], [0], [0], [0], [0]
#def
