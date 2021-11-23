# Import necessary packages

import numpy as np
import matplotlib.pyplot as plt
import pickle
from scipy.special import comb
import glob
import itertools
from point_source_cal.noisefunctions import readnoise
from point_source_cal.get_vars import get_vars
from astropy.io import fits
import seaborn as sns
import json

##################################################
# Set up plot style for consistency for paper!
from matplotlib.pyplot import rc
rc('text', usetex=True)
rc('font',family = 'serif')
rc('font',size = 10)
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')
rc('axes',labelsize='x-small')
##################################################


def plot_SDfamsize(eachregion,coadd_cat,wave,eachtargunc,date_cutoff,jsonfile='point_source_cal/bright_threshes.json',cat_dir='config/',goodmaps=False):
    '''
    This code is used to extract families of different sizes based on a target uncertianty in the final calibration.
    '''

###########
###########
# Preliminaries
###########
###########

    # Load the brightness thresholds that have been previously defined by running the make_FCF* code to observe
    # Ensemble uncertainties versus fluxes
    with open(jsonfile) as json_file:
        brightness_threshes = json.load(json_file)

    # Define a conversion factor to ensure that the 850 micron, low mass maps are in mJy/beam -- the 850 high mass regions and 450 micron data (all regions) is already in these units 
    if np.logical_and(wave == '850',eachregion in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']):  #-- Jy/beam! This is to match the coadd catalogue
        thresh_conversion_factor = 1000.0 # To multiply catalogue numbers by so it matches the sourceinfo file (Doug's catalogues were in Jy/beam)
    else: #-- mJy/beam! This is to match the coadd catalogues! 
        thresh_conversion_factor = 1.0 # To multiply catalogue numbers by so it matches the sourceinfo file -- no correction needed in this case (Steve's cats were in mJy/beam)


    # Collect the RMS noise from each observation of this region
    region_noises = {}
    date_scans,noises = readnoise(wave+'_noises.txt',eachregion) # This is always in mJy/beam for both 450 and 850
    region_noises[eachregion] = {}
    region_noises[eachregion]['noises'] = np.array(noises)
    region_noises[eachregion]['date_scans'] = np.array(date_scans)

    # Convert the target uncdertainty to the proper units and select the appropriate brightness threshold
    # for sources to include in the family of rhis region and wavelength.
    target_perc_unc = float(eachtargunc)/100.0
    brightnessthresh = brightness_threshes[wave][eachregion][eachtargunc]

###########
# Extract Source Information
###########

    # Get all of the source information and metadata corresponding to this region and wavelength
    if goodmaps:
        sourceinfo = pickle.load(open('pointsource_results/' + eachregion + '/' + eachregion + '_'+wave+'_Wcal_GoodMaps_sourcecat.bin', 'rb'))
        metadata = np.genfromtxt(sorted(list(glob.glob('pointsource_results/' + eachregion + '/*'+wave+'*_Wcal_GoodMaps_metadata.txt')))[-1], names=True,dtype=[('ID', '<f8'), ('Name', '<U28'), ('UT', '<f8'), ('JD', '<f8'), ('Obs', '<f8'), ('Elev', '<f8'),('Tau225', '<f8'), ('RMS', '<f8'), ('RMS_unit', '<f8')])
    else:
        sourceinfo = pickle.load(open('pointsource_results/' + eachregion + '/' + eachregion + '_'+wave+'_sourcecat.bin', 'rb'))
        metadata = np.genfromtxt(sorted(list(glob.glob('pointsource_results/' + eachregion + '/*'+wave+'*_Wcal_metadata.txt')))[-1], names=True,dtype=[('ID', '<f8'), ('Name', '<U28'), ('UT', '<f8'), ('JD', '<f8'), ('Obs', '<f8'), ('Elev', '<f8'),('Tau225', '<f8'), ('RMS', '<f8'), ('RMS_unit', '<f8')])


    # Find the number of sources above the given flux threshold: 
    # This is already in mJy/bm for 450 and 850 high mass --  Jy/bm for 850, low mass
    num_sources_above_flux = len(np.array(coadd_cat['PEAK_FLUX'])[np.where(np.array(coadd_cat['PEAK_FLUX']) > brightnessthresh)]) 

    #Find the lowest brightness source that is above the threshold -- multiply by 1000 if one of 8 original regions and 850 microns since coadd_cat is in Jy/beam, but maps are in mJy/beam - 2021-01-11
    lowest_flux_above_thresh = thresh_conversion_factor*min(np.array(coadd_cat['PEAK_FLUX'])[np.where(np.array(coadd_cat['PEAK_FLUX']) > brightnessthresh)])

    # The threshold RMS for maps to consider is: Target (RMS/(Ensemble Signal)) * Ensemble signal = Target RMS. 
    RMSthresh = target_perc_unc * lowest_flux_above_thresh * np.sqrt(num_sources_above_flux) # For both wavelengths this will be in mJy/beam now
    # Percentage of maps that meet this criteria:
    perc_maps = 100 * len(region_noises[eachregion]['noises'][np.where(region_noises[eachregion]['noises'] <= RMSthresh)]) / float(len(region_noises[eachregion]['noises']))

    # Pick out the variable sources
    var_list,var_names = get_vars(eachregion,wave)
    known_variables = []
    for eachknownvariable in var_list:
        known_variables.append(str(eachknownvariable).zfill(2))
    print('\n\t{} Known Variables:'.format(eachregion),known_variables,'\n\n')

############
############
# Keep the dates that survive the RMS threshold
# We will only take seriously data from maps whose RMS 
# is low enough to actually detect these sources and we'll ignore the other data
############
############

    # Extract RMS on each Date and get the DateScan
    # This is so we can compare, observation by observation, whether
    # The map meets the RMS criterion.
    metadata_datestring = []
    rms_on_date = []
    for i in range(len(metadata['Name'])):
        metadata_datestring.append(metadata['Name'][i].split('_')[1] + '_' + metadata['Name'][i].split('_')[2])
        rms_on_date.append(metadata['RMS'][i])

    metadata_datestring = np.array(metadata_datestring)
    rms_on_date = np.array(rms_on_date)

    # Now, for each source, we want to derive the average and normalised peakfluxes along with the 
    # list of dates for this RMS limit and add it to the sourceinfo dictionary for easy access 
    for eachsource in sourceinfo.keys():
        
        # First, order the datescans and the peak flux lists associated with each by date
        alldates = np.array([datedummy+'_'+scandummy for datedummy,scandummy in zip(np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))],np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))])])
        ordered_peakfluxes = np.array(sourceinfo[eachsource]['peakfluxes'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))]

        # Then, collect the dates and peakfluxlists that meet the RMS threshold 
        datelist_thisRMSlimit = []
        peakfluxlist_thisRMSlimit = []
        peakfluxes_to_normalise_date_cutoff = []
        for i in range(len(ordered_peakfluxes)):
            rms_thisdate = rms_on_date[np.where(metadata_datestring == np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i] + '_' +np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i])][0]
            if rms_thisdate <= RMSthresh:
                datelist_thisRMSlimit.append(alldates[i])
                peakfluxlist_thisRMSlimit.append(ordered_peakfluxes[i])
       
                # The date cutoff is for the normalisation -- so when we take the average, it doesn't change every time we get new data!
                if int(alldates[i].split('_')[0])<=date_cutoff:
                    peakfluxes_to_normalise_date_cutoff.append(ordered_peakfluxes[i])

        # Add results to sourceinfo dictionary
        sourceinfo[eachsource]['avg_peakflux'] = np.nanmean(peakfluxlist_thisRMSlimit)
        sourceinfo[eachsource]['norm_peakfluxes'] = np.array(peakfluxlist_thisRMSlimit) / np.nanmean(peakfluxes_to_normalise_date_cutoff)
        sourceinfo[eachsource]['dates_thisRMSlimit'] = np.array(datelist_thisRMSlimit)


    #print('\nDATES THAT SURVIVED RMS: ',datelist_thisRMSlimit,'\n')

#############
#############
# Pair up sources to find families!
#############
#############

    pairs_completed = []
    # Ensure we dont compare sources to themselves. Source pairs will be labelled 'Index1Index2', so  add e.g. '0202' to list of pairs to ignore
    for eachsource in sourceinfo.keys():
        pairs_completed.append(str(sourceinfo[eachsource]['index']).zfill(2) + str(sourceinfo[eachsource]['index']).zfill(2))

    # n choose k for number of unique pairs [ n!/(k!(n-k)!) ]
    n = len(sourceinfo.keys())
    k = 2
    pairnumber = comb(n, k)

    # Define a results dictionary that will hold pairnames and SD
    # measurements of lightcurve comparison: SD(norm_lightcurve1/norm_lightcurve2)
    results_dict = {}
    pairnames = []
    SDs = []
   
    # order the sourceinfo dictionary by brightness (brightest first)
    allindices = []
    allsources = []
    for eachsource in sourceinfo.keys(): # The long names of the sources: 'JCMTPP*'
        allindices.append(sourceinfo[eachsource]['index']) # An integer uniquely identifying the source in desending order of brightness
        allsources.append(eachsource) # The long names of the sources: 'JCMTPP*'

    sourceinfo_keys_brightness_order = np.array(allsources)[np.argsort(np.array(allindices))]

    # starting with the brightest sources:
    for eachsource in sourceinfo_keys_brightness_order:

        # construct the pairname: Index1Index2A
        # (second_ind+first_ind added to pairs completed at bottom of loop so no double counting)
        first_ind = str(sourceinfo[eachsource]['index']).zfill(2)
        for eachsource2nd in sourceinfo_keys_brightness_order:
            second_ind = str(sourceinfo[eachsource2nd]['index']).zfill(2)
            pairname = first_ind + second_ind 

            # If this pair hasn't been done yet 
            if pairname not in pairs_completed:
     
                # Are both sources above the brightness threshold for this family?
                # Take into consideration that the brightness threshold is in Jy/beam and sourceinfo is in mJy/beam in the case of 850 microns+original 8 maps -- otherwise the conversion does nothing- 2021-01-11
                if sourceinfo[eachsource]['avg_peakflux'] >= brightnessthresh*thresh_conversion_factor and sourceinfo[eachsource2nd]['avg_peakflux'] >= brightnessthresh*thresh_conversion_factor: 
                    # Find standard deviation of their ratio-ed normalised light curves (do they both go up and down at the same time, tracking the inherent JCMT calibration uncertainty?)
                    pairnames.append(pairname)
                    first_divided_by_second = sourceinfo[eachsource]['norm_peakfluxes'] / sourceinfo[eachsource2nd]['norm_peakfluxes']
                    SD = np.nanstd(first_divided_by_second,ddof=1) # Use sample SD, not Population SD (i.e. divide by n-1, ddof=1)
                    SDs.append(SD)
                    results_dict[pairname] = SD
                
                # place both orders of pairs in the "DO NOT REPEAT" list
                pairs_completed.append(pairname)
                pairs_completed.append(second_ind + first_ind)

    # Save SDs and pairnames as an array and sort them by SD! (lowest first)
    SDs = np.array(SDs)
    pairnames = np.array(pairnames)
    sorted_SDs = SDs[np.argsort(SDs)]
    sorted_pairnames = pairnames[np.argsort(SDs)]


    # Get a unique set of the individual individual source IDs within the family -- BUT LEAVE OUT THE VARIABLES!
    unsorted_sources = []
    for eachpair in pairnames:
        source1 = eachpair[0:2]
        source2 = eachpair[2:4]
        if source1 not in unsorted_sources:
            if source1 not in known_variables:
                unsorted_sources.append(source1)
        if source2 not in unsorted_sources:
            if source2 not in known_variables:
                unsorted_sources.append(source2)

    sources = sorted(unsorted_sources)
    source_set = set(sources)

##############
##############
# now gather information to give back to main program in order to plot FCF properties
##############
##############

    numcals_for_plot        = []  # The number of family members
    SD_of_best_fam_for_plot = []  # The lowest SD threshold that allows for this size of family
    fam_for_plot            = []  # The list of sources in the family
    FCF_dates_all           = []  # After removing RMS's that won't give robust data - these are the dates we are left with for each family size
    families_all            = []  # The families we derive for each family size
    FCFs_all                = []  # The FCFs we derive associated with the dates for each family size
    FCF_uncs_all            = []  # The FCF uncertainties we derive associated with the dates for each family size
    normfluxes_by_date_all  = []  # Normalised family fluxes by date for each family size
    med_FCF_unc_for_plot    = []  # Median FCF uncertainty for each family size
    med_FCF_cat_for_plot    = []  # Median FCF uncertainty categorical

    # Starting with a family size of 2 up to the maximum number of sources above the brightness threshold (that are not variables!) derive the above information.
    for eachfamilysize in np.arange(2, len(source_set) + 0.1, 1):

        eachfamilysize = int(eachfamilysize)
        numcals_for_plot.append(eachfamilysize)
        # The following line figures out how many combinations we can have between our sources in this family size
        subsets_of_this_size = list(itertools.combinations(source_set, eachfamilysize))

        # Find the highest SD in this combination of sources -- that's the SD threshold for all these sources to be in the same family
        highest_SD_in_each_combo = []
        list_of_combos = [] # We will get the optimal family using this parameter
        dummy = 0
        for eachcombo in subsets_of_this_size:
            dummy += 1
            SDs_in_this_combo_of_sources = []
            for eachpair, eachSD in zip(sorted_pairnames, sorted_SDs):
                if (eachpair[0:2] in eachcombo) and (eachpair[2:4] in eachcombo):
                    SDs_in_this_combo_of_sources.append(eachSD) # we have already derived the SDs between sources
            list_of_combos.append(eachcombo)
            highest_SD_in_each_combo.append(max(SDs_in_this_combo_of_sources))

        # Find the Lowest SD we can use for this number of sources and select that combination as the optimal family of this size 
        SD_of_best_fam_for_plot.append(min(highest_SD_in_each_combo))  # The lowest thres we could use to construct a fam of this number of sources
        fam_for_plot.append(np.array(list_of_combos)[np.argmin(np.array(highest_SD_in_each_combo))])

        # Now calculate median FCF unc using sourceinfo
        norm_flux_lists_RMS_corrected = []
        date_list_RMS_corrected = []
        # For each source in the family that has the lowest SD threshold
        for eachsource in np.array(list_of_combos)[np.argmin(np.array(highest_SD_in_each_combo))]:
            # loop through the sourceinfo catalogue and find the catalogue source name
            for eachcatsource in sourceinfo.keys():
                if sourceinfo[eachcatsource]['index'] == int(eachsource):
                    source_key = eachcatsource
            norm_flux_lists_RMS_corrected.append(sourceinfo[source_key]['norm_peakfluxes'])
            date_list_RMS_corrected.append(sourceinfo[source_key]['dates_thisRMSlimit'])

        # Extract the normalised fluxes by date
        # And the average = FCF
        # And the SD = FCF unc -- Reduce the FCF uncertainty by sqrt(N)
        normfluxes_by_date = []
        FCFs = []
        FCF_uncs = []
        for eachnormflux in np.array(norm_flux_lists_RMS_corrected).T:
            FCFs.append(np.average(eachnormflux))
            FCF_uncs.append(100 * np.std(eachnormflux, ddof=1) / np.sqrt(len(eachnormflux)))
            normfluxes_by_date.append(eachnormflux)

        # Compile all the information for this family size
        FCFs_all.append(FCFs)
        families_all.append(fam_for_plot)
        FCF_uncs_all.append(FCF_uncs)
        FCF_dates_all.append(date_list_RMS_corrected[0])
        normfluxes_by_date_all.append(np.array(normfluxes_by_date))

        med_FCF_unc_for_plot.append(np.median(FCF_uncs))
        if np.median(FCF_uncs)<1.0:
            med_FCF_cat_for_plot.append('<1.0\%')
        elif np.median(FCF_uncs)<2.0:
            med_FCF_cat_for_plot.append('1-2\%')
        elif np.median(FCF_uncs)<3.0:
            med_FCF_cat_for_plot.append('2-3\%')
        elif np.median(FCF_uncs)<5.0:
            med_FCF_cat_for_plot.append('3-5\%')
        elif np.median(FCF_uncs)<7.0:
            med_FCF_cat_for_plot.append('5-7\%')
        else:
            med_FCF_cat_for_plot.append('>7\%')

    family_results = {}
    family_results['SD_of_best_fam_for_plot'] = SD_of_best_fam_for_plot
    family_results['numcals_for_plot']        = numcals_for_plot
    family_results['med_FCF_cat_for_plot']    = med_FCF_cat_for_plot
    family_results['RMSthresh']               = RMSthresh
    family_results['FCF_dates_all']           = FCF_dates_all
    family_results['FCFs_all']                = FCFs_all
    family_results['FCF_uncs_all']            = FCF_uncs_all
    family_results['normfluxes_by_date_all']  = normfluxes_by_date_all
    family_results['families_all']            = families_all

    with open('config/Families_info_{}_{}_{}_{}.bin'.format(eachregion,wave,int(float(eachtargunc)),date_cutoff),'wb') as outfile:
        pickle.dump(family_results,outfile)

    return(SD_of_best_fam_for_plot, numcals_for_plot,med_FCF_cat_for_plot,RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,normfluxes_by_date_all,families_all)

