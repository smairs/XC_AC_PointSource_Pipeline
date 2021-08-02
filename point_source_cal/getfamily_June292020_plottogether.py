import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import pickle
from scipy.special import comb
import glob
import itertools
from point_source_cal.noisefunctions import readnoise
from point_source_cal.get_vars import get_vars
from astropy.io import fits
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns

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


def plot_SDfamsize(eachregion,wave,eachtargunc,date_cutoff,jsonfile='point_source_cal/bright_threshes.json'):


    brightness_threshes = json.load(jsonfile)

    if np.logical_and(wave == '850',region_to_run in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']):  #-- Jy/beam! This is to match the coadd catalogue
        thresh_conversion_factor = 1000.0 # To multiply catalogue numbers by so it matches the sourceinfo file (Doug's catalogues were in Jy/beam)
    else: #-- mJy/beam! This is to match the coadd catalogues! 
        thresh_conversion_factor = 1.0 # To multiply catalogue numbers by so it matches the sourceinfo file -- no correction needed in this case (Steve's cats were in mJy/beam)

    region_noises = {}
    date_scans,noises = readnoise(wave+'_noises.txt',eachregion) # This is always in mJy/beam for both 450 and 850
    region_noises[eachregion] = {}
    region_noises[eachregion]['noises'] = np.array(noises)
    region_noises[eachregion]['date_scans'] = np.array(date_scans)

    plotcolorbounds = [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 100.0]
    plotcolors      = sns.color_palette('colorblind',len(plotcolorbounds)-3)
    #plotcolors = list(mcolors.TABLEAU_COLORS)[0:len(plotcolorbounds)-3]
    plotcolors.insert(0,'k')
    plotcolors.append('0.75')

    SD_threshes_thistarg = []
    Numcals_thistarg = []
    Colors_for_plot_thistarg = []
    perc_maps_thistarg = []
    rmsthresh_thistarg = []
    targ_record_thistarg = []
    # Set up plot parameters
    fig, ax1 = plt.subplots()
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.05)
    # cax = fig.add_axes([0.9,0.27,0.05,0.5])
    cmap = mpl.colors.ListedColormap(plotcolors[1:-1])
    cmap.set_over('0.75')
    cmap.set_under('k')
    bounds = plotcolorbounds[1:-1]
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, extend='both', ticks=bounds, spacing='proportional',
                                   orientation='vertical', label='Median FCF Unc (\%)')
    props = dict(boxstyle='round', facecolor='#e9f8ff', alpha=1.0)
    #

    target_perc_unc = float(eachtargunc)/100.0
    brightnessthresh = brightness_threshes[wave][eachregion][eachtargunc]

    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('\nREGION = ' + eachregion + ', Flux Thresh: ', brightnessthresh, ', Targ Uncertainty: ', eachtargunc)
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    # coadd_cat is what we used to plot the ensemble signal plots

    if eachregion not in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']:
        coadd_cat_name = 'config/'+eachregion+'_'+wave+'_sourcecat_20201201.fits'
    elif wave == '450':
        coadd_cat_name = 'config/'+eachregion+'_'+wave+'_sourcecat_20200911.fits'
    elif wave == '850':
        coadd_cat_name = 'config/'+eachregion+'_'+wave+'_sourcecat_20170616.fits'

    coadd_cat = fits.getdata(coadd_cat_name)

    sourceinfo = pickle.load(open('pointsource_results/' + eachregion + '/' + eachregion + '_'+wave+'_sourcecat.bin', 'rb'))
    metadata = np.genfromtxt(
        sorted(list(glob.glob('pointsource_results/' + eachregion + '/*'+wave+'*_Wcal_metadata.txt')))[-1], names=True,
        dtype=[('ID', '<f8'), ('Name', '<U28'), ('UT', '<f8'), ('JD', '<f8'), ('Obs', '<f8'), ('Elev', '<f8'),
               ('Tau225', '<f8'), ('RMS', '<f8'), ('RMS_unit', '<f8')])


    num_sources_above_flux = len(
        np.array(coadd_cat['PEAK_FLUX'])[np.where(np.array(coadd_cat['PEAK_FLUX']) > brightnessthresh)]) # This is already in mJy/bm for 450 and Jy/bm for 850

    #Find the lowest brightness source that is above the threshold -- multiply by 1000 if one of 8 original regions and 850 microns since coadd_cat is in Jy/beam, but maps are in mJy/beam - 2021-01-11
    lowest_flux_above_thresh = thresh_conversion_factor*min(np.array(coadd_cat['PEAK_FLUX'])[
                                       np.where(np.array(coadd_cat['PEAK_FLUX']) > brightnessthresh)])

    RMSthresh = target_perc_unc * lowest_flux_above_thresh * np.sqrt(num_sources_above_flux) # For both wavelengths this will be in mJy/beam now
    perc_maps = 100 * len(region_noises[
                              eachregion]['noises'][
                              np.where(region_noises[eachregion]['noises'] <= RMSthresh)]) / float(
        len(region_noises[eachregion]['noises']))

    var_list,var_names = get_vars(eachregion,wave)
    known_variables = []
    for eachknownvariable in var_list:
        known_variables.append(str(eachknownvariable).zfill(2))
    print('\n\nKnown Variables:',known_variables,'\n\n')

    metadata_datestring = []
    rms_on_date = []
    tau225s = []
    for i in range(len(metadata['Name'])):
        metadata_datestring.append(metadata['Name'][i].split('_')[1] + '_' + metadata['Name'][i].split('_')[2])
        rms_on_date.append(metadata['RMS'][i])
        tau225s.append(metadata['Tau225'][i])

    metadata_datestring = np.array(metadata_datestring)
    rms_on_date = np.array(rms_on_date)
    tau225s = np.array(tau225s)

    totalnightnum = len(np.array(sourceinfo[list(sourceinfo.keys())[0]]['dates']))

    for eachsource in sourceinfo.keys():
        alldates2 = np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))]
        allscans2 = np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))]
        alldates1 = []
        for datedummy,scandummy in zip(alldates2,allscans2):
            alldates1.append(datedummy+'_'+scandummy)
        alldates1 = np.array(alldates1)
        ordered_peakfluxes1 = np.array(sourceinfo[eachsource]['peakfluxes'])[
            np.argsort(np.array(sourceinfo[eachsource]['dates']))]
        datelist_thistaulimit = []
        peakfluxlist_thistaulimit = []
        peakfluxes_to_normalise_date_cutoff = []
        for i in range(len(ordered_peakfluxes1)):
            #print(np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i] + '_' +np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i],metadata_datestring)
            print(metadata_datestring,np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i] + '_' +np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i])
            rms_thisdate = rms_on_date[np.where(metadata_datestring == np.array(sourceinfo[eachsource]['dates_reg'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i] + '_' +np.array(sourceinfo[eachsource]['scan'])[np.argsort(np.array(sourceinfo[eachsource]['dates']))][i])][0]
            if rms_thisdate <= RMSthresh:
                #tau225 = tau225s[np.where(metadata_datestring == np.array(sourceinfo[eachsource]['dates_reg'])[
                #    np.argsort(np.array(sourceinfo[eachsource]['dates']))][i] + '_' +
                #                          np.array(sourceinfo[eachsource]['scan'])[
                #                              np.argsort(np.array(sourceinfo[eachsource]['dates']))][i])][0]
                datelist_thistaulimit.append(alldates1[i])
                peakfluxlist_thistaulimit.append(ordered_peakfluxes1[i])
                #print(int(alldates1[i].split('_')[0]))
                if int(alldates1[i].split('_')[0])<=date_cutoff:
                    peakfluxes_to_normalise_date_cutoff.append(ordered_peakfluxes1[i])
        sourceinfo[eachsource]['avg_peakflux'] = np.mean(peakfluxlist_thistaulimit)
        sourceinfo[eachsource]['norm_peakfluxes'] = np.array(peakfluxlist_thistaulimit) / np.mean(
            peakfluxes_to_normalise_date_cutoff)
        sourceinfo[eachsource]['dates_thisRMSlimit'] = np.array(datelist_thistaulimit)


    print('DATES THAT SURVIVED RMS: ',datelist_thistaulimit)
    pairs_completed = []

    # Ensure we dont compare sources to themselves
    for eachsource in sourceinfo.keys():
        pairs_completed.append(
            str(sourceinfo[eachsource]['index']).zfill(2) + str(sourceinfo[eachsource]['index']).zfill(2))

    n = len(sourceinfo.keys())
    k = 2

    # n choose k for number of unique pairs
    pairnumber = comb(n, k)

    results_dict = {}

    pairnames = []
    SDs = []

    allindices = []
    allsources = []
    for eachsource in sourceinfo.keys():
        allindices.append(sourceinfo[eachsource]['index'])
        allsources.append(eachsource)

    sourceinfo_keys_brightness_order = np.array(allsources)[np.argsort(np.array(allindices))]

    for eachsource in sourceinfo_keys_brightness_order:
        first_ind = str(sourceinfo[eachsource]['index']).zfill(2)
        for eachsource2nd in sourceinfo_keys_brightness_order:
            second_ind = str(sourceinfo[eachsource2nd]['index']).zfill(2)
            pairname = first_ind + second_ind
            if pairname not in pairs_completed:
                if sourceinfo[eachsource]['avg_peakflux'] >= brightnessthresh*thresh_conversion_factor and sourceinfo[eachsource2nd][
                    'avg_peakflux'] >= brightnessthresh*thresh_conversion_factor: #take into consideration that the brightness threshold is in Jy/beam and sourceinfo is in mJy/beam in the case of 850 microns+original 8 maps -- otherwise the conversion does nothing- 2021-01-11
                    pairnames.append(pairname)
                    first_divided_by_second = sourceinfo[eachsource]['norm_peakfluxes'] / sourceinfo[eachsource2nd][
                        'norm_peakfluxes']
                    SD_num = []
                    for i in first_divided_by_second:
                        SD_num.append((i - np.mean(first_divided_by_second)) ** 2.0)
                    SD = np.sqrt(sum(SD_num) / (len(SD_num) - 1))
                    SDs.append(SD)

                    results_dict[pairname] = SD
                # place both orders of pairs in the "DO NOT REPEAT" list
                pairs_completed.append(pairname)
                pairs_completed.append(second_ind + first_ind)

    SDs = np.array(SDs)
    pairnames = np.array(pairnames)

    # Get the individual source names -- BUT LEAVE OUT THE VARIABLES!
    unsorted_sources = []
    for eachpair in pairnames:
        source1 = eachpair[0:2]
        source2 = eachpair[2:4]
        #print(source1,source2)
        if source1 not in unsorted_sources:
            if source1 not in known_variables:
                unsorted_sources.append(source1)
        if source2 not in unsorted_sources:
            if source2 not in known_variables:
                unsorted_sources.append(source2)
    sources = sorted(unsorted_sources)

    source_set = set(sources)


    sorted_SDs = SDs[np.argsort(SDs)]
    sorted_pairnames = pairnames[np.argsort(SDs)]


    numcals_for_plot = []
    SD_of_best_fam_for_plot = []
    fam_for_plot = []
    FCF_dates_all = []
    families_all = []
    FCFs_all = []
    FCF_uncs_all = []
    normfluxes_by_date_all = []
    med_FCF_unc_for_plot = []
    colors_for_plot = []
    print(sorted_pairnames,sorted_SDs)
    for eachpotentialcalnum in np.arange(2, len(source_set) + 0.1, 1):
        eachpotentialcalnum = int(eachpotentialcalnum)
        numcals_for_plot.append(eachpotentialcalnum)

        subsets_of_this_size = list(itertools.combinations(source_set, eachpotentialcalnum))

        highest_SD_in_each_combo = []
        list_of_combos = []
        dummy = 0
        for eachcombo in subsets_of_this_size:
            dummy += 1
            print('Working on Source Combo ', dummy, ' of ', len(subsets_of_this_size), ' -- Family Size: ',
                  eachpotentialcalnum, ' -- total cal num ', len(source_set))
            SDs_in_this_combo_of_sources = []
            for eachpair, eachSD in zip(sorted_pairnames, sorted_SDs):
                if (eachpair[0:2] in eachcombo) and (eachpair[2:4] in eachcombo):
                    SDs_in_this_combo_of_sources.append(eachSD)
            list_of_combos.append(eachcombo)
            highest_SD_in_each_combo.append(max(SDs_in_this_combo_of_sources))
        SD_of_best_fam_for_plot.append(min(
            highest_SD_in_each_combo))  # The lowest thres we could use to construct a fam of this number of sources
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

        normfluxes_by_date = []
        FCFs = []
        FCF_uncs = []
        for eachnormflux in np.array(norm_flux_lists_RMS_corrected).T:
            FCFs.append(np.average(eachnormflux))
            FCF_uncs.append(100 * np.std(eachnormflux, ddof=1) / np.sqrt(len(eachnormflux)))
            normfluxes_by_date.append(eachnormflux)

        FCFs_all.append(FCFs)
        families_all.append(fam_for_plot)
        FCF_uncs_all.append(FCF_uncs)
        FCF_dates_all.append(date_list_RMS_corrected[0])
        normfluxes_by_date_all.append(np.array(normfluxes_by_date))

        med_FCF_unc_for_plot.append(np.median(FCF_uncs))
        for eachplotcolorbound in range(len(plotcolorbounds) - 1):
            if np.logical_and(np.median(FCF_uncs) >= plotcolorbounds[eachplotcolorbound],
                              np.median(FCF_uncs) < plotcolorbounds[eachplotcolorbound + 1]):
                colors_for_plot.append(plotcolors[eachplotcolorbound])
            if np.median(FCF_uncs) > 14:
                colors_for_plot.append('0.75')

    return(SD_of_best_fam_for_plot, numcals_for_plot, colors_for_plot,RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,
           normfluxes_by_date_all,families_all)

