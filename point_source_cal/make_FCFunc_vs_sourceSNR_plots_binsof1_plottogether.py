# Import necessary Packages

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from point_source_cal.noisefunctions import readnoise
import glob
from point_source_cal.get_vars import get_vars
import seaborn as sns
from point_source_cal.getfamily_June292020_plottogether import plot_SDfamsize
from point_source_cal.get_previously_defined_family import get_previously_defined_family
from astropy.time import Time
import pickle
import os
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

def make_FCFunc_family_plots(region_to_run,wave,date_cutoff,find_new_family=False,target_uncertainties = [0.05],jsonfile='point_source_cal/bright_threshes.json',cat_dir = 'config/',goodmaps=False,gooddatescans=[]):
    '''
    This code plots the expected calibration uncertainty for each source (noise/[Flux*sqrt(N)]) 
    where N is the number of sources above that brightness threshold
    '''   

    # To scale the  Ensemble uncertainty plots 
    # (just a scaling factor to see the shape of the curve and select a proper brightness threshold)
    region_scaling = 100 #mJy/beam
    
    # Grab the source catalogue derived from the co-add for this region and wavelength
    if region_to_run not in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']:
        sourcecat = cat_dir+region_to_run+'_'+wave+'_sourcecat_20220303.fits'
    elif wave == '450':
        sourcecat = cat_dir+region_to_run+'_'+wave+'_sourcecat_20200911.fits'
    elif wave == '850':
        sourcecat = cat_dir+region_to_run+'_'+wave+'_sourcecat_20170616.fits'

    
    # Get the brightness thresholds for all regions and wavelengths organised by target uncertainty (previously compiled from manual testing) 
    with open(jsonfile) as json_file:
        brightness_threshes = json.load(json_file)

    # Make sure loop below will work - i.e. ensure there is information on brightness thresholds for each target uncertainty listed
    for i,eachtargunc in enumerate(target_uncertainties):
        try:
            brightness_threshes[wave][region_to_run][str(eachtargunc*100.0)]
        except KeyError:
            print('\tNO BRIGHTNESS THRESHOLD INFORMATION FOR A {}% target uncertainty!\n\n'.format(round(eachtargunc*100.0,2)))
            print('\n\tREDEFINING TARGET UNCERTAINTY TO 5%...')
            target_uncertainties[i] = 0.05

    # Eliminate any duplicates in target_uncertainty list
    target_uncertainties = sorted(list(set(target_uncertainties)))


    if np.logical_and(wave == '850',region_to_run in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']):  #-- Jy/beam! This is to match the coadd catalogue
        thresh_conversion_factor = 1000.0 # To multiply catalogue numbers by so it matches the sourceinfo file (Doug's catalogues were in Jy/beam)
    else: #-- mJy/beam! This is to match the coadd catalogues! 
        thresh_conversion_factor = 1.0 # To multiply catalogue numbers by so it matches the sourceinfo file -- no correction needed in this case (Steve's cats were in mJy/beam)

    
############
############
# Start loop
############
############

    # For each target uncertainty we want to get an ensemble signal uncertainty plot and a family!

    for eachtargunc in target_uncertainties:
        target_uncertainty_key = str(eachtargunc*100.0)
        region         = sourcecat.split('/')[-1].split('_')[0]

        #print('\n\n',region,'\n\n')

        # Read in the RMS noise of each observation of this region at this wavelength
        date_scans,noises = readnoise(wave+'_noises.txt',region)

        #####################
        # Make an index to indicate which sources are known variables
        # I will be rearranging the peak flux list later, so I want to
        # keep track of where the variables go.
        var_list,var_names = get_vars(region,wave)
        cat                = fits.getdata(sourcecat)
        var_index          = []
        var_name_ind       = []
        nameindex          = -1
        for i in range(len(cat['PEAK_FLUX'])):
            if i in var_list:
                nameindex += 1
                var_index.append(1)
                var_name_ind.append(var_names[nameindex])
            else:
                var_index.append(0)
                var_name_ind.append('NONE')
    
        var_index    = np.array(var_index)
        var_name_ind = np.array(var_name_ind)
        ###################

        # Sort everything by peak flux     
        sorting_index        = np.argsort(cat['PEAK_FLUX'])
        var_index_sorted     = var_index[sorting_index]
        var_names_sorted     = var_name_ind[sorting_index]
        indexes_sorted       = cat['ID'][sorting_index]
        peak_fluxes_sorted   = cat['PEAK_FLUX'][sorting_index]
   
        #######################################
        # Gather the plot information

        # Gather x, y, and hue data

        brightnesses_to_plot = []
        ensemble_uncertainty = []
        num_sources_to_plot  = []
        perc_maps_to_plot    = []
        for eachflux in peak_fluxes_sorted:
            num_sources_above_threshold = len(peak_fluxes_sorted[np.where(peak_fluxes_sorted>=eachflux)])
            # For x and y
            brightnesses_to_plot.append(eachflux*thresh_conversion_factor) # Convert from Jy/beam to mJy/beam if necessary - will do nothing if unnecessary
            ensemble_uncertainty.append(100.0*(region_scaling/(eachflux*np.sqrt(num_sources_above_threshold))))
            # For a proper index we can use for differentiating variable sources from regular sources
            num_sources_to_plot.append(len(peak_fluxes_sorted[np.where(peak_fluxes_sorted>=eachflux)]))
            # For colouring the plots by fraction of maps that could be used to achieve target uncertainty with this set of sources
            RMS_threshold_for_this_set_of_sources = eachtargunc*eachflux*np.sqrt(num_sources_above_threshold)
            # The fraction of plots with RMS noise low enough to use this full set of sources as a family
            perc_maps = 100*len(noises[np.where(noises<=RMS_threshold_for_this_set_of_sources)])/float(len(noises)) 
            if perc_maps <30:
                perc_maps_to_plot.append('<30\%')
            elif perc_maps <60:
                perc_maps_to_plot.append('30-60\%')
            elif perc_maps <90:
                perc_maps_to_plot.append('60-90\%')
            elif perc_maps <=100:
                perc_maps_to_plot.append('90-100\%')


        # Make everything arrays so they are easier to work with
        brightnesses_to_plot = np.array(brightnesses_to_plot)
        ensemble_uncertainty = np.array(ensemble_uncertainty)
        num_sources_to_plot  = np.array(num_sources_to_plot)
        perc_maps_to_plot    = np.array(perc_maps_to_plot)
    
####################
####################
# Begin Plotting
####################
####################

        sns.scatterplot(x=brightnesses_to_plot,y=ensemble_uncertainty,hue=perc_maps_to_plot)

        if len(np.where(var_index_sorted==1)[0])>0:
            sns.scatterplot(x=brightnesses_to_plot[np.where(var_index_sorted==1)],y=ensemble_uncertainty[np.where(var_index_sorted==1)],color='r',label='Variable',marker='*',s=120)
            for eachvar in range(len(brightnesses_to_plot[np.where(var_index_sorted==1)])):
                plt.text(brightnesses_to_plot[np.where(var_index_sorted==1)][eachvar],ensemble_uncertainty[np.where(var_index_sorted==1)][eachvar]+0.5,var_names_sorted[np.where(var_index_sorted==1)][eachvar],rotation=90)
        plt.xlim(xmin=200,xmax=3*1e5)
        plt.ylim(ymin=0,ymax=9)
        plt.grid(which='minor',alpha=0.2)
        plt.grid(which='major',alpha=0.5)
        plt.xlabel('Peak Flux (mJy beam$^{-1}$)')
        plt.ylabel('Ensemble Sig Unc For '+str(region_scaling)+' mJy RMS (\%)')
        plt.xscale('log')
    
        # Draw Vertical Threshold at cutoff flux for potential calibrators
        plt.axvline(x=brightness_threshes[wave][region][target_uncertainty_key],linestyle='dashed',color='k',linewidth=2,label='Cal Thresh for {}% uncertainty'.format(int(100*eachtargunc)))
   
        # Place region name on plot 
        props = dict(boxstyle='round', facecolor='#e9f8ff', alpha=1.0)
        if region in ['NGC2024','OMC23','OPHCORE','SERPS']:
            x_placement=0.25*max(brightnesses_to_plot)
        elif np.logical_and(region in ['NGC1333'],eachtargunc<0.05):
            x_placement=0.25*max(brightnesses_to_plot)
        else:
            x_placement = 0.75*max(brightnesses_to_plot)
        plt.text(x_placement, 0.8*max(ensemble_uncertainty), region,horizontalalignment='center',verticalalignment='top', fontsize=10, bbox=props)
    
        # Draw legend
        if region=='SERPM':
            plt.legend(loc='upper right')
        else:
            plt.legend(loc='upper left')

        plt.savefig('pointsource_results/Estimated_EnsembleUnc_versus_Flux_{}_targunc_{}.png'.format(region,int(100*eachtargunc)),dpi=300) 
        plt.clf()

###########
###########
# Get Families based on this target uncertainty
###########
########### 

        if find_new_family:

             SD_of_best_fam_for_plot,numcals_for_plot,FCF_uncs_cat,\
                 RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,\
                 normfluxes_by_date_all,families_all= plot_SDfamsize(region,cat,wave,str(100*eachtargunc),date_cutoff,goodmaps=goodmaps)

        else:
            try:
                SD_of_best_fam_for_plot,numcals_for_plot,FCF_uncs_cat,\
                     RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,\
                     normfluxes_by_date_all,families_all = get_previously_defined_family(region,wave,str(100*eachtargunc),date_cutoff,gooddatescans=gooddatescans) 
            except:
                print('\n\tNO PREVIOUS FAMILY FILE FOUND WITH PARAMETERS: {}, {} microns, {}% target unc, {} (date cutoff for normalisation)\n\n\t~~~~~~FINDING FAMILY FROM SCRATCH~~~~'.format(region,wave,str(100*eachtargunc),date_cutoff))
                SD_of_best_fam_for_plot,numcals_for_plot,FCF_uncs_cat,\
                    RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,\
                    normfluxes_by_date_all,families_all= plot_SDfamsize(region,cat,wave,str(100*eachtargunc),date_cutoff,goodmaps=goodmaps)


###########
###########
# Make Family Size versus Lowest SD for Family Size plots - color by median FCF uncertainty
###########
###########    
        props = dict(boxstyle='round', facecolor='#e9f8ff', alpha=1.0)
    
        sns.scatterplot(SD_of_best_fam_for_plot, numcals_for_plot, hue=FCF_uncs_cat,marker='o')
        plt.text(0.30*max(SD_of_best_fam_for_plot), 0.95*max(numcals_for_plot),region,horizontalalignment='center',verticalalignment='top', fontsize=10, bbox=props)
    
        plt.xlabel('Lowest SD Thresh For Family Size')
        plt.ylabel('Family Size')
        plt.ylim(ymin=1, ymax=14)
        plt.xlim(xmin=0.02, xmax=0.35)

        if goodmaps:
            plt.savefig('pointsource_results/FamilySize_versus_FamilySD_{}_targunc_{}_GoodMaps.png'.format(region,int(100*eachtargunc)),dpi=300)
        else:    
            plt.savefig('pointsource_results/FamilySize_versus_FamilySD_{}_targunc_{}.png'.format(region,int(100*eachtargunc)),dpi=300)
        plt.clf()
   
###########
###########
# Compile Family information and save it to a pickle file
###########
###########
 
        # The negative 1 index gives the highest cal number that meets the target uncertainty criterion
        cal_info_dict = {}
        if len(families_all)>0:
            if len(families_all[-1])>0:
                cal_info_dict['family'] = families_all[-1][-1]
            cal_info_dict['datescans'] = FCF_dates_all[-1]
            cal_info_dict['RelFCFs'] = FCFs_all[-1]
            cal_info_dict['RelFCFuncs'] = np.array(FCF_uncs_all[-1]) / 100

            if goodmaps:    
                pickle.dump(cal_info_dict, open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc'+str(int(float(target_uncertainty_key)))+'_'+wave+'_GoodMaps.pickle', 'wb'))
            else:
                pickle.dump(cal_info_dict, open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc'+str(int(float(target_uncertainty_key)))+'_'+wave+'.pickle', 'wb'))
###########
###########
# Plot FCFs by date for each size of family
###########
###########
    if len(families_all)>0:
        # Now, plot calfactors for each family number!
        for eachcalnum in range(len(FCFs_all)):
            if not os.path.exists('pointsource_results/FCFplots/'):
                os.system('mkdir pointsource_results/FCFplots/')
            FCF_dates_all_MJD = []
            for datescan in FCF_dates_all[eachcalnum]:
                datestr = datescan.split('_')[0]
                dateisot = datestr[0:4]+'-'+datestr[4:6]+'-'+datestr[6:8]+'T00:00:00.00'
                FCF_dates_all_MJD.append(Time(dateisot,format='isot').mjd)
            FCF_dates_all_MJD = np.array(FCF_dates_all_MJD)

            plt.errorbar(FCF_dates_all_MJD,FCFs_all[eachcalnum],yerr=np.array(FCF_uncs_all[eachcalnum])/100,color='blue')
            plt.suptitle(region+', Num Fam Members: '+str(numcals_for_plot[eachcalnum]))
            date_labels = []
            for i,normflux in enumerate(normfluxes_by_date_all[eachcalnum]):
                for j in normflux:
                    plt.scatter(FCF_dates_all_MJD[i],j,color='k')
                date_labels.append(str(Time(FCF_dates_all_MJD[i],format='jd').isot).split('T')[0])
            plt.ylabel('Rel FCF')
            plt.ylim(ymin=0.2,ymax=1.8)
        
            spaced_xticks = []
            spaced_xtick_labels = []
            for i,date,datelabel in zip(np.arange(0,len(FCF_dates_all_MJD),1),FCF_dates_all_MJD,date_labels):
                if i % 4 == 0:
                    spaced_xticks.append(date)
                    spaced_xtick_labels.append(datelabel)
            plt.xticks(spaced_xticks,spaced_xtick_labels,rotation=20)
            if goodmaps:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCF_with_time_'+wave+'_GoodMaps.png',format='png',dpi=300)
            else:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCF_with_time_'+wave+'.png',format='png',dpi=300)
            plt.clf()

############
############
# Plot FCF and FCF uncertainty hists for each size of family
############
############
        
        for eachcalnum in range(len(FCFs_all)):
            plt.hist(FCFs_all[eachcalnum],color='k')
            plt.suptitle(region+', Num Fam Members: '+str(numcals_for_plot[eachcalnum]))
            if goodmaps:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFhists_'+wave+'_GoodMaps.png',format='png',dpi=300)
            else:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFhists_'+wave+'.png',format='png',dpi=300)
            plt.clf()
        
        for eachcalnum in range(len(FCFs_all)):
            plt.hist(FCF_uncs_all[eachcalnum],color='k')
            plt.xlabel('Rel FCF Group Uncertainty (\%)')
            plt.ylabel('Epoch Count')
            plt.xlim(xmin=0,xmax=10)
            plt.ylim(ymin=0, ymax=12)
            plt.suptitle(region+', Num Fam Members: '+str(numcals_for_plot[eachcalnum]))
            if goodmaps:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFunchists_'+wave+'_GoodMaps.png',format='png',dpi=300)
            else:
                plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFunchists_'+wave+'.png',format='png',dpi=300)
            plt.clf()
