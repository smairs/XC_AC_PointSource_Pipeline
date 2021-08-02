# Import necessary Packages

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import numpy as np
from astropy.io import fits
from point_source_cal.noisefunctions import readnoise
import glob
from point_source_cal.get_vars import get_vars
import seaborn as sns
from point_source_cal.getfamily_June292020_plottogether import plot_SDfamsize
from point_source_cal.get_bright_threshes import get_bright_threshes
from mpl_toolkits.axes_grid1 import make_axes_locatable
from astropy.time import Time
import pickle
import os

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

#target_uncertainties = [0.05] # = 5%
def make_FCFunc_family_plots(region_to_run,wave,date_cutoff,target_uncertainties = [0.05],jsonfile='point_source_cal/bright_threshes.json'):
    
    regionstring = ''
    
    region_scaling = 100 #mJy/beam
    
    cat_dir  = 'config/'
    detection_significance_in_single_epoch = 5# sigma
    cat_list = sorted(list(glob.glob(cat_dir+'/'+region_to_run+'*'+wave+'*fits')))
    
    plotcolorbounds = [0,30,40,50,60,70,80,90,100]
    plotcolors      = sns.color_palette('colorblind',len(plotcolorbounds)-3)
    #plotcolors = list(mcolors.TABLEAU_COLORS)[0:len(plotcolorbounds)-3]
    plotcolors.insert(0,'0.75')
    plotcolors.append('k')
    
    plotdummy = -1
    
    brightness_threshes = json.load(jsonfile)

    # Make sure loop below will work - i.e. ensure there is information on brightness thresholds for each target uncertainty listed
    for i,eachtargunc in enumerate(target_uncertainties):
        try:
            brightness_threshes[wave][region_to_run][str(eachtargunc*100.0)]
        except KeyError:
            print('\n\nNO BRIGHTNESS THRESHOLD INFORMATION FOR A {}% target uncertainty!\n\n'.format(round(eachtargunc*100.0,2)))
            print('\n\nREDEFINING TARGET UNCERTAINTY TO 5%...\n\n')
            target_uncertainties[i] = 0.05

    # Eliminate any duplicates in target_uncertainty list
    target_uncertainties = sorted(list(set(target_uncertainties)))


    if np.logical_and(wave == '850',region_to_run in ['IC348','NGC1333','NGC2024','NGC2071','OMC23','OPHCORE','SERPM','SERPS']):  #-- Jy/beam! This is to match the coadd catalogue
        thresh_conversion_factor = 1000.0 # To multiply catalogue numbers by so it matches the sourceinfo file (Doug's catalogues were in Jy/beam)
    else: #-- mJy/beam! This is to match the coadd catalogues! 
        thresh_conversion_factor = 1.0 # To multiply catalogue numbers by so it matches the sourceinfo file -- no correction needed in this case (Steve's cats were in mJy/beam)

    fig,axs = plt.subplots(ncols=2,nrows=1,constrained_layout=False)
    
    for eachtargunc in target_uncertainties:
        target_uncertainty_key = str(eachtargunc*100.0)
        for eachcat in cat_list:
            region         = eachcat.split('/')[-1].split('_')[0]
            print('\n\n',region,'\n\n')
            var_list,var_names = get_vars(region,wave)
            date_scans,noises = readnoise(wave+'_noises.txt',region)
            if region == region_to_run:
                regionstring=regionstring+region+'_'
                # Make an array to indicate which sources are known variables
                # I will be rearranging the peak flux list later, so I want to
                # keep track of where the variables go.
                cat            = fits.getdata(eachcat)
                var_index      = []
                var_name_ind   = []
                nameindex      = -1
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
    
                sorting_index        = np.argsort(cat['PEAK_FLUX'])
                var_index_sorted     = var_index[sorting_index]
                var_names_sorted     = var_name_ind[sorting_index]
                indexes_sorted       = cat['ID'][sorting_index]
                peak_fluxes_sorted   = cat['PEAK_FLUX'][sorting_index]
    
                cal_bright_index = np.where(peak_fluxes_sorted>=brightness_threshes[wave][region][target_uncertainty_key])
    
                cal_IDs_no_vars =  indexes_sorted[cal_bright_index][np.where(var_index_sorted[cal_bright_index]==0)]
    
                brightnesses_to_plot = []
                num_sources_to_plot  = []
                #perc_maps_to_plot = []
                #colors_for_percs = []
                dummy = 0
                for eachflux in peak_fluxes_sorted:
                    dummy = dummy+1
                    brightnesses_to_plot.append(eachflux*thresh_conversion_factor) # Convert from Jy/beam to mJy/beam if necessary - will do nothing if unnecessary
                    num_sources_to_plot.append(dummy)
    
                 # Colouring the plots by source detection capability
                 #   RMS_threshold_for_this_set_of_sources = eachflux/detection_significance_in_single_epoch
                 #   perc_maps = 100*len(noises[np.where(noises<=RMS_threshold_for_this_set_of_sources)])/float(len(noises))
                 #   for eachplotcolorbound in range(len(plotcolorbounds)-1):
                 #       if np.logical_and(perc_maps>=plotcolorbounds[eachplotcolorbound],
                #       perc_maps<plotcolorbounds[eachplotcolorbound+1]):
                 #           colors_for_percs.append(plotcolors[eachplotcolorbound])
                 #   if perc_maps == 100.0:
                 #       colors_for_percs.append('k')
                 #   perc_maps_to_plot.append(100*len(noises[
                #   np.where(noises<=RMS_threshold_for_this_set_of_sources)])/float(len(noises)))
    
                num_sources_to_plot.reverse()
                brightnesses_to_plot = np.array(brightnesses_to_plot)
                num_sources_to_plot  = np.array(num_sources_to_plot)
                #perc_maps_to_plot    = np.array(perc_maps_to_plot)
                #colors_for_percs     = np.array(colors_for_percs)
    
                perc_maps_to_plot    = []
                colors_for_percs     = []
                print('Num sources, flux cutoff, RMS cutoff')
                for eachflux,eachnum in zip(brightnesses_to_plot,num_sources_to_plot):
                    if region in ['OMC23','SERPS']:
                        print(eachnum,eachflux,eachtargunc*eachflux*np.sqrt(eachnum))
                    # Ensemble Signal x Noise/Ensemble Signal
                    noise_required_for_this_sig_unc = eachtargunc*eachflux*np.sqrt(eachnum) # Noise should be in mJy/beam - it matches the "noises" file -- already converted, above
                    perc_maps = 100*len(noises[np.where(noises<=noise_required_for_this_sig_unc)])/float(len(noises))
                    for eachplotcolorbound in range(len(plotcolorbounds)-1):
                        if np.logical_and(perc_maps>=plotcolorbounds[eachplotcolorbound],
                                          perc_maps<plotcolorbounds[eachplotcolorbound+1]):
                            colors_for_percs.append(plotcolors[eachplotcolorbound])
                    if perc_maps == 100.0:
                        colors_for_percs.append('k')
                    perc_maps_to_plot.append(perc_maps)
    
                perc_maps_to_plot = np.array(perc_maps_to_plot)
                colors_for_percs  = np.array(colors_for_percs)
    
                # Design Custom Colourmap
                customcmap = mpl.colors.ListedColormap(plotcolors[1:-1])
                customcmap.set_over('k')
                customcmap.set_under('0.75')
                cmap_bounds = plotcolorbounds[1:-1]
                norm = mpl.colors.BoundaryNorm(cmap_bounds,customcmap.N)
                #
    
                #
                plotdummy=plotdummy+1
                #ax1 = axs[plotdummy,0]
                ax1 = axs[plotdummy]
                #divider = make_axes_locatable(ax1)
                #cax1 = divider.append_axes('right', size='5%', pad=0.05)
                #cb1 = mpl.colorbar.ColorbarBase(cax1, cmap=customcmap, norm=norm, extend='both', ticks=cmap_bounds,
                #                                spacing='proportional',
                #                                orientation='vertical', label='\% Maps w/ RMS that gives ' +
                #                                                              str(100 * eachtargunc) + '\% Unc')
                ax1.scatter(brightnesses_to_plot,100.0/(
                        np.sqrt(num_sources_to_plot)*brightnesses_to_plot/region_scaling),c=colors_for_percs,s=20)
                if len(np.where(var_index_sorted==1)[0])>0:
                    ax1.scatter(brightnesses_to_plot[np.where(var_index_sorted==1)],
                                100.0/(np.sqrt(num_sources_to_plot[np.where(var_index_sorted==1)])*brightnesses_to_plot[
                                    np.where(var_index_sorted==1)]/region_scaling),c=colors_for_percs[
                            np.where(var_index_sorted==1)],marker='*',s=120,label='Variable')
                    for eachvar in range(len(brightnesses_to_plot[np.where(var_index_sorted==1)])):
                        ax1.text(brightnesses_to_plot[np.where(var_index_sorted==1)][eachvar],
                                 100.0/(np.sqrt(num_sources_to_plot[np.where(var_index_sorted==1)][eachvar])*brightnesses_to_plot[
                                     np.where(var_index_sorted==1)][eachvar]/region_scaling)+0.5,var_names_sorted[
                                     np.where(var_index_sorted==1)][eachvar],rotation=90)
                ax1.set_xlim(xmin=200,xmax=3*1e5)
                ax1.set_ylim(ymin=0,ymax=9)
                ax1.grid(which='minor',alpha=0.2)
                ax1.grid(which='major',alpha=0.5)
                ax1.set_xlabel('Peak Flux (mJy beam$^{-1}$)')
                ax1.set_ylabel('Ensemble Sig Unc For '+str(region_scaling)+' mJy RMS (\%)')
                ax1.set_xscale('log')
    
                #fig.colorbar(mpl.cm.ScalarMappable(
                #    norm=norm,cmap=customcmap),extend='both',orientation='vertical',ticks=cmap_bounds,spacing='proportional',
                #    label='\% Maps w/ RMS that gives '+str(100*eachtargunc)+'\% Unc')
    
                # Draw Vertical Threshold at cutoff flux for potential calibrators
                ax1.axvline(x=brightness_threshes[wave][region][target_uncertainty_key],linestyle='dashed',color='k',linewidth=2,
                            label='Potential Cal Thresh')
    
                props = dict(boxstyle='round', facecolor='#e9f8ff', alpha=1.0)
                if region in ['NGC2024','OMC23','OPHCORE','SERPS']:
                    x_placement=0.25
                elif np.logical_and(region in ['NGC1333'],eachtargunc<0.05):
                    x_placement=0.25
                else:
                    x_placement = 0.75
    #            ax1.text(x_placement, 0.8, region+'\nTarget Uncertainty = '+target_uncertainty_key+'\%',
    #                     horizontalalignment='center',
    #                     verticalalignment='top', transform=ax1.transAxes, fontsize=4, bbox=props)
                ax1.text(x_placement, 0.8, region,
                         horizontalalignment='center',
                         verticalalignment='top', transform=ax1.transAxes, fontsize=10, bbox=props)
    
                if region=='SERPM':
                    ax1.legend(loc='upper right')
                elif region=='SERPS':
                    ax1.legend(loc='upper left')
                else:
                    ax1.legend(loc='upper left')
    
                # SD thresh plot
                plotcolorbounds = [0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 100.0]
                plotcolors = sns.color_palette('colorblind', len(plotcolorbounds) - 3)
                #plotcolors = list(mcolors.TABLEAU_COLORS)[0:len(plotcolorbounds)-3]
                plotcolors.insert(0, 'k')
                plotcolors.append('0.75')
    
                SD_of_best_fam_for_plot,numcals_for_plot,colors_for_plot,\
                    RMSthresh,FCF_dates_all,FCFs_all,FCF_uncs_all,\
                    normfluxes_by_date_all,families_all= plot_SDfamsize(region,wave,str(100*eachtargunc),date_cutoff)
    
                tick_index = np.round(np.linspace(0, len(brightnesses_to_plot) - 1, 3)).astype(int)
                source_num_ticks  = brightnesses_to_plot[tick_index]
                source_num_labels = []
                for i in num_sources_to_plot[tick_index]:
                    source_num_labels.append(str(i))
                ax2 = ax1.twiny()
                ax2.set_xlim(ax1.get_xlim())
                ax2.set_xscale('log')
                ax2.set_xticks(source_num_ticks)
                ax2.set_xticklabels(source_num_labels,fontsize=7)
                ax2.minorticks_off()
                ax2.set_xlabel('Number of Sources $\ge$ Peak Flux',fontsize=8)
                #divider = make_axes_locatable(ax2)
                #cax2 = divider.append_axes('right', size='5%', pad=0.1)
                #cb2 = mpl.colorbar.ColorbarBase(cax2, cmap=customcmap, norm=norm, extend='both', ticks=cmap_bounds,
                #                                spacing='proportional',
                #                                orientation='vertical', label='\% Maps w/ RMS that gives ' +
                #                                                              str(100 * eachtargunc) + '\% Unc')
    
                #ax3 = axs[plotdummy,1]
                ax3 = axs[plotdummy+1]
                divider = make_axes_locatable(ax3)
                cax = divider.append_axes('right', size='5%', pad=0.05)
                cmap = mpl.colors.ListedColormap(plotcolors[1:-1])
                cmap.set_over('0.75')
                cmap.set_under('k')
                bounds = plotcolorbounds[1:-1]
                norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
                cb = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, extend='both', ticks=bounds,
                                               spacing='proportional',
                                               orientation='vertical', label='Median FCF Unc (\%)')
    
                props = dict(boxstyle='round', facecolor='#e9f8ff', alpha=1.0)
    
                ax3.scatter(SD_of_best_fam_for_plot, numcals_for_plot, c=colors_for_plot,
                            cmap=cmap, norm=norm,
                            marker='o')
    #            ax3.text(0.30, 0.95,
    #                     region + '\nTarget Uncertainty = ' + str(100*eachtargunc) + '\%\n' + str(
    #                         round(perc_maps)) + '\% of maps w/ largest family\nRMS $<$ ' + str(
    #                         round(RMSthresh, 1)) + ' mJy/beam',
    #                     horizontalalignment='center',
    #                     verticalalignment='top', transform=ax3.transAxes, fontsize=4, bbox=props)
                ax3.text(0.30, 0.95,
                          region,horizontalalignment='center',
                         verticalalignment='top', transform=ax3.transAxes, fontsize=10, bbox=props)
    
                ax3.set_xlabel('Lowest SD Thresh For Family Size')
                ax3.set_ylabel('Family Size')
                ax3.set_ylim(ymin=1, ymax=14)
                ax3.set_xlim(xmin=0.02, xmax=0.35)
    
                fig.tight_layout(pad=0.45)
                fig.savefig('pointsource_results/'+region+'SignalUnc_versus_FaintSNR_'+wave+'_targunc'+str(int(float(target_uncertainty_key)))+
                            '.pdf',format='pdf')
                plt.clf()
                plt.close()
                plt.close(fig)
    
                # The negative 1 index gives the highest cal number that meets the target uncertainty criterion
                cal_info_dict = {}
                print(families_all)
                if len(families_all)>0:
                    if len(families_all[-1])>0:
                        cal_info_dict['family'] = families_all[-1][-1]
                    else:
                        print(families_all)
                    cal_info_dict['datescans'] = FCF_dates_all[-1]
                    cal_info_dict['RelFCFs'] = FCFs_all[-1]
                    cal_info_dict['RelFCFuncs'] = np.array(FCF_uncs_all[-1]) / 100
    
                    pickle.dump(cal_info_dict, open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc'+
                                                str(int(float(target_uncertainty_key)))+'_'+wave+'.pickle', 'wb'))
        #fig.tight_layout(pad=0.45)
        #fig.savefig(regionstring + 'SignalUnc_versus_FaintSNR_targunc' + str(int(float(target_uncertainty_key))) +
        #            '.pdf', format='pdf')
        #plt.clf()
        #plt.close()
        #plt.close(fig)
    if len(families_all)>0:
        # Now, plot calfactors for each family number!
        for eachcalnum in range(len(FCFs_all)):
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
            plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCF_with_time_'+wave+'.png',format='png',dpi=300)
            plt.clf()
        
        for eachcalnum in range(len(FCFs_all)):
            plt.hist(FCFs_all[eachcalnum],color='k')
            plt.suptitle(region+', Num Fam Members: '+str(numcals_for_plot[eachcalnum]))
            plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFhists_'+wave+'.png',format='png',dpi=300)
            plt.clf()
        
        for eachcalnum in range(len(FCFs_all)):
            plt.hist(FCF_uncs_all[eachcalnum],color='k')
            plt.xlabel('Rel FCF Group Uncertainty (\%)')
            plt.ylabel('Epoch Count')
            plt.xlim(xmin=0,xmax=10)
            plt.ylim(ymin=0, ymax=12)
            plt.suptitle(region+', Num Fam Members: '+str(numcals_for_plot[eachcalnum]))
            plt.savefig('pointsource_results/FCFplots/'+region+'_'+str(numcals_for_plot[eachcalnum])+'FamMems_FCFunchists_'+wave+'.png',format='png',dpi=300)
            plt.clf()
