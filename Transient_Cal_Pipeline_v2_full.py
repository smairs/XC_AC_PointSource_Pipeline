from jcmt_transient_alignment.HighMass_data_gen_EAO_newbox import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table
from jcmt_transient_alignment.create_makemap_script import make_pcor, makemap_infiles, create_makemap_script
from jcmt_transient_alignment.applycal import apply_relFCF_AC
from point_source_cal.make_coadds_metadata_tables import make_coadds_metadata
from point_source_cal.noisefunctions import make_noise_file
from point_source_cal.smooth_input_data import smooth_input_data
from point_source_cal.applyPScal import applyPScal
from point_source_cal.applyWeightedCal import applyWeightedCal
from point_source_cal.make_FCFunc_vs_sourceSNR_plots_binsof1_plottogether import make_FCFunc_family_plots
from point_source_cal.get_vars import get_vars
from weightedcal.weightedcal import weighted_cal
from transientclumps.pointing_check import pointing_check
from bookkeeping.generate_calfactors_file import generate_calfactors_file
from bookkeeping.make_final_lightcurves import make_final_lightcurves
from bookkeeping.family_members import family_members
from shift_centre.shift_centres_postprocess import shift_centres
from shift_centre.shift_maps import shift_maps
# Can install starlink py-wrapper from: https://github.com/Starlink/starlink-pywrapper ---- or use pip install starlink-pywrapper
from starlink import kappa
import subprocess
import os
import numpy as np
import glob


###############
###############
# Preliminaries/User Input
###############
###############

# The wavelengths to run. Can be ['850'], ['450'] (IF 850 micron data has already been produced! the 450 micron pipeline relies on that data for the pointing corrections!), or ['850','450']
waves = ['850']

# List of Regions To Run. Must select from: IC348, NGC1333, NGC2024, NGC2071, OMC23, OPHCORE, SERPM, SERPS, DR21C, DR21N, DR21S, M17, M17SWex, S255
regions_to_run = ['NGC1333','NGC2024','SERPM','SERPS','M17','M17SWex'] #['IC348', 'NGC1333', 'NGC2024', 'NGC2071', 'OMC23', 'OPHCORE', 'SERPM', 'SERPS', 'DR21C', 'DR21N', 'DR21S', 'M17', 'M17SWex', 'S255']

# The directories  micron, R3 image,'NGCS2071','OMC23','S255's you would like to process in the same order as the region names above
datadirs       = ['NGC1333','NGC2024','SERPM','SERPS','M17','M17SWex'] #['IC348', 'NGC1333', 'NGC2024', 'NGC2071', 'OMC23', 'OPHCORE', 'SERPM', 'SERPS', 'DR21C', 'DR21N', 'DR21S', 'M17', 'M17SWex', 'S255']

# Sets the number of XC/AC corrections to iteratively perform
Number_of_Alignment_Iterations = 1 

# Alternatively, you can set this to be a list, in the same order as regions_to_run and datadirs
# with a different integer for each region:

# Number_of_Alignment_Iterations = [1,1,1,1,1]

# Run Weighted and Source-Skewering Calibrations + Light Curves and Pointing Checks only?
pointsourceonly = False

# Run Cross-correlation and Auot-correlation methods only? (Involves reducing new maps with Starlink)
alignment_ACcalonly = True

# Has 450 micron data already been produced for this region? If so, you can skip re-runnning the first observation and just do the new ones
reduce_first_date_for_450 = False

# The date cutoff for flux normalisation when calibrating the data. We don't want the average flux to change over time for our calibrations!
PScal_cutoff_date = 20210410

# For the Point Source-Skewering calibration -- the theoretical target uncertainty of the selected family (determines how many 450 micron maps will be used in deriving a family)
target_uncertainties = [0.05]  # [0.05] means 5% target uncertainty -- keep it here for now, nothing else will work without manual updates
# For the Point Source-Skewering calibration -- do you want to derive a new family, or use the one that is already derived?
find_new_family = False

# For the weighted calibration 450 micron "Good Box"
TautimesAMCutoff     = 0.12
WeightedCalUncCutoff = 0.05

# Define default (nomincal) FCF corrections
mJy_arcsec2_850_FCF = 2340.0
mJy_beam_850_FCF    = 537000.0

mJy_arcsec2_450_FCF = 4710.0
mJy_beam_450_FCF    = 491000.0

# Top-Level directory in which point source results will be stored
results_dir = 'pointsource_results/'
#####################################################
#####################################################

###############
###############
# Begin Pipeline 
###############
###############

# Correct the format of the number of alignment iterations, if necessary
if isinstance(Number_of_Alignment_Iterations,int):
    XC_alignment_iterations_list = [Number_of_Alignment_Iterations]*len(regions_to_run)
else:
    XC_alignment_iterations_list = Number_of_Alignment_Iterations

##############
##############
# Part 1: Cross and Auto-Correlation Section.
##############
##############

if not pointsourceonly:
    
    for eachregion,eachdatadir,XC_alignment_iterations in zip(regions_to_run,datadirs,XC_alignment_iterations_list):

        # If there are no files to reduce in this region, we can't run the pipline on this region. The number of sdf files should be at least 2 because the first observation must be in each data directory as the comparison point for the relative alignment and AC calibration.        
        if len(sorted(glob.glob(eachdatadir+'/*sdf')))<2: 
            continue
        
        # If 450 microns was supplied first, change order since 850 microns needs to be run to derive pointing corrections 
        if waves == ['450','850']:
    
            # We may wish to run the cross correlation technique to align the images multiple times
            # so this pipeline works for 1 or more iterations of Colton's codes - but we only need to iterate
            # 850 microns! For 450 microns, we will just take the most recent pointing solution.
            # So 850 must be run BEFORE 450
    
            waves = ['850', '450']
    
        # Flag wavelength combinations that won't work
        if waves not in [['450'],['850'],['850','450'],['450','850']]:
            raise Exception("The list 'waves' must be one of: ['450'],['850'],['850','450'], or ['450','850']")
    
        for wave in waves:
            
            ######################################
            ######################################
            # 850 Micron Cross/Auto-correlation Methods
            ######################################
            ######################################
    
            if wave == '850':
                iteration_loop = np.arange(0,XC_alignment_iterations+1,1)
            
                for alignment_iteration in iteration_loop:

                    if alignment_iteration!=iteration_loop[-1]:
                        print('\n\nRunning the Cross-Correlation and Auto-Correlation, iteration {} - 850 microns....'.format(alignment_iteration))
                    else:
                        print('\n\nRunning the Cross-Correlation and Auto-Correlation for final convergence check - 850 microns....'.format(alignment_iteration))
                
                    ######################################
                    ######################################
                    ######################################
                    # Run XC and AC to get pointing offsets and flux cal factors
                    # This will happen 1 more time than the specified iterations 
                    # As a final check that the pointing and fluxcal looks good
                    ######################################
                    ######################################
                    ######################################
                    
                    # Run the Cross Correlation and Auto Correlation on the input data
                    # This will generate a binary file (pickle) with all the data on
                    # The XC and AC fits - giving us the first alignment and relative 
                    # flux cal values.
                
                    if alignment_iteration == 0:
                        make_data_dict(region=eachregion,datadir=eachdatadir,alignment_iteration=alignment_iteration,wavelength=wave)
                        make_table(eachregion,alignment_iteration=alignment_iteration,wavelength=wave)
                
                    # If the aligment has been run in a prior iteration - we need
                    # to point to the previously aligned dataset! So here, datadirs points
                    # to the CR3 files produced on the first iteration. We also
                    # Keep track of the tables by alignment number.
                
                    else:
                        make_data_dict(region=eachregion,datadir=eachregion + '_XCalign_'+str(alignment_iteration),alignment_iteration=alignment_iteration,wavelength=wave)
                        make_table(region=eachregion,alignment_iteration=alignment_iteration,wavelength=wave)
                
                    ######################################
                    ######################################
                    ######################################
                    # Perform Alignment and Rel Flux Cal
                    ######################################
                    ######################################
                    ######################################
                    
                    # Now perform the alignment and flux calibration based on the output table for each region
                    # This requires re-running the data reduction with newly created pointing
                    # correction files and then using kappa's cdiv to apply the Rel FCF
                    
                    
                    # First, create the pointing correction files based on the table for each region - the first file will be skipped
                    # Next, create files listing all the raw data for makemap to use - this creates one for each pointing correction file
                    # Then write a bash script using Python to perform the makemap call over all makemap "infiles"
                    # Run makemap
                    # Remove all the intermediate files

                    # The final iteration is only to check that the pointing corrections have converged - we don't want to reduce the maps again        
                    if alignment_iteration != iteration_loop[-1]:        
        
                        # Define a new directory for the cross/auto-correlation corrected data
                        new_direc_for_corrected_data = eachregion+'_XCalign_'+str(alignment_iteration+1)
                        if not os.path.exists(new_direc_for_corrected_data):
                            os.system('mkdir '+new_direc_for_corrected_data)

                        print('\nRunning Makemap...')

                        # Make Pointing Correction files to feed to makemap for all NEW (previously not corrected) observations
                        make_pcor("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table")

                        # Based on those pointing correction files, compile a list of raw data for each observation
                        makemap_infiles(eachregion,wave)
            
                        # Generate makemap script
                        create_makemap_script(wave)

                        # The *sh file below contains a loop through all supplied epochs
                        subprocess.call('sh ./makemaps.sh',shell=True)

                        # It is important to simply copy over the first epoch since that one is not being re-run with a new pointing correction (since the correction is [0,0])
                        if alignment_iteration == 0:
                            firstfile = sorted(glob.glob(eachdatadir+'/*'+wave+'_ER3.sdf'))[0]
                            newname = firstfile.split('/')[-1].split('_ER3.sdf')[0]+'_CR3.sdf'
                        else:
                            firstfile = sorted(glob.glob(eachregion+'_XCalign_'+str(alignment_iteration)+'/*'+wave+'_CR3_ACcal.sdf'))[0]
                            newname   = firstfile.split('/')[-1].split('_ACcal.sdf')[0]+'.sdf'
                        os.system('cp '+firstfile+' ./'+newname)


                        print('\nMaking the ACcal files - 850 microns....')

                        # Apply the relative FCF correction using kappa
                        apply_relFCF_AC("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table",wave)
                    
                        print('\nConverting the units to mJy/beam....')

                        newCR3files = sorted(list(glob.glob('*'+wave+'*CR3.sdf')))    
                        #Change to mJy/beam!
                        for eachCR3file in newCR3files:
                            firstconv  = mJy_arcsec2_850_FCF
                            convfactor = mJy_beam_850_FCF
                            pwfilename = eachCR3file.split('.sdf')[0]+'_pw.sdf'
                            kappa.cdiv(eachCR3file,firstconv,pwfilename)
                            mJybmfilename = eachCR3file.split('.sdf')[0]+'_mJybm.sdf'
                            kappa.cmult(pwfilename,convfactor,mJybmfilename)
                            kappa.setunits(mJybmfilename,'mJy/beam')
                            os.system('mv -f '+mJybmfilename+' '+eachCR3file)

                        print('\nCleaning Up....')
                        # Move the newly aligned and flux calibrated files to their own directory and cleanup
                        os.system('mv '+eachregion+'*_CR3.sdf '+eachregion+'*CR3_ACcal.sdf '+new_direc_for_corrected_data)
                        os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')
    
            ######################################
            ######################################
            # 450 Micron Cross/Auto-Correlation Methods
            ######################################
            ######################################
        
             ######################################
             # Perform Alignment and Rel Flux Cal based on *850* values
             ######################################
    
            if wave == '450':
                # Check if 850 micron tables already exist for this region and select the 2nd to most recent one! DO NOT CHANGE 850 IN THE LINE BELOW
                # (2nd to most since the most recent is simply a check of the final data - not the actual correction we used to produce the final data)
                if len(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table")) > 0:
                    best_850_correction_table = sorted(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table"))[-2]
                    alignment_iteration = int(best_850_correction_table.split('_')[-2])
                else:
                    raise Exception('This pipeline must first be run on 850 micron data in order to find the pointing offsets to apply to the 450 micron data. The pointing offsets are derived from whichever tables/Transient_'+eachregion+'_run_*_850.table file has the second highest "run" number since the solution converges at higher numbers of iterations, but the highest number is simply a check on the final maps, not the actual correction used to produce the final maps.')

                print('\nRunning Makemap - 450 microns....')
        
                # Throw the flag "reduce_firstepoch" to ensure that we make a 450 micron map for the first epoch if we need to! 
                # (850 always skips this step since no correction is applied anyway and we already have the map)
                make_pcor(best_850_correction_table, reduce_firstepoch=reduce_first_date_for_450)
                makemap_infiles(eachregion,wave)
                create_makemap_script(wave)
                # The *sh file below contains a loop through all supplied epochs
                subprocess.call('sh ./makemaps.sh',shell=True)

                print('\nConverting to units of mJy/beam - 450 microns....')

                newCR3files = sorted(list(glob.glob('*'+wave+'*CR3.sdf')))
                #Change to mJy/beam!
                for eachCR3file in newCR3files:
                    firstconv = mJy_arcsec2_450_FCF
                    convfactor = mJy_beam_450_FCF
                    pwfilename = eachCR3file.split('.sdf')[0]+'_pw.sdf'
                    kappa.cdiv(eachCR3file,firstconv,pwfilename)
                    mJybmfilename = eachCR3file.split('.sdf')[0]+'_mJybm.sdf'
                    kappa.cmult(pwfilename,convfactor,mJybmfilename)
                    kappa.setunits(mJybmfilename,'mJy/beam')
                    os.system('mv -f '+mJybmfilename+' '+eachCR3file)

                # Move the newly aligned and flux calibrated files to their own directory and cleanup
                os.system('mv *'+wave+'*CR3.sdf '+eachregion+'_XCalign_'+str(alignment_iteration+1))
   
                print('\nRunning Cross/Auto-correlation as a final check for convergence - 450 microns....') 
                # Perform a check on the final 450 micron data - we want to make sure the pointing corrections worked!
                make_data_dict(region=eachregion,datadir=eachregion + '_XCalign_'+str(alignment_iteration+1),alignment_iteration=alignment_iteration+1,wavelength=wave)
                make_table(region=eachregion,alignment_iteration=alignment_iteration+1,wavelength=wave)
    
                os.system('mv '+eachregion+'_XCalign_'+str(alignment_iteration+1)+'/*'+wave+'*CR3.sdf .')
    
                print('\nMaking the ACcal files - 450 microns....')
                # Apply the relative FCF correction using kappa
                apply_relFCF_AC("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration+1)+"_"+wave+".table",wave)
                 
                os.system('mv *_CR3.sdf *'+wave+'*CR3_ACcal.sdf '+eachregion+'_XCalign_'+str(alignment_iteration+1))
                os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')

#############################################################
#############################################################

############################
############################
# Part 2: Localised, point source Algorithm Start
############################
############################
   
if not alignment_ACcalonly:

    # Run on the highest iteration of the cross-correlation performed (presumably the most-aligned) 
    for eachregion,eachdatadir,XC_alignment_iterations in zip(regions_to_run,datadirs,XC_alignment_iterations_list):
    
        # Again, correct the wave order - it is possible to run the pipeline skipping Part 1 - so redefine here to make sure 850 microns comes first
        if waves == ['450','850']:
    
            # We may wish to run the cross correlation technique to align the images multiple times
            # so this pipeline works for 1 or more iterations of Colton's codes - but we only need to iterate
            # 850 microns! For 450 microns, we will just take the most recent pointing solution.
            # So 850 must be run BEFORE 450
    
            waves = ['850', '450']
    
        if waves not in [['450'],['850'],['850','450'],['450','850']]:
            raise Exception("The list 'waves' must be one of: ['450'],['850'],['850','450'], or ['450','850']")

        # Start the calibration pipeline    
        for wave in waves:        

            if wave == '450':
                mjypbmfactor     = mJy_beam_450_FCF 
                mjyparcsecfactor = mJy_arcsec2_450_FCF
            elif wave == '850':
                mjypbmfactor     = mJy_beam_850_FCF
                mjyparcsecfactor = mJy_arcsec2_850_FCF

            # Get most recent XC/AC iteration
            most_recent_XC = sorted(glob.glob(eachregion+'_XCalign_*'))[-1]
         
            # Now, with the new files, run the localised point source methods of
            # finding the relative flux calibration factors. 
            # At the same time, we can get the alignment factors based on the localised method
            # to test the cross correlation alignement went smoothly
            
            #################
            #################
            # Prepare for Calibrations
            #################
            #################

            # Smooth all the data (already in mJy/beam)! -- This now scales the flux upwards by ~16%, as well, since the beam size is changing!

            print('\nSmoothing Data....')           
 
            smooth_input_data(most_recent_XC,wave)   
           
            print('\nGenerating Initial Sourceinfo File....')

            # The following generates the sourceinfo file that the weighted calibration uses! 
            make_coadds_metadata(eachregion,most_recent_XC,wave,mjypbmfactor=mjypbmfactor,mjyparcsecfactor=mjyparcsecfactor)

            #################
            #################
            # The Weighted Calibration
            #################
            #################

            print('\nRunning The Weighted Calibration....')

            # Run Weighted Calibration program aand apply results
            # Get Variables and names to give to weighted_cal as sources to avoid (can be empty)!
            var_list,var_names = get_vars(eachregion,wave)
            sourceinfofile = sorted(glob.glob(results_dir+eachregion+'/'+'*_'+wave+'_sourceinfo.txt'))[-1]
            if wave == '450':
                weighted_cal(1000.0,var_list,sourceinfofile,results_dir+eachregion+'/',date_cutoff=PScal_cutoff_date,numiter=5,fid_450_perc=0.05) # Cutoff is in millijanskys!
            else:
                weighted_cal(100.0,var_list,sourceinfofile,results_dir+eachregion+'/',date_cutoff=PScal_cutoff_date,numiter=5,fid_450_perc=0.05) # Cutoff is in millijanskys!
            applyWeightedCal(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave)

           
            # Since we need to consider the weather at 450 microns, we define a box in Tau*Airmass versus WeightedCalUnc space to identify "Good Maps"
            if wave == '450':
                # Need to make sourceinfo and metadata file for GOOD BOX ONLY - so get a list of the datescans from the weightedcal function run, above
                weightedcallookup = np.genfromtxt(sorted(list(glob.glob(results_dir+eachregion+'/*'+wave+'*_calepoch_weightedmean.txt')))[-1],dtype=None,names=True)
                weighted_datescans1     = weightedcallookup['DateScan']
                weighted_datescans      = []
                for i in weighted_datescans1:
                    weighted_datescans.append(str(i)[0:8]+'_'+str(i)[8:])
                weighted_datescans = np.array(weighted_datescans)

                #get good maps based on a TauxAM cutoff and a Weighted Calibration Uncertainty Cutoff for each datescan
                goodmaps      = []
                gooddatescans = []
                for each450file in sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm_Wcal.sdf')):
                    datescan = each450file.split('_')[-6]+'_'+each450file.split('_')[-5]
                    tauvalue = (float(kappa.fitsval(each450file,'WVMTAUST').value)+float(kappa.fitsval(each450file,'WVMTAUEN').value))/2.0
                    AMvalue = (float(kappa.fitsval(each450file,'AMSTART').value)+float(kappa.fitsval(each450file,'AMEND').value))/2.0
                    if tauvalue*AMvalue<=TautimesAMCutoff:
                        if np.array(weightedcallookup['WeightedCalUnc'])[np.where(weighted_datescans==datescan)][0]<=WeightedCalUncCutoff:
                            goodmaps.append(each450file)
                            gooddatescans.append(each450file.split('_')[-6]+'_'+each450file.split('_')[-5])

                # If we have already identified good maps - make sure we skip the processing for those and only deal with new maps
                if len(list(glob.glob(results_dir+eachregion+'/*GoodMaps*metadata.txt')))>0:
                    previously_known_good_maps = []
                    prevmeta = np.genfromtxt(sorted(list(glob.glob(results_dir+eachregion+'/*GoodMaps*metadata.txt')))[-1],dtype=None,names=True)
                    for eachname in prevmeta['Name']:
                        strname = eachname.decode('utf-8')
                        previously_known_good_maps.append(strname.split('_')[-3]+'_'+strname.split('_')[-2])
                else:
                    previously_known_good_maps = []
                newgoodmaptally = 0
                for eachgooddatescan in gooddatescans:
                    if eachgooddatescan not in previously_known_good_maps:
                        newgoodmaptally+=1

                # Produce info files for the 450 micron "Good Maps"
                if newgoodmaptally>0:
                    print('\tNow Working With The Good Maps (sourceinfo file, metadata etc)....')
                    make_coadds_metadata(eachregion,most_recent_XC,wave,mjypbmfactor=mjypbmfactor,mjyparcsecfactor=mjyparcsecfactor,weighted=True,goodbox=True,goodboxfilelist=goodmaps)
      
            print('\nGenerating Souceinfo File Based On Weighted Calibration....')
            # The following generates the sourceinfo file and metadata file that uses the weighted calibration data! 
            make_coadds_metadata(eachregion,most_recent_XC,wave,mjypbmfactor=mjypbmfactor,mjyparcsecfactor=mjyparcsecfactor,weighted=True)
   
            ########################
            ########################
            # Point Source Skewer Method 
            ########################
            ########################
 
            # Generate noises_wave.txt -- To be used with the PScal - so this should refer to the NON ACcal, NON Wcal data -- 2021-05-14
            print('\nExtracting Noise Properties From Observations....')
            make_noise_file(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave) 
            
            # This will give all the files necessary to run make*plottogether*py which uses getfamily
            # and will generate FCF unc family plots as well as a dictionary of datescans and Rel FCFS
            # -- the dictionary will also include the family member IDs
            print('\nFinding Families And Applying The PScal....')
            make_FCFunc_family_plots(eachregion,wave,PScal_cutoff_date,find_new_family=False,target_uncertainties = [0.05],jsonfile='point_source_cal/bright_threshes.json',cat_dir = 'config/')
            
            # Apply to the data.
            
            applyPScal(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave,target_uncertainties)
            
            ###################
            ###################
            # Pointing Check
            ###################
            ###################

            # Run local point source alignment as well to check Colton's Results (850 microns only)
            if wave == '850':
                print('\nChecking The Pointing....')
                pointing_check(most_recent_XC,wave)

            ###################
            ###################
            # Bookkeeping - make nice files from which to extract information
            ###################
            ###################
            
            print('\nGenerating Final Information Files...')    
            # Get calfactors files        
            generate_calfactors_file(most_recent_XC,eachregion,wave)

            print('\nMaking Light Curves....')

            family_members(results_dir+eachregion+'/'+eachregion+'_PointSource_cal_info_dict_targunc5_'+wave+'.pickle',PScal_cutoff_date,wave)

            make_final_lightcurves(results_dir+eachregion+'/'+eachregion+'_'+wave+'_Wcal_sourcecat.bin',sorted(glob.glob(results_dir+eachregion+'/*'+wave+'_CalFactors.txt'))[-1],eachregion,wave)
            if wave == '450':
                print('\tNow making Good Maps Light Curves....')
                make_final_lightcurves(results_dir+eachregion+'/'+eachregion+'_'+wave+'_Wcal_GoodMaps_sourcecat.bin',sorted(glob.glob(results_dir+eachregion+'/*'+wave+'_CalFactors.txt'))[-1],eachregion,wave,GOODMAPS=True)


            #################
            #################
            # Shift All Coordinates in relevant files produced by point source algorithm to be the Centroid, not the first image
            #################
            #################

            print('\nShifting Centroids of sources and updating output files (creating *CoordFix.txt)...')
            sourceinfofiletoshift = sorted(glob.glob('{}/{}/{}*{}_Wcal_sourceinfo.txt'.format(results_dir,eachregion,eachregion,wave)))[-1]
            calfactorsfiletoshift = sorted(glob.glob('{}/{}/{}*{}_CalFactors.txt'.format(results_dir,eachregion,eachregion,wave)))[-1]
            try:
                variablesfiletoshift  = sorted(glob.glob('{}/{}/{}*{}_variables.txt'.format(results_dir,eachregion,eachregion,wave)))[-1]
            except:
                variablesfiletoshift  = 'NOFILE.txt'
            try:
                YSOcomparefiletoshift = sorted(glob.glob('{}/{}/{}_YSOcompare_{}.txt'.format(results_dir,eachregion,eachregion,wave)))[-1]
            except:
                YSOcomparefiletoshift  = 'NOFILE.txt'

            
            shift_centres(sourceinfofiletoshift,calfactorsfiletoshift,variablesfiletoshift,YSOcomparefiletoshift)


            # Now do the GoodMaps sourceinfo file!

            if wave == '450':
                sourceinfofiletoshift = sorted(glob.glob('{}/{}/{}*450_Wcal_GoodMaps_sourceinfo.txt'.format(results_dir,eachregion,eachregion)))[-1]
                calfactorsfiletoshift = 'NOFILE.txt' 
                variablesfiletoshift  = 'NOFILE.txt'
                YSOcomparefiletoshift = 'NOFILE.txt'
                shift_centres(sourceinfofiletoshift,calfactorsfiletoshift,variablesfiletoshift,YSOcomparefiletoshift)


            ### Then also shift the map to match!
            newmaps_to_shift = sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm_Wcal.sdf')) # We can include the first map here because Wcal is produced new every time!
            for eachnewmap in newmaps_to_shift:
                shift_maps(eachnewmap)

            ### Then also the coadds!
            newmaps_to_shift = sorted(glob.glob(results_dir+eachregion+'/*'+wave+'*Wcal*coadd.sdf')) # This includes GoodMaps!
            for eachnewmap in newmaps_to_shift:
                shift_maps(eachnewmap)


print("\n\n############\n############\nC'est Fini!\n############\n############\n\n") 
