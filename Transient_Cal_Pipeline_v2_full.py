from jcmt_transient_alignment.HighMass_data_gen_EAO_newbox import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table
from jcmt_transient_alignment.create_makemap_script import make_pcor, makemap_infiles, create_makemap_script
from jcmt_transient_alignment.applycal import apply_relFCF_AC
from point_source_cal.make_coadds_metadata_tables import make_coadds_metadata
from point_source_cal.noisefunctions import make_noise_file
from point_source_cal.smooth_input_data import smooth_input_data
from point_source_cal.changeunits import changeunits
from point_source_cal.applyPScal import applyPScal
from point_source_cal.applyWeightedCal import applyWeightedCal
from point_source_cal.make_FCFunc_vs_sourceSNR_plots_binsof1_plottogether import make_FCFunc_family_plots
from point_source_cal.get_vars import get_vars
from weightedcal.weightedcal import weighted_cal
from transientclumps.pointing_check import pointing_check
from transientclumps.GCfluxcal import GCfluxcal
from bookkeeping.generate_calfactors_file import generate_calfactors_file
from bookkeeping.make_final_lightcurves import make_final_lightcurves
from bookkeeping.family_members import family_members
# Can install starlink py-wrapper from: https://github.com/Starlink/starlink-pywrapper ---- or use pip install starlink-pywrapper
from starlink import kappa
import subprocess
import os
import numpy as np
import glob

regions_to_run = ['DR21C','DR21N','DR21S','M17','M17SWex']
datadirs       = ['DR21C','DR21N','DR21S','M17','M17SWex']
XC_alignment_iterations_list = [1,1,1,1,1,1,1] #This parameter only matters if running the 850 micron pipeline. - it sets the number of XC/AC corrections to iteratively perform
pointsourceonly = False
alignment_ACcalonly = True

reduce_first_date_for_450 = False

PScal_cutoff_date = 20210410

waves = ['850','450']
target_uncertainties = [0.05]  # [0.05] means 5% target uncertainty -- keep it here for now, nothing else will work without manual updates

if not pointsourceonly:
    
    for eachregion,eachdatadir,XC_alignment_iterations in zip(regions_to_run,datadirs,XC_alignment_iterations_list):
        
        if len(sorted(glob.glob(eachdatadir+'/*sdf')))<2: # -- If there are no files to reduce in this region, we can't run the pipline on this region. The number of sdf files should be at least 2 because the first observation must be in each data directory as the comparison point for the relative alignment and AC calibration.
            continue
        
        
        if waves == ['450','850']:
    
            # We may wish to run the cross correlation technique to align the images multiple times
            # so this pipeline works for 1 or more iterations of Colton's codes - but we only need to iterate
            # 850 microns! For 450 microns, we will just take the most recent pointing solution.
            # So 850 must be run BEFORE 450
    
            waves = ['850', '450']
    
        if waves not in [['450'],['850'],['850','450'],['450','850']]:
            raise Exception("The list 'waves' must be one of: ['450'],['850'],['850','450'], or ['450','850']")
    
        for wave in waves:
            
            #############################################################
    
            ############################
            ############################
            ############################
            ############################
            ############################
            ############################
            ############################
            # Global XC/AC Algorithm Start
            ############################
            ############################
            ############################
            ############################
            ############################
            ############################
            ############################
    
            ##############################################################
    
    
            ######################################
            ######################################
            ######################################
            ######################################
            ######################################
            # 850 Micron Pipeline
            ######################################
            ######################################
            ######################################
            ######################################
    
            if wave == '850':
                iteration_loop = np.arange(0,XC_alignment_iterations+1,1)
            
                for alignment_iteration in iteration_loop:
                
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
                    # The XC and AC fits - effectively giving us the first alignment and relative 
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
        
                    if alignment_iteration != iteration_loop[-1]:        
        
                        new_direc_for_corrected_data = eachregion+'_XCalign_'+str(alignment_iteration+1)
    
                        if not os.path.exists(new_direc_for_corrected_data):
                            os.system('mkdir '+new_direc_for_corrected_data)
                        make_pcor("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table")
                        makemap_infiles(eachregion,wave)
                        create_makemap_script(wave)
                        # The *sh file below contains a loop through all supplied epochs
                        subprocess.call('sh ./makemaps.sh',shell=True)
                        # It is important to simply copy over the first epoch since that one is not being re-run with a new pointgin correction (since the correction is [0,0])
                        if alignment_iteration == 0:
                            firstfile = sorted(glob.glob(eachdatadir+'/*'+wave+'_ER3.sdf'))[0]
                            newname = firstfile.split('/')[-1].split('_ER3.sdf')[0]+'_CR3.sdf'
                        else:
                            firstfile = sorted(glob.glob(eachregion+'_XCalign_'+str(alignment_iteration)+'/*'+wave+'_CR3_ACcal.sdf'))[0]
                            newname   = firstfile.split('/')[-1].split('_ACcal.sdf')[0]+'.sdf'
                        os.system('cp '+firstfile+' ./'+newname)
                        # Apply the relative FCF correction using kappa
                        apply_relFCF_AC("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table",wave)
                    
                        newCR3files = sorted(list(glob.glob('*'+wave+'*CR3.sdf')))    
                        #Change to mJy/beam!
                        for eachCR3file in newCR3files:
                            firstconv = 2340
                            convfactor = 537000
                            pwfilename = eachCR3file.split('.sdf')[0]+'_pw.sdf'
                            kappa.cdiv(eachCR3file,firstconv,pwfilename)
                            mJybmfilename = eachCR3file.split('.sdf')[0]+'_mJybm.sdf'
                            kappa.cmult(pwfilename,convfactor,mJybmfilename)
                            kappa.setunits(mJybmfilename,'mJy/beam')
                            os.system('mv -f '+mJybmfilename+' '+eachCR3file)

                        # Move the newly aligned and flux calibrated files to their own directory
                        os.system('mv '+eachregion+'*_CR3.sdf '+eachregion+'*CR3_ACcal.sdf '+new_direc_for_corrected_data)
                        os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')
    
            ######################################
            ######################################
            ######################################
            ######################################
            ######################################
            # 450 Micron Pipeline
            ######################################
            ######################################
            ######################################
            ######################################
        
             ######################################
             ######################################
             ######################################
             # Perform Alignment and Rel Flux Cal based on *850* values
             ######################################
             ######################################
             ######################################
    
    
            if wave == '450':
                # Check if 850 micron tables already exist for this region and select the 2nd to most recent one! DO NOT CHANGE 850 IN THE LINE BELOW
                # (2nd to most since the most recent is simply a check of the final data - not the actual correction we used to produce the final data)
                if len(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table")) > 0:
                    best_850_correction_table = sorted(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table"))[-2]
                    alignment_iteration = int(best_850_correction_table.split('_')[-2])
                else:
                    raise Exception('This pipeline must first be run on 850 micron data in order to find the pointing offsets to apply to the 450 micron data. The pointing offsets are derived from whichever tables/Transient_'+eachregion+'_run_*_850.table file has the second highest "run" number since the solution converges at higher numbers of iterations, but the highest number is simply a check on the final maps, not the actual correction used to produce the final maps.')
        
                # Throw the flag "reduce_firstepoch" to ensure that we make a 450 micron map for the first epoch if we need to! 
                # (850 always skips this step since no correction is applied anyway and we already have the map)
                make_pcor(best_850_correction_table, reduce_firstepoch=reduce_first_date_for_450)
                makemap_infiles(eachregion,wave)
                create_makemap_script(wave)
                # The *sh file below contains a loop through all supplied epochs
                subprocess.call('sh ./makemaps.sh',shell=True)
                newCR3files = sorted(list(glob.glob('*'+wave+'*CR3.sdf')))

                #Change to mJy/beam!
                for eachCR3file in newCR3files:
                    firstconv = 4710
                    convfactor = 491000
                    pwfilename = eachCR3file.split('.sdf')[0]+'_pw.sdf'
                    kappa.cdiv(eachCR3file,firstconv,pwfilename)
                    mJybmfilename = eachCR3file.split('.sdf')[0]+'_mJybm.sdf'
                    kappa.cmult(pwfilename,convfactor,mJybmfilename)
                    kappa.setunits(mJybmfilename,'mJy/beam')
                    os.system('mv -f '+mJybmfilename+' '+eachCR3file)
 
                os.system('mv *'+wave+'*CR3.sdf '+eachregion+'_XCalign_'+str(alignment_iteration+1))
    
                # Perform a check on the final 450 micron data - we want to make sure the pointing corrections worked!
                make_data_dict(region=eachregion,datadir=eachregion + '_XCalign_'+str(alignment_iteration+1),alignment_iteration=alignment_iteration+1,wavelength=wave)
                make_table(region=eachregion,alignment_iteration=alignment_iteration+1,wavelength=wave)
    
                os.system('mv '+eachregion+'_XCalign_'+str(alignment_iteration+1)+'/*'+wave+'*CR3.sdf .')
    
                # Apply the relative FCF correction using kappa
                apply_relFCF_AC("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration+1)+"_"+wave+".table",wave)
                 
                os.system('mv *_CR3.sdf *'+wave+'*CR3_ACcal.sdf '+eachregion+'_XCalign_'+str(alignment_iteration+1))
                os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')

#############################################################


if not alignment_ACcalonly:
    
    ############################
    ############################
    ############################
    ############################
    ############################
    ############################
    ############################
    # Localised, point source Algorithm Start
    ############################
    ############################
    ############################
    ############################
    ############################
    ############################
    ############################
    
    ##############################################################
    for eachregion,eachdatadir,XC_alignment_iterations in zip(regions_to_run,datadirs,XC_alignment_iterations_list):
    
        if waves == ['450','850']:
    
            # We may wish to run the cross correlation technique to align the images multiple times
            # so this pipeline works for 1 or more iterations of Colton's codes - but we only need to iterate
            # 850 microns! For 450 microns, we will just take the most recent pointing solution.
            # So 850 must be run BEFORE 450
    
            waves = ['850', '450']
    
        if waves not in [['450'],['850'],['850','450'],['450','850']]:
            raise Exception("The list 'waves' must be one of: ['450'],['850'],['850','450'], or ['450','850']")
    
        for wave in waves:        
            # Get most recent XC/AC iteration
            most_recent_XC = sorted(glob.glob(eachregion+'_XCalign_*'))[-1]
         
            # Now, with the new files, run the localised point source method of
            # finding the relative flux calibration factors. This will avoid any
            # variables and improve the original correction. At the same time,
            # we can get the alignment factors based on the localised method
            # to test the cross correlation alignement went smoothly
            
            # Change data units to mJy/bm ---- 2020-05-26 - DATA IS NOW ALREADY IN mJy/bm -- NO NEED TO CHANGE IT!
            
            #changeunits(most_recent_XC,wave) 
            
            # Smooth all the data!
            
            smooth_input_data(most_recent_XC,wave)   
           
            # The following generates the sourceinfo file that the weighted calibration uses! 
            make_coadds_metadata(eachregion,most_recent_XC,wave)

            # Run Doug's Calibration program
            var_list,var_names = get_vars(eachregion,wave)
            sourceinfofile = sorted(glob.glob('pointsource_results/'+eachregion+'/'+'*_'+wave+'_sourceinfo.txt'))[-1]
            print('SOURCEINFOFILE: ',sourceinfofile)
            if wave == '450':
                weighted_cal(1000.0,var_list,sourceinfofile,'pointsource_results/'+eachregion+'/',date_cutoff=PScal_cutoff_date,numiter=5,fid_450_perc=0.05) # Cutoff is in millijanskys!
            else:
                weighted_cal(100.0,var_list,sourceinfofile,'pointsource_results/'+eachregion+'/',date_cutoff=PScal_cutoff_date,numiter=5,fid_450_perc=0.05) # Cutoff is in millijanskys!


            applyWeightedCal(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave)

            if wave == '450':
                # Need to make sourceinfo and metadata file for GOOD BOX ONLY
                weightedcallookup = np.genfromtxt(sorted(list(glob.glob('pointsource_results/'+eachregion+'/*'+wave+'*_calepoch_weightedmean.txt')))[-1],dtype=None,names=True)
                weighted_datescans1     = weightedcallookup['DateScan']
                weighted_datescans      = []
                for i in weighted_datescans1:
                    weighted_datescans.append(str(i)[0:8]+'_'+str(i)[8:])
                    #print(str(i)[0:8]+'_'+str(i)[8:])
                weighted_datescans = np.array(weighted_datescans)

                goodmaps      = []
                gooddatescans = []
                #get good maps
                for each450file in sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm_Wcal.sdf')):
                    datescan = each450file.split('_')[-6]+'_'+each450file.split('_')[-5]
                    tauvalue = (float(kappa.fitsval(each450file,'WVMTAUST').value)+float(kappa.fitsval(each450file,'WVMTAUEN').value))/2.0
                    AMvalue = (float(kappa.fitsval(each450file,'AMSTART').value)+float(kappa.fitsval(each450file,'AMEND').value))/2.0
                    if tauvalue*AMvalue<=0.14:
                        if np.array(weightedcallookup['WeightedCalUnc'])[np.where(weighted_datescans==datescan)][0]<=0.05:
                            goodmaps.append(each450file)
                            gooddatescans.append(each450file.split('_')[-6]+'_'+each450file.split('_')[-5])

                previously_known_good_maps = []
                prevmeta = np.genfromtxt(sorted(list(glob.glob('pointsource_results/'+eachregion+'/*GoodMaps*metadata.txt')))[-1],dtype=None,names=True)
                for eachname in prevmeta['Name']:
                    strname = eachname.decode('utf-8')
                    previously_known_good_maps.append(strname.split('_')[-3]+'_'+strname.split('_')[-2])

                print(gooddatescans,previously_known_good_maps)

                newgoodmaptally = 0
                for eachgooddatescan in gooddatescans:
                    if eachgooddatescan not in previously_known_good_maps:
                        newgoodmaptally+=1

                if newgoodmaptally>0:
                    make_coadds_metadata(eachregion,most_recent_XC,wave,weighted=True,goodbox=True,goodboxfilelist=goodmaps)
      
            # The following generates the sourceinfo file and metadata file that uses the weighted calibration data! 
            make_coadds_metadata(eachregion,most_recent_XC,wave,weighted=True)
    
            # Generate noises_wave.txt -- To be used with the PScal - so this should refer to the NON ACcal, NON Wcal data -- 2021-05-14
            
            make_noise_file(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave) 
            
            # This will give all the files necessary to run make*plottogether*py which uses getfamily
            # and will generate FCF unc family plots as well as a dictionary of datescans and Rel FCFS
            # -- the dictionary will also include the family member IDs
            
            make_FCFunc_family_plots(eachregion,wave,target_uncertainties,PScal_cutoff_date) 
            
            # Apply to the data.
            
            applyPScal(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave,target_uncertainties)
            
            # Run local point source alignment as well to check Colton's Results (850 microns only)
            if wave == '850':
                pointing_check(most_recent_XC,wave)
                
            # Get calfactors files        
            generate_calfactors_file(most_recent_XC,eachregion,wave)
            family_members('pointsource_results/'+eachregion+'/'+eachregion+'_PointSource_cal_info_dict_targunc5_'+wave+'.pickle',PScal_cutoff_date,wave)
            make_final_lightcurves('pointsource_results/'+eachregion+'/'+eachregion+'_'+wave+'_Wcal_sourcecat.bin',sorted(glob.glob('pointsource_results/'+eachregion+'/*'+wave+'_CalFactors.txt'))[-1],eachregion,wave)
            if wave == '450':
                make_final_lightcurves('pointsource_results/'+eachregion+'/'+eachregion+'_'+wave+'_Wcal_GoodMaps_sourcecat.bin',sorted(glob.glob('pointsource_results/'+eachregion+'/*'+wave+'_CalFactors.txt'))[-1],eachregion,wave,GOODMAPS=True)
    
