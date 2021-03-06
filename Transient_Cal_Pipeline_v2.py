from jcmt_transient_alignment.HighMass_data_gen_EAO import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table
from jcmt_transient_alignment.create_makemap_script import make_pcor, makemap_infiles, create_makemap_script
from jcmt_transient_alignment.applycal import apply_relFCF_AC
from point_source_cal.make_coadds_metadata_tables import make_coadds_metadata
from point_source_cal.noisefunctions import make_noise_file
from point_source_cal.smooth_input_data import smooth_input_data
from point_source_cal.changeunits import changeunits
from point_source_cal.applyPScal import applyPScal
from point_source_cal.make_FCFunc_vs_sourceSNR_plots_binsof1_plottogether import make_FCFunc_family_plots
from transientclumps.pointing_check import pointing_check
from transientclumps.GCfluxcal import GCfluxcal
import subprocess
import os
import numpy as np
import glob

regions_to_run = ['DR21C']
datadirs       = ['DR21C']
wave = '450'
XC_alignment_iterations = 3 #This parameter only matters if running the 850 micron pipeline.
target_uncertainties = [0.05]  # [0.05] means 5% target uncertainty -- keep it here for now, nothing else will work without manual updates

# We may wish to run the cross correlation technique to align the images multiple times
# so this loop works for 1 or more iterations of Colton's codes - but we only need to iterate
# 850 microns! For 450 microns, we will just take the most recent pointing solution -- that is
# what the first part of this code does. So make sure 850 is run before 450!

if wave == '450':

    # Check if 850 micron tables already exist for this region
    for eachregion in regions_to_run:
        if len(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table")) > 0:
            best_850_correction_table = sorted(glob.glob("tables/Transient_"+eachregion+"_run_*_850.table"))[-1]
            iteration_loop = [int(best_850_correction_table.split('_')[-2])]
        else:
            raise Exception('This pipeline must first be run on 850 micron data in order to find the pointing offsets to apply to the 450 micron data. The pointing offsets are derived from whichever tables/Transient_'+eachregion+'_run_*_850.table file has the highest "run" number since the solution converges at higher numbers of iterations.')

if wave == '850':
    iteration_loop = np.arange(0,XC_alignment_iterations,1)

for alignment_iteration in iteration_loop:

    ######################################
    ######################################
    ######################################
    # STEP 1: Run XC and AC
    ######################################
    ######################################
    ######################################
    
    # Run the Cross Correlation and Auto Correlation on the input data
    # This will generate a binary file (pickle) with all the data on
    # The XC and AC fits - effectively giving us the first alignment and relative 
    # flux cal values.


    if wave == '850':
        if alignment_iteration == 0:
            make_data_dict(regions=regions_to_run,datadirs=datadirs,alignment_iteration=alignment_iteration,wavelength=wave)
            make_table(regions_to_run,alignment_iteration=alignment_iteration,wavelength=wave)

        # If the aligment has previously been run in iteration 0 - we need
        # to point to the already aligned data! So here, datadirs points
        # to the CR3 files produced on the first iteration. We also
        # Keep track of the tables by alignment number.

        else:
            make_data_dict(regions=regions_to_run,datadirs=[eachregion + '_XCalign_'+str(alignment_iteration) for eachregion in regions_to_run],alignment_iteration=alignment_iteration,wavelength=wave)
            make_table(regions=regions_to_run,alignment_iteration=alignment_iteration,wavelength=wave)

    if wave  == '450':
        make_data_dict(regions=regions_to_run,datadirs=datadirs,alignment_iteration=0,wavelength=wave)
        make_table(regions_to_run,alignment_iteration=0,wavelength=wave)

 
    ######################################
    ######################################
    ######################################
    # STEP 2: Perform Alignment and Rel Flux Cal
    ######################################
    ######################################
    ######################################
    
    # Now perform the alignment and flux calibration based on the output table for each region
    # This requires re-running the data reduction with newly created pointing
    # correction files and then using kappa's cdiv to apply the Rel FCF
    
    
    # First, create the pointing correction files based on the table for each region - the first file will be skipped
    # Next, create files listing all the raw data for makemap to use - this creates one for each pointing correction file
    # Then write a bash script using Python to perform the makemap call
    # Run makemap
    # Remove all the intermediate files

    for eachregion,eachdatadir in zip(regions_to_run,datadirs):
        if wave == '850':
            if not os.path.exists(eachregion+'_XCalign_'+str(alignment_iteration+1)):
                os.system('mkdir '+eachregion+'_XCalign_'+str(alignment_iteration+1))
        # In the line below, if a 450 catalogue is supplied to make_pcor, the make_pcor code will change it to 850 microns 
        # with the same iteration number
        make_pcor("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table")
        makemap_infiles(eachregion,wave)
        create_makemap_script(wave)
        subprocess.call('sh ./makemaps.sh',shell=True)
        if wave == '850':
            if alignment_iteration == 0:
                firstfile = sorted(glob.glob(eachdatadir+'/*'+wave+'_ER3.sdf'))[0]
                newname = firstfile.split('/')[-1].split('_ER3.sdf')[0]+'_CR3.sdf'
            else:
                firstfile = sorted(glob.glob(eachregion+'_XCalign_'+str(alignment_iteration-1)+'/*'+wave+'_CR3_ACcal.sdf'))[0]
                newname   = firstfile.split('/')[-1].split('_ACcal.sdf')[0]+'.sdf'
        elif wave == '450':
            firstfile = sorted(glob.glob(eachdatadir+'/*'+wave+'_ER3.sdf'))[0]
            newname = firstfile.split('/')[-1].split('_ER3.sdf')[0]+'_CR3.sdf'
        os.system('cp '+firstfile+' ./'+newname)
        # Apply the relative FCF correction using kappa
        apply_relFCF_AC("tables/Transient_"+eachregion+"_run_"+str(alignment_iteration)+"_"+wave+".table",wave)
        if wave == '850':
            # Move the newly aligned and flux calibrated files to their own directory
            os.system('mv *CR3_ACcal.sdf '+eachregion+'_XCalign_'+str(alignment_iteration+1))
        elif wave == '450':
            # Move the newly aligned and flux calibrated files to their own directory
            os.system('mv *CR3_ACcal.sdf '+eachregion+'_XCalign_'+str(alignment_iteration))
        # Remove the intermediate files
        os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt *_CR3.sdf')
     
#####
#####
#####   
#####
# Localised Algorithm Start
#####
#####
#####
#####

# Get most recent XC/AC iteration

for eachregion in regions_to_run:
    allfiles = os.listdir('./')
    XCalign_dirs = []
    for i in allfiles:
        if os.path.isdir(i):
            if len(i.split('_'))>1:
                if i.split('_')[-2] == 'XCalign':
                    if i.split('_')[0] == eachregion:
                        XCalign_dirs.append(i)
    align_numbers = []
    for i in XCalign_dirs:
        align_numbers.append(int(i.split('_')[-1]))
    most_recent_XC = eachregion+'_XCalign_'+str(max(align_numbers))
    
    # Now, with the new files, run the localised point source method of
    # finding the relative flux calibration factors. This will avoid any
    # variables and improve the original correction. At the same time,
    # we can get the alignment factors based on the localised method
    # to test the cross correlation alignement went smoothly
    
    # Change data units to mJy/bm

    changeunits(most_recent_XC,wave) 

    # Smooth all the data!

    smooth_input_data(most_recent_XC,wave)   
 
    make_coadds_metadata(eachregion,most_recent_XC,wave)

    # Generate noises_wave.txt
   
    make_noise_file(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave) 
    
    # This will give all the files necessary to run make*plottogether*py which uses getfamily
    # and will generate FCF unc family plots as well as a dictionary of datescans and Rel FCFS
    # -- the dictionary will also include the family member IDs

    make_FCFunc_family_plots([eachregion],wave,target_uncertainties) 

    # Apply to the data.

    applyPScal(sorted(glob.glob(most_recent_XC+'/*'+wave+'*sm.sdf')),wave,target_uncertainties)

    # Run local point source alignment as well to check Colton's Results

    if wave == '850':
        pointing_check(most_recent_XC,wave)
        GCfluxcal(most_recent_XC,wave)
    
