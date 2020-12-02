from jcmt_transient_alignment.HighMass_data_gen_EAO import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table
from jcmt_transient_alignment.create_makemap_script import make_pcor, makemap_infiles, create_makemap_script
import subprocess

regions_to_run = ['DR21C']
datadirs       = ['DR21C']
wave = '450'
XC_alignment_iterations = 1

# We may wish to run the cross correlation technique to align the images multiple times
# so this loop works for 1 or more iterations of Colton's codes.

for alignment_iteration in np.arange(0,XC_alignment_iterations,1):

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

    if alignment_iteration == 0:
        make_data_dict(regions=regions_to_run,datadirs=datadirs,alignment_iteration=alignment_iteration)
        make_table(regions_to_run,alignment_iteration=alignment_iteration)
    
    # If the aligment has previously been run in iteration 0 - we need
    # to point to the already aligned data! So here, datadirs points
    # to the CR3 files produced on the first iteration. We also
    # Keep track of the tables by alignment number.

    else:
        make_data_dict(regions=regions_to_run,datadirs=[eachregion + '_XCalign_'+str(alignment_iteration_number-1) for eachregion in regions_to_run],alignment_iteration=alignment_iteration)
        make_table(regions=regions_to_run,alignment_iteration=alignment_iteration)


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

    for eachregion in regions_to_run:
        make_pcor("tables/HM_"+eachregion+"_run_"+str(alignment_iteration)+".table")
        makemap_infiles("tables/HM_"+eachregion+"_run_"+str(alignment_iteration)+".table",eachregion,wave)
        create_makemap_script(wave)
        subprocess.call('makemaps.sh',shell=True)
        if not os.path.exists(eachregion+'_XCalign_'+str(alignment_iteration_number)):
            os.system('mkdir '+eachregion+'_XCalign_'+str(alignment_iteration_number))
        # Apply the relative FCF correction using kappa
        apply_relFCF_AC("tables/HM_"+eachregion+"_run_"+str(alignment_iteration)+".table",wave)
        # Move the newly aligned and flux calibrated files to their own directory
        os.system('mv *CR3_relcal.sdf '+eachregion+'_XCalign_'+str(alignment_iteration_number))
        # Remove the intermediate files
        os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')
        

# Now, with the new files, run the localised point source method of
# finding the relative flux calibration factors. This will avoid any
# variables and improve the original correction. At the same time,
# we can get the alignment factors based on the localised method
# to test the cross correlation alignement went smoothly

#Once config directory is set up and protocat/diskcat are int he right place....
#Run:
#make_coadds_metadata_tables.py

# This will give all the files necessary to run make*plottogether*py which uses getfamily
# Then I will get family members, FCFs, FCFuncs, etc. 
# Apply to the data.
# Find a way to run alignment as well to check Colton's Results

