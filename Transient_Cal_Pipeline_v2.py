from jcmt_transient_alignment.HighMass_data_gen_EAO import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table
from jcmt_transient_alignment.create_makemap_script import make_pcor, makemap_infiles, create_makemap_script
import subprocess

regions_to_run = ['DR21C']
wave = '450'

######################################
######################################
######################################
# STEP 1: Run XC and AC
######################################
######################################
######################################

# Run the Cross Correlation and Auto Correlation on the input data
# This will generate a binary file (pickle) with all the data on
# The XC and AC fits - effectively giving us alignement and relative 
# flux cal values.

make_data_dict(regions=regions_to_run)

# Harvest the data dictionary created above to make a simple txt table

make_table(regions_to_run)

######################################
######################################
######################################
# STEP 2: Perform Alignment and Rel Flux Cal
######################################
######################################
######################################

# Now perform the alignment and flux calibration based on those numbers
# This requires re-running the data reduction with newly created pointing
# correction files and then using kappa's cdiv to apply the Rel FCF


# First, create the pointing correction files based on the table for each region - the first file will be skipped
# Next, create files listing all the raw data for makemap to use - this creates one for each pointing correction file
# Then write a bash script using Python to perform the makemap call
# Run makemap
# Remove all the intermediate files

for i in regions_to_run:
    make_pcor("tables/HM_{:}.table".format(i))
    makemap_infiles("tables/HM_{:}.table".format(i),i,wave)
    create_makemap_script(wave)
    subprocess.call('makemaps.sh',shell=True)
    os.system('rm -f *pcor.txt  makemaps.sh *_'+wave+'.txt')


# Next, apply the relative FCF correction using kappa!


# ^^ May need to run the above steps multiple times - so the makemap "out" names might get confusing!


# Now, with the new files, run the localised point source method of
# finding the relative flux calibration factors. This will avoid any
# variables and improve the original correction. At the same time,
# we can get the alignment factors based on the localised method
# to test the cross correlation alignement went smoothly

