from jcmt_transient_alignment.HighMass_data_gen_EAO import make_data_dict
from jcmt_transient_alignment.HighMass_table_gen_EAO import make_table

regions_to_run = ['DR21C']

# Run the Cross Correlation and Auto Correlation on the input data
# This will generate a binary file (pickle) with all the data on
# The XC and AC fits - effectively giving us alignement and relative 
# flux cal values.

make_data_dict(regions=regions_to_run)

# Harvest the data dictionary created above to make a simple txt table

make_table(regions_to_run)

# Now perform the alignment and flux calibration based on those numbers
# This requires re-running the data reduction with newly created pointing
# correction files and then using kappa's cdiv to apply the Rel FCF




# ^^ May need to run the above section twice!


# Now, with the new files, run the localised point source method of
# finding the relative flux calibration factors. This will avoid any
# variables and improve the original correction. At the same time,
# we can get the alignment factors based on the localised method
# to test the cross correlation alignement went smoothly

