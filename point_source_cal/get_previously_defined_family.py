import pickle
import numpy as np

def get_previously_defined_family(region,wave,eachtargunc,date_cutoff):
    with open('config/Families_info_{}_{}_{}_{}.bin'.format(region,wave,int(float(eachtargunc)),date_cutoff),'rb') as infile:
        family_results = pickle.load(infile)
    return(family_results['SD_of_best_fam_for_plot'],family_results['numcals_for_plot'],family_results['med_FCF_cat_for_plot'],family_results['RMSthresh'],family_results['FCF_dates_all'],family_results['FCFs_all'],family_results['FCF_uncs_all'],family_results['normfluxes_by_date_all'],family_results['families_all'])

