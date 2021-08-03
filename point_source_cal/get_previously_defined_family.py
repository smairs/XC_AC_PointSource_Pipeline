import json

def get_previously_definied_family(region,wave,eachtargunc,date_cutoff):
    with open('config/Families_info_{}_{}_{}_{}.json'.format(region,wave,eachtargunc,date_cutoff),'w') as json_file:
        family_results = json.load(json_file)
        return(family_results['SD_of_best_fam_for_plot'],family_results['numcals_for_plot'],family_results['med_FCF_cat_for_plot'],family_results['RMSthresh'],family_results['FCF_dates_all'],family_results['FCFs_all'],family_results['FCF_uncs_all'],family_results['normfluxes_by_date_all'],family_results['families_all'])

