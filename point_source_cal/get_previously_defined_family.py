import pickle
import numpy as np

def get_previously_defined_family(region,wave,eachtargunc,date_cutoff,gooddatescans=[]):
    with open('config/Families_info_{}_{}_{}_{}.bin'.format(region,wave,int(float(eachtargunc)),date_cutoff),'rb') as infile:
        family_results = pickle.load(infile)
        if len(gooddatescans)==0:
            print('PS CAL BEING PERFORMED OVER ALL MAPS')
            return(family_results['SD_of_best_fam_for_plot'],family_results['numcals_for_plot'],family_results['med_FCF_cat_for_plot'],family_results['RMSthresh'],family_results['FCF_dates_all'],family_results['FCFs_all'],family_results['FCF_uncs_all'],family_results['normfluxes_by_date_all'],family_results['families_all'])

        else:
            print('PS CAL BEING PERFORMED OVER GOOD MAPS ONLY:')

            goodindex = [] 
            for i in family_results['FCF_dates_all'][0]:
                if i in gooddatescans:
                    goodindex.append(True)
                else:
                    goodindex.append(False)
            for eachind in range(len(family_results['FCF_dates_all'])):
                family_results['FCF_dates_all'][eachind]          = list(np.array(family_results['FCF_dates_all'][eachind])[goodindex])
                family_results['FCFs_all'][eachind]               = list(np.array(family_results['FCFs_all'][eachind])[goodindex])
                family_results['FCF_uncs_all'][eachind]           = list(np.array(family_results['FCF_uncs_all'][eachind])[goodindex])
                family_results['normfluxes_by_date_all'][eachind] = list(np.array(family_results['normfluxes_by_date_all'][eachind])[goodindex])

            print('PSCAL DATES: ',family_results['FCF_dates_all'][0])
            return(family_results['SD_of_best_fam_for_plot'],family_results['numcals_for_plot'],family_results['med_FCF_cat_for_plot'],family_results['RMSthresh'],family_results['FCF_dates_all'],family_results['FCFs_all'],family_results['FCF_uncs_all'],family_results['normfluxes_by_date_all'],family_results['families_all'])

