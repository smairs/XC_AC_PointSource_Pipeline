import numpy as np
import pickle
from starlink import kappa
import os

def applyPScal(imagelist,wave,target_uncertainties,goodmaps=False):

    region = imagelist[0].split('/')[-1].split('_')[0]
    target_uncertainty_key = str(int(100*target_uncertainties[0]))

    #Get calibration factors
    if goodmaps:
        calinfo = pickle.load(open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc'+target_uncertainty_key+'_'+wave+'_GoodMaps.pickle','rb'))
    else:
        calinfo = pickle.load(open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc'+target_uncertainty_key+'_'+wave+'.pickle','rb'))

    for eachimage in imagelist:
        if not goodmaps:
            if os.path.exists(eachimage.split('.sdf')[0]+'_PScal.sdf'):
                continue
        else:
            if os.path.exists(eachimage.split('.sdf')[0]+'_PScalGM.sdf'):
                continue
        date_scan = eachimage.split('/')[-1].split('_')[1]+'_'+eachimage.split('/')[-1].split('_')[2]
        if date_scan in calinfo['datescans']:
            calfactor = np.array(calinfo['RelFCFs'])[np.where(np.array(calinfo['datescans'])==date_scan)]
            print(calfactor) 
            if len(calfactor)>0:
                if not np.isnan(calfactor[0]):
                    if goodmaps:
                        kappa.cdiv(eachimage,calfactor[0],eachimage.split('.sdf')[0]+'_PScalGM.sdf')
                    else:
                        kappa.cdiv(eachimage,calfactor[0],eachimage.split('.sdf')[0]+'_PScal.sdf')
