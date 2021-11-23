import numpy as np
import pickle
from starlink import kappa
import glob
import os

def applyWeightedCal(imagelist,wave,goodmaps=False):

    region = imagelist[0].split('/')[-1].split('_')[0]

    #Get calibration factors
    if goodmaps:
        calinfo = np.genfromtxt(sorted(glob.glob('pointsource_results/'+region+'/*'+wave+'*_calepoch_weightedmean_GoodMaps.txt' ))[-1],dtype=None,names=True)        
    else:
        calinfo = np.genfromtxt(sorted(glob.glob('pointsource_results/'+region+'/*'+wave+'*_calepoch_weightedmean.txt' ))[-1],dtype=None,names=True)
    corrected_datescans = []
    for i in calinfo['DateScan']:
        corrected_datescans.append(str(i)[0:8]+'_'+str(i)[8:])

    for eachimage in imagelist:
        if not goodmaps:
            if os.path.exists(eachimage.split('.sdf')[0]+'_Wcal.sdf'):
                continue
        else:
            if os.path.exists(eachimage.split('.sdf')[0]+'_WcalGM.sdf'):
                continue
        date_scan = eachimage.split('/')[-1].split('_')[1]+'_'+eachimage.split('/')[-1].split('_')[2]
        if date_scan in corrected_datescans:
            calfactor = np.array(calinfo['Divisor'])[np.where(np.array(corrected_datescans)==date_scan)]
            #print(calfactor) 
            if len(calfactor)>0:
                if goodmaps:
                    kappa.cdiv(eachimage,calfactor[0],eachimage.split('.sdf')[0]+'_WcalGM.sdf')
                else:
                    kappa.cdiv(eachimage,calfactor[0],eachimage.split('.sdf')[0]+'_Wcal.sdf')
        
    
