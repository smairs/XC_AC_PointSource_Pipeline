import numpy as np
import pickle
from starlink import kappa
import glob
import os

def applyWeightedCal(imagelist,wave):

    region = imagelist[0].split('/')[-1].split('_')[0]

    #Get calibration factors
    calinfo = np.genfromtxt(sorted(glob.glob('pointsource_results/'+region+'/*'+wave+'*_calepoch_weightedmean.txt' ))[-1],dtype=None,names=True)
    corrected_datescans = []
    for i in calinfo['DateScan']:
        corrected_datescans.append(str(i)[0:8]+'_'+str(i)[8:])

    for eachimage in imagelist:
        if os.path.exists(eachimage.split('.sdf')[0]+'_Wcal.sdf'):
            continue
        date_scan = eachimage.split('/')[-1].split('_')[1]+'_'+eachimage.split('/')[-1].split('_')[2]
        if date_scan in corrected_datescans:
            calfactor = np.array(calinfo['Divisor'])[np.where(np.array(corrected_datescans)==date_scan)]
            print(calfactor) 
            if len(calfactor)>0:
                kappa.cdiv(eachimage,calfactor[0],eachimage.split('.sdf')[0]+'_Wcal.sdf')
        
    