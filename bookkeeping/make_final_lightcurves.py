import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle
from astropy.time import Time

def make_final_lightcurves(source_dict_file,calfactorfile,region,wave,GOODMAPS=False):
    # source_dict is the output_dir+'/'+region+'_'+wave+'_sourcecat.bin' file
    source_dict = pickle.load(open(source_dict_file,'rb'))

    calfactors = np.genfromtxt(calfactorfile,names=True,dtype=None)
    calfactor_dates = np.array(calfactors['Date_Scan'])

    Wcals  = np.array(calfactors['WeightedCal'])
    PScals = np.array(calfactors['PScalWeightedCal'])*np.array(calfactors['WeightedCal']) 

    # This next section is to ensure we are matching the correct Wcal to the correct date - it also selects the datapoints for GoodMaps
    for eachsource in source_dict.keys():
        good_ind = []
        dummy = -1
        for date1 in calfactor_dates:
            dummy=dummy+1
            date1 = str(date1)
            for date2,scan2 in zip(np.array(source_dict[eachsource]['dates_reg'])[np.argsort(np.array(source_dict[eachsource]['dates']))],np.array(source_dict[eachsource]['scan'])[np.argsort(np.array(source_dict[eachsource]['dates']))]):
                if date2 == date1[0:8]:
                    if scan2 == date1[8:]:
                        good_ind.append(dummy)
        

        dates_to_plot    = np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]
        uncal_peakfluxes = Wcals[good_ind]*np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]
        cal_peakfluxes   = np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]

        # UNCAL Plot
        sns.scatterplot(x=dates_to_plot,y=uncal_peakfluxes,color='grey',alpha=0.3,marker='s',label='Uncal') 

        # WCAL PLOTS
        sns.scatterplot(x=dates_to_plot,y=cal_peakfluxes,color='b',marker='o',label='Wcal')

        # Lines, Titles and Labels
        plt.suptitle('Index '+str(source_dict[eachsource]['index'])+': '+eachsource.replace('_',''))
        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])+np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])-np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
        plt.legend(loc='lower left')
        plt.ylabel('Flux (mJy/beam)')
        mindate = min(source_dict[eachsource]['dates'])
        maxdate = max(source_dict[eachsource]['dates'])
        xticksforlabels = np.linspace(mindate,maxdate,5)
        xticklabels = []
        for i in xticksforlabels:
            xticklabels.append(str(Time(i,format='jd').isot).split('T')[0])
        plt.xticks(xticksforlabels,xticklabels,rotation=20)
        if GOODMAPS:
            plt.savefig('pointsource_results/light_curves/'+region+'/'+'index_'+str(source_dict[eachsource]['index']).zfill(4)+'_'+eachsource+'_'+wave+'_Wcal_GoodMaps_lightcurve.pdf',format='pdf')
        else:
            plt.savefig('pointsource_results/light_curves/'+region+'/'+'index_'+str(source_dict[eachsource]['index']).zfill(4)+'_'+eachsource+'_'+wave+'_Wcal_lightcurve.pdf',format='pdf')
        plt.clf()
 
