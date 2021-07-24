import numpy as np
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
        dummy = -1                

# 2021-05-14 --- UPDATING SO WE JUST HAVE WCAL IN THE PLOTS! Note that this commented out section assumed the source_dict_file referrred to the ACcal data. THIS IS NO LONGER THE CASE -- now the source_dict_cile refers to the Wcal data.

#        #ACCAL PLOTS
#        plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label=str(source_dict[eachsource]['index'])+':'+eachsource.replace('_','')+' - ACcal',color='b')
#        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])+np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
#        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])-np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
#
#        #PSCAL PLOTS
#        plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/PScals,label='PScal',color='darkgoldenrod',marker='s')
#        plt.axhline(y=np.nanmean(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/PScals)+np.nanstd(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/PScals,ddof=1),color='darkgoldenrod',linestyle='dashed')
#        plt.axhline(y=np.nanmean(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/PScals)-np.nanstd(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/PScals,ddof=1),color='darkgoldenrod',linestyle='dashed')
#
#
#        #WCAL PLOTS
#        plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/Wcals,label='Wcal',color='g',marker='*')
#        plt.axhline(y=np.nanmean(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/Wcals)+np.nanstd(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/Wcals,ddof=1),color='g',linestyle='dashed')
#        plt.axhline(y=np.nanmean(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/Wcals)-np.nanstd(np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))]/Wcals,ddof=1),color='g',linestyle='dashed')
#
#        #FIDUCIAL
#        #plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes']),color='k',linestyle='solid')
#        #plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])+source_dict[eachsource]['sd_fiducial'],color='k',linestyle='dotted')
#        #plt.axhline(y=np.average(source_dict[eachsource]['peakfluxes'])-source_dict[eachsource]['sd_fiducial'],color='k',linestyle='dotted')
#        plt.legend(loc='lower left')
#        plt.savefig('pointsource_results/light_curves/'+region+'/'+'index_'+str(source_dict[eachsource]['index']).zfill(4)+'_'+eachsource+'_'+wave+'_lightcurve.pdf',format='pdf')
#        plt.clf()

        #UNCAL PLOTS
        #plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,Wcals*np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label=str(source_dict[eachsource]['index'])+':'+eachsource.replace('_','')+' - Uncal',color='grey',alpha=0.3,marker='s')
        plt.scatter(np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))],Wcals[good_ind]*np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label='Uncal',color='grey',alpha=0.3,marker='s')
        #plt.axhline(y=Wcals*np.nanmean(source_dict[eachsource]['peakfluxes'])+Wcals*np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='grey',alpha=0.3,linestyle='dashed')
        #plt.axhline(y=Wcals*np.nanmean(source_dict[eachsource]['peakfluxes'])-Wcals*np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='grey',alpha=0.3,linestyle='dashed')

        #WCAL PLOTS
        #plt.scatter((np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))]-np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))][0])/365.24,np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label=str(source_dict[eachsource]['index'])+':'+eachsource.replace('_','')+' - Wcal',color='b')
        plt.scatter(np.array(source_dict[eachsource]['dates'])[np.argsort(np.array(source_dict[eachsource]['dates']))],np.array(source_dict[eachsource]['peakfluxes'])[np.argsort(np.array(source_dict[eachsource]['dates']))],label='Wcal',color='b')
        plt.suptitle('Index '+str(source_dict[eachsource]['index'])+': '+eachsource.replace('_',''))
        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])+np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
        plt.axhline(y=np.nanmean(source_dict[eachsource]['peakfluxes'])-np.nanstd(source_dict[eachsource]['peakfluxes'],ddof=1),color='b',linestyle='dashed')
        plt.legend(loc='lower left')
        plt.ylabel('Flux (mJy/beam)')
        mindate = min(source_dict[eachsource]['dates'])
        maxdate = max(source_dict[eachsource]['dates'])
        xticksforlabels = np.arange(mindate,maxdate+((maxdate-mindate)/5.0),(maxdate-mindate)/5.0)
        xticklabels = []
        for i in xticksforlabels:
            xticklabels.append(str(Time(i,format='jd').isot).split('T')[0])
        plt.xticks(xticksforlabels,xticklabels,rotation=20)
        if GOODMAPS:
            plt.savefig('pointsource_results/light_curves/'+region+'/'+'index_'+str(source_dict[eachsource]['index']).zfill(4)+'_'+eachsource+'_'+wave+'_Wcal_GoodMaps_lightcurve.pdf',format='pdf')
        else:
            plt.savefig('pointsource_results/light_curves/'+region+'/'+'index_'+str(source_dict[eachsource]['index']).zfill(4)+'_'+eachsource+'_'+wave+'_Wcal_lightcurve.pdf',format='pdf')
        plt.clf()
 
