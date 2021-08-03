import numpy as np
import glob
import pickle
import os

def generate_calfactors_file(datadir,region,wave):

    # Check for a previous Calfactors file!
    previous_file = len(sorted(list(glob.glob(datadir+'/'+region+'_*'+wave+'_CalFactors.txt'))))
    if previous_file>0:
        previous_filename = sorted(list(glob.glob(datadir+'/'+region+'_*'+wave+'_CalFactors.txt')))[-1]
    
    # In the following line we should ALWAYS USE 850 -- since the 450 files were generated using 850 pointing corrections
    pointing_file = np.genfromtxt(sorted(glob.glob('tables/Transient_'+region+'_run_*_850.table'))[-2],names=True,dtype=None) # the second-to-last file, since the last file is the final check, not the pointing corrections we used
    # In the next lines, we need to use whatever wavelength we are interested in
    if wave == '450':
        # only 1 is generated for 450 -- and this is the ACcal correction we use
        AC_file = np.genfromtxt(sorted(glob.glob('tables/Transient_'+region+'_run_*_'+wave+'.table'))[-1],names=True,dtype=None)
        previous_AC_files = sorted(glob.glob('tables/Transient_'+region+'_run_*_'+wave+'.table'))[0:-1] #Previous Pointing corrections already taken into consideration in tables!
    if wave == '850':
        # again, use the second-to-last file, since the last file is the final check, not the ACcal corrections we used
        AC_file = np.genfromtxt(sorted(glob.glob('tables/Transient_'+region+'_run_*_'+wave+'.table'))[-2],names=True,dtype=None)
        previous_AC_files = sorted(glob.glob('tables/Transient_'+region+'_run_*_'+wave+'.table'))[0:-2] #Previous Pointing corrections already taken into consideration in tables!


    datescans_pointing = pointing_file['Key']
    datescans_ACcal    = AC_file['Key']
    
    datescans_pointing_converted = []
    datescans_ACcal_converted = []
    for i in datescans_pointing:
         scan = i.decode('utf-8').split('-')[-1].zfill(5)
         datescans_pointing_converted.append(i.decode('utf-8').split('-')[0]+'_'+scan)
    for i in datescans_ACcal:
         scan = i.decode('utf-8').split('-')[-1].zfill(5)
         datescans_ACcal_converted.append(i.decode('utf-8').split('-')[0]+'_'+scan)


    datescans_pointing_converted = np.array(datescans_pointing_converted)
    datescans_ACcal_converted = np.array(datescans_ACcal_converted)
    dx = pointing_file['dx']
    dy = pointing_file['dy'] 
    if len(previous_AC_files) == 0: #Previous Pointing corrections already taken into consideration in tables!
        ACcals = AC_file['CalF_gauss']
    else:
        prev_ACcals_combined = []
        for eachdate in datescans_ACcal_converted:
            prev_ACcal_by_date = []
            for eachprevAC in previous_AC_files:
                prevAC = np.genfromtxt(eachprevAC,names=True,dtype=None)
                prevACcals_alldates = np.array(prevAC['CalF_gauss'])
                prevdates = prevAC['Key']
                prevdates_converted = []
                for i in prevdates:
                    prevdates_converted.append(i.decode('utf-8').split('-')[0]+'_'+i.decode('utf-8').split('-')[-1].zfill(5))
                thisdate_ind = np.where(np.array(prevdates_converted)==eachdate)
                prev_ACcal_by_date.append(prevACcals_alldates[thisdate_ind])
            most_recent_ACcal = AC_file['CalF_gauss'][np.where(datescans_ACcal_converted==eachdate)]
            prev_ACcal_by_date.append(most_recent_ACcal)
            prev_ACcals_combined.append(np.prod(prev_ACcal_by_date))
        ACcals = np.array(prev_ACcals_combined)

    weighted_file           = np.genfromtxt(sorted(glob.glob('pointsource_results/'+region+'/*'+wave+'*_calepoch_weightedmean.txt'))[-1],dtype=None,names=True)
    weighted_datescans1     = np.array(weighted_file['DateScan'])
    weighted_datescans      = []
    for i in weighted_datescans1:
        weighted_datescans.append(str(i)[0:8]+'_'+str(i)[8:])
    weighted_calfactors     = np.array(weighted_file['Divisor'])[np.argsort(weighted_datescans)]
    weighted_calfactoruncs  = np.array(weighted_file['WeightedCalUnc'])[np.argsort(weighted_datescans)]
    weighted_datescans      = np.array(weighted_datescans)[np.argsort(np.array(weighted_datescans))]
    
    try:
        PS_file          = pickle.load(open('pointsource_results/'+region+'/'+region+'_PointSource_cal_info_dict_targunc5_'+wave+'.pickle','rb'))
        PS_datescans     = np.array(PS_file['datescans'])
        PS_calfactors    = np.array(PS_file['RelFCFs'])[np.argsort(PS_datescans)]
        PS_calfactoruncs = np.array(PS_file['RelFCFuncs'])[np.argsort(PS_datescans)]
        PS_datescans     = np.array(PS_file['datescans'])[np.argsort(PS_datescans)]
    except:
        PS_datescans     = np.zeros(len(weighted_datescans))*np.nan
        PS_calfactors    = np.zeros(len(weighted_datescans))*np.nan
        PS_calfactoruncs = np.zeros(len(weighted_datescans))*np.nan
        PS_datescans     = np.zeros(len(weighted_datescans))*np.nan


    #Always 850, below -- no pointing check done at 450 becasue it will be the same!
    pointing_check_file = np.genfromtxt('pointsource_results/'+region+'/'+region+'_850_local_pointing_corrections.txt',dtype=None,names=True)
    pointingcheckdatescans = []
    for i in pointing_check_file['Datescan']:
        pointingcheckdatescans.append(str(i)[0:8]+'_'+str(i)[8:])
    pointingcheckdatescans = np.array(pointingcheckdatescans)
    pointingcheckxoff = pointing_check_file['Xoff'][np.argsort(pointingcheckdatescans)]
    pointingcheckyoff = pointing_check_file['Yoff'][np.argsort(pointingcheckdatescans)]
    pointingcheckdatescans = pointingcheckdatescans[np.argsort(pointingcheckdatescans)]

    output = open('pointsource_results/'+region+'/'+region+'_'+str(datescans_pointing_converted[-1])+'_'+wave+'_CalFactors.txt','w') 
    # 2021-05-14 --- Updated to reflect that now the Wcal and PScal are done on the NON ACcal data. The ACcal and PScal are treated as "checks" to the Wcal
    #output.write('#Date_Scan\tdx\tdy\tdx_check\tdy_check\tACcal\tWeightedCal\tWeightedCal_unc\tPScal/WeightedCal\tPScal_unc\n')
    output.write('#Date_Scan\tdx\tdy\tdx_check\tdy_check\tWeightedCal\tWeightedCal_unc\tACcal/WeightedCal\tPScal/WeightedCal\tPScal_unc\n')

    ##if previous_file>0:
    ##    prev_file = open(previous_filename,'r')
    ##    header = prev_file.readline()
    ##    line = prev_file.readline()
    ##    while line:
    ##        output.write(line)
    ##        line = prev_file.readline()

    for ds_point,dxi,dyi in zip(datescans_pointing_converted,dx,dy):
        output_ds = ds_point
        output_dx = dxi
        output_dy = dyi

        Weighted_index = np.where(weighted_datescans == ds_point)
        if len(weighted_calfactors[Weighted_index])<1:
            output_weighted = 'nan'
            output_weightedunc = 'nan'
        else:
            output_weighted = round(weighted_calfactors[Weighted_index][0],3)
            output_weightedunc = round(weighted_calfactoruncs[Weighted_index][0],4)

        ACcal_index = np.where(datescans_ACcal_converted==ds_point)
        if len(ACcals[ACcal_index])<1:
            output_ac = 'nan'
        else:
            output_ac = round(ACcals[ACcal_index][0]/weighted_calfactors[Weighted_index][0],3)

        PScal_index = np.where(PS_datescans==ds_point)
        if len(PS_calfactors[PScal_index])<1:
            output_ps = 'nan'
            output_psunc = 'nan'
        else:
            output_ps = round(PS_calfactors[PScal_index][0]/weighted_calfactors[Weighted_index][0],3)
            output_psunc = round(PS_calfactoruncs[PScal_index][0],4)

        PointingCheck_index = np.where(pointingcheckdatescans==ds_point)
        if len(pointingcheckxoff[PointingCheck_index])<1:
            dx_check = 'nan'
            dy_check = 'nan'
        else:
            dx_check = str(round(pointingcheckxoff[PointingCheck_index][0],2))
            dy_check = str(round(pointingcheckyoff[PointingCheck_index][0],2))

        output.write(output_ds+'\t'+str(round(output_dx,2))+'\t'+str(round(output_dy,2))+'\t'+dx_check+'\t'+dy_check+'\t'+str(output_weighted)+'\t'+str(output_weightedunc)+'\t'+str(output_ac)+'\t'+str(output_ps)+'\t'+str(output_psunc)+'\n')
    output.close()

