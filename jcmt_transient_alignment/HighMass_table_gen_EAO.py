import pickle
import numpy as np
import os

def make_table(region='DR21C',alignment_iteration=0,wavelength='850'):
    """
    :param regions: a string representing the region to run
    :param alignment_iteration: there will be multiple iterations of the alignment run - this 0-based integer designates which alignment iteration the output table describes
    """
    previous_data = []
    if os.path.exists('tables/Transient_'+region+'_run_'+str(alignment_iteration)+'_'+wavelength+'.table'):
        with open('tables/Transient_'+region+'_run_'+str(alignment_iteration)+'_'+wavelength+'.table') as previoustable:
            content = previoustable.readlines()
        content = [x.strip() for x in content]
        for i in content:
            thisline = []
            if i[0] == '#':
                continue
            for j in i.split():
                thisline.append(j)
            previous_data.append(thisline)
        previous_data = np.array(previous_data)

    with open("data/data_Transient_"+region+"_run_"+str(alignment_iteration)+"_"+wavelength+".pickle", 'rb') as data:
        data = pickle.load(data)
    data=data[region]
    px2fwhm = 2.355 / np.sqrt(2)


    if os.path.exists('tables/Transient_'+region+'_run_'+str(alignment_iteration-1)+'_'+wavelength+'.table'):
        num_new_files = len(np.array(list(data[wavelength]['header']['airmass'].keys())))
        indices_to_grab_from_last_table = np.arange(-1*num_new_files+1,0,1) # Uf we are running one new file, num_new_files will be 2: One for the new file and one for the reference
        oldtable = np.genfromtxt('tables/Transient_'+region+'_run_'+str(alignment_iteration-1)+'_'+wavelength+'.table',dtype=None,names=True)
        lastdx = [oldtable['dx'][0]]
        lastdy = [oldtable['dy'][0]]
        for eachind in indices_to_grab_from_last_table:
            lastdx.append(oldtable['dx'][eachind])
            lastdy.append(oldtable['dy'][eachind])
        lastdx = np.array(lastdx)
        lastdy = np.array(lastdy)
    else:
        lastdx = 0.0
        lastdy = 0.0

    keyarray = np.array(list(data[wavelength]['header']['airmass'].keys()))
    
    JD = np.array(list(data[wavelength]['header']['julian_date'].values()),dtype=str)[np.argsort(keyarray)].T
    #key = np.array(list(data[wavelength]['header']['airmass'].keys()), dtype=str)[np.argsort(keyarray)].T
    airmass = np.array(list(data[wavelength]['header']['airmass'].values()), dtype=float)[np.argsort(keyarray)].T
    tau225 = np.array(list(data[wavelength]['header']['t225'].values()), dtype=float)[np.argsort(keyarray)].T
    
    elev = np.array(list(data[wavelength]['header']['elevation'].values()), dtype=float)[np.argsort(keyarray)].T
   
    if wavelength  == '850': 
        M_850 = np.array(list(data['850']['linear']['m'].values()), dtype=float)[np.argsort(keyarray)].T
        M_err_850 = np.array(list(data['850']['linear']['m_err'].values()), dtype=float)[np.argsort(keyarray)].T
    
        B_850 = np.array(list(data['850']['linear']['b'].values()), dtype=float)[np.argsort(keyarray)].T
        B_err_850 = np.array(list(data['850']['linear']['b_err'].values()), dtype=float)[np.argsort(keyarray)].T
    
        AMP_850 = np.array(list(data['850']['AC']['amp'].values()), dtype=float)[np.argsort(keyarray)].T
        AMP_err_850 = np.array(list(data['850']['AC']['amp_err'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGX_850 = np.array(list(data['850']['AC']['sig_x'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGX_fwhm_850 = SIGX_850 * px2fwhm
        SIGX_err_850 = np.array(list(data['850']['AC']['sig_x_err'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGY_850 = np.array(list(data['850']['AC']['sig_y'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGY_fwhm_850 = SIGY_850 * px2fwhm
        SIGY_err_850 = np.array(list(data['850']['AC']['sig_y_err'].values()), dtype=float)[np.argsort(keyarray)].T
        THETA_850 = np.array(list(data['850']['AC']['theta'].values()), dtype=float)[np.argsort(keyarray)].T
        THETA_err_850 = np.array(list(data['850']['AC']['theta_err'].values()), dtype=float)[np.argsort(keyarray)].T
        #dx = np.array(list(data['850']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[0] * 3
 
        # 2020-01-05 -- Flip the X axis based on tests looking at the pointing files and the generated images
        # I really thought it was the Y-axis that should be reversed - but tests indicate it is the X-axis -- opposite to common sense
        dx = lastdx+np.array(list(data['850']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[0] * -3
        dy = lastdy+np.array(list(data['850']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[1] * 3
        key = np.array(sorted(list(data[wavelength]['header']['airmass'].keys())), dtype=str).T

        cal_measure_850 = np.sqrt(-M_850)
        #cal_frac_linear_850 = cal_measure_850/cal_measure_850.mean()
        # 2021-04-07 Change normalisation from mean to first epoch! This is essential so these numbers stay the same as we get more data!
        cal_frac_linear_850 = cal_measure_850/cal_measure_850[0]
    
        beam_linear_850 = np.sqrt(-B_850/M_850)
        beam_linear_fwhm_850 = beam_linear_850 * px2fwhm
        beam_gaussian_850 = np.sqrt(SIGY_850*SIGX_850)
        beam_gaussian_fwhm_850 = beam_gaussian_850 * px2fwhm
        #cal_frac_gauss_850 = np.sqrt(AMP_850/AMP_850.mean())/(beam_gaussian_850/beam_gaussian_850.mean())
        # 2021-04-07 Change normalisation from mean to first epoch! This is essential so these numbers stay the same as we get more data!
        cal_frac_gauss_850 = np.sqrt(AMP_850/AMP_850[0])/(beam_gaussian_850/beam_gaussian_850[0])

        
        hdr = "Key JulianDate ELEV Airmass Tau225 " \
              "CalF_lin CalF_gauss Beam_lin Beam_gauss " \
              "Beam_lin_fwhm Beam_gauss_fwhm Sigx_fwhm Sigy_fwhm " \
              "dx dy " \
              "B B_err M M_err " \
              "Amp Amp_err Sigx Sigx_err Sigy Sigy_err Theta Theta_err "
    
        arr = [key, JD, elev, airmass, tau225,
               cal_frac_linear_850, cal_frac_gauss_850, beam_linear_850, beam_gaussian_850,
               beam_linear_fwhm_850, beam_gaussian_fwhm_850, SIGX_fwhm_850, SIGY_fwhm_850,
               dx, dy,
               B_850, B_err_850, M_850, M_err_850,
               AMP_850, AMP_err_850, SIGX_850, SIGX_err_850, SIGY_850, SIGY_err_850, THETA_850, THETA_err_850,
               ]

    elif wavelength == '450':
        M_450 = np.array(list(data['450']['linear']['m'].values()), dtype=float)[np.argsort(keyarray)].T
        M_err_450 = np.array(list(data['450']['linear']['m_err'].values()), dtype=float)[np.argsort(keyarray)].T
    
        B_450 = np.array(list(data['450']['linear']['b'].values()), dtype=float)[np.argsort(keyarray)].T
        B_err_450 = np.array(list(data['450']['linear']['b_err'].values()), dtype=float)[np.argsort(keyarray)].T
    
        AMP_450 = np.array(list(data['450']['AC']['amp'].values()), dtype=float)[np.argsort(keyarray)].T
        AMP_err_450 = np.array(list(data['450']['AC']['amp_err'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGX_450 = np.array(list(data['450']['AC']['sig_x'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGX_fwhm_450 = SIGX_450 * px2fwhm
        SIGX_err_450 = np.array(list(data['450']['AC']['sig_x_err'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGY_450 = np.array(list(data['450']['AC']['sig_y'].values()), dtype=float)[np.argsort(keyarray)].T
        SIGY_fwhm_450 = SIGY_450 * px2fwhm
        SIGY_err_450 = np.array(list(data['450']['AC']['sig_y_err'].values()), dtype=float)[np.argsort(keyarray)].T
        THETA_450 = np.array(list(data['450']['AC']['theta'].values()), dtype=float)[np.argsort(keyarray)].T
        THETA_err_450 = np.array(list(data['450']['AC']['theta_err'].values()), dtype=float)[np.argsort(keyarray)].T
        #dx450 = np.array(list(data['450']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[0] * 2
        # 2020-01-05 -- Flip the X axis based on tests looking at the pointing files and the generated images
        # I really thought it was the Y-axis that should be reversed - but tests indicate it is the X-axis -- opposite to common sense
        dx450 = lastdx+np.array(list(data['450']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[0] * -2
        dy450 = lastdy+np.array(list(data['450']['XC']['alignment'].values()), dtype=float)[np.argsort(keyarray)].T[1] * 2
        key = np.array(sorted(list(data[wavelength]['header']['airmass'].keys())), dtype=str).T

        M_450[M_450 > 0] = -0.0001
    
        cal_measure_450 = np.sqrt(-M_450)
        #cal_frac_linear_450 = cal_measure_450/cal_measure_450.mean()
        # 2021-04-07 Change normalisation from mean to first epoch! This is essential so these numbers stay the same as we get more data!
        cal_frac_linear_450 = cal_measure_450/cal_measure_450[0]
    
        beam_linear_450 = np.sqrt(-B_450/M_450) / np.sqrt(2)
        beam_linear_fwhm_450 = beam_linear_450 * px2fwhm
        beam_gaussian_450 = np.sqrt(SIGY_450*SIGX_450)
        beam_gaussian_fwhm_450 = beam_gaussian_450 * px2fwhm
        #cal_frac_gauss_450 = np.sqrt(AMP_450/AMP_450.mean())/(beam_gaussian_450/beam_gaussian_450.mean())
        # 2021-04-07 Change normalisation from mean to first epoch! This is essential so these numbers stay the same as we get more data!
        cal_frac_gauss_450 = np.sqrt(AMP_450/AMP_450[0])/(beam_gaussian_450/beam_gaussian_450[0])

        hdr = "Key JulianDate ELEV Airmass Tau225 " \
              "CalF_lin CalF_gauss Beam_lin Beam_gauss " \
              "Beam_lin_fwhm Beam_gauss_fwhm Sigx_fwhm Sigy_fwhm " \
              "dx dy " \
              "B B_err M M_err " \
              "Amp Amp_err Sigx Sigx_err Sigy Sigy_err Theta Theta_err "

        arr = [key, JD, elev, airmass, tau225,
               cal_frac_linear_450, cal_frac_gauss_450, beam_linear_450, beam_gaussian_450,
               beam_linear_fwhm_450, beam_gaussian_fwhm_450, SIGX_fwhm_450, SIGY_fwhm_450,
               dx450, dy450,
               B_450, B_err_450, M_450, M_err_450,
               AMP_450, AMP_err_450, SIGX_450, SIGX_err_450, SIGY_450, SIGY_err_450, THETA_450, THETA_err_450,
               ]

    if not os.path.exists('tables'):
        os.system('mkdir tables')

    if len(previous_data)>0:
        combinedsizelist = []
        for eachelement in previous_data.T[0]:
            combinedsizelist.append(eachelement)
        for eachelement in arr[0][1:]: #Make sure to skip first date since it will already be in table!
            combinedsizelist.append(eachelement)
        table = np.array(np.zeros(len(combinedsizelist)))
        for i,j in zip(previous_data.T,arr): 
            combined_data = []
            for eachelementini in i:
                combined_data.append(eachelementini)
            for eachelementinj in j[1:]: #Make sure to skip first date since it will already be in table!
                combined_data.append(eachelementinj)
            table = np.vstack((table,np.array(combined_data,dtype=str)))
    else:
        table = np.array(np.zeros(len(key)))
        for val in arr:
            table = np.vstack((np.array(table), np.array(val, dtype=str)))
    np.savetxt("tables/Transient_"+region+"_run_"+str(alignment_iteration)+"_"+wavelength+".table", np.array(table)[1:].T, fmt="%s", header=hdr)
