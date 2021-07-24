import numpy as np
from astropy.io import fits
import json
import glob
from starlink import kappa

def GCfluxcal(direc,wave):

    sdf_files = sorted(glob.glob(direc+'/*'+wave+'*mJybmsm.sdf'))

    region = sdf_files[0].split('_')[0]
    jsonregion = region
    if region == 'OPHCORE':
        jsonregion = 'OPH_CORE'
    if region == 'SERPENSSOUTH':
        jsonregion = 'SERPENS_SOUTH'
    if region == 'SERPENSMAIN':
        jsonregion = 'SERPENS_MAIN'

    if wave == '450':
        print('Skipping GC Catalogue Fluxcal since that is only appropriate for 850 microns')
    elif wave == '850':
        date_cutoff = '20170301'

    #Get PIDENT refs of known calibrators -- needs to be updated for high mass regions
    family_data = json.load(open('config/family_{}.json'.format(wave)))
    good_sources = family_data[jsonregion]['A3']

    calinfo = {}
    for eachfile in sdf_files:
        datescan = eachfile.split('_')[-6]+'_'+eachfile.split('_')[-5]
        calinfo[datescan] = {}
        calinfo[datescan]['PIDENT_ref'] = []
        calinfo[datescan]['Flux'] = []

        culled_cat = eachfile.split('_mJybmsm.sdf')[0]+'_Jypbmsm_crop'+'_log_cull.FIT'
        cullcat = fits.getdata(culled_cat)
        for eachgoodsource in good_sources:
            catalogue_ind = np.where(cullcat['PIDENT_ref'] == eachgoodsource)
            calinfo[datescan]['PIDENT_ref'].append(eachgoodsource)
            calinfo[datescan]['Flux'].append(cullcat['Peak_'+datescan.split('_')[0]][catalogue_ind])

    cal_norm_fluxes_by_date = []
    for eachdatescan in sorted(list(calinfo.keys())):
        if int(eachdatescan.split('_')[0])<int(date_cutoff):
            cal_norm_fluxes_by_date.append(calinfo[eachdatescan]['Flux'])

    cal_norm_fluxes = np.array(cal_norm_fluxes_by_date).mean(axis=0)
    
    for eachfile in sdf_files:
        datescan = eachfile.split('_')[-6]+'_'+eachfile.split('_')[-5]
        calfactor = np.mean(calinfo[datescan]['Flux']/cal_norm_fluxes)
        calfactor_unc = np.std(calinfo[datescan]['Flux']/cal_norm_fluxes,ddof=1)/np.sqrt(len(cal_norm_fluxes))
        kappa.cdiv(eachfile,calfactor,eachfile.split('.sdf')[0]+'_GCcal.sdf') 

