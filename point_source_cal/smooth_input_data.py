from starlink import kappa
import glob
import os

def smooth_input_data(datadir,wave):
    if wave == '850':
        beamsize    = 14.5 # arcsec - 14.5 is what was assumed in A3 in 2015 and used throughout
        smoothsize  = 6  # arcsec
    if wave == '450':
        beamsize    = 9.8 # arcsec
        smoothsize  = 4  # arcsec
    beam_factor = (beamsize**2+smoothsize**2)/(beamsize**2)
    inputfiles = sorted(glob.glob(datadir+'/*'+wave+'*CR3.sdf')) # Already in mJy/beam
    for eachfile in inputfiles:
        kappa.gausmooth(eachfile,eachfile.split('.sdf')[0]+'_mJybmsm.sdf',2)
        kappa.cmult(eachfile.split('.sdf')[0]+'_mJybmsm.sdf',beam_factor,eachfile.split('.sdf')[0]+'_mJybmsm_cor.sdf')
        os.system('mv -f '+eachfile.split('.sdf')[0]+'_mJybmsm_cor.sdf '+eachfile.split('.sdf')[0]+'_mJybmsm.sdf')
