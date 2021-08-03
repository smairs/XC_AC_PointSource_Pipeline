from starlink import kappa
import glob

def smooth_input_data(datadir,wave):
    inputfiles = sorted(glob.glob(datadir+'/*'+wave+'*CR3.sdf')) # Already in mJy/beam
    for eachfile in inputfiles:
        kappa.gausmooth(eachfile,eachfile.split('.sdf')[0]+'_mJybmsm.sdf',2)
