import subprocess
import glob
import os

def changeunits(datadir,wave):
    # 2021-05-14 -- Changing from ACcal-based results to WrightedCal-based results
    #inputfiles = sorted(glob.glob(datadir+'/*'+wave+'*ACcal.sdf'))
    inputfiles = sorted(glob.glob(datadir+'/*'+wave+'*CR3.sdf'))
    for eachfile in inputfiles:
        subprocess.call('/stardev/bin/oracdr/src/etc/picard_start.sh -nodisplay -log sf UNCALIBRATE_SCUBA2_DATA '+eachfile,shell=True)
        if os.path.exists(eachfile.split('/')[-1].split('.sdf')[0]+'_uncal.sdf'):
            subprocess.call('/stardev/bin/oracdr/src/etc/picard_start.sh -nodisplay -log sf CALIBRATE_SCUBA2_DATA '+eachfile.split('/')[-1].split('.sdf')[0]+'_uncal.sdf',shell=True)
            os.system('mv '+eachfile.split('/')[-1].split('.sdf')[0]+'_uncal_cal.sdf '+eachfile.split('.sdf')[0]+'_mJybm.sdf')
            os.system('rm -f '+eachfile.split('/')[-1].split('.sdf')[0]+'_uncal.sdf')
        else:
            subprocess.call('/stardev/bin/oracdr/src/etc/picard_start.sh -nodisplay -log sf CALIBRATE_SCUBA2_DATA '+eachfile,shell=True)
            os.system('mv '+eachfile.split('/')[-1].split('.sdf')[0]+'_cal.sdf '+eachfile.split('.sdf')[0]+'_mJybm.sdf')
