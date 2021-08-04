import os
import numpy as np
import glob

def make_pcor(tablefile,reduce_firstepoch=False):

    table = np.genfromtxt(tablefile,dtype=None,names=True)

    if reduce_firstepoch:
        for eachdate,eachdx,eachdy in zip(table['Key'],table['dx'],table['dy']):
            eachdatestr = eachdate.decode('utf-8')
            if len(glob.glob(tablefile.split('_')[1]+'/*'+eachdatestr.split('-')[0]+'_'+eachdatestr.split('-')[1].zfill(5)+'*'))>0:
                pcorfile = open(tablefile.split('_')[1]+'_'+eachdatestr.split('-')[0]+'_'+eachdatestr.split('-')[1].zfill(5)+'_pcor.txt','w')
                pcorfile.write('#SYSTEM=TRACKING\n')
                pcorfile.write('#TAI DLON DLAT\n')
                pcorfile.write('1 '+str(eachdx)+' '+str(eachdy)+'\n')
                pcorfile.write('10000000 '+str(eachdx)+' '+str(eachdy))
                pcorfile.close()
                print('POINTING FOR '+eachdatestr+': '+str(eachdx)+', '+str(eachdy))

    else:
        for eachdate,eachdx,eachdy in zip(table['Key'],table['dx'],table['dy']):
            if not np.logical_and(str(eachdx)==str(0.0),str(eachdy)==str(0.0)):
                if not np.logical_and(str(eachdx)==str(-0.0),str(eachdy)==str(0.0)):
                    eachdatestr = eachdate.decode('utf-8')
                    if len(glob.glob(tablefile.split('_')[1]+'/*'+eachdatestr.split('-')[0]+'_'+eachdatestr.split('-')[1].zfill(5)+'*'))>0:
                        pcorfile = open(tablefile.split('_')[1]+'_'+eachdatestr.split('-')[0]+'_'+eachdatestr.split('-')[1].zfill(5)+'_pcor.txt','w')
                        pcorfile.write('#SYSTEM=TRACKING\n')
                        pcorfile.write('#TAI DLON DLAT\n')
                        pcorfile.write('1 '+str(eachdx)+' '+str(eachdy)+'\n')
                        pcorfile.write('10000000 '+str(eachdx)+' '+str(eachdy))
                        pcorfile.close()
                        print('POINTING FOR '+eachdate.decode('utf-8')+': '+str(eachdx)+', '+str(eachdy))
    if not os.path.exists('pointing_files_posterity'):
        os.system('mkdir pointing_files_posterity')
    os.system('cp *_pcor.txt pointing_files_posterity')

def makemap_infiles(region,wave):

    for eachpcorfile in sorted(glob.glob('*pcor.txt')):
        date = eachpcorfile.split('_')[-3]
        scan = eachpcorfile.split('_')[-2]
        if len(glob.glob(region+'/*'+date+'_'+scan+'*'))>0:
            os.system('ls -1d /jcmtdata/raw/scuba2/s'+wave[0]+'*/'+date+'/'+scan+'/*sdf > '+region+'_'+date+'_'+scan+'_'+wave+'.txt')

def create_makemap_script(wave):
    with open('makemaps.sh','w') as makemapsscript:
        makemapsscript.write('#!/bin/bash\n')
        makemapsscript.write('# Run makemap with pointing corrections - calibrate them and crop them\n')
        makemapsscript.write('# User-defined variables on which maps to reduce\n')
        makemapsscript.write('WAVES='+wave+'\n')
        makemapsscript.write('STARLINK_DIR=/stardev\n')
        makemapsscript.write('ORAC_DIR=/stardev/bin/oracdr/src\n')
        makemapsscript.write('for WAVE in $WAVES\n')
        makemapsscript.write('  do\n')
        makemapsscript.write('    INFILES=*$WAVE.txt\n')
        makemapsscript.write('    # Select proper pixel size\n')
        makemapsscript.write('    if [[ "$WAVE" == 850 ]]; then\n')
        makemapsscript.write('      pixsize=3\n')
        makemapsscript.write('    elif [[ "$WAVE" == 450 ]]; then\n')
        makemapsscript.write('      pixsize=2\n')
        makemapsscript.write('    else\n')
        makemapsscript.write('    echo Pick a wavelength of either 850 or 450\n')
        makemapsscript.write('    fi\n')
        makemapsscript.write('    # Run a loop to reduce R1, R2, R3, and R4 maps\n')
        makemapsscript.write('    for f in $INFILES\n')
        makemapsscript.write('      do\n')
        makemapsscript.write('        if [[ $f =~ (.*)_(.*)_(.*)_(.*)\.txt$ ]]; then\n')
        makemapsscript.write('          region=${BASH_REMATCH[1]}\n')
        makemapsscript.write('          thisdate=${BASH_REMATCH[2]}\n')
        makemapsscript.write('          scan=${BASH_REMATCH[3]}\n')
        makemapsscript.write('          thiswave=${BASH_REMATCH[4]}\n')
        makemapsscript.write('        fi\n')
        makemapsscript.write('    ##\n')
        makemapsscript.write('    # Run Makemap With The Pointing Correction File\n')
        makemapsscript.write('    ##\n')
        makemapsscript.write('        cmdR1map="$STARLINK_DIR/bin/smurf/makemap in=^${region}_${thisdate}_${scan}_${thiswave}.txt out=${region}_${thisdate}_${scan}_${thiswave}_CR3.sdf ref=config/${region}_R3_extmask_${thiswave}.sdf config=^config/dimmconfig_R3.lis pixsize=$pixsize pointing=${region}_${thisdate}_${scan}_pcor.txt"\n')
        makemapsscript.write('        echo $cmdR1map\n')
        makemapsscript.write('        $cmdR1map\n')
        makemapsscript.write('    ##\n')
        makemapsscript.write('    # Run Picard to calibrate the data\n')
        makemapsscript.write('    ##\n')
        makemapsscript.write('        if [[ "$thiswave" == 850 ]]; then\n')
        makemapsscript.write('          cmdR1cal="${ORAC_DIR}/etc/picard_start.sh -log sf -nodisplay CALIBRATE_SCUBA2_DATA -recpars=\'FCF_CALTYPE=ARCSEC,FCF=2.34\' ${region}_${thisdate}_${scan}_850_CR3.sdf" \n')
        makemapsscript.write('        else\n')
        makemapsscript.write('          cmdR1cal="${ORAC_DIR}/etc/picard_start.sh -log sf -nodisplay CALIBRATE_SCUBA2_DATA -recpars=\'FCF_CALTYPE=ARCSEC,FCF=4.71\' ${region}_${thisdate}_${scan}_450_CR3.sdf"\n')
        makemapsscript.write('        fi\n')
        makemapsscript.write('        echo $cmdR1cal\n')
        makemapsscript.write('        $cmdR1cal\n')
        makemapsscript.write('    ##\n')
        #makemapsscript.write('    # Crop the calibrated maps\n')
        #makemapsscript.write('    ##\n')
        #makemapsscript.write('        cmdR1crop="${ORAC_DIR}/etc/picard_start.sh -nodisplay -log sf -recpars crop.ini CROP_SCUBA2_IMAGES ${region}_${thisdate}_${scan}_${thiswave}_CR3_cal.sdf"\n')
        #makemapsscript.write('        echo $cmdR1crop\n')
        #makemapsscript.write('        $cmdR1crop\n')
        makemapsscript.write('    ##\n')
        makemapsscript.write('    # Rename\n')
        makemapsscript.write('    ##\n')
        #makemapsscript.write('        cmdR1mv="mv ${region}_${thisdate}_${scan}_${thiswave}_CR3_cal_crop.sdf ${region}_${thisdate}_${scan}_${thiswave}_CR3_crop.sdf"\n')
        makemapsscript.write('        cmdR1del="rm -f ${region}_${thisdate}_${scan}_${thiswave}_CR3.sdf"\n')
        makemapsscript.write('        cmdR1mv="mv ${region}_${thisdate}_${scan}_${thiswave}_CR3_cal.sdf ${region}_${thisdate}_${scan}_${thiswave}_CR3.sdf"\n')
        makemapsscript.write('        echo $cmdR1mv\n')
        makemapsscript.write('        $cmdR1del\n')
        makemapsscript.write('        $cmdR1mv\n')
        makemapsscript.write('  done\n')
        makemapsscript.write('done')
        makemapsscript.close()
