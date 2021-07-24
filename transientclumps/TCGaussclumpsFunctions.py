import subprocess
import math
import numpy as np
import astropy.io.fits as apfits

### get_noise ###
# Takes an image in .sdf format and calculates the noise in the central region
# of the map.
#
# Keywords:
# img_name - name of the target image in .sdf format (string) -- required
# Output:
# noise - the noise in the central region of the map (float).

def get_noise(img_name):

    #Convert the .sdf to .fits.
    fits_name = img_name[:-4]+'.fits'
    sdf_conv_command = '$CONVERT_DIR/ndf2fits '+img_name+' !'+fits_name
    subprocess.call(sdf_conv_command, shell=True)

    #Read in the variance extension of the file.
    vardata = apfits.getdata(fits_name, 1)

    #Find the noise in the central portion of the map.
    vardata_rad = vardata.shape[1]/2
    #Pick a region in the middle of the map to find the median noise. The X and
    # Y ranges will be the same.
    vardata_lo = int((vardata_rad-(vardata_rad/2)))
    vardata_hi = int((vardata_rad+(vardata_rad/2)))
    vardata_mid = vardata[0, vardata_lo:vardata_hi, vardata_lo:vardata_hi]
    #2021-05-26 -- added the factor of 1000, since the variance array did not seem to update when we converted to Jy/beam
    noise = 0.01 # 0.01 Jy/beam = 10 mJy/beam
    #print('\n\n\n',math.sqrt(np.median(vardata_mid)),'\n\n\n')
    #if math.sqrt(np.median(vardata_mid))>1:
    #    noise = math.sqrt(np.median(vardata_mid))/1000
    #else:
    #    noise = math.sqrt(np.median(vardata_mid))*537 # This puts it in Jy/beam!!
    return noise

#def


### run_gaussclumps ###
# Takes an image in .sdf format and runs gaussclumps.
#
# Keywords:
# img_name - name of the target image in .sdf format (string). -- required
# param_file - name of the gaussclumps parameter file (string). -- required

def run_gaussclumps(img_name, param_file):

    #Default outputs are "img_name"+suffix.extension
    outputfile = img_name[:-4]+'_clumps.sdf'
    outcatfile = img_name[:-4]+'_log.FIT'
    logfile = img_name[:-4]+'_log.txt'
    method = 'GaussClumps'
    noise_str = str( ('{:.5e}').format( get_noise(img_name) ) )

    #Write out the .sh script that will be executed.
    shellscript = open('./findclumpscript.sh', 'w')
    shellscript.write('#!/bin/bash\n')
    shellscript.write('{ . $CUPID_DIR/cupid.sh; }\n')
    shellscript.write("findclumps '"
    +img_name+"' '"+outputfile+"' '"+outcatfile+"' '"+method+"' CONFIG=^'"
    +param_file+"' LOGFILE='"+logfile+"' REPCONF=TRUE RMS="+noise_str+
    " DECONV=FALSE MSG_FILTER=3 SHAPE=Ellipse WCSPAR=TRUE BACKOFF=TRUE")
    shellscript.close()
    subprocess.call('source ./findclumpscript.sh', shell='TRUE')

#def
