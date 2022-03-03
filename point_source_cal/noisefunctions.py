import subprocess
import math
import numpy as np
import astropy.io.fits as apfits
import os
import glob

def get_noise(fname,wave):
    '''
    This Function Determines the Noise of a given map
    '''
    ### get_noise ###

    #2021-05-26 --- The variance is not multiplied by the new mJy/beam factor for some reason
    if wave == '450':
        mJybeamfactor = 491000 
    elif wave == '850':
        mJybeamfactor = 537000

    #Convert the .sdf to .fits.
    #print('Converting '+fname+' to fits file: '+fname[:-4]+'.fits')
    fits_name = fname[:-4]+'.fits'
    sdf_conv_command = '$CONVERT_DIR/ndf2fits '+fname+' !'+fits_name
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
    #Fudge factor - some variance arrays don't update for some reason when we perform mJy/beam conversion - so the noise measurment is off by that factor!
    noise = mJybeamfactor*math.sqrt(np.nanmedian(vardata_mid))
    if noise>10000:
        noise = math.sqrt(np.nanmedian(vardata_mid))
    os.system('rm '+fits_name)
    return noise


def make_noise_file(inputfiles,wave):
    
    noises  = []
    images  = []
    
    noisefile = open(wave+'_noises.txt','w')
    noisefile.write('#Image                      Noise\n')
    for eachimage in inputfiles:
        noise  = get_noise(eachimage,wave)
        noisefile.write(eachimage.split('/')[-1]+' '+str(round(noise,4))+'\n')
    noisefile.close() 


def readnoise(filename,region):
    catalogue = np.genfromtxt(filename,dtype=None,names=True)
    date_scan = []
    noise     = []
    try:
        for eachfile,eachnoise in zip(catalogue['Image'],catalogue['Noise']):
            eachfile = eachfile.decode('utf-8')
            if eachfile.split('_')[0] == region:
                date_scan.append(eachfile.split('_')[-6]+'_'+eachfile.split('_')[-5])
                noise.append(eachnoise)
    except TypeError:
        try:
            for eachfile,eachnoise in zip([catalogue['Image']],[catalogue['Noise']]):
                eachfile = eachfile.decode('utf-8')
                if eachfile.split('_')[0] == region:
                    date_scan.append(eachfile.split('_')[-6]+'_'+eachfile.split('_')[-5])
                    noise.append(eachnoise)
        except IndexError:
           for eachfile,eachnoise in zip([catalogue['Image'].item()],[catalogue['Noise'].item()]):
                eachfile = eachfile.decode('utf-8')
                if eachfile.split('_')[0] == region:
                    date_scan.append(eachfile.split('_')[-6]+'_'+eachfile.split('_')[-5])
                    noise.append(eachnoise)

    return(np.array(date_scan),np.array(noise))
