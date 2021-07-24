import subprocess
import pdb
from starlink import kappa, picard
import os

##### Image Preparation Functions #####

### crop_image ###
# Takes an image in .sdf form and crops it to a given radius. The standard
# output is "img_name"_crop.sdf
#
# Keywords:
# img_name - name of the image (string). -- required
# crop_radius - radius, in arcseconds, of the crop (int).
# crop_method - method used for the crop (string).
# See SUN

def crop_image(img_name, crop_radius=1200, crop_method='CIRCLE'):

    #Make the cropping parameter file.
    crop_parms = open('crop.ini', 'w')
    crop_parms.write('[CROP_SCUBA2_IMAGES]\n')
    crop_parms.write('CROP_METHOD = '+str(crop_method)+'\n')
    crop_parms.write('MAP_RADIUS = '+str(crop_radius)+'\n')
    crop_parms.close()

    #Perform the cropping.
    crop_command = '${ORAC_DIR}/etc/picard_start.sh CROP_SCUBA2_IMAGES '
    crop_command += '-log f -recpars crop.ini '+img_name+' ${1+"$@"};'
    subprocess.call(crop_command, shell=True)
    print('\nCROP = DONE\n')

#def

#### smooth_image ###
## Takes an image in .sdf form and smooths it with a kernal in .sdf form.
## output is "img_name"_smooth.sdf
##
## Keywords:
## img_name - name of the image (string). -- required
## kern_name - name of the kernal .sdf to smooth with (string). -- required
#
#def smooth_image(img_name, kern_name):
#
#    #Determine the size of the kernal file.
#    naxis1_command = '$KAPPA_DIR/fitsval ndf='+kern_name+' keyword=NAXIS1'
#    naxis2_command = '$KAPPA_DIR/fitsval ndf='+kern_name+' keyword=NAXIS2'
#
#    naxis1_proc = subprocess.Popen(naxis1_command, shell=True,
#                                    stdout=subprocess.PIPE)
#    naxis1_stdout = naxis1_proc.communicate()
#    naxis1_size = int(str(naxis1_stdout[0])[2:-3])
#
#    naxis2_proc = subprocess.Popen(naxis2_command, shell=True,
#                                    stdout=subprocess.PIPE)
#    naxis2_stdout = naxis2_proc.communicate()
#    naxis2_size = int(str(naxis2_stdout[0])[2:-3])
#
#    #Assume that the center of the kernal is at the center of the .sdf
#    xcenter = round(naxis1_size/2)
#    ycenter = round(naxis2_size/2)
#
#    #Write the smoothing command.
#    smooth_command = '$KAPPA_DIR/convolve in='+img_name+' psf='+kern_name+' '
#    smooth_command += 'out='+img_name[:-4]+'_smooth.sdf '
#    smooth_command += 'xcentre='+str(xcenter)+' ycentre='+str(ycenter)
#    subprocess.call(smooth_command, shell=True)
#    print('\nSMOOTH = DONE\n')
#
##def

#### smooth_image_kappa ###
## Takes an image in .sdf form and smooths it with a kernal in .sdf form.
## output is "img_name"_smooth.sdf
##
## Keywords:
## img_name - name of the image (string). -- required
## pix_scale - number of pixels to smooth by (Gauss FWHM) (int). -- required
#
##def smooth_image_kappa(img_name, pix_scale):
#
#    kappa.gausmooth(img_name,img_name.split('.sdf')[0]+'_smooth.sdf',pix_scale)
#
#    print('\nSMOOTH = DONE\n')
#
##def

### unitconv_image ###
# Takes an image in .sdf form and converts the units to Jy/Beam. Default output
# is "img_name"_jypbm.sdf
#
# Keywords:
# img_name - name of the image (string). -- required

def unitconv_image(img_name):
    kappa.cdiv(img_name,1000,img_name.split('_mJybmsm')[0]+'_Jypbmsm_crop.sdf')
    kappa.setunits(img_name.split('_mJybmsm')[0]+'_Jypbmsm_crop.sdf','Jy/beam')

#def

### prepare_image ###
# Takes an image in .sdf form and prepares it by performing cropping, smoothing
# and unit conversions.
#
# Keywords:
# img_name - name of the image (string). --required
# kern_name - smoothing kernal name (string). --required
# kern_fwhm - size of the smoothing kernal FWHM (float). --required
# jypbm_conv - Conversion factor to turn image units into Jy/Beam (float).
# Assumes maps in mJy/sq.arcsecs. Conversion factor from:
# http://pipelinesandarchives.blogspot.ca/2012/02/scuba-2-calibration-redux.html
# crop_radius - radius of the crop (int).
# crop_method - method for cropping (string).

def prepare_image(img_name, crop_radius=1200, crop_method='CIRCLE'):

    #Crop the image.
    crop_image(img_name, crop_radius, crop_method)

    #Smooth the image.
    img_name = img_name.split('/')[-1][:-4]+'_crop.sdf'
    unitconv_image(img_name)

#def
