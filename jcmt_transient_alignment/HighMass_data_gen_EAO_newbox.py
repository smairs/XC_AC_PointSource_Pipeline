from datetime import datetime
import numpy as np
import numpy.ma as ma
#import seaborn as sns
import os as os
import operator as op
import pickle
import starlink
from astropy.time import Time
from itertools import product
from astropy.io import fits
from astropy.modeling.models import Gaussian2D
from astropy.modeling.fitting import LevMarLSQFitter
from astropy.convolution import interpolate_replace_nans, Gaussian2DKernel, convolve
from scipy.optimize import curve_fit
from collections import defaultdict

from starlink import kappa, convert
starlink.wrapper.change_starpath("/stardev")
#starlink.wrapper.change_starpath("/home/cobr/star-2018A/")

#sns.set_style('whitegrid')
#sns.color_palette('colorblind')


def default_to_regular(d):
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def colourbar(mappable):
    """
    :param mappable: a map axes object taken as input to apply a colourbar to

    :return: Edits the figure and subplot to include a colourbar which is scaled to the map correctly
    """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    figure_one = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return figure_one.colorbar(mappable, cax=cax, format='%g')


def fourier_gaussian_function(axis_one, axis_two, scale=1.0, sigma_x=1.0, sigma_y=1, theta=0):
    xo = axis_one.shape[0] // 2
    yo = axis_two.shape[1] // 2
    if sigma_x == 0:
        sigma_x = 1
    if sigma_y == 0:
        sigma_y = 1
    sigma_x = 1 / sigma_x
    sigma_y = 1 / sigma_y

    a = np.cos(theta) ** 2 / (2 * sigma_x ** 2) + np.sin(theta) ** 2 / (2 * sigma_y ** 2)
    b = -np.sin(2 * theta) / (4 * sigma_x ** 2) + np.sin(2 * theta) / (4 * sigma_y ** 2)
    c = np.sin(theta) ** 2 / (2 * sigma_x ** 2) + np.cos(theta) ** 2 / (2 * sigma_y ** 2)

    fourier_gaussian = scale * np.exp(
        (-4 * np.pi ** 2 / (axis_one.shape[0] ** 2)) *
        (a * (axis_one - xo) ** 2 +
         2 * b * (axis_one - xo) * (axis_two - yo) +
         c * (axis_two - yo) ** 2))
    return fourier_gaussian


def gaussian_fit_ac(auto_correlation):
    # figuring out where I need to clip to, realistically, this SHOULD be at the physical centre (200,200)
    width = 7
    y_max, x_max = np.where(auto_correlation == auto_correlation.max())
    y_max, x_max = int(np.amax(y_max)), int(np.amax(x_max))
    # Setting the middle auto_correlation point to be our estimated value of B for a better fit.

    mask = np.zeros(auto_correlation.shape)
    mask[y_max, x_max] = 1
    ac_masked = ma.masked_array(auto_correlation, mask=mask)

    # clipping map further to better fit a gaussian profile to it
    auto_correlation = ac_masked[y_max - width:y_max + width + 1, x_max - width:x_max + width + 1]

    # generating the gaussian to fit
    x_mesh, y_mesh = np.meshgrid(np.arange(auto_correlation.shape[0]), np.arange(auto_correlation.shape[1]))
    gauss_init = Gaussian2D(
        amplitude=auto_correlation.max(),
        x_mean=auto_correlation.shape[1] // 2,  # location to start fitting gaussian
        y_mean=auto_correlation.shape[0] // 2,  # location to start fitting gaussian
    )
    fitting_gauss = LevMarLSQFitter()  # Fitting method; Levenberg-Marquardt Least Squares algorithm
    best_fit_gauss = fitting_gauss(gauss_init, x_mesh, y_mesh, auto_correlation)  # The best fit for the map
    gauss_model = best_fit_gauss(x_mesh, y_mesh)  # the model itself (if we want to plot it
    try:
        ac_error = np.sqrt(np.diag(fitting_gauss.fit_info['param_cov']))
    except ValueError:
        ac_error = np.ones(10) * -5
    amplitude = float(best_fit_gauss.amplitude.value)
    amplitude_error = ac_error[0]
    sigma_x = float(best_fit_gauss.x_stddev.value)
    sigma_x_error = ac_error[3]
    sigma_y = float(best_fit_gauss.y_stddev.value)
    sigma_y_error = ac_error[4]
    theta = float(best_fit_gauss.theta.value)
    theta_error = ac_error[5]

    return [[amplitude, sigma_x, sigma_y, theta],
            [amplitude_error, sigma_x_error, sigma_y_error, theta_error]], gauss_model


def gaussian_fit_xc(x_correlation):
    # import numpy as np
    # from astropy.modeling.models import Gaussian2D
    # from astropy.modeling.fitting import LevMarLSQFitter
    # figuring out where i need to clip to
    y_center = x_correlation.shape[0] // 2
    x_center = x_correlation.shape[1] // 2  # centre of the Cross-Corr maps default: (200,200)
    width = 7
    y_max, x_max = np.where(x_correlation == x_correlation.max())
    y_max = int(y_max)
    x_max = int(x_max)

    print(x_correlation.shape[0],x_correlation.shape[1])
    print(y_max,x_max) 

    # clipping map further to better fit a gaussian profile to it
    x_correlation = x_correlation[y_max - width:y_max + width + 1, x_max - width:x_max + width + 1]
    if x_correlation.shape[0]==x_correlation.shape[1]:
        # subtracting half the side to then add the mean values after
        x_max -= x_correlation.shape[1] // 2
        y_max -= x_correlation.shape[0] // 2
        # generating the gaussian to fit.

        x_mesh, y_mesh = np.meshgrid(np.arange(x_correlation.shape[0]), np.arange(x_correlation.shape[1]))
        gauss_init = Gaussian2D(
            amplitude=x_correlation.max(),
            x_mean=np.where(x_correlation == x_correlation.max())[1],  # location to start fitting gaussian
            y_mean=np.where(x_correlation == x_correlation.max())[0],  # location to start fitting gaussian
            # fixed={},  # any fixed parameters
            bounds={
                # 'amplitude': (x_correlation.max() * 0.90, x_correlation.max() * 1.10),
                'x_mean': (int(np.where(x_correlation == x_correlation.max())[1]) - 1,
                           int(np.where(x_correlation == x_correlation.max())[1]) + 1),
                'y_mean': (int(np.where(x_correlation == x_correlation.max())[0]) - 1,
                           int(np.where(x_correlation == x_correlation.max())[0]) + 1)
            },  # allowing var in amplitude to better fit gauss
        )
        fitting_gauss = LevMarLSQFitter()  # Fitting method; Levenberg-Marquardt Least Squares algorithm
        print(x_mesh,y_mesh)
        best_fit_gauss = fitting_gauss(gauss_init, x_mesh, y_mesh, x_correlation)  # The best fit for the map

        # now we can get the location of our peak fitted gaussian and add them back to get a total offset
        y_max += best_fit_gauss.y_mean.value  # Finding the distance from 0,0 to the centre gaussian
        x_max += best_fit_gauss.x_mean.value  # and y.
        try:
            x_correlation_error = np.sqrt(np.diag(fitting_gauss.fit_info['param_cov']))
        except ValueError:
            x_correlation_error = np.ones(10) * -5
        offset = (x_center - x_max, y_center - y_max)
        offset_err = (x_correlation_error[1], x_correlation_error[2])
        return offset, offset_err
    else:
        return np.nan,np.nan


def correlate(epoch_1=None, epoch_2=None, highmass_offset_x=0, highmass_offset_y=0, clipped_side=400, clip_only=False, psd=False):
    """
    :param epoch_1:
        2-Dimensional numpy array. Default: None
        When only epoch_1 is passed it is auto correlated with itself
    :param epoch_2:
        2-Dimensional numpy array. Default: None
        When both epoch_1 and epoch_2 are passed the two arrays are cross correlated
    :param clipped_side:
        Integer. Default: 400.
        The length of one side of the clipped array.
    :param clip_only:
        Boolean. Default: False
        When True is passed to clip_only it will only clip epoch_1
    :param psd:
        Boolean. Default: False
        When true is passed the power spectrum is returned
    :return:
    """
    from numpy.fft import fft2, ifft2, fftshift
    if clip_only:
        mid_map_x, mid_map_y = (epoch_1.shape[1] // 2) + highmass_offset_x, (epoch_1.shape[0] // 2) + highmass_offset_y
        clipped_epoch = epoch_1[
                        mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                        mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                        ]
        return clipped_epoch

    elif psd:
        # Clip_only is called first- so offsets have already been applied! we don't need to add them here
        mid_map_x, mid_map_y = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        clipped_epoch = epoch_1[
                        mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                        mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                        ]
        psd = fft2(clipped_epoch) * fft2(clipped_epoch).conj()
        return fftshift(psd)

    elif epoch_1 is None:
        raise Exception('You need to pass a 2D map for this function to work')

    elif epoch_2 is None:
        # Clip_only is called first- so offsets have already been applied! we don't need to add them here
        mid_map_x, mid_map_y = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        clipped_epoch = epoch_1[
                        mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                        mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                        ]
        ac = ifft2(fft2(clipped_epoch) * fft2(clipped_epoch).conj())
        return fftshift(ac)

    else:
        # Clip_only is called first- so offsets have already been applied! we don't need to add them here
        mid_map_x_1, mid_map_y_1 = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        mid_map_x_2, mid_map_y_2 = epoch_2.shape[1] // 2, epoch_2.shape[0] // 2
        clipped_epoch_1 = epoch_1[
                          mid_map_y_1 - clipped_side // 2:mid_map_y_1 + clipped_side // 2 + 1,
                          mid_map_x_1 - clipped_side // 2:mid_map_x_1 + clipped_side // 2 + 1
                          ]
        clipped_epoch_2 = epoch_2[
                          mid_map_y_2 - clipped_side // 2:mid_map_y_2 + clipped_side // 2 + 1,
                          mid_map_x_2 - clipped_side // 2:mid_map_x_2 + clipped_side // 2 + 1
                          ]
        x_correlation = ifft2(fft2(clipped_epoch_1) * fft2(clipped_epoch_2).conj())
        return fftshift(x_correlation)


def f(independent, m, b):
    """
    :param independent: independent variable
    :param m: slope
    :param b: intercept
    :return: y: a quadratic
    """
    dependent = m * independent ** 2 + b
    return dependent


def f_linear(p, independent):
    """
    :param independent: independent variable
    :param p: fitting parameters
    :return: y: a linear monomial
    """

    dependent = p[0] * independent + p[1]
    return dependent


def amp(epoch):
    from numpy import sqrt
    return sqrt(epoch.real ** 2 + epoch.imag ** 2)


def beam_fit(sigma, power_spectrum, required_length_scale):
    from numpy import meshgrid, arange, sqrt
    from numpy.fft import ifft2, fftshift

    axis_1_size = axis_2_size = power_spectrum.shape[0]
    axis_1, axis_2 = meshgrid(arange(axis_1_size), arange(axis_2_size))
    numeric_gaussian = fourier_gaussian_function(axis_1, axis_2, sigma_x=sigma, sigma_y=sigma)  # guess!
    output_gaussian = amp(fftshift(ifft2(numeric_gaussian * numeric_gaussian * power_spectrum)))  # guess amplitude
    [[_, sigma_x, sigma_y, _], _], _ = gaussian_fit_ac(output_gaussian)

    length_scale = sqrt(sigma_x * sigma_y)

    dif = required_length_scale - length_scale
    return abs(dif)


def make_data_dict(region='DR21C',datadir='DR21C',alignment_iteration=0,DIST=7,length=200,kernel_sigma=6,wavelength = '850'):
    """
    :param region: The region to run
    :param datadir: The directory that holds the data for the region listed in parameter region.
    :param alignment_iteration: there will be multiple iterations of the alignment run - this 0-based integer designates which alignment iteration the output file describes
    :param DIST: the distance used for linear fitting and gaussian fitting (use width = RADIUS*2 + 1)
    :param length: the distance used for linear fitting and gaussian fitting (use width = RADIUS*2 + 1)
    :param kernel_sigma: the smoothing kernel (in pixels) to subtract largescale structure (high pass filter)
    """
    # + ===================== +
    # | Global parameters     |
    # + ===================== +
    tol = 0.05
    
    align_smooth_kernel = Gaussian2DKernel(x_stddev=kernel_sigma, y_stddev=kernel_sigma)
    
    data = defaultdict(dict)
    
    data[region] = defaultdict(dict)
    Dates850 = []
    Dates450 = []
    DataRoot = datadir + "/"  # where all the data is stored
    files = []
    for eachfile in os.listdir(DataRoot):
        if os.path.isfile(os.path.join(DataRoot, eachfile)):
            if np.logical_or(eachfile.split('_')[-1] == 'CR3.sdf',eachfile.split('_')[-1] == 'ER3.sdf'):
                if wavelength in eachfile:
                    files.append(eachfile)
    #files = [f for f in os.listdir(DataRoot) if (os.path.isfile(os.path.join(DataRoot, f)) and
    #                                             os.path.join(DataRoot, f)[-4:] ==".sdf")]  # all the files in dir for this wavelength
    files = sorted(files)  # sorting to ensure we select the correct first region
    
    if wavelength == '450':
        scale = 2
    elif wavelength == '850':
        scale = 3
    else:
        scale = 0
    data[region][wavelength] = defaultdict(dict)
    data[region][wavelength]['epoch'] = defaultdict(list)
    
    data[region][wavelength]['dates'] = list()  # to collect all of the dates in the data[region] set
    data[region][wavelength]['JCMT_offset'] = defaultdict(str)  # to use the date as the index
    data[region][wavelength]['header'] = defaultdict(dict)
    
    data[region][wavelength]['XC'] = defaultdict(dict)
    data[region][wavelength]['XC']['offset'] = defaultdict(list)
    data[region][wavelength]['XC']['offset_err'] = defaultdict(list)
    data[region][wavelength]['XC']['alignment'] = defaultdict(list)
    
    data[region][wavelength]['linear'] = defaultdict(dict)
    data[region][wavelength]['linear']['m'] = defaultdict(dict)
    data[region][wavelength]['linear']['m_err'] = defaultdict(dict)
    data[region][wavelength]['linear']['b'] = defaultdict(dict)
    data[region][wavelength]['linear']['b_err'] = defaultdict(dict)
    
    data[region][wavelength]['AC'] = defaultdict(dict)
    data[region][wavelength]['AC']['beam'] = defaultdict(list)
    data[region][wavelength]['AC']['amp'] = defaultdict(list)
    data[region][wavelength]['AC']['amp_err'] = defaultdict(list)
    data[region][wavelength]['AC']['sig_x'] = defaultdict(list)
    data[region][wavelength]['AC']['sig_x_err'] = defaultdict(list)
    data[region][wavelength]['AC']['sig_y'] = defaultdict(list)
    data[region][wavelength]['AC']['sig_y_err'] = defaultdict(list)
    data[region][wavelength]['AC']['theta'] = defaultdict(list)
    data[region][wavelength]['AC']['theta_err'] = defaultdict(list)
    
    
    FEN = files[0]
    FilePath = datadir + "/" + FEN
    OutPath = datadir + "/" + FEN.split('.sdf')[0] + ".fit"
    if os.path.isfile(OutPath):
        pass
    else:
        convert.ndf2fits(FilePath, OutPath)
    FilePath = OutPath
    print('\n\nFIRST EPOCH: '+FilePath+'\n\n')
    FirstEpoch = fits.open(FilePath)  # opening the file in astropy
    FirstEpochData = FirstEpoch[0].data[0]  # Numpy data array for the first epoch
    FirstEpochCentre = np.array([FirstEpoch[0].header['CRPIX1'], FirstEpoch[0].header['CRPIX2']])
   
    if region == 'DR21C':
        highmass_offset_x = -50
        highmass_offset_y = 0
    elif region == 'DR21N':
        highmass_offset_x = -100
        highmass_offset_y = 0
    else:
        highmass_offset_x = 0
        highmass_offset_y = 0
 
    # middle of the map of the first epoch
    FED_MidMapX = (FirstEpochData.shape[1] // 2) + highmass_offset_x
    FED_MidMapY = (FirstEpochData.shape[0] // 2) + highmass_offset_y
    FirstEpochVec = np.array([FirstEpochCentre[0] - FED_MidMapX,
                              FirstEpochCentre[1] - FED_MidMapY])
    FirstEpochData = FirstEpochData[
                     FED_MidMapY - length:FED_MidMapY + length + 1,
                     FED_MidMapX - length:FED_MidMapX + length + 1]
    FirstEpochData_smooth = convolve(FirstEpochData, align_smooth_kernel, normalize_kernel=False)
    FirstEpochData -= FirstEpochData_smooth
    for fn in files:
        if wavelength in fn:
            FilePath = datadir + "/" + fn
    
            tau225_start = float(kappa.fitsval(FilePath, 'WVMTAUST').value)
            tau225_end = float(kappa.fitsval(FilePath, 'WVMTAUEN').value)
            tau225 = sum([tau225_start, tau225_end]) / 2
    
            AirMass_start = float(kappa.fitsval(FilePath, 'AMSTART').value)
            AirMass_end = float(kappa.fitsval(FilePath, 'AMEND').value)
            AirMass = sum([AirMass_start, AirMass_end]) / 2
    
            elev_start = float(kappa.fitsval(FilePath, 'ELSTART').value)
            elev_end = float(kappa.fitsval(FilePath, 'ELEND').value)
            elev = int(round(sum([elev_start, elev_end]) / 2, 0))
    
            OutPath = datadir + "/" + fn[:-4] + ".fit"
    
            if os.path.isfile(OutPath):
                pass
            else:
                convert.ndf2fits(FilePath, OutPath)
    
            FilePath = OutPath
            hdul = fits.open(FilePath)  # opening the file in astropy
            date = ''.join(str(hdul[0].header['DATE-OBS']).split('T')[0].split('-'))  # extract date from the header
            date += '-' + str(hdul[0].header['OBSNUM'])
            JulianDate = str(float(hdul[0].header['MJD-OBS']) + 2400000.5)
            print('Epoch: {:14}'.format(date))
            data[region][wavelength]['header']['airmass'][date] = AirMass
            data[region][wavelength]['header']['t225'][date] = tau225
            data[region][wavelength]['header']['julian_date'][date] = JulianDate
            data[region][wavelength]['header']['elevation'][date] = elev
            data[region][wavelength]['dates'].append(date)
            centre = (hdul[0].header['CRPIX1'], hdul[0].header['CRPIX2'])  # JCMT's alleged centre is
            hdu = hdul[0]  # a nice compact way to store the data for later.
    
            # data[region][wavelength]['epoch'][date].append(hdu)o
            Epoch = hdu.data[0]  # map of the region
            Map_of_Region = interpolate_replace_nans(correlate(Epoch,highmass_offset_x=highmass_offset_x,highmass_offset_y=highmass_offset_y, clip_only=True),
                                                     Gaussian2DKernel(5))
            Map_of_Region_smooth = convolve(Map_of_Region, align_smooth_kernel, normalize_kernel=False)
            Map_of_RegionXC = Map_of_Region - Map_of_Region_smooth
    
            XC = correlate(epoch_1=Map_of_RegionXC, epoch_2=FirstEpochData).real
            PS = correlate(Map_of_Region, psd=True)
            AC = correlate(Map_of_Region).real  # auto correlation of the map
            centre = (hdul[0].header['CRPIX1'], hdul[0].header['CRPIX2'])  # JCMT's alleged centre is
            Vec = np.array([centre[0] - ((hdul[0].shape[2] // 2) + highmass_offset_x),
                            centre[1] - ((hdul[0].shape[1] // 2) + highmass_offset_y)])
            JCMT_offset = FirstEpochVec - Vec  # JCMT offset from headers
            data[region][wavelength]['JCMT_offset'][date] = JCMT_offset  # used for accessing data later.
    
            [[AMP, SIGX, SIGY, THETA], [AMP_ERR, SIGX_ERR, SIGY_ERR, THETA_ERR]], _ = gaussian_fit_ac(AC)
            offset, offset_err = gaussian_fit_xc(XC)
            alignment = JCMT_offset - offset
            Length_Scale = np.sqrt(SIGX * SIGY)
    
            data[region][wavelength]['XC']['offset'][date] = offset * scale
            data[region][wavelength]['XC']['offset_err'][date] = offset_err
            data[region][wavelength]['XC']['alignment'][date] = alignment
    
            data[region][wavelength]['AC']['beam'][date] = Length_Scale
            data[region][wavelength]['AC']['amp'][date] = AMP
            data[region][wavelength]['AC']['amp_err'][date] = AMP_ERR
            data[region][wavelength]['AC']['sig_x'][date] = SIGX
            data[region][wavelength]['AC']['sig_x_err'][date] = SIGX_ERR
            data[region][wavelength]['AC']['sig_y'][date] = SIGY
            data[region][wavelength]['AC']['sig_y_err'][date] = SIGY_ERR
            data[region][wavelength]['AC']['theta'][date] = THETA
            data[region][wavelength]['AC']['theta_err'][date] = THETA_ERR
    
            Clipped_Map_of_Region_LENGTH = np.arange(0, Map_of_Region.shape[0])
            loc = list(product(Clipped_Map_of_Region_LENGTH, Clipped_Map_of_Region_LENGTH))
            MidMapX = AC.shape[1] // 2  # middle of the map x
            MidMapY = AC.shape[0] // 2  # and y
            radius, AC_pows = [], []
            for idx in loc:  # Determining the power at a certain radius
                r = ((idx[0] - MidMapX) ** 2 + (idx[1] - MidMapY) ** 2) ** (1 / 2)
                AC_pow = AC[idx[0], idx[1]].real
                radius.append(r)
                AC_pows.append(AC_pow)
            radius, AC_pows = zip(*sorted(list(zip(radius, AC_pows)), key=op.itemgetter(0)))
            radius = np.array(radius)
            AC_pows = np.array(AC_pows)
    
            num = len(radius[np.where(radius <= DIST)])
            opt_fit_AC, cov_mat_AC = curve_fit(f, radius[1:num], AC_pows[1:num])
            err = np.sqrt(np.diag(cov_mat_AC))
    
            M = opt_fit_AC[0]
            M_err = err[0]
            B = opt_fit_AC[1]
            B_err = err[1]
    
            data[region][wavelength]['linear']['m'][date] = M
            data[region][wavelength]['linear']['m_err'][date] = M_err
            data[region][wavelength]['linear']['b'][date] = B
            data[region][wavelength]['linear']['b_err'][date] = B_err
    
    
    data = default_to_regular(data)
    if not os.path.exists('data'):
        os.system('mkdir data')
    with open("data/data_Transient_"+region+"_run_"+str(alignment_iteration)+"_"+wavelength+".pickle", 'wb') as OUT:
        pickle.dump(data, OUT)
