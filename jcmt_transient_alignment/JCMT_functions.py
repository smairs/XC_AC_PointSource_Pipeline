def default_to_regular(d):
    """
    converts a default dictionary (and any nested defaultdicts)
    into a python standard dictionary
    """
    from collections import defaultdict
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def defaultify(d):
    """
    converts a standard python dictionary to a defaultdict
    """
    from collections import defaultdict
    if not isinstance(d, dict):
        return d
    return defaultdict(lambda: None, {k: defaultify(v) for k, v in d.items()})


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
    """
    This function was made to generate, in fourier space, a gaussian kernel used to convolve
    an epoch up to a common beam size. It will only create a circular gaussian and as is will
    not deal with elliptical gaussians.
    """
    import numpy as np
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
    """
    using astropy's modelling and fitting modules
    to fit a 2d gaussian to an auto-correlation
    of an epoch from the JCMT-Transient survey
    """
    import numpy as np
    import numpy.ma as ma
    from astropy.modeling.models import Gaussian2D
    from astropy.modeling.fitting import LevMarLSQFitter
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
    except:
        ac_error = np.ones(10) * -5
    amplitude = float(best_fit_gauss.amplitude.value)
    amplitude_error = ac_error[0]
    sigma_x = float(best_fit_gauss.x_stddev.value)
    sigma_x_error = ac_error[3]
    sigma_y = float(best_fit_gauss.y_stddev.value)
    sigma_y_error = ac_error[4]
    theta = float(best_fit_gauss.theta.value)
    theta_error = ac_error[5]

    return [
               [amplitude, sigma_x, sigma_y, theta],
                [amplitude_error, sigma_x_error, sigma_y_error, theta_error]
           ], gauss_model


def gaussian_fit_xc(x_correlation):
    """
    Using astropypy's modelling and fitting modules to fit a gaussian
    to the cross-correlation of two epochs in the JCMT transient survey.
    
    """
    import numpy as np
    from astropy.modeling.models import Gaussian2D
    from astropy.modeling.fitting import LevMarLSQFitter
    # figuring out where i need to clip to
    y_center = x_correlation.shape[0] // 2
    x_center = x_correlation.shape[1] // 2  # centre of the Cross-Corr maps default: (200,200)
    width = 7
    y_max, x_max = np.where(x_correlation == x_correlation.max())
    y_max = int(y_max)
    x_max = int(x_max)

    # clipping map further to better fit a gaussian profile to it
    x_correlation = x_correlation[y_max - width:y_max + width + 1, x_max - width:x_max + width + 1]
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
    best_fit_gauss = fitting_gauss(gauss_init, x_mesh, y_mesh, x_correlation)  # The best fit for the map
    gauss_model = best_fit_gauss(x_mesh, y_mesh)  # the model itself (if we want to plot it

    # now we can get the location of our peak fitted gaussian and add them back to get a total offset
    y_max += best_fit_gauss.y_mean.value  # Finding the distance from 0,0 to the centre gaussian
    x_max += best_fit_gauss.x_mean.value  # and y.
    try:
        x_correlation_error = np.sqrt(np.diag(fitting_gauss.fit_info['param_cov']))
    except:
        x_correlation_error = np.ones(10) * -5
    offset = (x_center - x_max, y_center - y_max)
    offset_err = (x_correlation_error[1], x_correlation_error[2])
    return offset, offset_err


def correlate(epoch_1=None, epoch_2=None, clipped_side=400, clip_only=False, psd=False):
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
        mid_map_x, mid_map_y = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        clipped_epoch = epoch_1[mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                                mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                                ]
        return clipped_epoch

    elif psd:
        mid_map_x, mid_map_y = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        clipped_epoch = epoch_1[mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                                mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                                ]
        psd = fft2(clipped_epoch) * fft2(clipped_epoch).conj()
        return fftshift(psd)

    elif epoch_1 is None:
        raise Exception('You need to pass a 2D map for this function to work')

    elif epoch_2 is None:
        mid_map_x, mid_map_y = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        clipped_epoch = epoch_1[mid_map_y - clipped_side // 2:mid_map_y + clipped_side // 2 + 1,
                                mid_map_x - clipped_side // 2:mid_map_x + clipped_side // 2 + 1
                                ]
        ac = ifft2(fft2(clipped_epoch) * fft2(clipped_epoch).conj())
        return fftshift(ac)

    else:
        mid_map_x_1, mid_map_y_1 = epoch_1.shape[1] // 2, epoch_1.shape[0] // 2
        mid_map_x_2, mid_map_y_2 = epoch_2.shape[1] // 2, epoch_2.shape[0] // 2
        clipped_epoch_1 = epoch_1[mid_map_y_1 - clipped_side // 2:mid_map_y_1 + clipped_side // 2 + 1,
                                  mid_map_x_1 - clipped_side // 2:mid_map_x_1 + clipped_side // 2 + 1
                                  ]
        clipped_epoch_2 = epoch_2[mid_map_y_2 - clipped_side // 2:mid_map_y_2 + clipped_side // 2 + 1,
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
    """
    Computers the amplitude of a complex number
    """

    from numpy import sqrt
    return sqrt(epoch.real ** 2 + epoch.imag ** 2)


def beam_fit(sigma, PS, OUT_LENSCALE):
    """
    Creating a function using for a Brent-Dekker minimization routine implemented by scipy.
    
    """
    from numpy import meshgrid, arange, sqrt
    from numpy.fft import ifft2, fftshift

    x_size = y_size = PS.shape[0]
    x, y = meshgrid(arange(x_size), arange(y_size))
    LOW_G2d = fourier_gaussian_function(x, y, sigma_x=sigma, sigma_y=sigma)  # guess!
    TEST_low = amp(fftshift(ifft2(LOW_G2d * LOW_G2d * PS)))  # guess amplitude
    [[_, sigma_x, sigma_y, _], _], _ = gaussian_fit_ac((TEST_low))

    TEST_LenScale = sqrt(sigma_x * sigma_y)

    dif = OUT_LENSCALE - TEST_LenScale
    return abs(dif)
    
