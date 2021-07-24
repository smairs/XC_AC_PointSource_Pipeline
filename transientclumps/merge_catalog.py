import numpy as np
from astropy.io import fits

def merge_catalog(
        cat_file, ref_file, date, obsnum, filename, ref_index, cat_index):
    cat_hdus = fits.open(cat_file)
    ref_hdus = fits.open(ref_file)

    cat = cat_hdus[1].data
    ref = ref_hdus[1].data

    # If this is a 1-element tuple - find the length of its single entry.
    n_match = len(ref_index[0]) if isinstance(ref_index, tuple) else len(ref_index)

    dx = np.zeros(n_match)
    dy = np.zeros(n_match)
    peak_ratio = np.zeros(n_match)
    radius_ratio = np.zeros(n_match)
    date_scan_ref = np.zeros(n_match)
    date_scan_cat = np.zeros(n_match)

    for i in range(n_match):
        # The "delta x". This source's x position minus the reference x position
        # - properly including the declination argument
        dx[i] = -1.0 * (
            (ref[ref_index][i][3] - cat[cat_index][i][3]) * 3600.0
            * np.cos(ref[ref_index][i][4] * np.pi / 180.0))

        # The "delta y". This source's y position minus the reference x position
        dy[i] = -1.0 * (ref[ref_index][i][4] - cat[cat_index][i][4]) * 3600.0

        # The "peak ratio". This source's peak brightness divided by the
        # reference peak brightness
        peak_ratio[i] = 1.0 / ((ref[ref_index][i][8] / cat[cat_index][i][8]))

        # The "radius ratio". This source's effective radius divided by the
        # reference effective radius.
        radius_ratio[i] = 1.0 / (
            (np.sqrt(np.multiply(ref[ref_index][i][13], ref[ref_index][i][14])) * 1.5) /
            (np.sqrt(np.multiply(cat[cat_index][i][13], cat[cat_index][i][14])) * 1.5))

        # The date.scan value for the reference map
        date_scan_ref[i] = '20000000.999'  # ref_date + '.' + str(ref_obsnum)

        date_scan_cat[i] = date + '.' + str(obsnum)

    columns = fits.ColDefs([
        fits.Column(name='dx', format='D', array=dx),
        fits.Column(name='dy', format='D', array=dy),
        fits.Column(name='peak_ratio', format='D', array=peak_ratio),
        fits.Column(name='radius_ratio', format='D', array=radius_ratio),
        fits.Column(name='date.scan_ref', format='11A', array=date_scan_ref),
    ])

    # Now merge the new columns above with the reference catalogue and give all
    # the headings "_ref" suffixes
    for i in range(len(ref.columns)):
        name = ref.columns.names[i]

        columns.add_col(fits.Column(
            name=name + '_ref',
            format=ref.columns.formats[i],
            array=ref[name][ref_index]))

    columns.add_col(
        fits.Column(name='date.scan_' + date, format='11A', array=date_scan_cat))

    # Now merge the catalog columns, adding the date as a suffix
    for i in range(len(cat.columns)):
        name = cat.columns.names[i]

        columns.add_col(fits.Column(
            name=name + '_' + date,
            format=cat.columns.formats[i],
            array=cat[name][cat_index]))

    fits.BinTableHDU.from_columns(columns).writeto(filename)
