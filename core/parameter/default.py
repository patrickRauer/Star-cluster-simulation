from . import parameter, cols


def get_all(path):
    path0 = '{}{}/{}'
    file_path = '00{}0M.track.eep.fits'
    sdss = parameter.Survey('SDSS', path0.format(path, '{}/sdss', file_path),
                            cols.sdss_cols, cols.sdss_cols_use)
    panstarrs = parameter.Survey('PAN-STARRS', path0.format(path, '{}/panstarrs', file_path),
                                 cols.pan_cols, cols.pan_cols_use)
    two_mass = parameter.Survey('2MASS', path0.format(path, '{}/2mass', file_path),
                                cols.two_mass_cols, cols.two_mass_cols_use)
    wise = parameter.Survey('WISE', path0.format(path, '{}/wise', file_path),
                            cols.wise_cols, cols.wise_cols_use)
    galex = parameter.Survey('GALEX', path0.format(path, '{}/galex', file_path),
                             cols.galex_cols, cols.galex_cols_use)
    mega = parameter.Survey('MegaCam', path0.format(path, '{}/megacam', file_path),
                            cols.megacam_cols, cols.megacam_cols_use)
    decam = parameter.Survey('DECam', path0.format(path, '{}/decam', file_path),
                             cols.decam_cols, cols.decam_cols_use)
    return parameter.Parameter([sdss, panstarrs, two_mass, wise, galex, mega, decam])
