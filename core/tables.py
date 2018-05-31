from astropy.table import Table
from astropy.io import fits
from scipy.interpolate import interp1d
import numpy as np
import os

pan_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf',
            'PS_g', 'PS_r', 'PS_i', 'PS_z', 'PS_y', 'PS_w',
            'PS_open', 'phase']

pan_cols_use = ['PS_g', 'PS_r', 'PS_i', 'PS_z', 'PS_y', 'PS_w']

wise_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf',
             'WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4', 'phase']
wise_cols_use = ['WISE_W1', 'WISE_W2', 'WISE_W3', 'WISE_W4']

two_mass_cols = ['star_age', 'log_Teff', 'log_g', 'log_L',
                 'Z_surf', 'Bessell_U',
                 'Bessell_B', 'Bessell_V', 'Bessell_R', 'Bessell_I',
                 '2MASS_J', '2MASS_H', '2MASS_Ks', 'Kepler_Kp', 'Kepler_D51',
                 'Hipparcos_Hp', 'Tycho_B', 'Tycho_V', 'Gaia_G', 'Gaia_BP',
                 'Gaia_RP']

two_mass_cols_use = ['2MASS_J', '2MASS_H', '2MASS_Ks', 'Kepler_Kp', 'Kepler_D51',
                     'Gaia_G', 'Gaia_BP', 'Gaia_RP']

galex_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf', 'GALEX_FUV',
              'GALEX_NUV', 'phase']
galex_cols_use = ['GALEX_FUV', 'GALEX_NUV']

sdss_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf', 'SDSS_u',
             'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', 'phase']
sdss_cols_use = ['SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z']


main_path = ''


def all_cols():
    return [pan_cols, wise_cols, two_mass_cols, galex_cols,
            sdss_cols]


def all_use_cols():
    return [pan_cols_use, wise_cols_use, two_mass_cols_use,
            galex_cols_use, sdss_cols_use]


def all_path():
    path0 = '{}{}/{}'
    pan_path = path0.format(main_path, 'panstarrs', '{:07d}M.track.eep.fits')
    wise_path = path0.format(main_path, 'wise', '{:07d}M.track.eep.fits')
    two_mass_path = path0.format(main_path, '2mass', '{:07d}M.track.eep.fits')
    sdss_path = path0.format(main_path, 'sdss', '{:07d}M.track.eep.fits')
    galex_path = path0.format(main_path, 'galex', '{:07d}M.track.eep.fits')
    paths = [pan_path, wise_path, two_mass_path, galex_path,
             sdss_path]
    return paths


def create_zero_table(use_cols, times):
    cols = ['star_age']
    cols.extend(use_cols)
    tab0 = Table(rows=np.zeros((times, len(cols))), names=cols)
    for k in np.linspace(0.1, 2, num=78 // 2):
        tab0[str(k)] = 0
    tab0 = np.array(tab0)
    return tab0


def rename_columns(tab, cols):
    """
    Rename columns
    :param tab: The table with columns name to change
    :type tab: astropy.table.Table
    :param cols: The new names of the columns
    :type cols: list
    :return: The input table with the new column names
    :rtype: astropy.table.Table
    """
    for i, c in enumerate(cols):
        tab.rename_column('col{}'.format(i + 1), c)
    return tab


def read_data(path, all_columns):
    """
    Read the isochrones from a file
    :param path: The path to the file
    :type path: str
    :param all_columns: Name of th columns
    :type all_columns: list
    :return: A table with the input data
    :rtype: astropy.table.Table
    """
    # if the file has a fits equvilent
    if os.path.exists(path.split('.cmd')[0] + '.fits'):
        with fits.open(path.split('.cmd')[0] + '.fits') as fi:
            tab = np.array(fi[1].data)

    else:
        if 'fits' in path:
            try:
                with fits.open(path) as fi:
                    return np.array(fi[1].data)
            except OSError as e:
                print(path)
                raise OSError(e)

        tab = Table.read(path, format='ascii')
        tab = rename_columns(tab, all_columns)
        tab.write(path.split('.cmd')[0] + '.fits')

    return tab


def load_inter(path, cols, all_columns):
    """
    Interpolate the magnitudes over the star age

    :param path: Path to the file
    :type path: str
    :param cols: List with used column names
    :type cols: list
    :param all_columns: List with all column names
    :type all_columns: list
    :return: Interpolated magnitude over the star age
    :rtype: interp1d
    """
    out = {}
    tab = read_data(path, all_columns)

    for c in cols:
        out[c] = interp1d(tab['star_age'], tab[c], bounds_error=False,
                          fill_value=9999)
    return out


def interpolate(path, m_s, cols, all_columns, inter=None):
    """

    :param path:
    :param m_s:
    :param cols:
    :param all_columns:
    :param inter:
    :return:
    """
    if inter is None:
        inter = {}
    inter[m_s] = load_inter(path, cols, all_columns)
    return inter
