import numpy as np


def calc_colors(out, ages):
    rs = []
    for a in ages:
        dt = out['star_age']-a
        dt = np.abs(dt)
        p = np.where(dt == np.min(dt))[0][0]
        rs.append((out['PS_g'][p]-out['PS_r'][p],
                   out['PS_r'][p]-out['PS_i'][p],
                   out['PS_i'][p]-out['PS_z'][p],
                   out['PS_z'][p]-out['PS_y'][p],
                   out['2MASS_J'][p]-out['2MASS_H'][p],
                   out['2MASS_H'][p]-out['2MASS_Ks'][p],
                   out['WISE_W1'][p]-out['2MASS_J'][p],
                   out['GALEX_NUV'][p]-out['SDSS_u'][p],
                   round(out['star_age'][p]/1e9, 2)))
    return rs
