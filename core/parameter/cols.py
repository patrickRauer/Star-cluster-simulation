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
                     'Gaia_G', 'Gaia_BP', 'Gaia_RP', 'Bessell_B', 'Bessell_V', 'Bessell_I']

galex_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf', 'GALEX_FUV',
              'GALEX_NUV', 'phase']
galex_cols_use = ['GALEX_FUV', 'GALEX_NUV']

sdss_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf', 'SDSS_u',
             'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z', 'phase']
sdss_cols_use = ['SDSS_u', 'SDSS_g', 'SDSS_r', 'SDSS_i', 'SDSS_z']

megacam_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf',
                'CFHT_u', 'CFHT_g', 'CFHT_r', 'CFHT_i_new', 'CFHT_i_old',
                'CFHT_z', 'phase']

megacam_cols_use = ['CFHT_u', 'CFHT_g', 'CFHT_r', 'CFHT_i_new', 'CFHT_i_old',
                    'CFHT_z']

decam_cols = ['star_age', 'log_Teff', 'log_g', 'log_L', 'Z_surf',
              'DECam_u', 'DECam_g', 'DECam_r', 'DECam_i', 'DECam_z',
              'DECam_Y', 'phase']
decam_cols_use = ['DECam_u', 'DECam_g', 'DECam_r', 'DECam_i', 'DECam_z',
                  'DECam_Y']
