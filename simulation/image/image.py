from . import positions
import numpy as np
import os
import shutil
from astropy.table import Table


def skymaker(path, save_path, band):
    """
    Runs SkyMaker
    :param path: The path to the input list
    :type path: str
    :param save_path: The path to the output file
    :type save_path: str
    :param band:
        The used band (to identify to right configuration file), if there is no such
        configuration file doesn't exists, then the default configuration are used.
    :type band: str
    :return:
    """
    if os.path.exists('./skymaker_configs/{}.conf'.format(band)):
        os.system('sky {} -c ./skymaker_configs/{}.conf'.format(path, band))
    else:
        os.system('sky {}'.format(path))
    shutil.copy2('{}sky.fits'.format(path.split('list')[0]), save_path)


def scale_factor(diameter, distance_module):
    """
    estimates the scale factor for the tangential component based on the distance module
    :param diameter: The diameter of the source
    :param distance_module: The distance module in magnitude
    :return: The scale factor and the distance derived from the distance module
    """
    distance = 10*pow(10., distance_module/5)
    scale = np.arctan(diameter/distance)
    scale = np.rad2deg(scale)
    scale *= 3600
    return scale, distance


def mass2str(mass):
    """
    Convert the mass value into the MIST mass format
    :param mass: The mass in solar masses
    :return: The MIST mass format
    """
    mass = str(mass)
    if len(mass) < 3:
        mass = '0'+mass
    return mass


class ImageSimulation:

    pos2d = None
    masses = None
    magnitude = None
    time_id = 60
    distance_module = 0
    # pixel scale of different surveys/cameras
    pixel_scale = {'PS_g': 0.25,
                   'PS_r': 0.25,
                   'PS_i': 0.25,
                   'PS_z': 0.25,
                   'PS_y': 0.25,
                   'SDSS_u': 0.396,
                   'SDSS_g': 0.396,
                   'SDSS_r': 0.396,
                   'SDSS_i': 0.396,
                   'SDSS_z': 0.396,
                   '2MASS_J': 1.0,
                   '2MASS_H': 1.0,
                   '2MASS_Ks': 1.0,
                   'WISE_W1': 1.375,
                   'WISE_W2': 1.375,
                   'WISE_W3': 1.375,
                   'WISE_W4': 1.375}

    def __init__(self, n, scale, magnitude, distance_module=15):
        self.pos2d = positions.positions_2d(n, scale)
        self.magnitude = magnitude

        self.distance_module = distance_module

    def __create_skymaker_band__(self, path, band, center):
        out_string = '100\t{: 9.4f}\t{: 9.4f}\t{: 9.4f}\n'
        with open('{}list_{}.dat'.format(path, band), 'w') as f:
            out = []
            for i, pn in enumerate(self.pos2d):
                m = self.magnitude[band][0][i][self.time_id]
                p = pn / self.pixel_scale[band]
                f.write(out_string.format(p[0] + center,
                                          p[1] + center,
                                          m + self.distance_module))
                out.append((i,
                            p[0] + center,
                            p[1] + center,
                            m + 20))
        skymaker('{}list_{}.dat'.format(path, band),
                 '{}img_{}.fits'.format(path, band),
                 band)
        out = Table(rows=out, names=['id', 'x', 'y', 'm'])
        out.write('{}tab_{}.fits'.format(path, band), overwrite=True)

    def create_skymaker_list(self, path, bands):
        if type(bands) != list:
            bands = [bands]
        if path[-1] != '/':
            path += '/'

        center = 1000
        center = center / 2
        for b in bands:
            self.__create_skymaker_band__(path, b, center)


def perform_image_creation(path, band):
    os.system('sky {} -c ./skymaker_configs/{}.conf'.format(path, band))


def perform_sextractor():
    os.system('sex sky.fits')


class ImageSim:
    """
    Class to simulate images based on the HCSS synthetic photometry
    """

    distance_module = 15
    fl = 0
    phot_sim = None
    pos2d = None
    time_id = 10
    stars = []
    tot_mag = 0
    scale = 0
    # pixel scale of different surveys/cameras
    pixel_scale = {'PS_g': 0.25,
                   'PS_r': 0.25,
                   'PS_i': 0.25,
                   'PS_z': 0.25,
                   'PS_y': 0.25,
                   'SDSS_u': 0.396,
                   'SDSS_g': 0.396,
                   'SDSS_r': 0.396,
                   'SDSS_i': 0.396,
                   'SDSS_z': 0.396,
                   '2MASS_J': 1.0,
                   '2MASS_H': 1.0,
                   '2MASS_Ks': 1.0,
                   'WISE_W1': 1.375,
                   'WISE_W2': 1.375,
                   'WISE_W3': 1.375,
                   'WISE_W4': 1.375,
                   'VST_u': 0.21,
                   'VST_g': 0.21,
                   'VST_r': 0.21,
                   'VST_i': 0.21,
                   'VST_z': 0.21}

    def __init__(self, phot_sim, diameter, distance_module=15):
        self.phot_sim = phot_sim
        self.distance_module = distance_module
        self.scale, distance = scale_factor(diameter, distance_module)
        self.calc_pos()

    def calc_pos(self):
        """
        Calculates the position of the different sources in a 2D projection
        :return:
        """
        self.pos2d = positions.positions_2d(self.phot_sim.parameters.number_stars,
                                            self.scale)

    def convert_flux2mag(self, mass, survey_band):
        """
        Converts flux of a star to magnitude. If the flux is
        too small the magnitude is set to 99.

        :param mass: The mass of the star
        :param survey_band: The bandpass name
        :return: The magnitude of the star in the band
        """
        # take the flux of the stars with the mass at the time_id
        fl = self.phot_sim.flux[mass][survey_band][self.time_id]
        # add the flux of this source to the total flux
        self.fl += fl
        mag = -2.5*np.log10(fl)
        # if the magnitude is too small
        if mag > 99:
            mag = 99
        return mag

    def get_masses(self):
        masses = self.phot_sim.mass_sample(1)[0]
        m = []
        for i, k in enumerate(masses):
            s = self.phot_sim.parameters.masses[i]
            m.extend(np.linspace(s, s, k))
        m = np.array(m)
        m *= 100
        m = np.int16(m)
        return m

    def get_magnitude(self, mass, band):
        """
        Returns the magnitude of a star in a certain distance

        :param mass: The mass of the star
        :type mass: float
        :param band: The name of the band
        :type band: str
        :return: The magnitude
        :rtype: float
        """
        mass = mass2str(mass)
        mag = self.convert_flux2mag(mass, band)
        mag += self.distance_module
        return mag

    def create_skymaker_output(self, path, survey, band):
        """
        Creates the source input file for SkyMaker
        :param path: The path where the file should be stored
        :type path: str
        :param survey: A survey object with all the information of the survey
        :type survey: HCSC_Simulation.core.parameter.parameter.Survey
        :param band: The name of the required band
        :type band: str
        :return:
        """
        # reset the stars variable
        self.stars = []
        # select the right pixel scale for the band
        pixel_scale = self.pixel_scale[band]
        # formatted string for the SkyMaker input list
        out_string = '100\t{: 9.4f}\t{: 9.4f}\t{: 9.4f}\n'
        # center of the image
        center = 500

        # path to the SkyMaker input list
        file_path = '{}{}'.format(path, survey.name)

        # open the SkyMaker input list and overwrite the old data
        with open(file_path, 'w') as f:
            self.phot_sim.interpolate_survey(survey)
            m = self.get_masses()

            # add the sources of the cluster to the input list
            for mass, p in zip(m, self.pos2d):
                mag = self.get_magnitude(mass, band)
                # ignore to faint objects (increase performance of the image simulation)
                if mag > 30:
                    continue

                # convert equa-coordinates to pixel coordinates
                c = p / pixel_scale
                c += center
                f.write(out_string.format(c[0], c[1], mag))

            # convert flux to the magnitude
            self.tot_mag = -2.5 * np.log10(self.fl)

            # add the used stars to the input list, every star mass gets one star
            for i, mass in enumerate(np.unique(m)):
                mag = self.get_magnitude(mass, band)
                # calculate the position of the next star on the image
                x = 50 * i
                y = 50 * (int(x / 950) + 1)
                x = x % 950 + 50
                self.stars.append((x, y, mag))
                f.write(out_string.format(x, y, mag))

        # create a astropy table from the star list
        self.stars = Table(rows=self.stars, names=['x', 'y', 'mag'])
        # store the star list as a fits file
        self.stars.write('{}star_list.fits'.format(path), overwrite=True)

        perform_image_creation(file_path, band)
        perform_sextractor()
