from core import tables, mass_function
from simulation.photometric import colors
from core.parameter.default import get_all
from astropy.table import vstack, Table
import numpy as np
import os
import time


def convert_mass_to_string(m):
    """
    Converts the mass to the mass format of the MIST data files

    :param m: The mass of the required star
    :type m: float
    :return: Mass string in the mist format (:math:`10^{-4} M_\odot`)
    :rtype: str
    """
    m = round(m, 3)
    m_s = str(m)
    if '.' in m_s:
        m_s = '{}{}'.format(*m_s.split('.'))
    while len(m_s) < 4:
        m_s = m_s + '0'
    return m_s


def interp_flux(inter, m_s, flux, ages):
    """
    Converts the interpolated magnitude to flux

    :param inter: A dict with the interpolated magnitudes
    :type inter: dict
    :param m_s: The mass of the star as a string
    :type m_s: str
    :param flux: A dict with all fluxes
    :type flux: dict
    :param ages: A array with the to simulate ages
    :type ages: numpy.ndarray
    :return: The updated flux dict
    """
    for c in inter[m_s]:
        flux[m_s][c] = np.power(10, -0.4 * inter[m_s][c](ages))
    return flux


def matrix_calc(flux, mass_sample):
    """
    Calculate the total magnitude of a simulated mass sample

    :param flux: The flux of the stars with the different masses
    :type flux: numpy.ndarray
    :param mass_sample: The simulated mass sample(s)
    :type mass_sample: numpy.ndarray
    :return: The total magnitude
    :rtype: numpy.ndarray
    """
    m = np.dot(flux, mass_sample)
    return -2.5*np.log10(m)


class SimPointSource:

    parameters = None
    colors = []
    magnitudes = []
    flux = {}
    inter = {}
    tab0 = None
    cluster_mass = []
    directory = ''
    time0 = 0
    last_number_stars = 0

    def __init__(self, parameters):
        if type(parameters) == str:
            path0 = parameters
            parameters = get_all(parameters)
            parameters.path0 = path0
        self.parameters = parameters

    def interpolate(self, survey, m, m_s):
        """
        Interpolate the fluxes

        :param survey: The survey object
        :param m: The mass of the star
        :param m_s: The mass of the star as a string
        :return:
        """
        self.inter = tables.interpolate(survey.path.format(self.directory, int(m * 1000)),
                                        m_s, survey.use_cols, survey.cols, self.inter)
        self.flux = interp_flux(self.inter, m_s, self.flux, self.parameters.ages)

    def add_fluxes(self, mass_sel, m, out, survey):
        """
        Add up the fluxes of the single stars inside the cluster.

        :param mass_sel: The array with the different masses of the stars inside the cluster.
        :type mass_sel: numpy.ndarray
        :param m: The chosen mass
        :type m: float
        :param out: The array to store all simulated data
        :type out: numpy.ndarray
        :param survey: The survey object
        :type survey: parameter.parameter.Survey
        :return: The updated output array and the number of sources with these mass
        """
        # todo interpolate all masses
        p = np.where(mass_sel == m)[0]
        sources = float(len(p))
        if sources > 0:
            m_s = convert_mass_to_string(m)

            for c in survey.use_cols:
                out[c] += sources * self.flux[m_s][c]
        return out, sources

    def get_data(self, mass_sel, survey):
        """
        Simulate the color of a cluster with stars which have the given masses.

        :param mass_sel: The chosen masses of stars inside the cluster
        :type mass_sel: numpy.ndarray
        :param survey: The survey object with the used columns
        :type: parameter.parameter.Survey
        :return: The photometric data and masses of this cluster
        """
        out = self.tab0.copy()

        out['star_age'] = self.parameters.ages
        masses = []
        mass = 0
        for m in np.unique(mass_sel):
            out, star_mass = self.add_fluxes(mass_sel, m, out, survey)
            masses.append(star_mass)
            mass += star_mass
            out[str(m)] = star_mass
        # out['masses'] = np.array(masses)
        self.cluster_mass.append(mass)
        for c in survey.use_cols:
            if c != 'star_age' and c != 'masses':
                out[c] = -2.5 * np.log10(out[c])

        return out

    def combine_data(self):
        """
        Simulate the photometry of one cluster for different surveys.
        The cluster itself is treated as a point source.
        :return: A table with the photometric data for all bands of the different surveys
        """
        mass_sel = mass_function.get_masses(self.parameters.number_stars,
                                            self.parameters.masses,
                                            mass_function.get_mass_function(self.parameters.masses))
        out = self.tab0.copy()
        d = None

        for survey in self.parameters.surveys:
            d = self.get_data(mass_sel, survey)
            for c in survey.use_cols:
                out[c] = d[c]

        try:
            out['star_age'] = d['star_age']

            for k in np.linspace(0.1, 2, num=78 // 2):
                out[str(k)] = d[str(k)]
            return out
        except TypeError:
            return None

    def stats(self, n=25):
        """
        Simulate n different clusters

        :param n: The number of the simulation runs
        :type n: int
        :return:
            A table with a selection of colors of the simulated clusters and
            a table with the complete photometry of all surveys of all clusters.
        """
        colors_out = []
        total = []
        times = []

        self.tab0 = tables.create_zero_table(self.parameters.get_all_cols(), 14)
        self.time0 = time.time()

        for i in range(n):

            out = self.combine_data()
            colors_out.extend(colors.calc_colors(out, self.parameters.ages))
            total.append(out)
            if i % 100 == 0:
                times.append(time.time()-self.time0)
                print('run: {}'.format(i))
                print('{:6.2f} sec left'.format((n-i-1)*np.median(times)/100))
                self.time0 = time.time()

        total = np.vstack(total)
        return Table(rows=colors_out, names=['gr', 'ri', 'iz', 'zy',
                                             'jh', 'hk', 'w1j', 'NUV_u', 'age']), total

    def metal(self, m, n=1000):
        """
        Simulate a bunch of clusters for a given metallicity.

        :param m:
            The metallicity which is also the name of the directory with the synthetic photometry
        :type m: str
        :param n: The number of simulation runs
        :type n: int
        :return:
            A table with a selection of different colors of all simulated clusters and
            the complete photometry of all surveys of all simulated clusters.
        """
        self.last_number_stars = n
        self.directory = m
        self.interpolate_all()
        st, total = self.stats(n=n)

        fe_h = float('{}.{}'.format(*m.split('_')))
        st['fe/h'] = fe_h
        total = Table(total)
        total['fe/h'] = fe_h
        for k in np.linspace(0.1, 2, num=78 // 2):
            k = str(k)
            m = total[k]
            del total[k]
            total[k] = m[:, 0]
        return st, total

    def save(self, path):
        """
        Saves the simulated data as a fits file.

        :param path: The path to the directory
        :type path: str
        :return:
        """
        if path[-1] != '/':
            path += '/'
        self.magnitudes.meta['RUNS'] = self.last_number_stars
        self.magnitudes.write('{}simulated_magnitudes_time.fits'.format(path),
                              overwrite=True)

        self.colors.write('{}simulated_colors_time.fits'.format(path),
                          overwrite=True)

    def interpolate_survey(self, survey):
        """
        Interpolate the magnitudes over the star age for all masses for this survey

        :param survey: Survey object of the required survey
        :type survey: HCSC_Simulation.core.parameter.parameter.Survey
        :return:
        """
        for m in self.parameters.masses:

            k = convert_mass_to_string(round(m, 3))

            if k in self.inter.keys():
                if survey.use_cols[0] in self.inter[k].keys():
                    continue
            self.inter = tables.interpolate(survey.path.format(self.directory, k),
                                            k, survey.use_cols, survey.cols, self.inter)
            if k not in self.flux.keys():
                self.flux[k] = {}
            self.flux = interp_flux(self.inter, k, self.flux, self.parameters.ages)

    def interpolate_all(self):
        for survey in self.parameters.surveys:
            self.interpolate_survey(survey)

    def start(self, n=100):
        """
        Starts the simulation of the clusters

        :param n: The number of simulation runs
        :type n: int
        :return:
        """
        for f in os.listdir(self.parameters.path0):
            if os.path.isfile(self.parameters.path0 + f):
                continue

            self.inter = {}
            self.flux = {}
            color_data, mags = self.metal(f, n)
            self.colors.append(color_data)
            self.magnitudes.append(mags)

        self.colors = vstack(self.colors)
        self.magnitudes = vstack(self.magnitudes)
        print('simulation done')

    def to_matrix(self):
        """
        Converts the flux vector/matrix to a tensor.
        The dimensions are length of the used bands, number of used time steps and
        the number of used masses.
        :return:
        """
        # load the names of the used bands
        bands = list(self.flux[list(self.flux.keys())[0]].keys())
        flux = np.zeros((len(bands),
                         len(self.flux[list(self.flux.keys())[0]][bands[0]]),
                         len(self.parameters.masses)))

        for i, b in enumerate(bands):
            for k, m in enumerate(self.parameters.masses):
                f = self.flux[convert_mass_to_string(m)][b]
                flux[i, :, k] = f

        return flux

    def mass_counter(self, prop):
        mass_sel = mass_function.get_masses(self.parameters.number_stars,
                                            self.parameters.masses,
                                            prop)
        mass_sel = np.append(mass_sel, self.parameters.masses)
        return np.unique(mass_sel, return_counts=True)[-1]-1

    def mass_sample(self, n):
        prob = mass_function.get_mass_function(self.parameters.masses)
        out = np.zeros((n, len(self.parameters.masses)))
        t0 = time.time()
        stars = self.parameters.number_stars

        masses = self.parameters.masses
        for i in range(n):
            if i % 100 == 0:
                print(i)
                print(time.time()-t0)
                t0 = time.time()
                print(masses.shape)
            out[i] = mass_function.get_mass_sample(stars, masses, prob)
        return out

    def get_magnitudes(self, n, path):
        try:
            if path[-1] != '/':
                path += '/'
        except IndexError:
            pass

        number = self.mass_sample(n)
        for d in self.parameters.metallicities:
            self.directory = d
            self.interpolate_all()
            m = self.to_matrix()

            m = matrix_calc(m, np.transpose(number))
            bands = list(self.flux[list(self.flux.keys())[0]].keys())
            tab = Table()
            for i, b in enumerate(bands):
                k = m[i]
                tab[b] = np.transpose(k)
            tab['age'] = np.tile(self.parameters.ages, (len(tab), 1))
            tab['mass'] = np.int32(number)

            if path == '':
                return tab

            tab.write(path+'sim_mags_{}.fits'.format(d), overwrite=True)
            print(np.median(np.sum(number, axis=-1)))
        # m = m.reshape((m.shape[0]*m.shape[-1], m.shape[1]))
        # print(m.shape)
        # m = Table(rows=m)
        return tab


if __name__ == '__main__':
    sps = SimPointSource('/Volumes/UNTITLED/hcss_colors/mist/new2/')
    # sps.start(10000)
    # print(sps.magnitudes)
    # sps.directory = '0_5'
    # sps.interpolate_all()
    magnitudes = sps.get_magnitudes(10000, '/Volumes/UNTITLED/hcss_colors/mist/')
    magnitudes.write('/Volumes/UNTITLED/hcss_colors/mist/test_new2_0_4_s_0_3_ms.fits',
                     overwrite=True)
    # print(m.shape)
    print(len(magnitudes.colnames))
    # m.show_in_browser(jsviewer=True)
