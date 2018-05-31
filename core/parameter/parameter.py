from datetime import datetime
import numpy as np


class Survey:
    """
    Class to provide the information of surveys which are used to
    for the photometric simulation
    """

    name = ''
    path = ''
    cols = ''
    use_cols = ''

    def __init__(self, name, path, cols, use_cols):
        self.name = name
        self.path = path
        self.cols = cols
        self.use_cols = use_cols


class Parameter:
    """
    Parameter class for HCSS photometric simulation
    """
    surveys = []

    ages = None
    number_stars = 100000
    masses = None
    metallicities = None

    path0 = ''

    def __init__(self, surveys, ages=None, number_stars=10000, path0='', masses=None,
                 metallicities=None):
        """
        Basic parameter class which provides the parameters for the simulation

        :param surveys: A list with the required surveys
        :type surveys: list
        :param ages:
            The times of the simulated magnitudes, default is None which will be converted to
            84 time steps between 1Gyr and 13Gyr
        :type ages: numpy.ndarray
        :param number_stars: The amount of stars per simulated cluster
        :type number_stars: int
        :param path0: The root path to the stellar tracks
        :param path0: str
        :param masses:
            The provides masses of the stellar tracks, default is None which will convert to
            0.01 solar mass sampling rate between 0.1 and 0.8 solar masses and between
            0.8 and 2 solar masses a sampling rate of 0.001 solar mass.
        :type masses: numpy.ndarray
        :param metallicities: List with the metallicities
        :type metallicities: list
        """
        self.surveys = surveys

        # if no ages are given
        if ages is None:
            ages = np.linspace(1e9, 13e9, num=6*84)
        self.ages = ages
        self.number_stars = number_stars
        self.path0 = path0

        # if no masses are given
        if masses is None:
            # combine the asymmetric mass sample
            masses = np.linspace(0.1, 2, num=191)
            masses = np.append(masses, np.linspace(0.8, 1.89, 1091))
            # take every only onetime
            masses = np.unique(masses)
        self.masses = masses

        if metallicities is None:
            metallicities = ['-0_2', '0_0', '0_2', '0_4', '0_5']
        self.metallicities = metallicities

    def get_all_cols(self):
        """
        Returns all used columns
        :return: All the used columns for the output file
        :rtype: list
        """
        out = []
        for s in self.surveys:
            out.extend(s.use_cols)
        return out

    def save(self, path):
        """
        Writes the parameter data to a text file

        :param path: Save place
        :type path: str
        :return:
        """
        if path[-1] != '/':
            path += '/'
        with open('{}parameters.rst'.format(path), 'w') as f:
            f.write('=====================\n')
            f.write('Simulation parameters\n')
            f.write('=====================\n\n')

            f.write('date: {}\n\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

            f.write('Number of stars\n===============\n{}\n\n'.format(self.number_stars))
            f.write('Ages (yr)\n=========\n')
            for a in self.ages:
                f.write('{}\n'.format(str(a)))
            f.write('\n')

            f.write('Masses (solar mass)\n===================\n')
            for m in self.masses:
                f.write('{:05.2f}\n'.format(m))
            f.write('\n')

            f.write('Surveys\n=======\n')
            for s in self.surveys:
                f.write('{}: '.format(s.name))
                for i, c in enumerate(s.use_cols):
                    if i == 0:
                        f.write('{}'.format(c))
                    else:
                        f.write(', {}'.format(c))
                f.write('\n')

    def save_html(self, path):
        """
        Saves the output as a HTML file
        :param path:
        :return:
        """
        if path[-1] != '/':
            path += '/'
        with open('{}parameters.html'.format(path), 'w') as f:
            f.write('<html><body>')
            f.write('<h3>Simulation parameters</h3><br>\n')
            f.write('date: {}<br><br>\n'.format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))

            f.write('<h4>Number of stars</h4>{}\n\n'.format(self.number_stars))

            f.write('<ul><h4>Ages (yr)</h4>\n')
            for a in self.ages:
                f.write('<li>{}</li>\n'.format(str(a)))
            f.write('</ul><br>\n')

            f.write('<ul><h4>Masses (solar mass)</h4>\n')
            for m in self.masses:
                f.write('<li>{:05.2f}</li>\n'.format(m))
            f.write('</ul><br>\n')

            f.write('<h3>Surveys</h3>\n')
            f.write('<table>')
            for s in self.surveys:
                f.write('<tr>')
                f.write('<th>{}</th>'.format(s.name))
                for i, c in enumerate(s.use_cols):
                    f.write('<td style="padding-left:10px">{}</td>'.format(c))
                f.write('</tr>\n')
            f.write('</table>')
            f.write('</body></html>')


class PlotParameter:

    def __init__(self):
        pass
