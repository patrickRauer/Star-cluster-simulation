from astropy.table import Table, join
from sklearn.cluster import DBSCAN
import numpy as np
import pylab as pl


def get_identification(x1, y1, x2, y2, eps=1.4):
    """
    Creates a cross match identification between two samples based on
    a DBSCAN algorithm.

    :param x1: The x-values of the first sample
    :type x1: numpy.ndarray
    :param y1: The y-values of the first sample
    :type y1: numpy.ndarray
    :param x2: The x-values of the second sample
    :type x2: numpy.ndarray
    :param y2: The y-sample of the second sample
    :type y2: numpy.ndarray
    :param eps:
        The maximal distance between two points,
        default is 1.4 which is a proper value for pixel coordinates
    :type eps: float
    :return:
        The labels (identification numbers) for both samples in the order of the input
    :rtype: numpy.ndarray, numpy.ndarray
    """
    l1 = len(x1)
    # prepare the input values for the DBSCAN
    x = np.zeros((l1, 2))
    x[:l1, 0] = x1
    x[:l1, 1] = y1
    x[l1:, 0] = x2
    x[l1:, 1] = y2

    # create a DBSCAN object with a min sample size of 2 for the cross matching
    db = DBSCAN(eps=eps, min_samples=2)
    # fit the coordinates
    db.fit(x)

    # return the labels of the samples
    return db.labels_[:l1], db.labels_[l1:]


class ImgSimSExtractor:
    """
    This class provides the interaction with the SExtractor results and the
    image simulation input list
    """
    sextractor_detection = None
    simulated_image = None
    cross_match = None

    def __init__(self, image_sim, path='./test.fits'):
        self.sextractor_detection = Table.read(path)
        self.simulated_image = image_sim

    def __cross_match__(self):
        """
        Performs the cross-match between the SExtractor and image simulation list

        :return:
        """
        # copy the original data to avoid problems later
        sext = self.sextractor_detection.copy()
        stars = self.simulated_image.stars.copy()

        # get the identifications between the simulation list and the results of SExtractor
        l1, l2 = get_identification(sext['X_IMAGE'],
                                    sext['Y_IMAGE'],
                                    stars['x'],
                                    stars['y'])
        # set the label to the tables
        sext['id'] = l1
        stars['id'] = l2

        # exclude such sources which doesn't have a cross-match
        sext = sext[sext['id'] >= 0]
        stars = stars[stars['id'] >= 0]

        # join both tables with id as key
        self.cross_match = join(sext, stars, keys='id')

        # calculate all distances to the center of the image
        # and set the total simulated magnitude to the central source
        r = np.hypot(self.cross_match['X_IMAGE']-500,
                     self.cross_match['Y_IMAGE']-500)
        self.cross_match['mag'][r == np.min(r)] = self.simulated_image.tot_mag

    def plot_magnitudes(self):
        """
        Plot the synthetic magnitude vs the extracted magnitude
        :return:
        """
        fig = pl.figure()
        fig.clf()
        sp = fig.add_subplot(111)
        sp.scatter(self.cross_match['mag'],
                   self.cross_match['MAG_BEST'],
                   marker='.')
        sp.set_xlabel('synthetic magnitude')
        sp.set_ylabel('extracted magnitude')
        pl.show()
