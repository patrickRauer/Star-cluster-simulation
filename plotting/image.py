#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 13:30:06 2018

@author: patrickr
"""

from astropy.io import fits
from matplotlib.colors import LogNorm
import pylab as pl
import numpy as np


path = '/Users/patrickr/Documents/GitHub/HCSC_Simulation/sky.fits'
with fits.open(path) as fi:
    pl.clf()
    sp = pl.subplot()
    sp.imshow(fi[0].data[450:550,
              450:550], norm=LogNorm(vmin=1e3, vmax=1e5), cmap='gray_r')
    sp.contour(np.linspace(0, 100, num=100), 
               np.linspace(0, 100, num=100), fi[0].data[450:550,
              450:550], [3.43e3])
    sp.text(0, 5, 'stars$_b$=10000')
    sp.text(0, 10, 'r$_{eff}$=0.1ages = np.linspace(1e9, 13e9, num=6*14)pc')
    sp.text(0, 15, '$\\Delta m$=15')
    sp.text(0, 20, 'Filter: SDSS g')
    sp.set_xlabel('X [pixel]')
    sp.set_ylabel('Y [pixel]')
    pl.show()