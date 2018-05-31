#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 15 12:32:07 2018

@author: patrickr
"""

from astropy.io import fits
from astropy.table import Table
import numpy as np
import pylab as pl


def load_sdss_classification():
    path = '/Users/patrickr/Downloads/sdss_wise_mass.fits'
    with fits.open(path) as fi:
        d = fi[1].data
        p = np.where((d['u'] > 0) &
                     (d['g'] > 0) &
                     (d['r'] > 0) &
                     (d['i'] > 0) &
                     (d['z'] > 0))[0]
        d = d[p]
        for c in np.unique(d['subClass']):
            p = np.where(d['subClass'] == c)[0]
            d['subClass'][p] = c[:2]
        p = np.where(d['subClass'] != '')[0]
        d = d[p]
#        print(np.unique(d['subClass']))
        d = Table(d)
        d['classID'] = 0
        for i, c in enumerate(np.unique(d['subClass'])):
            p = np.where(d['subClass'] == c)[0]
            d['classID'][p] = i
            
        d['ug'] = d['u']-d['g']
        d['gr'] = d['g']-d['r']
        d['ri'] = d['r']-d['i']
        d['iz'] = d['i']-d['z']
        d['jh'] = d['j_m_2mass']-d['h_m_2mass']
        d['hk'] = d['h_m_2mass']-d['k_m_2mass']
        d['w1j'] = d['w1mpro']-d['j_m_2mass']
        p = np.where((d['ug'] > 0.5) &
                     (d['ug'] < 3.5) &
                     (d['gr'] > 0) &
                     (d['gr'] < 2) &
                     (d['ri'] > 0) &
                     (d['ri'] < 1.5) &
                     (d['iz'] > -0.05) &
                     (d['iz'] < 0.9) &
                     (d['jh'] > 0) &
                     (d['jh'] < 1.2) &
                     (d['hk'] > -0.2) &
                     (d['hk'] < 0.75))[0]
        d = d[p]
        return np.array(d)
    
def load_sim(met):
    path = '/Volumes/UNTITLED/hcss_colors/mist/results/full/{}/'
    file_name = 'sim_mags_{}.fits'
    out = []
    for i in ['1', '2', '3', '4']:
        
        with fits.open(path.format(i)+file_name.format(met)) as fi:
            out.append(np.array(fi[1].data))
    out = np.vstack(out)
    return out


def transform_hist_x(x):
    return x[:-1] + (x[1]-x[0])/2

def hist1d(x, min_value, max_value):
    hist, e = np.histogram(x, bins=100, 
                           range=[min_value, max_value])
    hist = hist/np.sum(hist)
    e = transform_hist_x(e)
    return hist, e
    
    
def plot_color_color(c1, c2, c3, c4):
    sdss = load_sdss_classification()
    
    conv = {'u': 'SDSS_u', 'g': 'SDSS_g', 'r': 'SDSS_r',
            'i': 'SDSS_i', 'z': 'SDSS_z', 'j_m_2mass': '2MASS_J',
            'h_m_2mass': '2MASS_H', 'k_m_2mass': '2MASS_Ks',
            'w1mpro': 'WISE_W1'}
    label_conv = {'u': 'u$_{SDSS}$', 'g': 'g$_{SDSS}$', 'r': 'r$_{SDSS}$',
                  'i': 'i$_{SDSS}$', 'z': 'z$_{SDSS}$', 
                  'j_m_2mass': 'J$_{2MASS}$',
                  'h_m_2mass': 'H$_{2MASS}$', 'k_m_2mass': 'K$_{s}_{2MASS}$',
                  'w1mpro': 'W1$_{WISE}$'}
    sim = load_sim('0_0')
    
    sim_shape = sim['PS_r'].shape
    sim_shape = sim_shape[0]*sim_shape[1]*sim_shape[2]
    c1 = 'j_m_2mass'
    c2 = 'h_m_2mass'
    c3 = 'h_m_2mass'
    c4 = 'k_m_2mass'
    
    pl.clf()
    sp = pl.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    
    cc1 = sdss[c1]-sdss[c2]
    cc2 = sdss[c3]-sdss[c4]
    sp.hist2d(cc1, cc2, bins=200,
              cmap='gray_r')
    
    sim_cc1 = (sim[conv[c1]]-sim[conv[c2]]).reshape((sim_shape, ))
    sim_cc2 = (sim[conv[c3]]-sim[conv[c4]]).reshape((sim_shape, ))
    hist, x, y = np.histogram2d(sim_cc1,
                                sim_cc2,
                                bins=100,
                                range=[[np.min(cc1), np.max(cc1)],
                                        [np.min(cc2), np.max(cc2)]])
    hist = hist/np.sum(hist)
    
    x = transform_hist_x(x)
    y = transform_hist_x(y)
    
    ct = sp.contour(x, y, np.transpose(hist), [0.001, 0.005, 0.01, 0.05],
                    cmap='jet')
    
    sp.set_xlabel('{} - {}'.format(label_conv[c1], label_conv[c2]))
    sp.set_ylabel('{} - {}'.format(label_conv[c3], label_conv[c4]))
    
    sp = pl.subplot2grid((3, 3), (0, 0), colspan=2)
    hist, e = hist1d(sim_cc1, np.min(cc1), np.max(cc1))
    sp.step(e, hist, color='k')
    
    sp.set_xlim(np.min(cc1), np.max(cc1))
    sp.set_ylim(0, 0.1)
    pl.yticks([0.05, 0.1], [0.05, 0.1], rotation=30)
    
    sp.get_xaxis().set_visible(False)
    
    sp = pl.subplot2grid((3, 3), (1, 2), rowspan=2)
    hist, e = hist1d(sim_cc2, np.min(cc2), np.max(cc2))
    sp.step(hist, e, color='k')
    
    sp.set_ylim(np.min(cc2), np.max(cc2))
    sp.set_xlim(0, 0.1)
    
    sp.get_yaxis().set_visible(False)
    pl.xticks([0.05, 0.1, 0.2], [0.05, 0.1, 0.2], rotation=-30)
    
    pl.subplots_adjust(wspace=0, hspace=0)
