#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  9 13:18:04 2018

@author: patrickr
"""

from astropy.io import fits
import numpy as np
import pylab as pl

def panstarrs_color_curve():
    time_id = 60
    
    with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['PS_g']-d['PS_r'], axis=0), 
                    np.std(d['PS_g']-d['PS_r'], axis=0),label='g-r',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['PS_r']-d['PS_i'], axis=0), 
                    np.std(d['PS_r']-d['PS_i'], axis=0), label='r-i',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['PS_i']-d['PS_z'], axis=0), 
                    np.std(d['PS_i']-d['PS_z'], axis=0), label='i-z',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['PS_z']-d['PS_y'], axis=0), 
                    np.std(d['PS_z']-d['PS_y'], axis=0), label='z-y',
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('color')
        pl.legend(loc='best')

def sdss_color_curve():
    with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_g']-d['SDSS_r'], axis=0), 
                    np.std(d['SDSS_g']-d['SDSS_r'], axis=0),label='g-r',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_r']-d['SDSS_i'], axis=0), 
                    np.std(d['SDSS_r']-d['SDSS_i'], axis=0), label='r-i',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_i']-d['SDSS_z'], axis=0), 
                    np.std(d['SDSS_i']-d['SDSS_z'], axis=0), label='i-z',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_u']-d['SDSS_g'], axis=0), 
                    np.std(d['SDSS_u']-d['SDSS_g'], axis=0), label='u-g',
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('color')
        pl.legend(loc='best')

def twomass_color_curve():
    

     with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['2MASS_J']-d['2MASS_H'], axis=0), 
                    np.std(d['2MASS_J']-d['2MASS_H'], axis=0),label='J-H',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['2MASS_H']-d['2MASS_Ks'], axis=0), 
                    np.std(d['2MASS_H']-d['2MASS_Ks'], axis=0), label='H-K',
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('color')
        pl.legend(loc='best')
def twomass_color_color():
    

     with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(np.median(d['2MASS_J']-d['2MASS_H'], axis=0), 
                    np.median(d['2MASS_H']-d['2MASS_Ks'], axis=0),
                    xerr=np.std(d['2MASS_J']-d['2MASS_H'], axis=0),  
                    yerr=np.std(d['2MASS_H']-d['2MASS_Ks'], axis=0),
                    capsize=2)
#        sp.errorbar(d['age'][0]/1e9, label='H-K',
#                    capsize=2)
        
        sp.set_xlabel('J-H')
        sp.set_ylabel('H-K$_s$')
#        pl.legend(loc='best')
        
def twomass_color_h_ks_compare():
    

     with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['2MASS_H']-d['2MASS_Ks'], axis=0), 
                    np.std(d['2MASS_H']-d['2MASS_Ks'], axis=0),
                    label='min. 0.1 $M_\odot$',
                    capsize=2)
        fi = fits.open('/Volumes/UNTITLED/hcss_colors/mist/sim_mags_0_4_0_3ms.fits')
        d = fi[1].data
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['2MASS_H']-d['2MASS_Ks'], axis=0), 
                    np.std(d['2MASS_H']-d['2MASS_Ks'], axis=0),
                    label='min. 0.3 $M_\odot$',
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('H-K$_s$')
        pl.legend(loc='best')
        
        

def sdss_panstarrs_dif_curve():
    
    with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_g']-d['PS_g'], axis=0), 
                    np.std(d['SDSS_g']-d['PS_g'], axis=0),label='g',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_r']-d['PS_r'], axis=0), 
                    np.std(d['SDSS_r']-d['PS_r'], axis=0),label='r',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_i']-d['PS_i'], axis=0), 
                    np.std(d['SDSS_i']-d['PS_i'], axis=0),label='i',
                    capsize=2)
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['SDSS_z']-d['PS_z'], axis=0), 
                    np.std(d['SDSS_z']-d['PS_z'], axis=0),label='z',
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('SDSS-PANSTARRS')
        pl.legend(loc='best')

def sdss_color_curve2():
    with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.scatter(d['mass'][:, 161], d['mass'][:, 160],c=d['SDSS_i'][:, 55])
#        sp.hist(d['SDSS_i'][:, 55], bins=100, 
#                histtype='step')
#        sp.hist(d['2MASS_H'][:, 55], bins=100, 
#                histtype='step')
#        sp.hist(d['SDSS_i'][:, 44], bins=100, 
#                histtype='step')
        
        sp.set_xlabel('i')
        sp.set_ylabel('counts')
#        pl.legend(loc='best')
        
def sdss_color_curve3():
    with fits.open(path) as fi:
        d = fi[1].data
        print(d['mass'][0].shape)
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(range(84), 
                    np.median(d['SDSS_r']-d['SDSS_i'], axis=0), 
                    np.std(d['SDSS_r']-d['SDSS_i'], axis=0), label='r-i',
                    capsize=2)
#        sp.hist(d['SDSS_i'][:, 44], bins=100, 
#                histtype='step')
        
        sp.set_xlabel('i')
        sp.set_ylabel('counts')
#        pl.legend(loc='best')
        
def wise_2mass_curve():
    
     with fits.open(path) as fi:
        d = fi[1].data
        
        pl.clf()
        sp = pl.subplot()
        sp.errorbar(d['age'][0]/1e9, 
                    np.median(d['WISE_W1']-d['2MASS_J'], axis=0), 
                    np.std(d['WISE_W1']-d['2MASS_J'], axis=0),
                    capsize=2)
        
        sp.set_xlabel('age [Gyr]')
        sp.set_ylabel('W1-J$')
        
def contour_sdss_2mass_wise():
     with fits.open(path) as fi:
        d = fi[1].data
        
        back = fits.open('/Users/patrickr/Downloads/sdss_wise_mass.fits')
        back = back[1].data
        
        pl.clf()
        sp = pl.subplot()
        
        sp.hist2d(back['u']-back['g'],
                  back['w1mpro']-back['r'],
                  bins=200, cmap='gray_r',
                  range=[[0.5, 2.9], [-4.5, 2.5]])
        
        gr = d['SDSS_u']-d['SDSS_g']
        print(np.min(back['g']-back['r']))
        w1j = d['WISE_W1']-d['SDSS_r']
        print(np.min(gr), np.max(gr))
        print(np.min(w1j), np.max(w1j))
        gr = gr.reshape((gr.shape[0]*gr.shape[1], ))
        w1j = w1j.reshape((w1j.shape[0]*w1j.shape[1], ))
        hist, x, y = np.histogram2d(gr,
                                    w1j,
                                    bins=50)
        x = x[:-1]+(x[1]-x[0])/2
        y = y[:-1]+(y[1]-y[0])/2
        print(x[0], x[-1])
        print(y[0], y[-1])
        hist = hist/np.sum(hist)
        print(np.max(hist))
        sp.contour(x, y, np.transpose(hist), [0.001, 0.005, 0.01, 0.05])
        
        sp.set_xlabel('g-r')
        sp.set_ylabel('W1-J')
        
#        sp.set_xlim(np.min(x), np.max(x))
#        sp.set_ylim(np.min(y), np.max(y))
    
        
path = '/Volumes/UNTITLED/hcss_colors/mist/results/full/1/sim_mags_0_5.fits'