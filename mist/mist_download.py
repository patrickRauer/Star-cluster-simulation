#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 10 07:27:20 2018

Downloads the mist interpolated isochronoes
@author: patrickr
"""
from splinter import Browser
from astropy.table import Table
from threading import Thread
from selenium.common.exceptions import StaleElementReferenceException, JavascriptException
import numpy as np
import shutil
import urllib
import time
import zipfile
import os
from . import cols

# converter form standard names to the MIST names
survey_conv = {'SDSS': 'SDSSugriz',
               'MegaCam': 'CFHTugriz',
               'DECam': 'DECam',
               'GALEX': 'GALEX',
               'PanSTARRS': 'PanSTARRS',
               '2MASS': 'UBVRIplus',
               'UKIDSS': 'UKIDSS',
               'WISE': 'WISE'}


def remove_files(path):
    """
    Deletes not needed files

    :param path: Path to the directory which is to check
    :type path: str
    :return:
    """
    for f in os.listdir(path):
        if 'cmd' not in f:
            os.remove(path+f)


def rename_columns(tab, survey):
    """
    Renames the table column

    :param tab: The table with the columns
    :type tab: astropy.table.Table
    :param survey: The name of the survey
    :type survey: str
    :return: The table with the modified column names
    :rtype: astropy.table.Table
    """
    for i, c in enumerate(tab.colnames):
        tab.rename_column(c, cols.survey_cols[survey][i])
    return tab


def convert_files(survey, path):
    """
    Converts the MIST data from an ascii file to a fits file

    :param survey: THe name of the used survey
    :type survey: str
    :param path: The path to the directory
    :type path: str
    :return:
    """
    for f in os.listdir(path):
        if 'cmd' in f and '._' not in f and 'fits' not in f:
            tab = Table.read(path+f, format='ascii')
            tab = rename_columns(tab, survey)
            tab.write(path+f.split('cmd')[0]+'fits', overwrite=True)
            os.remove(path+f)


def copy_files(directory, temp_path):
    """
    Copies the data from one directory to the new one and removes to old data

    :param directory: The final directory
    :type directory: str
    :param temp_path: The original directory
    :type temp_path: str
    :return:
    """
    path = temp_path
    if not os.path.exists(directory):
        os.makedirs(directory)
    for f in os.listdir(path):
        if 'fits' in f and '._' not in f:
            shutil.copy2(path+f, directory+f)
            os.remove(path+f)


def get_mist_data(browser, fe_h, min_mass, max_mass, mass_step, survey,
                  directory):
    """
    Download a MIST data set
    :param browser: Browser object which is used as the interface to the MIST web-page
    :param fe_h: The metallicity of the stars in Fe/H
    :type fe_h: float
    :param min_mass:
        The minimal mass of a star in solar masses.
        Is the mass smaller than 0.1 it will set to 0.1.
    :type min_mass: float
    :param max_mass: The maximal mass of a star in solar masses.
    :type max_mass: float
    :param mass_step:
        The mass resolution in solar masses, the minimal resolution is :math:`10^{-5} M_\odot`
    :type mass_step: float
    :param survey: The name of the survey
    :type survey: str
    :param directory: The path to the directory where the data will be stored
    :type directory: str
    :return:
    """
    # Visit URL
    url = "http://waps.cfa.harvard.edu/MIST/interp_tracks.html"
    browser.visit(url)
    # fill the form on the web-page
    browser.fill('mass_range_low', str(min_mass))
    browser.fill('mass_range_high', str(max_mass))
    browser.fill('mass_range_delta', str(mass_step))
    browser.fill('new_met_value', str(fe_h))
    browser.choose('mass_type', 'range')
    browser.choose('output_option', 'photometry')
    browser.select('output', survey_conv[survey])

    # click the button the start the server side process
    browser.find_by_tag('button')[0].click()
    link = browser.find_by_tag('a')[0]

    link = link['href']
    wait = True

    time_start = 0
    # wait until the server is done
    while wait:
        try:
            # extract the download link
            link = browser.find_by_tag('a')[0]
            link = link['href']
            # if the link includes zip the server is done and then stop the loop
            if 'zip' in link:
                wait = False
                break

        except StaleElementReferenceException:
            pass
        except JavascriptException:
            print('wait')
        time.sleep(1)
        time_start += 1
        # if the loop runs too long
        if time_start > 1500:
            # raise an error
            raise TimeoutError('Something went wrong! Stop the waiting loop after 1500secs.')

    temp_name = 'temp_{}'.format(str(fe_h))

    # download the data
    urllib.urlretrieve(str(link), temp_name+".zip")

    # extract the data
    zip_ref = zipfile.ZipFile(temp_name+".zip", 'r')
    if not os.path.exists('./{}/'.format(temp_name)):
        os.makedirs(os.path.abspath('./{}/'.format(temp_name)))
    zip_ref.extractall('./{}/'.format(temp_name))
    zip_ref.close()
    remove_files('./{}/'.format(temp_name))

    # convert the data to fits
    convert_files(survey, './{}/'.format(temp_name))

    fe_h_str = str(fe_h).split('.')
    # copy the data to the final place
    copy_files(directory + '{}_{}/'.format(fe_h_str[0], fe_h_str[1]) +
               survey + '/', './{}/'.format(temp_name))


def download_survey_metallicity(fe_h,
                                min_mass, max_mass, mass_step,
                                survey, directory):
    """
    Downloads the MIST data for a specific metallicity and a specific survey
    :param fe_h: The metallicity of the stars in Fe/H
    :type fe_h: float
    :param min_mass:
        The minimal mass of a star in solar masses.
        Is the mass smaller than 0.1 it will set to 0.1.
    :type min_mass: float
    :param max_mass: The maximal mass of a star in solar masses.
    :type max_mass: float
    :param mass_step:
        The mass resolution in solar masses, the minimal resolution is :math:`10^{-5} M_\odot`
    :type mass_step: float
    :param survey: The name of the survey
    :type survey: str
    :param directory: The path to the directory where the data will be stored
    :type directory: str
    :return:
    """
    # generate the mass steps for the download
    masses = np.arange(min_mass, max_mass, mass_step*10)
    print_str = 'download: Fe/H={} min_mass={} max_mass={} survey={}'

    # open a new browser
    browser = Browser(executable_path='/Users/patrickr/Documents/add/geckodriver')
    for m1, m2 in zip(masses[:-1], masses[1:]):
        print(print_str.format(fe_h, m1, m2, survey))
        try:
            # download the isochronos for the current mass range
            get_mist_data(browser, fe_h, m1, m2, mass_step, survey,
                          directory)
        except ValueError as e:
            print(e)
    browser.quit()


def download_survey(min_mass, max_mass, mass_step, metallicities,
                    survey, directory):
    """
    Downloads all the interpolated isochrones for a specific survey.
    Every metallicity will be downloaded in separate process to increase the
    total speed of the download because the slowest part here is the server response.

    :param min_mass:
        The minimal mass of a star in solar masses.
        Is the mass smaller than 0.1 it will set to 0.1.
    :type min_mass: float
    :param max_mass: The maximal mass of a star in solar masses.
    :type max_mass: float
    :param mass_step:
        The mass resolution in solar masses, the minimal resolution is :math:`10^{-5} M_\odot`
    :type mass_step: float
    :param metallicities: The metallicity of the stars in Fe/H
    :type metallicities: list
    :param survey: The name of the survey
    :type survey: str
    :param directory: The path to the directory where the data will be stored
    :type directory: str
    :return:
    """
    th = []
    for fe_h in metallicities:
        t = Thread(target=download_survey_metallicity,
                   args=(fe_h,
                         min_mass, max_mass, mass_step,
                         survey, directory))
        t.start()
        th.append(t)

    for t in th:
        t.join()


def download_all(min_mass, max_mass, mass_step, metallicities,
                 directory):
    """
    Downloads the MIST data for all surveys

    :param min_mass:
        The minimal mass of a star in solar masses.
        Is the mass smaller than 0.1 it will set to 0.1.
    :type min_mass: float
    :param max_mass: The maximal mass of a star in solar masses.
    :type max_mass: float
    :param mass_step:
        The mass resolution in solar masses, the minimal resolution is :math:`10^{-5} M_\odot`
    :type mass_step: float
    :param metallicities: The metallicity of the stars in Fe/H
    :type metallicities: list
    :param directory: The path to the directory where the data will be stored
    :type directory: str
    :return:
    """
    if min_mass < 0.1:
        min_mass = 0.1
    for survey in list(survey_conv.keys()):
        download_survey(min_mass, max_mass, mass_step, metallicities,
                        survey, directory)
