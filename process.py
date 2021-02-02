#!/usr/bin/env python
# coding: utf-8

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os
from scipy import optimize
from matplotlib.ticker import FormatStrFormatter


def load_fits(r, galaxy_title, line_name):

    galaxy_file = galaxy_title.replace(' ', '_')

    line_measurement = 0
    addendum = ''

    if line_name == '[OIII]':
        line_measurement = 52
        addendum = 'sw_1'

    if line_name == '[NIII]':
        line_measurement = 57
        addendum = 'sw_2'


    titlename = galaxy_title + ' ' + line_name + ' ' + str(line_measurement) + ' ' + r'$\mu$m, r = {}'.format(str(r))
    filename = galaxy_file + '_' + line_name + '_r={}'.format(str(r))
    
    loaddir = '/Users/thepoetoftwilight/Documents/SOFIA_FIFI_Cycle-8/Data/' + galaxy_file + '_' + addendum + '/'


    savedir = loaddir + 'r=' + str(r) + '/'

    fits_name = galaxy_file + '_' + addendum + '.fits'
    fits_spec_name = galaxy_file + '_' + addendum + '_r={}.fits'.format(str(r))
    
    hdulist = fits.open(loaddir + fits_name)
    hdulist_spec = fits.open(savedir + fits_spec_name)
    
    return hdulist, hdulist_spec


def renormalize_spectra(central_fluxes_unfiltered):
    
    central_fluxes_unfiltered = [flux*(10**(-26)) for flux in central_fluxes_unfiltered]

    c = 3*10**8

    for i in range(0, len(wavelengths)):

        wavelength = wavelengths[i]

        central_fluxes_unfiltered[i] = central_fluxes_unfiltered[i]*(c/(wavelength*10**(-6))**2)*(10**(-6))

    central_fluxes_unfiltered = np.array(central_fluxes_unfiltered)
    
    return central_fluxes_unfiltered


def filter_spectra(lower_ind, upper_ind, wavelengths, central_fluxes_unfiltered): 

    wavelengths_filtered = wavelengths[lower_ind:upper_ind+1]
    central_fluxes_filtered = central_fluxes_unfiltered[lower_ind:upper_ind+1]
    
    return wavelengths_filtered, central_fluxes_filtered


def identify_continuum(left_cut, right_cut, wavelengths_filtered, central_fluxes_filtered):
    
    for i in range(0, len(wavelengths_filtered)):

        wavelength = wavelengths_filtered[i]

        if(wavelength < left_cut or wavelength > right_cut):
            central_fluxes_filtered_continuum.append(central_fluxes_filtered[i])

        else:
            central_fluxes_filtered_continuum.append(float("Nan"))

    central_fluxes_filtered_continuum = np.array(central_fluxes_filtered_continuum)

    idcont = np.isfinite(central_fluxes_filtered_continuum)
    cont_params = np.polyfit(wavelengths_filtered[idcont], central_fluxes_filtered_continuum[idcont], 1)
    
    return central_fluxes_filtered_continuum, cont_params


def subtract_continuum(central_fluxes_filtered_continuum, cont_params, wavelengths_filtered, central_fluxes_filtered):
    
    cont_line_filtered = cont_params[0]*wavelengths_filtered + cont_params[1]
    
    id_act = np.isnan(central_fluxes_filtered_continuum)

    wavelengths_act = wavelengths_filtered[id_act]
    central_fluxes_filtered_act = central_fluxes_filtered[id_act]
    cont_line_filtered_act = cont_line_filtered[id_act]

    central_fluxes_continuum_subtracted_act = central_fluxes_filtered_act - cont_line_filtered_act
    
    return wavelengths_act, central_fluxes_continuum_subtracted_act


def fit_gaussian(wavelengths_act, central_fluxes_continuum_subtracted_act):
    
    heights = central_fluxes_continuum_subtracted_act/np.max(central_fluxes_continuum_subtracted_act)

    centers = wavelengths_act

    # Information about the peak in the numerical PDF
    peak_ind = np.where(heights == np.max(heights))[0][0]
    peak_height = np.max(heights)

    # mu is where the numerical PDF peaks
    mu = centers[peak_ind]

    # Estimating sigma using FWHM
    sigma = 0

    for i in range(0, peak_ind):
        if(heights[i] >= peak_height/2):
            sigma = (mu - centers[i])/np.sqrt(2*np.log(2))
            break

    # First fit a Gaussian

    guess_params = np.array([peak_height, sigma])
    fit_params, fit_covar = optimize.curve_fit(lambda centers, peak_height, sigma: 
                                               fit_func_1(centers, peak_height, mu, sigma), 
                                               centers, heights, p0=guess_params)

    fit_params = [fit_params[0], mu, fit_params[1]]

    fit_params[0] *= np.max(central_fluxes_continuum_subtracted_act)
    
    return(fit_params)


def overlay(wavelengths_filtered, fit_params):
    
    wavelength_range = np.arange(np.min(wavelengths_filtered), np.max(wavelengths_filtered), 0.001)
    cont_line_filtered = cont_params[0]*wavelength_range + cont_params[1]

    fitted_central_fluxes_act = fit_func_1(wavelength_range, *fit_params) + cont_line_filtered
    
    return wavelength_range, cont_line_filtered, fitted_central_fluxes_act


def compute_fluxes(wavelengths_filtered, wavelength_range, cont_line_filtered, central_fluxes_filtered, 
                   fitted_central_fluxes_act):
    
    flux_continuum = np.trapz(cont_line_filtered, x = wavelength_range, dx = wav_diff)

    flux_obs = np.trapz(central_fluxes_filtered, x = wavelengths_filtered, dx = wav_diff)
    flux_obs_sub = flux_obs - flux_continuum

    flux_fit = np.trapz(fitted_central_fluxes_act, x = wavelength_range, dx = wav_diff)
    flux_fit_sub = flux_fit - flux_continuum
    
    return flux_continuum, flux_obs_sub, flux_fit_sub