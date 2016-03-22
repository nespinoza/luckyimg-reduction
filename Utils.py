# -*- coding: utf-8 -*-
import numpy as np

def get_sigma_mad(x):
    med = np.median(x)
    return 1.4826*np.median(np.abs(x-med))

from scipy.ndimage.filters import median_filter, gaussian_filter
def guess_gaussian_parameters(d):
    """
    This image guesses the maximum intensity of the 
    image by obtaining a smoothed version of the image 
    via median + gaussian filtering. Then, finds the 
    maximum of the smoothed image. Using the median-filtered 
    image, the algorithm also estimates the width of the PSF and 
    with this value and an estimate of the volume, the amplitude of 
    such gaussian.

    Input

    d       Numpy array that contains the values of the pixels.

    Output

    x0      x coordinate of the image maximum
    
    y0      y coordinate of the image maximum

    sigma   Width of the PSF

    A       Estimated amplitude of the gaussian function
    """

    # First, smooth the image with a median filter. For this, find 
    # the optimal window as the square-root of the geometric mean 
    # between the sizes of the array. This step is good to kill outliers:
    window = int(np.sqrt(np.sqrt(d.shape[0]*d.shape[1])))
    if window % 2 == 0:
        window = window + 1 
    d_mfiltered = median_filter(d,size = window)

    # Next, smooth it with a gaussian filter:
    d_gfiltered = gaussian_filter(d_mfiltered,sigma = window)

    # Now, find the maximum of this image:
    y0, x0 = np.where(d_gfiltered == np.max(d_gfiltered))

    # Take the first element. This helps in two cases: (1) only one maximum has 
    # been found, the outputs are numpy arrays and you want to extract the numbers 
    # and (2) in case there are several maximums (this has not happened so 
    # far but could in principle), pick the first of them:
    y0 = y0[0]
    x0 = x0[0]

    # Now estimate the width of the PSF by taking a "cross-plot" using the 
    # maximum values found:
    x_cut = d[:,int(x0)]
    sigma_x = (np.sum(x_cut*(np.abs(np.arange(len(x_cut))-y0)))/np.sum(x_cut))/3.
    y_cut = d[int(y0),:]
    sigma_y = (np.sum(y_cut*(np.abs(np.arange(len(y_cut))-x0)))/np.sum(y_cut))/3.
    sigma = np.sqrt(sigma_x*sigma_y)

    # (Under) estimate amplitude assuming a gaussian function:
    A = (np.sum(d-np.median(d))/(2.*np.pi*sigma**2))

    return x0,y0,sigma,2.*A

def moffat(x,y,A,x0,y0,sigma,beta):
    first_term = ((x-x0)**2 + (y-y0)**2)/sigma**2
    return A*(1. + first_term)**(-beta)

def assymetric_gaussian(x, y, A, x0, y0, sigma_x, sigma_y, theta):
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    return A*np.exp( - (a*((x-x0)**2) + 2*b*(x-x0)*(y-y0) +\
                        c*((y-y0)**2)))

def gaussian(x, y, A, x0, y0, sigma):
    a = ((x-x0)**2 + (y-y0)**2)/(2.*(sigma**2))
    return A*np.exp(-a)

import lmfit

def modelPSF(params,mesh):
    W = params['W'].value
    ag = (1.-W)*assymetric_gaussian(mesh[0],mesh[1],params['Ag'].value,params['x0'].value,\
                        params['y0'].value,params['sigma_x'].value,params['sigma_y'].value,\
                        params['theta'].value)
    mof = W*moffat(mesh[0],mesh[1],params['Am'].value,params['x0'].value,\
                        params['y0'].value,params['sigma_m'].value,params['beta'].value)
    return ag+mof+params['bkg'].value

def fitPSF(d,x0,y0,sigma,A):
    """
    This function fits a sum of a rotated gaussian plus a 
    moffat profile.
    """
    mesh = np.meshgrid(np.arange(d.shape[0]),np.arange(d.shape[1]))
    flattened_d = d.flatten()

    def residuals(params):
        return flattened_d - (modelPSF(params,mesh)).flatten()

    prms = lmfit.Parameters()
    prms.add('x0', value = x0, min = 0, max = d.shape[0], vary = True)
    prms.add('y0', value = y0, min = 0, max = d.shape[1], vary = True)
    prms.add('W' , value = 0.5, min = 0, max = 1., vary = True)
    prms.add('Ag', value = A, min = 0, max = np.sum(d-np.median(d)), vary = True)
    prms.add('Am', value = A, min = 0, max = np.sum(d-np.median(d)), vary = True)
    prms.add('sigma_x', value = sigma, min = 0, max = d.shape[0]/3., vary = True)
    prms.add('sigma_y', value = sigma, min = 0, max = d.shape[1]/3., vary = True)
    prms.add('sigma_m', value = sigma, min = 0, max = d.shape[1]/3., vary = True)
    prms.add('beta', value = 1., min = 0, max = 10.)
    prms.add('bkg', value = np.median(d), min = np.median(d)-10*get_sigma_mad(d), \
                    max = np.median(d)+10*get_sigma_mad(d), vary = True)
    prms.add('theta', value = np.pi/4., min = 0, max = np.pi/2.)
    result = lmfit.minimize(residuals, prms)
    return result.params

from photutils import CircularAperture,CircularAnnulus,aperture_photometry
def getApertureFluxes(subimg,x_cen,y_cen,Radius,sky_sigma,GAIN):
    apertures = CircularAperture([(x_cen,y_cen)],r=Radius)
    rawflux_table = aperture_photometry(subimg, apertures, \
            error=sky_sigma, effective_gain = GAIN)
    return rawflux_table['aperture_sum'][0],rawflux_table['aperture_sum_err'][0]    
