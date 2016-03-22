import numpy as np

def get_mad(x):
    med = np.median(x)
    return np.median(np.abs(x-med))

from scipy.ndimage.filters import median_filter, gaussian_filter
def guess_gaussian_parameters(d):
    """
    This image guesses the maximum intensity of the 
    image by obtaining a smoothed version of the image 
    via median + gaussian filtering. Then, finds the 
    maximum of the smoothed image. Using the median-filtered 
    image, the algorithm also estimates the width of the PSF.

    Input

    d       Numpy array that contains the values of the pixels.

    Output

    x0      x coordinate of the image maximum
    
    y0      y coordinate of the image maximum

    sigma   Width of the PSF
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

    return x0,y0,np.sqrt(sigma_x*sigma_y)
