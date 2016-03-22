import numpy as np

from scipy.ndimage.filters import median_filter, gaussian_filter
def guess_maximum(d):
    """
    This image guesses the maximum intensity of the 
    image by smoothing it with a median filter, smoothing 
    with a gaussian filter and then finding the maximum of 
    the smoothed image

    Input

    d       Numpy array that contains the values of the pixels.

    Output

    x0      x coordinate of the image maximum
    
    y0      y coordinate of the image maximum
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
    y0, x0 = np.where(d_g == np.max(d_g))

    # Take the first element. This helps in two cases: (1) only one maximum has 
    # been found, the outputs are numpy arrays and you want to extract the numbers 
    # and (2) in case there are several maximums (this has not happened so 
    # far but could in principle), pick the first of them:
    y0 = y0[0]
    x0 = x0[0]

    # Now return the values found:
    return x0,y0
