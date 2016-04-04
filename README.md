# luckyimg-reduction

Lucky Imaging Reduction (contrast curve generator)
--------------------------------------------------
This code takes lucky imaging observations and generates 
contrast curves around the brightest star on the image 
(useful for exoplanet studies).

The steps are:

    1. Model the brightest star as a weighted sum of a 
       rotated, assymetric gaussian and a Moffat profile.

    2. Generate a residual image between this model and the 
       original image which will be used to estimate the 
       noise at each pixel (in order to include systematic and 
       statistical uncertainties).

    3. With this noise estimation, extract photometry at the 
       center and at different radii in a N pixel radius and 
       use the flux at the center and the flux at each radii 
       plus 5-sigma in order to generate a threshold function, 
       which will aid in defining the flux ratio at which one 
       considers a detection.

    4. Inject fake signals and recover them using the threshold 
       function. The largest contrast at which one can consider 
       a detection is the 5-sigma contrast limit.

The code performs the modelling, and generates the threshold 
function and contrast curve at the end. So far, it supports 
lucky imaging observations done with AstraLux Sur, but it should 
be easy to port it to other instruments.

Dependencies
------------
The code relies on the following open source libraries:

    - Numpy
    - Scipy
    - Matplotlib
    - lmfit (https://lmfit.github.io/lmfit-py/)
    - photutils (https://photutils.readthedocs.org/en/latest/)
