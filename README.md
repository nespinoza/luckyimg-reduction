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
       original image.

    3. Use this image by, at each position r and angle theta 
       in the image, estimate on a nxn box the noise (sigma). 

    4. Inject sources with the same but scaled (i.e., with different 
       magnitudes) PSFs at each position. Detect them by trying to detect 
       5 pixels above 5-sigma. The magnitude at which we start getting 4 or 
       less pixels above 5-sigma defines our contrast at that position and  
       angle.

    5. For each radius, average all the contrasts in the angles. This 
       defines our contrast curve.

The radii at which the contrast curve is sampled is defined by the effective 
FWHM, which is the median fwhm in all directions calculated numerically using the 
modelled PSF.

The code performs the modelling, and generates the contrast curves at 
the end. So far, it supports lucky imaging observations done with AstraLux Sur, 
but it should be easy to port it to other instruments.

Dependencies
------------
The code relies on the following open source libraries:

    - Numpy
    - Scipy
    - Matplotlib
    - lmfit (https://lmfit.github.io/lmfit-py/)
