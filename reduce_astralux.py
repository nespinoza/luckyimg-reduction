# -*- coding: utf-8 -*-
import os
import numpy as np
import pyfits
import pickle
import Utils

#################### USER DEFINITIONS #######################
data_folder = '/Volumes/SeagateEHD/data/AstraLux/22122015/'
filename = 'TDRIZZLE_0010_006_HATS754005_SDSSz__000'
# Scale of the image in arcsecs/pixel:
scale = 23*1e-3 
# Size of the box that will be used at each point to estimate 
# the noise. Also is the radius of the circle used to get the 
# photometry at that point:
N = 5 
#############################################################

# Create output directory if non-existent for the current image:
out_dir = filename.split('.')[0]
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# If not already done, model the input image. If already done, 
# obtain saved data:
if not os.path.exists(out_dir+'/model_image.fits'):
    # First, extract image data:
    d,h = pyfits.getdata(data_folder+filename, header=True)

    # Guess centroid by maximum intensity; also estimate approximate 
    # width of the PSF by weighted median-absolute distance from the 
    # estimated center and use it to estimate amplitude:
    x0,y0,sigma,A = Utils.guess_gaussian_parameters(d)

    # Estimate model of moffat + rotated gaussian:
    out_params = Utils.fitPSF(d,x0,y0,sigma,A)

    # Save output parameters:
    fout = open(out_dir+'/out_params.pkl','wb')
    pickle.dump(out_params,fout)
    fout.close()

    # Generate model image:
    model = Utils.modelPSF(out_params,\
                           np.meshgrid(np.arange(d.shape[0]),np.arange(d.shape[1])))

    # Generate residual image:
    res = model - d

    # Save images:
    pyfits.PrimaryHDU(model).writeto(out_dir+'/model_image.fits')
    pyfits.PrimaryHDU(d).writeto(out_dir+'/original_image.fits')
    pyfits.PrimaryHDU(res).writeto(out_dir+'/residual_image.fits')
else:
    # If everything already done, read data:
    model = pyfits.getdata(out_dir+'/model_image.fits')
    d = pyfits.getdata(out_dir+'/original_image.fits')
    res = pyfits.getdata(out_dir+'/residual_image.fits')
    par = open(out_dir+'/out_params.pkl','r')
    out_params = pickle.load(par)
    par.close()

# Get centroids:
x0,y0 = out_params['x0'].value,out_params['y0'].value

# Remove estimated background: 
d = d - out_params['bkg']

# Now generate 5-sigma contrast curves. For this, first find 
# closest distance to edges of the image:
right_dist = int(np.floor(np.abs(x0 - d.shape[0])))-N
left_dist = int(np.ceil(x0))-N
up_dist = int(np.floor(np.abs(y0 - d.shape[1])))-N
down_dist = int(np.ceil(y0))-N
max_radius = np.min([right_dist,left_dist,up_dist,down_dist])

# And extract photometry around an N pixel radius around the center 
# (i.e., the target star):
flux_target,error_target = Utils.getApertureFluxes(d,x0,y0,N,0.0,1.0)

# Generate the contrast curve by first generating the Threshold Function.
# This is generated at each radius by extracting 
# photometry at 30 angles around an n-pixel radius of the 
# original image. The noise of the pixels is estimated with the 
# residual image, and the 5-sigma flux at that position is obtained 
# as the flux + 5sigma. Then, compare the total flux at the center 
# with the total (5-sigma) flux at that point.

# First, define the radii that will be explored:
radii = np.arange(0,max_radius,10)

# Initialize arrays that will save the threshold functions:
T = np.zeros(len(radii))
T_err = np.zeros(len(radii))

# Initialize arrays that will save the contrast curves:
contrast = np.zeros(len(radii))
contrast_err = np.zeros(len(radii))

# Initialize magnitude contrasts to be explored:
possible_contrasts = np.linspace(0,10,100)

# Initialize the angles:
thetas = np.linspace(0,2*np.pi,27)
for i in range(len(radii)):
    print i,'of',len(radii),'(',radii[i],')'
    # Generate vector that saves the threshold functions 
    # at a given angle:
    c_T = np.zeros(len(thetas))
    for j in range (len(thetas)):
        # Get current pixel to use as center around which we will
        # extract the photometry:
        c_x = x0 + int(radii[i]*np.cos(thetas[j]))
        c_y = y0 + int(radii[i]*np.sin(thetas[j]))

        # Get nxn sub-image at the current pixel:
        c_subimg = res[c_x-(N/2)-1:c_x+(N/2),\
                       c_y-(N/2)-1:c_y+(N/2)]

        # Estimate the (empirical) standard-deviation of the pixels
        # in the box:
        sigma = np.sqrt(np.var(c_subimg))

        # Extract the total (5-sigma) flux around the current center 
        # in a N pixel radius:
        flux_obj,error_obj = Utils.getApertureFluxes(d+5.*sigma,c_x,c_y,N,0.0,1.0)

        # Get the flux ratio between the flux of the 
        # star and flux at current position + 5-sigma the measured 
        # standard-deviation. This defines our threshold function:
        c_T[j] = flux_obj/flux_target

    # Take out nans:
    idx = np.where(~np.isnan(c_T))[0]
    # Get threshold functions as the median of all the thresholds at all angles;
    # estimate errors empirically from this distribution:
    T[i] = np.median(c_T[idx])
    T_err[i] = np.sqrt(np.var(c_T[idx])*(len(idx))/np.double(len(idx)-1.))

    # Now inject fake source on the images at each position, and see when we 
    # recover it. First, set background of the model to zero:
    out_params['bkg'].value = 0.0
    
    # Initialize the array that will save the contrasts at all angles:
    c_contrast = np.zeros(len(thetas))
    for j in range(len(thetas)):
        # Get current position:    
        c_x = x0 + int(radii[i]*np.cos(thetas[j]))
        c_y = y0 + int(radii[i]*np.sin(thetas[j]))
        # Set the model PSF at the center of the current position:
        out_params['x0'].value = c_x
        out_params['y0'].value = c_y
        # Generate the fake source. We will scale it below to match 
        # different contrasts:
        fake_signal = Utils.modelPSF(out_params,\
                      np.meshgrid(np.arange(d.shape[0]),np.arange(d.shape[1])))
        for k in range(len(possible_contrasts)):
            # Generate the scaling factor:
            scaling_factor = 10**(possible_contrasts[k]/2.51)
            # Construct fake image:
            fake_image = d + (fake_signal/scaling_factor)
            # Extract aperture photometry:
            flux_obj,error_obj = Utils.getApertureFluxes(fake_image,c_x,c_y,N,0.0,1.0)
            # If the contrast of the fake image is lower than the threshold limit, then 
            # the contrast before this one was the detection limit. Save it, break out of 
            # this loop, and go for the next position (next angle):
            if (flux_obj/flux_target < T[i]):
                if k != 0:
                    c_contrast[j] = possible_contrasts[k-1]      
                else:
                    c_contrast[j] = 0.0
                break
    idx = np.where(~np.isnan(c_contrast))[0]
    contrast[i] = np.median(c_contrast[idx])
    contrast_err[i] = np.sqrt(np.var(c_contrast[idx])*len(idx)/np.double(len(idx)-1.))

# Convert radii in pixels to arseconds:
radii = radii*scale

# Save results:
fout = open(out_dir+'/threshold_function_'+filename+'.dat','w')
fout.write('# Radius ('') \t T_5s(Radius) \t Error\n')
for i in range(len(radii)):
    fout.write('{0: 3.3f} \t {1: 3.3f} \t {2: 3.3f} \n'.format(radii[i],T[i],T_err[i]))
fout.close()

# Save results:
fout = open(out_dir+'/contrast_curve_'+filename+'.dat','w')
fout.write('# Radius ('') \t Magnitude Contrast \t Error\n')
for i in range(len(radii)):
            fout.write('{0: 3.3f} \t {1: 3.3f} \t {2: 3.3f} \n'.format(radii[i],\
                                                    contrast[i],contrast_err[i]))
fout.close()

# Plot final results to the user
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.errorbar(radii,T,yerr=T_err)
plt.title('Threshold function')
plt.xlabel(r'$r$ (arcsec)')
plt.ylabel(r'$T_{5\sigma}(r)$')
plt.show()

plt.errorbar(radii,contrast,yerr=contrast_err)
plt.title('Magnitude contrast')
plt.xlabel(r'$r$ (arcsec)')
plt.ylabel(r'$\Delta m$')
plt.show()
