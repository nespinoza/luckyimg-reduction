# -*- coding: utf-8 -*-
import os
import numpy as np
import pyfits
import pickle
import Utils

#################### USER DEFINITIONS #######################
data_folder = '/Volumes/SeagateEHD/data/AstraLux/22122015/'
filename = 'TDRIZZLE_0010_029_HATS606007_SDSSz__000.fits'
scale = 23*1e-3 # arcsecs/pixel
N = 5 # size of the box to perform the contrast profiles
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
right_dist = int(np.floor((x0 - d.shape[0])**2 ))
left_dist = int(np.ceil(x0))
up_dist = int(np.floor((y0 - d.shape[1])**2 ))
down_dist = int(np.ceil(y0))
max_radius = np.min([right_dist,left_dist,up_dist,down_dist])

# And extract photometry around an N pixel radius around the center 
# (i.e., the target star):
flux_target,error_target = Utils.getApertureFluxes(d,x0,y0,N,0.0,1.0)

# Generate the contrast curve by, at each radius, extracting 
# photometry at 30 angles around an N-pixel radius of the 
# original image, estimate the noise of the pixels with the 
# residual image, and summing 5 times that noise to the 
# image. Then, compare the flux at the center with the flux 
# at that point.

# First, define the radii that will be explored:
radii = np.arange(0,max_radius,5)

# Initialize arrays that will save the contrast curve:
contrasts = np.zeros(len(radii))
contrasts_err = np.zeros(len(radii))

# Initialize the angles:
thetas = np.linspace(0,2*np.pi,30)
for i in range(len(radii)):
    # Generate vector that saves the contrast at a given 
    # angle:
    c_contrasts = np.zeros(len(thetas))
    for j in range (len(thetas)):
        # Get current pixel to use as center around which we will
        # extract the photometry:
        c_x = x0 + int(radii[i]*np.cos(thetas[j]))
        c_y = y0 + int(radii[i]*np.sin(thetas[j]))

        # Get NxN sub-image at the current pixel:
        c_subimg = res[c_x-(N/2)-1:c_x+(N/2),\
                       c_y-(N/2)-1:c_y+(N/2)]

        # Estimate the standard-deviation of the pixels
        # in the box:
        sigma = np.sqrt(np.var(c_subimg))

        # Extract photometry around the current center in a NxN box by adding 5-times 
        # the noise:
        flux_obj,error_obj = Utils.getApertureFluxes(d+5.*sigma,c_x,c_y,N,0.0,1.0)

        # Get the magnitude contrast between the flux of the 
        # star and flux at current position + 5-sigma the measured 
        # standard-deviation. This defines our 5-sigma contrast at this 
        # radius/angle:
        c_contrasts[j] = -2.51*np.log10(flux_target/flux_obj)

    # Take out nans:
    idx = np.where(~np.isnan(c_contrasts))[0]
    # Get contrasts as the median of all the contrasts at all angles;
    # estimate errors empirically from this distribution:
    contrasts[i] = np.median(c_contrasts[idx])
    contrasts_err[i] = np.sqrt(np.var(c_contrasts[idx]))

# Convert radii in pixels to arseconds:
radii = radii*scale

# Save results:
fout = open(out_dir+'/contrast_curve_'+filename+'.dat','w')
fout.write('# Radius ('') \t Delta mag \t Error on Delta Mag\n')
for i in range(len(radii)):
    fout.write('{0: 3.3f} \t {1: 3.3f} \t {2: 3.3f} \n'.format(radii[i],contrasts[i],contrasts_err[i]))
fout.close()

# Plot final results to the user
import matplotlib.pyplot as plt
plt.style.use('ggplot')
plt.errorbar(radii,contrasts,yerr=contrasts_err)
xlabel('Radius (arcsec)')
ylabel('Delta mag')
show()
