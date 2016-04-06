# -*- coding: utf-8 -*-
import os
import numpy as np
import pyfits
import pickle
import Utils

#################### USER DEFINITIONS #######################

data_folder = '/Volumes/SeagateEHD/data/AstraLux/22122015/'
filename = 'TDRIZZLE_0010_029_HATS606007_SDSSz__000'
# Minimum magnitude contrast to be explored:
min_m = 0
# Maximum magnitude contrast to be explored:
max_m = 10
# Contrast steps:
contrast_steps = 100
# Scale of the image in arcsecs/pixel:
scale = 23*1e-3 

#############################################################
print '\n\t     AstraLux contrast curve generator'
print '\t-----------------------------------------------'
print '\tAuthors: Nestor Espinoza (nespino@astro.puc.cl)'
print '\t         Andres Jordan (ajordan@astro.puc.cl)\n'
# Create output directory if non-existent for the current image:
out_dir = filename.split('.')[0]
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

# If not already done, model the input image. If already done, 
# obtain saved data:
if not os.path.exists(out_dir+'/model_image.fits'):
    print '\t > Modelling the PSF...'
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
    print '\t > Saving results...'
    pyfits.PrimaryHDU(model).writeto(out_dir+'/model_image.fits')
    pyfits.PrimaryHDU(d).writeto(out_dir+'/original_image.fits')
    pyfits.PrimaryHDU(res).writeto(out_dir+'/residual_image.fits')
else:
    print '\t > PSF already modelled. Extracting data...'
    # If everything already done, read data:
    model = pyfits.getdata(out_dir+'/model_image.fits')
    d = pyfits.getdata(out_dir+'/original_image.fits')
    res = pyfits.getdata(out_dir+'/residual_image.fits')
    par = open(out_dir+'/out_params.pkl','r')
    out_params = pickle.load(par)
    par.close()

# Define the step in radius at which we will calculate the contrasts. This is 
# calculated in terms of the "effective FWHM", which we calculate numerically from 
# the model PSF, by trying different radii and angles and finding the positions at which 
# the flux is half the peak flux.
max_flux_model = np.max(model) 
radii = np.linspace(0,5.*((out_params['sigma_x'].value+out_params['sigma_y'].value)/2.),100)
thetas = np.linspace(0,2*np.pi,100)
fwhms = np.zeros(len(thetas))
for j in range(len(thetas)):
    for i in range(len(radii)):
        c_x = out_params['x0'].value + radii[i]*np.cos(thetas[j])
        c_y = out_params['y0'].value + radii[i]*np.sin(thetas[j])
        if model[int(c_y),int(c_x)]<max_flux_model/2.:
           fwhms[j] = np.sqrt((out_params['x0'].value-int(c_x))**2 + (out_params['y0'].value-int(c_y))**2)
           break

print '\t Effective FWHM:',np.median(fwhms),'+-',np.sqrt(np.var(fwhms)),' pix (',np.median(fwhms)*scale,'+-',np.sqrt(np.var(fwhms))*scale,' arcsecs)'
radii_step = np.median(fwhms)
N = np.median(fwhms)

# Convert the radius step to int:
radii_step = int(radii_step)

# Get centroids:
x0,y0 = out_params['x0'].value,out_params['y0'].value

# Remove estimated background from original image: 
d = d - out_params['bkg'].value

# Now generate 5-sigma contrast curves. For this, first find 
# closest distance to edges of the image:
right_dist = int(np.floor(np.abs(x0 - d.shape[0])))-N
left_dist = int(np.ceil(x0))-N
up_dist = int(np.floor(np.abs(y0 - d.shape[1])))-N
down_dist = int(np.ceil(y0))-N
max_radius = np.min([right_dist,left_dist,up_dist,down_dist])

# Generate the contrast curve by, for a given radius and using the 
# residual image, injecting fake sources with the same parameters as
# the fitted PSF but scaled in order to see at what scale (i.e., magnitude)
# we detect the injected signal at 5-sigma). A detection is defined if 
# more than 5 pixels are above the 5-sigma noise level of this residual 
# image at that position.

# First, define the radii that will be explored:
radii = np.arange(radii_step,max_radius,radii_step)

# Initialize arrays that will save the contrast curves:
contrast = np.zeros(len(radii))
contrast_err = np.zeros(len(radii))

# Initialize magnitude contrasts to be explored:
possible_contrasts = np.linspace(min_m,max_m,contrast_steps)

# Now inject fake source on the images at each position, and see when we 
# recover it. First, set background of the model to zero:
out_params['bkg'].value = 0.0

print '\t > Generating contrast curves...'
for i in range(len(radii)):
    # Define the number of angles that, given the radius, have 
    # independant information:
    if radii[i] != 0:
        n_thetas = np.min([int((2.*np.pi)/((2.*N)/radii[i])),30])
        thetas = np.linspace(0,2*np.pi,n_thetas)
    # Generate vector that saves the threshold functions 
    # at a given angle:
    c_contrast = np.zeros(len(thetas))
    for j in range (len(thetas)):
        # Get current pixel to use as center around which we will
        # extract the photometry:
        c_x = x0 + int(np.round(radii[i]*np.cos(thetas[j])))
        c_y = y0 + int(np.round(radii[i]*np.sin(thetas[j])))

        # Get nxn sub-image at the current pixel:
        c_subimg = res[c_y-(N/2)-1:c_y+(N/2),\
                       c_x-(N/2)-1:c_x+(N/2)]

        # Estimate the (empirical) standard-deviation of the pixels
        # in the box:
        sigma = np.sqrt(np.var(c_subimg))

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
            # Construct fake image on top of the residual image, cut the portion under
            # analysis:
            fake_image = (res + (fake_signal/scaling_factor))[c_y-(N/2)-1:c_y+(N/2),\
                                                              c_x-(N/2)-1:c_x+(N/2)]
            # If our detection limit (i.e., 5 pixels or more are above 5-sigma) is not accomplished,
            # then the source cannot be detected and this defines our 5-sigma contrast:
            if (len(np.where(fake_image>5*sigma)[0])<5):
                if k != 0:
                    c_contrast[j] = possible_contrasts[k-1]      
                else:
                    c_contrast[j] = 0.0
                break

    idx = np.where((~np.isnan(c_contrast))&(c_contrast!=0.0))[0]
    contrast[i] = np.median(c_contrast[idx])
    contrast_err[i] = np.sqrt(np.var(c_contrast[idx])*len(idx)/np.double(len(idx)-1.))

# Convert radii in pixels to arseconds:
radii = radii*scale

# Save results:
fout = open(out_dir+'/contrast_curve_'+filename+'.dat','w')
fout.write('# Radius ('') \t Magnitude Contrast \t Error\n')
for i in range(len(radii)):
            fout.write('{0: 3.3f} \t {1: 3.3f} \t {2: 3.3f} \n'.format(radii[i],\
                                                    contrast[i],contrast_err[i]))
fout.close()

# Plot final results to the user
import matplotlib.pyplot as plt
plt.errorbar(radii,contrast,yerr=contrast_err)
plt.title('Magnitude contrast')
plt.xlabel(r'$r$ (arcsec)')
plt.ylabel(r'$\Delta m$')
plt.gca().invert_yaxis()
plt.show()
