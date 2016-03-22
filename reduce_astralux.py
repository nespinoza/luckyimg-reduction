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

x0,y0 = out_params['x0'].value,out_params['y0'].value

# Now generate 5-sigma contrast curves. For this, first find 
# closest distance to edges of the image:
right_dist = int(np.floor((x0 - d.shape[0])**2 ))
left_dist = int(np.ceil(x0))
up_dist = int(np.floor((y0 - d.shape[1])**2 ))
down_dist = int(np.ceil(y0))
max_radius = np.min([right_dist,left_dist,up_dist,down_dist])

# Generate N-pixel subimg around the center of the original image:
original_subimg = d[int(x0)-(N/2)-1:int(x0)+(N/2),\
                    int(y0)-(N/2)-1:int(y0)+(N/2)]

# Get the total flux in this box:
stellar_flux = np.sum(original_subimg)

# Generate the contrast curve by, at each radius, generating 30 
# NxN subimg (one at each of the 30 sampled angles around the 
# estimated centroid of the image) using the residual from the 
# model fit, and comparing it with the original_subimg above, 
# generated with the original image:
radii = np.arange(1,max_radius)
contrasts = np.zeros(len(radii))
contrasts_err = np.zeros(len(radii))
thetas = np.linspace(0,2*np.pi,30)
for i in range(len(radii)):
    c_contrasts = np.zeros(len(thetas))
    for j in range (len(thetas)):
        # Get current pixel to use as center:    
        c_x = x0 + int(radii[i]*np.cos(thetas[j]))
        c_y = y0 + int(radii[i]*np.sin(thetas[j]))
        # Get NxN sub-image at the current pixel:
        c_subimg = res[c_x-(N/2)-1:c_x+(N/2),\
                       c_y-(N/2)-1:c_y+(N/2)]

        # Estimate the standard-deviation of the pixels
        # in the box:
        sigma = np.sqrt(np.var(c_subimg))

        # Get the magnitude contrast between the flux of the 
        # star and 5-sigma the measured standard-deviation. This 
        # defines our 5-sigma contrast at this pixel:
        c_contrasts[j] = -2.51*np.log10(stellar_flux/(5.*sigma))
    contrasts[i] = np.median(c_contrasts)
    contrasts_err[i] = np.sqrt(np.var(c_contrasts))

# Convert radii in pixels to arseconds:
radii = radii*scale

# Save results:
fout = open(out_dir+'/contrast_curve_'+filename+'.dat','w')
fout.write('# Radius ('') \t Delta mag \t Error on Delta Mag\n')
print len(contrasts_err),len(radii)
for i in range(len(radii)):
    print len(contrasts_err)
    print i,radii[i],contrasts[i],contrasts_err[i]
    fout.write('{0: 3.3f} \t {1: 3.3f} \t {2: 3.3f} \n'.format(radii[i],contrasts[i],contrasts_err[i]))
fout.close()
