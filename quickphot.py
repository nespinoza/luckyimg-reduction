# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pyfits
import glob
import os
from photutils import CircularAperture,CircularAnnulus,aperture_photometry

EGAIN = 1.0 # electrons/ADU
ERON = 3.0 # electrons

def getApertureFluxes(subimg,x_cen,y_cen,Radius,sky_sigma,GAIN):
    apertures = CircularAperture([(x_cen,y_cen)],r=Radius)
    rawflux_table = aperture_photometry(subimg, apertures, \
            error=sky_sigma, effective_gain = GAIN)
    return rawflux_table['aperture_sum'][0],rawflux_table['aperture_sum_err'][0]

from scipy.ndimage.filters import gaussian_filter
def get_centroid(original_subimg):
    subimg = gaussian_filter(original_subimg,10)
    total = np.sum(subimg)
    x0 = np.sum(np.sum(subimg,axis=0)*np.arange(subimg.shape[0]))/np.sum(np.sum(subimg,axis=0))
    y0 = np.sum(np.sum(subimg,axis=1)*np.arange(subimg.shape[1]))/np.sum(np.sum(subimg,axis=1))
    return x0,y0

def read_dfile(fname):
    ff = open(fname,'r')
    im = []
    t = []
    f = []
    ferr = []
    while True:
        line = ff.readline()
        if line != '':
            d1,d2,d3,d4 = line.split()
            im.append(d1)
            t.append(np.double(d2))
            f.append(np.double(d3))
            ferr.append(np.double(d4))
        else:
            break
    ff.close()
    return im,t,f,ferr

def bin_data(t,rf_mag,rf_mag_err,n_bin = 10):
    times_bins = []
    mags_bins = []
    errors_bins = []
    for i in range(0,len(t),n_bin):
        times_bins.append(np.median(t[i:i+n_bin-1]))
        mags_bins.append(np.median(rf_mag[i:i+n_bin-1]))
        errors_bins.append(np.sqrt(np.sum(rf_mag_err[i:i+n_bin-1]**2))/np.double(n_bin))
    return times_bins,mags_bins,errors_bins

#############################################
filename = 'out_qphot.dat'
images = ['original_image.fits']
target_center = [251.,236.]
comp_center = [199.976,415.976]#[1080,660]
first_time = True
hs = 80
ap_rad = 78
user_filt = 'i'
#############################################

print 'Distance (pixels):',np.sqrt((target_center[0]-comp_center[0])**2+\
                                   (target_center[1]-comp_center[1])**2)

print 'Distance (arcsec):',23*1e-3*np.sqrt((target_center[0]-comp_center[0])**2+\
                                                   (target_center[1]-comp_center[1])**2)
if not first_time:
   all_images,times,rel_fluxes,rel_fluxes_err = read_dfile(filename)
else:
    all_images = []
    times = []
    rel_fluxes = []
    rel_fluxes_err = []

if not os.path.exists('phot_images_target'):
   os.mkdir('phot_images_target')

if not os.path.exists('phot_images_comp'):
   os.mkdir('phot_images_comp')

for im in images:
    d,h = pyfits.getdata(im,header=True)
    if True:#(h['FILTER'] == user_filt) and (im not in all_images):
        all_images.append(im)
        time = 0.00#(hh + (mm/60.) + (ss/3600.))-24.
        subimg = d[int(target_center[0]-hs):int(target_center[0]+hs),int(target_center[1]-hs):int(target_center[1]+hs)]
        sky = np.median(subimg)
        subimg = subimg - sky
        imm = plt.imshow(subimg)
        imm.set_clim(0,10)
        x0,y0 = hs,hs#get_centroid(subimg)
        plt.plot(x0,y0,'wx',markersize=15,alpha=0.5)
        circle = plt.Circle((x0,y0),ap_rad,color='white',fill=False)
        plt.gca().add_artist(circle)
        plt.savefig('phot_images_target/'+im+'.png')
        plt.close()
        f,ferr = getApertureFluxes(subimg,x0,y0,ap_rad,np.sqrt(sky),EGAIN)
        subimg = d[int(comp_center[0]-hs):int(comp_center[0]+hs),int(comp_center[1]-hs):int(comp_center[1]+hs)]
        sky = np.median(subimg)
        subimg = subimg - sky
        imm = plt.imshow(subimg)
        imm.set_clim(0,10)
        x0,y0 = hs,hs#get_centroid(subimg)
        plt.plot(x0,y0,'wx',markersize=15,alpha=0.5)
        circle = plt.Circle((x0,y0),ap_rad,color='white',fill=False)
        plt.gca().add_artist(circle)
        plt.savefig('phot_images_comp/'+im+'.png')
        plt.close()
        fc,ferrc = getApertureFluxes(subimg,x0,y0,ap_rad,np.sqrt(sky),EGAIN)
        times.append(time)
        rel_fluxes.append(f/fc)
        err = np.sqrt((ferr/fc)**2 + (f*ferrc/(fc**2))**2)
        rel_fluxes_err.append(err)      
plt.close()
print 'Writing...'
ff = open(filename,'w')
times = np.array(times)
rel_fluxes = np.array(rel_fluxes)
rel_fluxes_err = np.array(rel_fluxes_err)
for i in range(len(all_images)):
    ff.write(all_images[i]+'  '+str(times[i])+'  '+str(rel_fluxes[i])+'  '+str(rel_fluxes_err[i])+'\n')
ff.close()
plt.style.use('ggplot')
diff_mags = -2.5*np.log10(rel_fluxes)
err = (1.08574/rel_fluxes)*rel_fluxes_err
print diff_mags
print err
idx = np.argsort(times)
