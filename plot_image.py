# -*- coding: utf-8 -*-
import lmfit
import pickle
from astropy.io import fits as pyfits
import numpy as np
from pyx import *

####################################################
fname = 'TDRIZZLE_0010_010_HATS755004_SDSSz__000/'
band = "z'"
plot_title = 'HATS-49'
date = '2015-12-23'
min_contrast = -4.
####################################################

image,h = pyfits.getdata(fname+'/original_image.fits',header=True)
image = image - np.median(image)
image = (image - np.min(image))+0.01
image = np.log(image/np.max(image))

# Parameter fit data:
thefile = open(fname+'/out_params.pkl','rb')
params = pickle.load(thefile,encoding='latin1')
xcen,ycen,sigma_x,sigma_y = params['x0'].value,params['y0'].value,\
                            params['sigma_x'].value,params['sigma_y'].value

sigma = (sigma_x+sigma_y)/2.
xpixels = np.arange(image.shape[0])+0.5
ypixels = np.arange(image.shape[1])+0.5
scale_x = 15.20*1e-3
scale_y = 15.20*1e-3
xarcsec = (xpixels - (image.shape[0]/2.))*scale_x
yarcsec = (ypixels - (image.shape[1]/2.))*scale_y
# Convert centroid to arcsec:
xcen = (xcen + 0.5 - (image.shape[0]/2.))*scale_x
ycen = (ycen + 0.5 - (image.shape[1]/2.))*scale_y
print('centroid (ra,dec):',xcen,ycen)
sigma = sigma*(scale_x)
print('sigma = ',sigma)
print('fwhm  = ',sigma*2.355)
x,y = np.meshgrid(xarcsec-xcen,yarcsec-ycen)
data = list(zip(x.flat,y.flat,image.flat))
# Plot:
unit.set(xscale=2.5)
text.set(mode="latex")
text.preamble(r"\usepackage{color}")
text.preamble(r"\usepackage{wasysym}")
c = canvas.canvas()
g = c.insert(graph.graphxy(height=20, width=20, key=graph.key.key(pos='br'),
                  x=graph.axis.linear(min=-3.6, max=4., title='$\Delta$ R.A. (arcsec)'),
                  y=graph.axis.linear(min=-4.0, max=3.2, title='$\Delta$ DEC (arcsec)')))

g.plot(graph.data.points(data, x=1, y=2, color=3, title = None),
       [graph.style.density(gradient=color.gradient.YellowBlack,
                            coloraxis=graph.axis.linear(min=min_contrast, max=0.0,
                                                        title="Relative flux (log-scale)"))])

g.plot(graph.data.values(x=[-60,60],y=[0,0], title = None),\
                       styles = [graph.style.line([color.cmyk.Grey,\
                                 style.linestyle.dashed,\
                                 style.linewidth.thin])])

g.plot(graph.data.values(x=[0,0],y=[-60,60], title = None),\
                       styles = [graph.style.line([color.cmyk.Grey,\
                                 style.linestyle.dashed,\
                                 style.linewidth.thin])])

# Plot circles:
xc,yc = g.pos(0.,0.)#g.pos(xcen,ycen)
xcr,ycr = g.pos(3.,ycen)
rc = np.sqrt((xc-xcr)**2 + (yc-ycr)**2)
c.stroke(path.circle(xc,yc,np.abs(rc)),[style.linewidth.thin, style.linestyle.solid, color.cmyk.Black])

# Plot circles:
xc,yc = g.pos(0.,0.)#g.pos(xcen,ycen)
xcr,ycr = g.pos(1.,ycen)
rc = np.sqrt((xc-xcr)**2 + (yc-ycr)**2)
c.stroke(path.circle(xc,yc,np.abs(rc)),[style.linewidth.thin, style.linestyle.dashed, color.cmyk.Black])

# Plot circles:
############ IF A SOURCE, UNCOMMENT THIS #####################
#xc,yc = g.pos(-1.52613,-3.48374)
#rc = 1#np.sqrt((xc-xcr)**2 + (yc-ycr)**2)
#c.stroke(path.circle(xc,yc,np.abs(rc)),[style.linewidth.thin, style.linestyle.dashed, color.cmyk.Blue])
##############################################################
#x0,y0 = g.pos(0.,0.)
#xr,yr = g.pos(3.98*3.,0.)
#r = np.sqrt((x0-xr)**2+(y0-yr)**2)
#c.stroke(path.circle(x0,y0,np.abs(r)),[style.linewidth.Thick, style.linestyle.solid, color.cmyk.Black])

# Plot text:
xtext,ytext = g.pos(-50,50)
#g.text(xtext,ytext,r"(b)")
g.text(g.width/2, g.height + 0.3, plot_title+' '+date+" (AstraLux Sur+"+band+")", 
                  [text.halign.center, text.valign.bottom])
#kg = graph.graphx(-5,5,length=10,direction="horizontal")
#g.plot(graph.data.points(zero_data, x=1, y=2, color=3),
#       styles = [graph.style.density(gradient=color.gradient.Gray,
#                            keygraph=kg)])

#g.insert(kg)
c.writeEPSfile(plot_title+'.eps')#,write_mesh_as_bitmap = True,write_mesh_as_bitmap_resolution=5)
c.writePDFfile(plot_title+'.pdf')#,write_mesh_as_bitmap = True,write_mesh_as_bitmap_resolution=5)
