import numpy as np
from pyx import text as pyx_text
from pyx import trafo
from pyx import *

def text_pyx(g, x_coord, y_coord, text_input, text_size = -2, color = None, rotation = 0.):
    """
    Function that draws text in a given plot
    INPUTS:
        g       (Object) A graph-type object to which you want to add the text
        x_coord     (Double) x-coordinate (in plot units) at which you want to place
                the text
        y_coord     (Double) y-coordinate (in plot units) at which you want to place
                the text
        text_input  (String) Text that you want to add.
        text_size   (int, optional) Text size of the text added to the plot. Default is -2.
        color       (instance) Color instance that defines the color that you want the 
                text to have. Default is black.
    """

    # First define the text attributes:
    textattrs = [pyx_text.size(text_size),pyx_text.halign.center, pyx_text.vshift.middlezero, trafo.rotate(rotation)]

    # Now convert plot positions to pyx's:
    x0,y0 = g.pos(x_coord, y_coord)

    # If no color is given, draw black text. If color is given, draw text with the input color:
    if color is None:
        g.text(x0,y0,text_input,textattrs)
    else:
        # First, check which was the input color palette:
        color_dict = color
        try:
            color_string = str(color_dict.c)+','+str(color_dict.m)+','+str(color_dict.y)+','+str(color_dict.k)
            color_palette = 'cmyk'
        except:
                        color_string = str(color_dict.r)+','+str(color_dict.g)+','+str(color_dict.b)
                        color_palette = 'rgb'
        # Now draw the text:
        g.text(x0, y0, r"\textcolor["+color_palette+"]{"+color_string+"}{"+text_input+"}",textattrs)

fname = 'TDRIZZLE_0010_010_HATS755004_SDSSz__000'
band = "z'"
title = 'HATS-49'
max_contrast = 4.0
r26,m26,dm26 = np.loadtxt(fname+'/contrast_curve_'+fname+'.fits.dat',unpack=True)
m26 = np.abs(m26)
unit.set(xscale=1.5)
legend_text_size = -2
legend_pos = 'tr'
the_color1 = color.cmyk.MidnightBlue
text.set(mode="latex")
text.preamble(r"\usepackage{color}")
text.preamble(r"\usepackage{wasysym}")
c = canvas.canvas()
g = c.insert(graph.graphxy(height=13, width=15,\
                  x = graph.axis.linear(min=0, max=np.max(r26), title='Radial distance (arcsec)'),
                  y = graph.axis.linear(min=max_contrast+0.5, max=0.0, title="$\Delta "+band+"$")))

g.plot(graph.data.values(x=r26,y=m26, title = None),\
                        styles = [graph.style.line([the_color1,\
                                                    style.linestyle.solid,\
                                                    style.linewidth.thick])])


### Bands around data:
h = graph.graphxy(height=13,width=15, x=graph.axis.linkedaxis(g.axes["x"]), y=graph.axis.linkedaxis(g.axes["y"]))
dupp = h.plot(graph.data.values(x=r26, y=m26+dm26), [graph.style.line()])
dlow = h.plot(graph.data.values(x=r26, y=m26-dm26), [graph.style.line()])

h.doplot()
upp = dupp.path.reversed()
low = dlow.path
x0, y0 = low.atend()
x1, y1 = upp.atbegin()
connect1 = path.line(x0, y0, x1, y1)

area = low << connect1 << upp
area.append(path.closepath())

g.layer("filldata").draw(area, [deco.filled([color.gray(0.8)])])

text_pyx(g, 8, 2.25, title, text_size = -2, color = the_color1, \
                                   rotation = 0.0)

####### UNCOMMENT THIS TO PLOT A SOURCE ############
#distance = 2.85 # arcsecs
#delta_mag = 2.457
#delta_mag_err = 0.013
#planet_color = color.cmyk.Red
#g.plot(graph.data.values(x=[distance], \
#                         y=[delta_mag], \
#                         ymin = [delta_mag+delta_mag_err], \
#                         ymax = [delta_mag-delta_mag_err],title = None),\
#            styles = [graph.style.symbol(graph.style.symbol.circle,\
#            symbolattrs = [deco.filled([planet_color]),
#                        deco.stroked([planet_color])],\
#                        size = 0.15),\
#                        graph.style.errorbar(errorbarattrs = [planet_color])])
####################################################

print(dir(c))
c.writeEPSfile('contrast_curve_'+title)#,meshasbitmap = True,meshasbitmap_resolution=5)
c.writePDFfile('contrast_curve_'+title)#,meshasbitmap = True,meshasbitmap_resolution=5)
