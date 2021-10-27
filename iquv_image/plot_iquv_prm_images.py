#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
from astropy.io import fits

# Sans-serif fonts for Nature
plt.rcParams.update({
    "text.usetex": False,
    "font.family": "sans-serif",
    "font.size": 7,
    "font.sans-serif": ["Helvetica"]})

cm = 1/2.54  # centimeters in inches

def format_inset(a, framecolor="white"):
    lon = a.coords[0]
    lat = a.coords[1]
    lon.set_ticks_visible(False)
    lon.set_ticklabel_visible(False)
    lat.set_ticks_visible(False)
    lat.set_ticklabel_visible(False)
    lon.set_axislabel('')
    lat.set_axislabel('')
    a.coords.frame.set_color(framecolor)
    return a

tran_ra = 246.99796
tran_dec = -52.58453

i_map = fits.open('1200354592-MFS-pbcorr-image-I.fits')
q_map = fits.open('1200354592-MFS-pbcorr-image-Q.fits')
u_map = fits.open('1200354592-MFS-pbcorr-image-U.fits')
v_map = fits.open('1200354592-MFS-pbcorr-image-V.fits')
p_map = fits.open('FDF_maxPI.fits')
rm_map = fits.open('FDF_peakRM.fits')

wcs = WCS(i_map[0].header, naxis=2)

# select polarised source region from the RM map
p_data = p_map[0].data[0]
rm_data = rm_map[0].data[0]

# Flag the <7-sigma regions
p_mean = np.mean(p_data)
p_std = np.std(p_data)
rm_data[p_data<p_mean+7*p_std] = 0

# Set up constants to be used across all images
boxwidth = 50 # pixel width of inset box
cent = i_map[0].header["NAXIS1"]/2 # center of inset
# image settings
iargs = {"vmin" : "-0.5", "vmax" : "1.4", "origin" : "lower", "cmap" : "cubehelix"}
# box settings
bargs = {"width" : boxwidth, "height" : boxwidth, "facecolor" : 'none', "edgecolor" : 'white', "lw" : 0.5}

fig = plt.figure(figsize=[18*cm,18*cm])

ax = fig.add_axes([0.1, 0.7, 0.3, 0.3], projection = wcs)
im = ax.imshow(i_map[0].data[0,0], **iargs)
ax.text(150, 1400, 'Stokes I', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('')
lat.set_axislabel('Declination (J2000)')
lon.set_ticklabel_visible(False)
cutout = Cutout2D(i_map[0].data[0,0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.11, 0.71, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, **iargs)
cutax = format_inset(cutax)

ax = fig.add_axes([0.4, 0.7, 0.3, 0.3], projection = wcs)
ax.imshow(q_map[0].data[0,0], **iargs)
ax.text(150, 1400, 'Stokes Q', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_ticklabel_visible(False)
lat.set_ticklabel_visible(False)
lon.set_axislabel(None)
lat.set_axislabel('')
cutout = Cutout2D(q_map[0].data[0,0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.41, 0.71, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, **iargs)
cutax = format_inset(cutax)

ax = fig.add_axes([0.1, 0.4, 0.3, 0.3], projection = wcs)
ax.imshow(u_map[0].data[0,0], **iargs)
ax.text(150, 1400, 'Stokes U', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('')
lat.set_axislabel('Declination (J2000)')
lon.set_ticklabel_visible(False)
cutout = Cutout2D(u_map[0].data[0,0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.11, 0.41, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, **iargs)
cutax = format_inset(cutax)

ax = fig.add_axes([0.4, 0.4, 0.3, 0.3], projection = wcs)
ax.imshow(v_map[0].data[0,0], vmin=-0.5, vmax=1.5, origin='lower',cmap='cubehelix')
ax.text(150, 1400, 'Stokes V', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
cax = fig.add_axes([0.71, 0.41, 0.015, 0.28])
char = plt.colorbar(im, cax = cax)
char.set_label('Jy/beam')#, rotation=270, fontsize=11)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_ticklabel_visible(False)
lat.set_ticklabel_visible(False)
lon.set_axislabel('')
lat.set_axislabel('')
cutout = Cutout2D(v_map[0].data[0,0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.41, 0.41, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, **iargs)
cutax = format_inset(cutax)

ax = fig.add_axes([0.1, 0.1, 0.3, 0.3], projection = wcs)
im = ax.imshow(p_map[0].data[0], vmin=-0.5, vmax=1.5, origin='lower',cmap='cubehelix')
ax.text(150, 1400, 'Polarised Intensity', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('Right Ascension (J2000)')
lat.set_axislabel('Declination (J2000)')
cutout = Cutout2D(p_map[0].data[0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.11, 0.11, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, **iargs)
cutax = format_inset(cutax)

ax = fig.add_axes([0.4, 0.1, 0.3, 0.3], projection = wcs)
im = ax.imshow(rm_data, vmin=-75, vmax=75, origin='lower',cmap='RdBu')
ax.text(150, 1400, 'Rotation Measure', backgroundcolor = "white")
ax.add_patch(Rectangle((cent-boxwidth/2, cent-boxwidth/2), **bargs))
cax = fig.add_axes([0.71, 0.11, 0.015, 0.28])
char = plt.colorbar(im, cax = cax)
char.set_label('rad m$^{-2}$')#, rotation=270, fontsize=11)
lon = ax.coords[0]
lat = ax.coords[1]
lon.set_axislabel('Right Ascension (J2000)')
lat.set_ticklabel_visible(False)
cutout = Cutout2D(rm_map[0].data[0], (cent, cent), (boxwidth, boxwidth), wcs = wcs)
cutax = fig.add_axes([0.41, 0.11, 0.1, 0.1], projection = cutout.wcs)
cutim = cutax.imshow(cutout.data, vmin=-75, vmax=75, origin='lower',cmap='RdBu')
cutax = format_inset(cutax, "black")

plt.subplots_adjust(hspace=0.0, wspace=0.0)
plt.savefig('iquvprm.png', facecolor='w', dpi=300, bbox_inches='tight')
plt.savefig('iquvprm.pdf', facecolor='w', dpi=300, bbox_inches='tight')
