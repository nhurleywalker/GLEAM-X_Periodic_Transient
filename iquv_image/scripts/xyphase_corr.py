#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:04:50 2019

@author: zha292

correct xy phase of U/V chan images
"""

import numpy as np
from astropy.io import fits
import glob, os, sys

obsid = sys.argv[1]

# load the xyphase coefficients

def func(x, a, b):
	return a*x+b

popt = np.genfromtxt('../xyphase_coe_20180204.txt')

# start correcting

# xylist = glob.glob(obsid+'-0*XY-image_b.fits')
xylist = glob.glob(obsid+'-*-XY-image_b.fits')
xylist.sort()

# xymfs = fits.open(obsid + '-MFS-XY-image_b.fits')
# xymfs_d = xymfs[0].data*0.
# yxmfs = fits.open(obsid + '-MFS-XYi-image_b.fits')
# yxmfs_d = yxmfs[0].data*0.

i = 0

for xy_name in xylist:
	xy_file = fits.open(xy_name)
	freq = xy_file[0].header['CRVAL3']
	xy_data = xy_file[0].data

	yx_name = xy_name.replace('-XY-image_b', '-XYi-image_b')
	yx_file = fits.open(yx_name)
	yx_data = yx_file[0].data

	phase = func(freq, *popt)/180.*np.pi

	xy_data_corr = xy_data*np.cos(phase) - yx_data*np.sin(phase)
	yx_data_corr = xy_data*np.sin(phase) + yx_data*np.cos(phase)
	# u_data_corr = u_data*np.cos(phase) - v_data*np.sin(phase)
	# v_data_corr = u_data*np.sin(phase) + v_data*np.cos(phase)

	xy_file[0].data = xy_data_corr
	yx_file[0].data = yx_data_corr

	xy_corr_name = xy_name.replace('-XY-image_b', '-XY-image')
	yx_corr_name = yx_name.replace('-XYi-image_b', '-XYi-image')

	xy_file.writeto(xy_corr_name)
	yx_file.writeto(yx_corr_name)	
	xy_file.close()
	yx_file.close()

# 	if np.all(xy_data_corr==0):
# 		print(xy_name)
# 		continue
# 	else:
# 		xymfs_d = xymfs_d+xy_data_corr
# 		yxmfs_d = yxmfs_d+yx_data_corr
# 		i = i+1

# xymfs[0].data = xymfs_d/i
# yxmfs[0].data = yxmfs_d/i

# xymfs.writeto(obsid + '-MFS-XY-image.fits')
# yxmfs.writeto(obsid + '-MFS-XYi-image.fits')

