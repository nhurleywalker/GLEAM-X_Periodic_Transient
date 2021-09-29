#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 14:04:50 2019

@author: zha292

Make Q,U spectrum from channel images
"""

import numpy as np
from astropy.io import fits
from astropy import wcs
from astropy.coordinates import SkyCoord
import astropy.units as u
import glob, os, sys


def getnoise(data):
	# hdu = fits.open(img_name)
	# data = hdu[0].data
	rms_initial = np.std(data)
	rms = np.std(data[np.logical_and(data>-2.5*rms_initial, data<2.5*rms_initial)])
	# hdu.close()
	return rms

# flist = np.genfromtxt('quocka_select.csv', dtype=str)
s = sys.argv[1]
obsid = sys.argv[2]


# rah = s[1:3]
# ram = s[3:5]
# ras = s[5:7]
# decd = s[7:10]
# decm = s[10:12]
# decs = s[12:14]

# coor0 = SkyCoord(rah+':'+ram+':'+ras+' '+decd+':'+decm+':'+decs, unit=(u.hourangle, u.deg))
# coor0 = SkyCoord('17:37:41.0', '-56:32:38', unit=(u.hourangle, u.deg))

# coor = [coor0.ra.deg, coor0.dec.deg]

coor = [246.99796, -52.58453]

c_light = 299792458.0

# get the peak flux pixel from mfs images
qlist = glob.glob('./'+obsid+'/*Q.fits')
qlist.sort()
ulist = glob.glob('./'+obsid+'/*U.fits')
ulist.sort()
ilist = glob.glob('./'+obsid+'/*I.fits')
ilist.sort()

# print qlist

frame = qlist[-1]
img_f = fits.open(frame)
wcs_f = wcs.WCS(img_f[0].header).dropaxis(3).dropaxis(2)
pix_coor = wcs_f.wcs_world2pix(coor[0], coor[1], 0)
img_f.close()

stokes_file = open(s+'_'+obsid+'.txt', 'w')
rmsyn_file = open(s+'_'+obsid+'.rmsyn.txt', 'w')

for q_name in qlist[0:-1]:

	peak_x = int(np.round(pix_coor[0]))
	peak_y = int(np.round(pix_coor[1]))

	q_img = fits.open(q_name)
	chan = q_img[0].header['CRVAL3']
	data_q = q_img[0].data[0,0]
	peak_q = data_q[peak_y, peak_x]
	# box_q = data_q[img_size-box_r:img_size-box_l, img_size-box_t:img_size-box_b]
	noise_q = getnoise(data_q)
	q_img.close()

	u_name = q_name.replace('-Q.', '-U.')
	u_img = fits.open(u_name)
	data_u = u_img[0].data[0,0]
	peak_u = data_u[peak_y, peak_x]
	# box_u = data_u[img_size-box_r:img_size-box_l, img_size-box_t:img_size-box_b]
	noise_u = getnoise(data_u)
	u_img.close()

	i_name = q_name.replace('-Q.', '-I.')
	i_img = fits.open(i_name)
	data_i = i_img[0].data[0,0]
	peak_i = data_i[peak_y, peak_x]
	# box_u = data_u[img_size-box_r:img_size-box_l, img_size-box_t:img_size-box_b]
	noise_i = getnoise(data_i)
	i_img.close()

	v_name = q_name.replace('-Q.', '-V.')
	v_img = fits.open(v_name)
	data_v = v_img[0].data[0,0]
	peak_v = data_v[peak_y, peak_x]
	# box_u = data_u[img_size-box_r:img_size-box_l, img_size-box_t:img_size-box_b]
	noise_v = getnoise(data_v)
	v_img.close()

	# print q_name

	stokes_file.write(str(chan)+' '+str(peak_i)+' '+str(peak_q)+' '+str(peak_u)+' '+str(noise_i)+' '
		+str(noise_q)+' '+str(noise_u)+' '+str(peak_v)+' '+str(noise_v)+'\n')
	rmsyn_file.write(str(chan)+' '+str(peak_i)+' '+str(peak_q)+' '+str(peak_u)+' '+str(noise_i)+' '
                +str(noise_q)+' '+str(noise_u)+'\n')
