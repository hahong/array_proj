#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by labuser on 2010-10-22.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
#import os
import cPickle as pk
import numpy as np
import pylab as pl

def main():
	cdict1 = {'red':   ((0.0, 0.0, 0.0),
	                   (1.0, 0.0, 0.0)),

	         'green': ((0.0, 0.0, 0.0),
	                   (1.0, 1.0, 1.0)),

	         'blue':  ((0.0, 0.0, 0.0),
	                   (1.0, 0.0, 0.0))
	        }
	my_cmap = pl.matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict1,256)
	imdir = '../../../Desktop/ChaboSetup/NATURAL300/'
	mov_mat = pk.load(open('./data_postproc/Chabo_20100818_RSVPNicole_A_001.movie.100ms.smooth.pk'))
	t_fr = mov_mat['frate']
	for chind in range(96):
	    if np.max(t_fr[chind,:]) != 0.0:
	        t_fr[chind,:] = t_fr[chind,:]/np.max(t_fr[chind,:])
	t_fr = t_fr**2
	fr_indices = np.flatnonzero(mov_mat['stm_index'] == 42)
	fr_list = list(fr_indices[np.flatnonzero(np.diff(fr_indices) > 10)+1])
	fr_list.append(fr_indices[0])
	pl.figure(figsize = (4,3),dpi = 100)
	for st_ind in sorted(fr_list)[:10]:
		for fr_ind in range(st_ind-10,st_ind+40):
			print fr_ind
			ax = pl.subplot(1,2,1)
			if np.max(t_fr[:,fr_ind]) != 0.0:
				pl.imshow(np.reshape(t_fr[:,fr_ind],(8,12)),cmap = my_cmap,interpolation='nearest')
				imname = 'Nat300_' + str(mov_mat['stm_index'][fr_ind])
				im = pl.imread(imdir+imname+'.png')
				ay = pl.subplot(1,2,2)
				pl.imshow(im,cmap = 'gray')
				ax.set_axis_off()
				ay.set_axis_off()
				pl.savefig('../movies/array_movie/fr_%04d.png' % fr_ind)
				pl.clf()
	pass


if __name__ == '__main__':
	main()

