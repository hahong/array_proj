#!/usr/bin/env python
# encoding: utf-8
"""
summarize.py

Created by najiblocal on 2010-08-25.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import glob
import pylab as pl
import numpy as np

def get_dp(fname):
	for ind,x in enumerate(fname):
		#subplot(2,7,ind+1)
		temp = np.loadtxt(x,delimiter=',',skiprows=1)[:,1:]
		tempmu = np.mean(temp,1)
		# old: tempmax = np.mean(np.sort(temp,1)[:,-10:],1)
		tempmax = np.empty(tempmu.shape)
		# calculate the mean of top-absolute-10% for each channel
		for i_row, row in enumerate(temp):
			nsel = int(len(row) * 0.1)              # number of 10%
			rabs = [(abs(r), r) for r in row]
			tempmax[i_row] = np.mean( \
							[r for _, r in sorted(rabs, reverse=True)[:nsel]])
		if ind == 0:
			dprimemu = tempmu
			dprimemx = tempmax
		else:
			dprimemu = np.c_[dprimemu,tempmu]
			dprimemx = np.c_[dprimemx,tempmax]
	return dprimemu,dprimemx
	
	
def show_best(a_sortind):
	i_A = np.nonzero((a_sortind >= 0) & (a_sortind <= 95))[0]
	i_M = np.nonzero((a_sortind >= 96) & (a_sortind <= 96*2 - 1))[0]
	i_P = np.nonzero((a_sortind >= 96*2) & (a_sortind <= 96*3 - 1))[0]
	print 'A:', a_sortind[i_A]
	print 'M:', a_sortind[i_M]
	print 'P:', a_sortind[i_P]
	print '#A, #M, #P:', len(i_A), len(i_M), len(i_P)
	
	
def plot_cumul(arr):   # arr should be sorted
	y = []
	for x in arr:
		n = np.sum(arr > x)
		y.append(n)
	pl.plot(arr, y)
	
def main(TARGET, LASTN=0):
	if len(sys.argv) < 2:
		print 'summarize.py <dir>'
		return
		
	dirname = sys.argv[1]
	fname = glob.glob(dirname+'/*blank_bstrp.csv')[-LASTN:]
	#fname = glob.glob(dirname+'/*blank.csv')[-LASTN:]
	print fname
	
	fnameA = [x for x in fname if TARGET + '_A_' in x]
	dpmuA,dpmxA = get_dp(fnameA)
	fnameM = [x for x in fname if TARGET + '_M_' in x]
	dpmuM,dpmxM = get_dp(fnameM)
	fnameP = [x for x in fname if TARGET + '_P_' in x]
	print fnameP
	dpmuP,dpmxP = get_dp(fnameP)
	
	dpmu = np.r_[np.mean(dpmuA,1),np.mean(dpmuM,1),np.mean(dpmuP,1)]
	dpmx = np.r_[np.mean(dpmxA,1),np.mean(dpmxM,1),np.mean(dpmxP,1)]	
	#dpmu = np.r_[dpmuA,dpmuM,dpmuP]
	#dpmx = np.r_[dpmxA,dpmxM,dpmxP)]	
	# old: sortindmx = set(np.argsort(dpmx)[-128:])
	# old: sortindmu = set(np.argsort(dpmu)[-128:])
	sortindmx = set(np.argsort(np.abs(dpmx))[-128:])
	sortindmu = set(np.argsort(np.abs(dpmu))[-128:])
	arr_sortindmu = np.array(list(sortindmu))
	arr_sortindmx = np.array(list(sortindmx))
	
	print 'mean:', sortindmu
	print 'max10:', sortindmx
	print 'common:', sortindmu & sortindmx
	print '#common:', len(sortindmu&sortindmx)
	print 'thr:', dpmu[np.argsort(dpmu)[128]], dpmx[np.argsort(dpmx)[128]]
	show_best(arr_sortindmx)
	
	pl.figure()
	pl.plot(dpmu,'k')
	pl.figure()
	pl.plot(dpmx,'k')
	pl.figure()
	pl.subplot(3,2,1)
	pl.imshow(dpmuA,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	pl.subplot(3,2,2)
	pl.imshow(dpmxA,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	pl.subplot(3,2,3)
	pl.imshow(dpmuM,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	pl.subplot(3,2,4)
	pl.imshow(dpmxM,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	pl.subplot(3,2,5)
	pl.imshow(dpmuP,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	pl.subplot(3,2,6)
	pl.imshow(dpmxP,vmin=-2,vmax=2,aspect='auto',interpolation='nearest')
	
	pl.figure()
	pl.subplot(3,2,1)
	pl.hist(dpmu[:96], bins=40, range=(-1, 3))
	pl.subplot(3,2,2)
	pl.hist(dpmx[:96], bins=40, range=(-1, 3))
	pl.subplot(3,2,3)
	pl.hist(dpmu[96:96*2], bins=40, range=(-1, 3))
	pl.subplot(3,2,4)
	pl.hist(dpmx[96:96*2], bins=40, range=(-1, 3))
	pl.subplot(3,2,5)
	pl.hist(dpmu[96*2:96*3], bins=40, range=(-1, 3))
	pl.subplot(3,2,6)
	pl.hist(dpmx[96*2:96*3], bins=40, range=(-1, 3))
	
	pl.figure()
	pl.subplot(121)
	pl.hist(dpmu[arr_sortindmu], bins=40)
	pl.subplot(122)
	pl.hist(dpmx[arr_sortindmx], bins=40)

	pl.figure()
	pl.subplot(121)
	plot_cumul(sorted(dpmu))
	pl.subplot(122)
	plot_cumul(sorted(dpmx))
	
	
	#pl.colorbar()


if __name__ == '__main__':
	#main('Nicole')
	#main('', LASTN=27)   # last 3 days
	main('')   # all days
	pl.show()

