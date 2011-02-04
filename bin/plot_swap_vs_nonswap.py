#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by labuser on 2010-11-15.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import pylab as pl
import cPickle as pk
import glob
import numpy as np



def main():
	fns = glob.glob('./out_njm/*out*')
	fig = pl.figure(1)
	s_p = []
	[s_p.append(fig.add_subplot(2,3,s_p_in)) for s_p_in in range(1,7,1)]
	pair_a_ns = []
	pair_a_s = []
	pair_a_fov = []
	pair_b_ns = []
	pair_b_s = []
	pair_b_fov = []
	for fn in fns:
		fn_data = pk.load(open(fn))
		if fn.split('/')[-1].split('_')[1]<='20101001':
			pair_a_ns.extend(fn_data['ns_dPN'])
			pair_a_s.extend(fn_data['s_dPN'])
			pair_a_fov.extend(fn_data['fov_dPN'])
		else:
			pair_b_ns.extend(fn_data['ns_dPN'])
			pair_b_s.extend(fn_data['s_dPN'])
			pair_b_fov.extend(fn_data['fov_dPN'])
	s_p[0].plot(pair_a_ns,pair_a_s,'k.',alpha = 0.25)
	s_p[0].errorbar(np.mean(pair_a_ns),np.mean(pair_a_s),xerr = np.std(pair_a_ns),yerr=np.std(pair_a_s),fmt = 'o-',color = 'k',capsize = 0,ms = 6,elinewidth = 2)
	s_p[1].plot(pair_a_fov,pair_a_s,'k.',alpha = 0.25)
	s_p[1].errorbar(np.mean(pair_a_fov),np.mean(pair_a_s),xerr = np.std(pair_a_fov),yerr=np.std(pair_a_s),fmt = 'o-',color = 'k',capsize = 0,ms = 6,elinewidth = 2)
	s_p[2].plot(pair_a_fov,pair_a_ns,'k.',alpha = 0.25)
	s_p[2].errorbar(np.mean(pair_a_fov),np.mean(pair_a_ns),xerr = np.std(pair_a_fov),yerr=np.std(pair_a_ns),fmt = 'o-',color = 'k',capsize = 0,ms = 6,elinewidth = 2)
	s_p[3].plot(pair_b_ns,pair_b_s,'r.',alpha = 0.25)
	s_p[3].errorbar(np.mean(pair_b_ns),np.mean(pair_b_s),xerr = np.std(pair_b_ns),yerr=np.std(pair_b_s),fmt = 'o-',color = 'r',mec ='r',mfc = 'r',capsize = 0,ms = 6,elinewidth = 2)
	s_p[4].plot(pair_b_fov,pair_b_s,'r.',alpha = 0.25)
	s_p[4].errorbar(np.mean(pair_b_fov),np.mean(pair_b_s),xerr = np.std(pair_b_fov),yerr=np.std(pair_b_s),fmt = 'o-',color = 'r',mec ='r',mfc = 'r',capsize = 0,ms = 6,elinewidth = 2)
	s_p[5].plot(pair_b_fov,pair_b_ns,'r.',alpha = 0.25)
	s_p[5].errorbar(np.mean(pair_b_fov),np.mean(pair_b_ns),xerr = np.std(pair_b_fov),yerr=np.std(pair_b_ns),fmt = 'o-',color = 'r',mec ='r',mfc = 'r',capsize = 0,ms = 6,elinewidth = 2)
	
	#print 'non swap position = ' + str(np.median(pair_b_ns))
	#print 'swap position = ' + str(np.median(pair_b_s))
	#print 'fovea position = ' + str(np.median(pair_b_fov))
	[t_s_p.plot([-50,50],[0,0],'k--',alpha = 0.125) for t_s_p in s_p]
	[t_s_p.plot([0,0],[-50,50],'k--',alpha = 0.125) for t_s_p in s_p]
	[t_s_p.plot([-50,50],[-50,50],'k--',alpha = 0.125) for t_s_p in s_p]
	[t_s_p.plot([50,-50],[-50,50],'k--',alpha = 0.125) for t_s_p in s_p]
	[t_s_p.set_aspect('equal') for t_s_p in s_p]
	for ax in s_p:
		ax.spines['left'].set_position(('outward',10))
		ax.spines['bottom'].set_position(('outward',10))
		ax.set_xticks([-50,-25,0,25,50])
		ax.set_xticklabels([])
		ax.set_yticks([-50,-25,0,25,50])
		ax.set_yticklabels([])
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.spines['top'].set_color('none')
		ax.spines['right'].set_color('none')
	s_p[3].set_xticklabels([-50,-25,0,25,50])
	s_p[3].set_yticklabels([-50,-25,0,25,50])
	s_p[3].set_xlabel('ns_dPN')
	s_p[3].set_ylabel('s_dPN')
	s_p[0].set_ylabel('s_dPN')
	s_p[4].set_xlabel('fov_dPN')
	s_p[4].set_ylabel('s_dPN')
	s_p[1].set_ylabel('s_dPN')
	s_p[5].set_xlabel('fov_dPN')
	s_p[5].set_ylabel('ns_dPN')
	s_p[2].set_ylabel('ns_dPN')
	
	#pl.show()
	fig.savefig('./plots/dPN_summary_corr085_dFR40.pdf')
	pass


if __name__ == '__main__':
	main()

