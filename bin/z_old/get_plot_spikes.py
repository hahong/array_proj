#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by labuser on 2010-10-18.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
sys.path.append('lib')
import pylab as pl
import numpy as np
import collections as cl
import numpy.linalg #numpy's linear algebra code to do svd and pca.
import cPickle as pk
from scipy.cluster.vq import vq, whiten, kmeans2
from joblib import Memory
from common_fn import *

mem = Memory(cachedir='.z_tmp_get_plot_spikes.py', verbose=0)

# -----------------------------------------------------------------------------
DEF_HIGHLIGHT = 'r'

def get_color_code(centroids, threshold=1.8, 
		rules_default=['#aaaaaa','#aaaaaa','#aaaaaa','#aaaaaa', '#aaaaaa'], 
		highlight=DEF_HIGHLIGHT):
    com = centroids.mean(axis=0)
    #print com
    index = [(fastnorm(com - centroids[i]), i) for i in range(len(centroids))]
    index = sorted(index, reverse=True)
    print index
    rule = list(rules_default)
	# first biggest deviation
    if fastnorm(centroids[index[0][1]] - com) < threshold: return rule
    rule[index[0][1]] = highlight
    #if fastnorm(centroids[index[1][1]] - centroids[index[2][1]]) > threshold * 0.85:
    #    rule[index[1][1]] = 'g-'
    #print index
    return rule

def draw_color_coded(x, waves, label, col, invisible=True, highlight=DEF_HIGHLIGHT):
	for i, wav in enumerate(waves):
		if col[label[i]] == highlight or not invisible:
			pl.plot(x, wav, color=col[label[i]])

	
# -----------------------------------------------------------------------------
VOLT_CONV = 0.249
TIME_CONV = 0.033
N_PRETHRESHOLD = 10
SEL_ELEC_A = [17,24, 61, 73, 69, 30, 34, 65, 18, 42, 50, 51,53, 71]   # selected electrodes to process
SEL_ELEC_A = range(1,97,1)   # selected electrodes to process
#SEL_ELEC_A = [17,53]
N_SNIPPET_MAX = 1000
N_SNIPPET_PTS = 48

def _get_waveforms(fn_mwk, fn_nev, sel=SEL_ELEC_A, nmax=N_SNIPPET_MAX, npts=N_SNIPPET_PTS):
	nelec = cl.defaultdict(int)
	wavs = np.empty((len(sel), nmax, npts))    # elec #, snippet #, data point #
	iids = cl.defaultdict(list)
	iimgs = np.empty((len(sel), nmax)) 
	trels = np.empty((len(sel), nmax))

	# -- collect snippets for each image presentation
	# override_elecs: process only these electrodes. 1-based.
	for arg in getspk(fn_mwk, fn_nev, override_elecs=sel):
		# -- preps
		wav = np.array(arg['wav']['unsorted'])      # unsorted raw snippet
		ch = arg['ch']; t_rel = arg['t_rel']
		t0 = arg['t_imgonset']; iid = arg['imgid']; iimg = arg['iimg']

		# -- main
		
		iarr = sel.index(ch)
		isnp = nelec[ch]
		if isnp >= nmax:
			if all([nelec[c] >= nmax for c in nelec]): break
			continue

		wavs[iarr,isnp,:] = wav
		iids[iarr].append(iid)
		iimgs[iarr,isnp] = iimg
		trels[iarr,isnp] = t_rel
		nelec[ch] += 1
	
	return nelec, wavs, iids, iimgs, trels
get_waveforms = mem.cache(_get_waveforms)

def get_cluster_centroids(wavs):
	pc_sd = []
	k_cent = []
	k_label = []
	T_all = []
	for wfind,wf in enumerate(wavs):
		T = pca_eigvec(wf)
		pc_sd.append(np.std(np.dot(wf,T),0))
		pc = whiten(np.dot(wf, T))
		n_labels = 3
		cent, w_label = kmeans2(pc, n_labels)
		k_cent.append(cent)
		k_label.append(w_label)
		T_all.append(T)
	
	return k_cent,pc_sd,k_label,T_all


def sort_on_centroids(fn_mwk, fn_nev,k_cent_all,pc_sd_all,T_all,sel=SEL_ELEC_A, nmax=N_SNIPPET_MAX, npts=N_SNIPPET_PTS,c_v=VOLT_CONV):
	#this spikes sorts based on how close a spike is to the farthest centroid from center of mass of all the centroids
	sorted_data = {}
	for arg in getspk(fn_mwk, fn_nev, override_elecs=sel):
			# -- preps
			w = np.array(arg['wav']['unsorted'])
			w *= c_v
			ch = arg['ch']
			if ch not in sorted_data.keys():
				sorted_data[ch] = {'good_sp':[],'bad_sp':[],'g_sp_t_abs':[],'g_sp_t_rel':[],'g_sp_t_imgonset':[],'g_sp_im_id':[],'g_sp_im_ind':[],'b_sp_t_abs':[],'b_sp_t_rel':[],'b_sp_t_imgonset':[],'b_sp_im_id':[],'b_sp_im_ind':[]}
			ch_ind = sel.index(ch)
			#pick the centroid and the spikes sd for that channel
			pc_sd = pc_sd_all[ch_ind]
			k_cent = k_cent_all[ch_ind]
			k_dist = np.array([fastnorm(k-k_cent.mean(axis=0)) for k in k_cent])
			T = T_all[ch_ind]
			b_c = np.flatnonzero(k_dist == k_dist.max())
			#compute the distance of the waveform from the the centroids and determine whether it is a good spike or a bad spike.
			w_dist = np.array([fastnorm((np.dot(w,T)/pc_sd)-k_cent[k]) for k in range(len(k_cent))])
			if  np.flatnonzero(w_dist == w_dist.min())==b_c and k_dist.max() > 1.75:
				sorted_data[ch]['good_sp'].append(w)
				sorted_data[ch]['g_sp_t_abs'].append(arg['t_abs'])
				sorted_data[ch]['g_sp_t_rel'].append(arg['t_rel'])
				sorted_data[ch]['g_sp_t_imgonset'].append(arg['t_imgonset'])
				sorted_data[ch]['g_sp_im_id'].append(arg['imgid'])
				sorted_data[ch]['g_sp_im_ind'].append(arg['iimg'])
			else:
				sorted_data[ch]['b_sp_t_abs'].append(arg['t_abs'])
				sorted_data[ch]['b_sp_t_rel'].append(arg['t_rel'])
				sorted_data[ch]['b_sp_t_imgonset'].append(arg['t_imgonset'])
				sorted_data[ch]['b_sp_im_id'].append(arg['imgid'])
				sorted_data[ch]['b_sp_im_ind'].append(arg['iimg'])
				
	return sorted_data

def main(fn_mwk, fn_nev, fn_out, c_v=VOLT_CONV, c_t=TIME_CONV, n_prethr=N_PRETHRESHOLD,):
	#fn_mwk = '../analysis/data_merged/Chabo_20100818_RSVPNicole_A_001.mwk'
	#fn_nev = '../../Blackrock/Data/Chabo_20100818_RSVPNicole_A_001.nev'
	nelec, wavs, iids, iimgs, trels = get_waveforms(fn_mwk, fn_nev, sel=SEL_ELEC_A, nmax=N_SNIPPET_MAX)
	wavs *= c_v    # convert to uV
	t = (np.array(range(wavs.shape[-1])) - n_prethr) * c_t
	k_cent_all,pc_sd_all,k_label_all,T_all = get_cluster_centroids(wavs)
	sorted_data = sort_on_centroids(fn_mwk, fn_nev,k_cent_all,pc_sd_all,T_all)
	f = open(fn_out,'wb')
	pk.dump(sorted_data,f)
	f.close()
	b_c = []
	for k_ind,k_cent in enumerate(k_cent_all):
		k_dist = np.array([fastnorm(k-k_cent.mean(axis=0)) for k in k_cent])
		if k_dist.max() > 1.75:
			b_c.append((np.flatnonzero(k_dist == k_dist.max()),k_ind+1))
	print b_c
	print len(b_c)
	# -- spike sorting
	pl.figure()
	good_sp = []
	bad_sp = []
	#np.flatnonzero(k_dist == k_dist.max()) == np.flatnonzero(w_dist == w_dist.min()):
	#for k in k_cent:
	#	print fastnorm(k-k_cent.mean(axis=0))
	#	print fastnorm(k-)
	#cols = get_color_code(k_cent, threshold=2.)
	#pl.subplot(10,10,SEL_ELEC_A[wfind])
	#draw_color_coded(t, wf, k_label, cols)
	#best_c =  np.flatnonzero(k_label == b_c)
	foo =  [a for a in good_sp if a not in best_c]
	pl.subplot(151)
	#pl.plot(wf[best_c,:].T,'g')
	#pl.ylim([-200,200])
	#pl.subplot(152)
	#pl.plot(np.array(good_sp).T,'b')
	#pl.ylim([-200,200])
	#pl.subplot(155)
	#pl.plot(np.array(bad_sp).T,color = [0.5,0.5,0.5])
	#pl.ylim([-200,200])
	print 'Done.'
	#pl.show()

if __name__ == '__main__':
	if len(sys.argv) != 4:
		print 'get_plot_spikes.py <mwk> <nev> <output file name>'
		sys.exit(1)
	fn_mwk, fn_nev, fn_out = sys.argv[1:4]
	print '* Processing:', fn_mwk, fn_nev, fn_out
	main(fn_mwk, fn_nev, fn_out)

