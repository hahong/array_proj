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
import common_fn as ss
import cPickle as pk

VOLT_CONV = 0.249
N_SNIPPET_PTS = 48

def get_waveforms(fn_mwk, fn_nev, npts=N_SNIPPET_PTS, c_v=VOLT_CONV,
                  verbose=True):
    wv = {}
    wvsq = {}
    nelec = cl.defaultdict(int)
    iimg0 = None

    # -- collect snippets for each image presentation
    # override_elecs: process only these electrodes. 1-based.
    for arg in ss.getspk(fn_mwk, fn_nev):
        # -- preps
        wav = np.array(arg['wav']['unsorted']) * c_v      # unsorted raw snippet
        wavsq = wav * wav
        ch = arg['ch']; t_rel = arg['t_rel']
        t0 = arg['t_imgonset']; iid = arg['imgid']; iimg = arg['iimg']
        if verbose:
            if iimg0 == None: iimg0 = iimg
            if iimg != iimg0:
                n = nelec[1]
                print 'At: %d ([1]: %f)' % (iimg, np.mean(np.sqrt(wvsq[1]/n - wv[1]*wv[1]/n/n)))
                iimg0 = iimg
        # -- main
        nelec[ch] += 1
        if not wv.has_key(ch): wv[ch] = wav
        else: wv[ch] += wav
        if not wvsq.has_key(ch): wvsq[ch] = wavsq
        else: wvsq[ch] += wavsq
        
    for ch in nelec:
        wv[ch] /= nelec[ch]
        wvsq[ch] /= nelec[ch]
     
    return wv, wvsq


def main(fn_mwk, fn_nev, fn_out):
    wv, wvsq = get_waveforms(fn_mwk, fn_nev)
    dev = {}
    for ch in wv:
        dev[ch] = np.mean(np.sqrt(wvsq[ch] - wv[ch]*wv[ch]))
    fp = open(fn_out, 'wb')
    pk.dump({'wv':wv, 'wvsq':wvsq, 'dev':dev}, fp)
    fp.close()

    # print wv[1], wvsq[1], np.mean(np.sqrt(wvsq - wv*wv))




if __name__ == '__main__':
	if len(sys.argv) != 4:
		print 'get_plot_spikes.py <mwk> <nev> <output file name>'
		sys.exit(1)
	fn_mwk, fn_nev, fn_out = sys.argv[1:4]
	print '* Processing:', fn_mwk, fn_nev, fn_out
	main(fn_mwk, fn_nev, fn_out)

