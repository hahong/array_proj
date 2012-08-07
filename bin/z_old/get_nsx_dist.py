#!/usr/bin/env python
# encoding: utf-8

import sys
import os
sys.path.append('lib')
import numpy as np
import collections as cl
import cPickle as pk
from mworks.data import MWKFile
from mergeutil import Merge, BRReader

C_STIM = '#announceStimulus'
I_STIM_ID = 2
T_START =  75000
T_STOP = 175000
C_SUCCESS = 'number_of_stm_shown'   # one visual stimulus should have one "success" event within
T_SUCCESS = 250000                  # 250ms time window in order to be considered as valid one.
OVERRIDE_DELAY_US = 300
VOLT_CONV = 0.249

# ----------------------------------------------------------------------------
def get_stim_info(mf, c_stim=C_STIM, extinfo=False):
    stims = mf.get_events(codes=[c_stim])
    stims = [x for x in stims if x.value['type'] == 'image']
    # when the stimulus was shown? (in us).
    img_onset = [x.time for x in stims]
    # ..and each corresponding image id
    # img_id = [int(x.value['name'].split('_')[1]) for x in stims]
    if not extinfo:
        img_id = [x.value['name'] for x in stims]
    else:
        img_id = [(x.value['name'], x.value['pos_x'], x.value['pos_y'], \
                x.value['rotation'], x.value['size_x'], x.value['size_y']) for x in stims]
    assert (len(img_id) == len(img_onset))
    return img_onset, img_id

# ----------------------------------------------------------------------------
def nsxload(fn_mwk, fn_nev, fn_nsx, override_delay_us=OVERRIDE_DELAY_US,
            verbose=False, extinfo=False, c_success=C_SUCCESS, data_only=True,
            list_form=False):
    mf = MWKFile(fn_mwk)
    mf.open()
    br = BRReader(fn_nev, fn_nsx)
    br.open()

    # read TOC info from the "merged" mwk file
    toc = mf.get_events(codes=[Merge.C_MAGIC])[0].value
    c_spikes = toc[Merge.K_SPIKE]                  # get the code name from the toc

    # when the visual stimuli presented is valid?
    t_success = [ev.time for ev in mf.get_events(codes=[c_success])]
    t_success = np.array(t_success)

    img_onset, img_id = get_stim_info(mf, extinfo=extinfo)
    n_stim = len(img_onset)

    # MAC-NSP time translation
    if override_delay_us != None: 
        t_delay = toc['align_info']['delay']
        t_adjust = int(np.round(override_delay_us - t_delay))
    else: 
        t_adjust = 0
    a, b = toc['align_info']['params']
    f = lambda t_mwk: float(t_mwk - b) / a
    t_start = T_START - t_adjust
    t_stop = T_STOP - t_adjust

    # actual calculation -------------------------------
    for i in range(n_stim):
        t0 = img_onset[i]; iid = img_id[i]
        # check if this presentation is successful. if it's not ignore this.
        if np.sum((t_success > t0) & (t_success < (t0 + T_SUCCESS))) < 1: continue

        if verbose: print 'At', (i + 1), 'out of', n_stim
   
        t_nsx_start = f(t0 + t_start)
        t_nsx_stop = f(t0 + t_stop)
        yield [data for data in br.nsx_read_range_ts(t_nsx_start, t_nsx_stop, \
                    data_only=data_only, list_form=list_form)], br.nsx_chn_order
# ----------------------------------------------------------------------------
def count(dst, dat, rkey):
    assert dst.shape[1] == len(rkey.keys())
    assert dst.shape[0] == dat.shape[0]

    for r in range(dat.shape[0]):
        for v in dat[r]:
            try:
                c = rkey[v]
                dst[r, c] += 1
            except KeyError:
                pass

def get_waveforms(fn_mwk, fn_nev, fn_nsx, c_v=VOLT_CONV, verbose=True, rng=3000):
    lbl = [l for l in range(-rng, rng+1)]  # true value of each column
    dst = None
    rkey = {}
    for c, k in enumerate(lbl): rkey[k] = c

    for dat0, chn_order in nsxload(fn_mwk, fn_nev, fn_nsx, list_form=True, verbose=True):
        dat = np.array(dat0).T
        if dst == None:
            nch = dat.shape[0]
            dst = np.zeros((nch, rng*2 + 1)) 
            dst_diff = np.zeros((nch, rng*2 + 1))
        count(dst, dat, rkey)
        count(dst_diff, np.diff(dat), rkey)

    return dst, dst_diff, lbl, chn_order


def main(fn_mwk, fn_nev, fn_nsx, fn_out):
    dst, dst_diff, lbl, chn_order = get_waveforms(fn_mwk, fn_nev, fn_nsx)

    fp = open(fn_out, 'wb')
    pk.dump({'dst':dst, 'dst_diff':dst_diff, 'lbl':lbl, 'chn_order':chn_order}, fp)
    fp.close()
    """
    pl.figure()
    for i, cnt in enumerate(dst):
        pl.subplot(10, 10, i+1)
        pl.plot(lbl, cnt)

    pl.figure()
    for i, cnt in enumerate(dst_diff):
        pl.subplot(10, 10, i+1)
        pl.plot(lbl, cnt)
        pl.xlim(-1000, 1000)

    pl.show()
    """

if __name__ == '__main__':
	if len(sys.argv) != 5:
		print 'nsxcalc.py <mwk> <nev> <nsx> <output file name>'
		sys.exit(1)
	fn_mwk, fn_nev, fn_nsx, fn_out = sys.argv[1:5]
	print '* Processing:', fn_mwk, fn_nev, fn_nsx, fn_out
	main(fn_mwk, fn_nev, fn_nsx, fn_out)

