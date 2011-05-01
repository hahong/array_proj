#!/usr/bin/env python

import numpy as np
import cPickle as pk
import struct
import sys
import os
sys.path.append('lib')
import numpy.linalg                  # numpy's linear algebra code to do svd and pca.
from multiprocessing import Process, Queue
from collections import defaultdict, namedtuple
from mworks.data import MWKFile
from mergeutil import Merge, BRReader, PLXReader

C_STIM = '#announceStimulus'
I_STIM_ID = 2
T_START = -100000
T_STOP = 250000
C_SUCCESS = 'number_of_stm_shown'   # one visual stimulus should have one "success" event within
T_SUCCESS = 250000                  # 250ms time window in order to be considered as valid one.
OVERRIDE_DELAY_US = 300
T_REJECT = 10
N_REJECT = 50
DEFAULT_N_PCA = 3
N_PRE_PT = 11
SEARCH_RNG = [6, 15]

# ----------------------------------------------------------------------------
def xget_events(mf, **kwargs):
    """Memory saving trick"""
    def xloader(q, mf, kwargs):
        evs = mf.get_events(**kwargs)
        q.put([(ev.value, ev.time, ev.empty) for ev in evs])

    q = Queue()
    p = Process(target=xloader, args=(q, mf, kwargs))
    p.start()
    evs = q.get()
    p.join()
    
    Event = namedtuple('Event', 'value time empty')
    return [Event(*ev) for ev in evs]

# ----------------------------------------------------------------------------
def get_stim_info(mf, c_stim=C_STIM, extinfo=False, \
        dynstim='ds_9999', rotrmn=180, blank='zb_99999',\
        exclude_img=None):
    stims0 = xget_events(mf, codes=[c_stim])
    stims = []
    for x in stims0:
        if type(x.value) == int:
            continue
        if x.value['type'] == 'image':
            if exclude_img != None:
                ignore = False
                for patt in exclude_img:
                    # if it's in the excluded list, don't put that one!!
                    if patt in x.value['name']: 
                        ignore = True
                        break
                if ignore: continue
            stims.append(x)
        elif x.value['type'] == 'dynamic_stimulus' or x.value['type'] == 'blankscreen':
            stims.append(x)
        elif x.value['type'] == 'image_directory_movie':
            if type(x.value['current_stimulus']) == int: continue
            stims.append(x)
        # otherwise, ignore

    # when the stimulus was shown? (in us).
    img_onset = [x.time for x in stims]
    # ..and each corresponding image id
    # img_id = [int(x.value['name'].split('_')[1]) for x in stims]
    img_id = []
    for x in stims:
        if x.value['type'] == 'image':
            if not extinfo: iid = x.value['name']
            else: iid = (x.value['name'], x.value['pos_x'], x.value['pos_y'], \
                x.value['rotation'], x.value['size_x'], x.value['size_y'])
        elif x.value['type'] == 'dynamic_stimulus':
            if not extinfo: iid = dynstim
            else: iid = (dynstim, x.value['xoffset'], x.value['yoffset'], \
                x.value['rotation'] % rotrmn, x.value['width'], x.value['height'])
        elif x.value['type'] == 'image_directory_movie':
            if not extinfo: iid = x.value['current_stimulus']['name']
            else: iid = (x.value['current_stimulus']['name'], \
                    x.value['current_stimulus']['pos_x'], \
                    x.value['current_stimulus']['pos_y'], \
                    x.value['current_stimulus']['rotation'], \
                    x.value['current_stimulus']['size_x'], \
                    x.value['current_stimulus']['size_y'])
        elif x.value['type'] == 'blankscreen':
            if not extinfo: iid = blank
            else: iid = (blank, 0, 0, 0, 0, 0)
        img_id.append(iid)
    
    assert (len(img_id) == len(img_onset))
    return img_onset, img_id

# ----------------------------------------------------------------------------
def load_spike_data(neu_filename):
    ext = os.path.splitext(neu_filename)[1]
    if ext.lower() == '.nev':
        nf = BRReader(neu_filename)
    else:
        nf = PLXReader(neu_filename)
    return nf

# ----------------------------------------------------------------------------
def seq_search(iterable, target):
    """do sequential search"""
    for i, e in enumerate(iterable):
        if e != target: continue
        return i
    return None


# ----------------------------------------------------------------------------
def sort_uniq(base, *args):
    """sort and remove duplicates based on `base` and apply on to `args`"""
    if len(args) == 0: return None
    res = []
    # sort
    si = np.argsort(base)
    base = np.array(base[si])
    for arg in args:
        res.append(np.array(arg[si]))
    # remove duplicates
    di = np.nonzero(np.diff(base) == 0)[0]
    si = list(set(range(len(base))) - set(list(di)))
    for i in xrange(len(res)):
        res[i] = np.array(res[i][si])
    return res

# ----------------------------------------------------------------------------
def getspk(fn_mwk, fn_nev, override_elecs=None, \
        override_delay_us=OVERRIDE_DELAY_US, verbose=False, extinfo=False, \
        c_success=C_SUCCESS, t_start0=T_START, t_stop0=T_STOP, full=True, new_thr=None):
    mf = MWKFile(fn_mwk)
    mf.open()
    br = load_spike_data(fn_nev)
    br.open()
    valid_elec_ids = br.spike_id

    # read TOC info from the "merged" mwk file
    toc = xget_events(mf, codes=[Merge.C_MAGIC])[0].value
    c_spikes = toc[Merge.K_SPIKE]                  # get the code name from the toc

    # when the visual stimuli presented is valid?
    t_success = [ev.time for ev in xget_events(mf, codes=[c_success])]
    t_success = np.array(t_success)

    img_onset, img_id = get_stim_info(mf, extinfo=extinfo)
    n_stim = len(img_onset)

    # MAC-NSP time translation
    if override_delay_us != None: 
        t_delay = toc['align_info']['delay']
        t_adjust = int(np.round(override_delay_us - t_delay))
    else: 
        t_adjust = 0
    t_start = t_start0 - t_adjust
    t_stop = t_stop0 - t_adjust

    # memory saving trick ------------------------------
    def xloader(q, mf, code, time_range):
        spikes = mf.get_events(codes=[code], time_range=time_range)
        q.put([(s.value['id'], s.value['foffset'], s.time) for s in spikes])

    # actual calculation -------------------------------
    for i in range(n_stim):
        t0 = img_onset[i]; iid = img_id[i]
        # check if this presentation is successful. if it's not ignore this.
        if np.sum((t_success > t0) & (t_success < (t0 + T_SUCCESS))) < 1: continue

        if verbose: print 'At', (i + 1), 'out of', n_stim
   
        # TODO: unused due to the memory management
        # spikes = mf.get_events(codes=[c_spikes], time_range=[t0+t_start, t0+t_stop])
        # for s in spikes:

        q = Queue()
        p = Process(target=xloader, args=(q, mf, c_spikes, [t0 + t_start, t0 + t_stop]))
        p.start()
        spk_info = q.get()
        p.join()

        for (ch, pos, t_abs) in spk_info:
            if ch not in valid_elec_ids: continue
            if override_elecs != None and ch not in override_elecs: continue 
            # waveform
            try:
                wav_info = br.read_once(pos=pos, proc_wav=True)
            except Exception, e:
                print '*** Exception:', e
                continue

            # -- apply new threshold if requested
            if new_thr != None:
                if new_thr.has_key('mult'):
                    lthr = br.chn_info[ch]['low_thr']
                    hthr = br.chn_info[ch]['high_thr']
                    if lthr == 0: thr0 = hthr
                    else: thr0 = lthr
                    thr = thr0 * new_thr['mult']
                elif new_thr.has_key('abs'):
                    thr = new_thr['abs']
            
                wf0 = wav_info['waveform']
                wf = set_new_threshold(wf0, thr)
                
                if wf == None: continue   # `wf0` is smaller than `thr`
                wav_info['waveform'] = wf

            # -- done. yield the results
            if full:
                # relative time
                t_rel = int(t_abs + t_adjust - t0)
                wav = {'unsorted':wav_info['waveform']}
                rtn = {'t_abs':t_abs, 't_rel':t_rel, 'ch':ch, 
                       'wav_info':wav_info, 'wav':wav, 
                       't_imgonset':t0, 'imgid':iid, 'iimg':i, 'pos':pos}
                yield rtn
            else:
                yield wav_info, i



# -----------------------------------------------------------------------------
# if there are more than `N_REJET` spikes within `T_REJECT`us window,
# invalidate all of them.
def invalidate_artifacts(buf0, t_reject=T_REJECT, n_reject=N_REJECT, verbose=True):
    ti_all = [(b['timestamp'], i) for i, b in enumerate(buf0)]
    ti_all = sorted(ti_all)
    t_all = np.array([t[0] for t in ti_all])
    i_all = [t[1] for t in ti_all]

    nb = len(buf0)
    ri = range(nb)
    i = 0
    while i < nb - 1:
        ii = []
        t0 = t_all[i]
        for j in xrange(i + 1, nb):
            if t_all[j] < t0 + t_reject: ii.append(j)
            else: break
        i = j

        if len(ii) < n_reject: continue
        for ix in ii:
            try:
                ri.remove(i_all[ix])
            except ValueError:
                pass

    buf = [buf0[i] for i in ri]
    if verbose and len(buf) != nb:
        print '* Rejecting', nb - len(buf), 'spikes.'
    return buf


# -----------------------------------------------------------------------------
# Parse the options in the command line.
def parse_opts(opts0):
    opts = {}
    # parse the stuff in "opts"
    for opt in opts0:
        parsed = opt.split('=')
        key = parsed[0].strip()
        if len(parsed) > 1:
            cmd = parsed[1].strip()
        else:
            cmd = ''
        opts[key] = cmd

    return opts



def prep_files(flist, sep=',', extchk=True):
    flist = flist.split(sep)
    if flist[0][0] == '+':
        flist = [f.strip() for f in open(flist[0][1:]).readlines()]
    if extchk: assert all([os.path.exists(f) for f in flist])

    return flist

# -----------------------------------------------------------------------------
# Set new threshold `thr`. If the `waveform` cannot pass `thr` returns None.
# The new waveform is re-aligned based on the steepest point.
# The returned new waveform has `n_pre` points before the alignment point.
def set_new_threshold(wavform, thr, n_pre=N_PRE_PT, rng=SEARCH_RNG, i_chg=20):
    wav = np.array(wavform)
    sgn = np.sign(thr)
    if np.max(wav[rng[0]:rng[1]] * sgn) < np.abs(thr): return None   # reject

    """ NOT USED -- GIVES IMPRECISE RESULT
    # -- align: find the steepest point having the same sign as `sgn`
    df = np.diff(wav)
    si = np.argsort(-sgn * df)   # reverse sorted
    for i in si:
        if np.sign(wav[i]) == sgn: break
    """
    # -- align: find the point where waveform crosses `thr`
    n = len(wav)
    for i in range(n - 1):
        if sgn*wav[i] <= sgn*thr and sgn*thr <= sgn*wav[i+1]: break
    if i == n - 2: return None   # although i could be n - 2, it's highly likely an artifact
    n_shift = n_pre - i - 1      # > 0: right shift, < 0: left shift
    if n_shift == 0: return wav
    
    wavnew = np.empty(wav.shape)
    wavnew[n_shift:] = wav[:-n_shift]   # PBC shifting
    wavnew[:n_shift] = wav[-n_shift:]

    # -- done: but if the spike doesn't change its sign within `i_chg`, reject.
    if np.max(-sgn * wavnew[n_pre:i_chg]) < 0: return None

    """ DEBUG
    if np.abs(n_shift) > 3:
        print '!!!', n_shift, '/', i, '/', n
        print '---',  np.max(-sgn * wavnew[n_pre:i_chg])
        print list(wav)
        print list(wavnew)
    """

    return wavnew

# -----------------------------------------------------------------------------
# fastnorm: from Nicolas' code
def fastnorm(x):
    xv = x.ravel()
    return np.dot(xv, xv)**(1/2.)


# fastsvd: from Nicolas' code
def fastsvd(M):
    h, w = M.shape
    # -- thin matrix
    if h >= w:
        # subspace of M'M
        U, S, V = np.linalg.svd(np.dot(M.T, M))
        U = np.dot(M, V.T)
        # normalize
        for i in xrange(w):
            S[i] = fastnorm(U[:,i])
            U[:,i] = U[:,i] / S[i]
    # -- fat matrix
    else:
        # subspace of MM'
        U, S, V = np.linalg.svd(np.dot(M, M.T))
        V = np.dot(U.T, M)
        # normalize
        for i in xrange(h):
            S[i] = fastnorm(V[i])
            V[i,:] = V[i] / S[i]
    return U, S, V


def pca_eigvec(M, pca_threshold = DEFAULT_N_PCA):
    U,S,V = fastsvd(M)
    eigvectors = V.T
    eigvectors = eigvectors[:, :pca_threshold]
    # this gives PCA: 
    # M = np.dot(M, eigvectors)
    return eigvectors

# ----------------------------------------------------------------------------
if __name__ == '__main__':
    print 'This script is not supposed to be excuted directly.'
