#!/usr/bin/env python
# encoding: utf-8

import sys
import os
sys.path.append('lib')
import numpy as np
import pywt as wt
import collections as cl
import cPickle as pk
import shutil as sh
import tempfile
from common_fn import *
from scipy import stats
from scipy.stats import kstest
from joblib import Parallel, delayed

N_SNIPPET_PER_IMG = 2
I_THRESHOLD = 10
N_FEATDIM = 12
IIMG_HALT = 500
SUFF_FET = '.fet.'
SUFF_INF = '.inf.'
FEAT_METD = 'wav'
SEP = ','
WAVEDEC_LEV = 5
EPS = 1e-5

APPLY_NEW_THR = True
M_REJECT = 5.5/3.5
NCPU = -1

PARRUN = '++parrun'

# -----------------------------------------------------------------------------
def get_few_waveforms(fn_mwks, fn_nevs, nmax=N_SNIPPET_PER_IMG, ihalt=IIMG_HALT, \
        apply_new_thr=APPLY_NEW_THR, reject=M_REJECT, fill_empty=True, wavedec_lev=WAVEDEC_LEV):
    if apply_new_thr: new_thr = {'mult':reject}
    else: new_thr = None
    wavs = cl.defaultdict(list)   # wavs[ch] = [[snippet 1], [snippet 2], ...]

    # -- collect snippets for each image presentation
    # override_elecs: process only these electrodes. 1-based.
    for fn_mwk, fn_nev in zip(fn_mwks, fn_nevs):
        prev_iimg = -1

        for arg in getspk(fn_mwk, fn_nev,  
                t_start0=75000, t_stop0=175000, new_thr=new_thr):
            iimg = arg['iimg']
            ch = arg['ch']
            wav = arg['wav']['unsorted']      # unsorted raw snippet

            # if this is the beginning of the new image presentation...
            if prev_iimg != iimg:
                if iimg >= ihalt: break
                print 'At:', iimg, '  \r',
                sys.stdout.flush()
                prev_iimg = iimg
                nwavs = cl.defaultdict(int)
            # if we already collected `n_max` snippets for this image, skip this one.
            if nwavs[ch] >= nmax: continue

            wavs[ch].append(wav)
            nwavs[ch] += 1

    if fill_empty:
        for ch in range(min(wavs), max(wavs)):
            if wavs.has_key(ch): continue
            print '* Empty channel', ch
            wavs[ch].append(np.zeros(np.shape(wav)))

    print '\nSampling Done.'
    return wavs


# -----------------------------------------------------------------------------
def get_feat_transf(wavs, n_featdim=N_FEATDIM, metd=FEAT_METD):
    feat_proj = {}
    for ch in wavs:
        if len(wavs[ch]) > 1:
            ndim = len(wavs[ch][0])
            break
    
    for ch in wavs:
        # -- prep/clean-up
        n = len(wavs[ch])
        if n < 2:
            print '* Really bad channel', ch, 'of length', n
            feat_proj[ch] = np.zeros((ndim, n_featdim))
            continue
        elif n < n_featdim:
            print '* Bootstrapping samples for the (bad) channel', ch
            M = [wavs[ch][i % n] for i in range(n_featdim)]
            M = np.array(M)
        else:
            M = np.array(wavs[ch])

        # -- actual computation
        # PCA-based --------------------------------------------- 
        if metd == 'pca':
            feat_proj[ch] = pca_eigvec(M, pca_threshold=n_featdim)
            #pc = np.dot(wf, T)

        # Wavelet decomposition (from Wave_clus) ----------------
        elif metd == 'wav':
            # -- wavelet decomposition
            W = np.empty(M.shape)
            npts = len(M)
            for i in range(npts):
                dec = wt.wavedec(M[i], 'haar', level=wavedec_lev)
                W[i] = np.concatenate(dec)[:ndim]

            # -- K-S test
            dev = []   # deviation from norm
            for i in range(ndim): 
                w = np.array(W[:,i])
                # if the stddv is too small, ignore this entry
                if w.std() < EPS:
                    print '* Ch %d: ignore %d-th entry' % (ch, i)
                    dev.append(-1)
                    continue
                w = (w - w.mean()) / w.std()
                D, p = kstest(w, 'norm') 
                dev.append(D)

            dev = np.array(dev)
            # index to the `n_featdim` biggest features
            ii = np.argsort(-dev)[:n_featdim]
            T = np.zeros((ndim, n_featdim))
            for c, r in enumerate(ii): T[r,c] = 1
            feat_proj[ch] = T
            # DEBUG:
            # print '*M', ch, M[1], M[2]
            # print '*W', ch, W[1], W[2]
            # print '**', ch, W[:10,0], W[:10,3], W[:10,4]
            # print '**', ch, T.T.tolist()

        # Error! God why? whyyyyy? ------------------------------- 
        else:
            raise ValueError, 'wrong `metd` specified!'

    return feat_proj



# -----------------------------------------------------------------------------
def write_features(T, fn_mwk, fn_nev, fn_out, full=False, n_featdim=N_FEATDIM, metd=FEAT_METD, \
        apply_new_thr=APPLY_NEW_THR, reject=M_REJECT, wavedec_lev=WAVEDEC_LEV):
    ff = {}; fi = {}
    for ch in T:
        ff[ch] = open(fn_out + SUFF_FET + str(ch) + '.tmp', 'wt')
        fi[ch] = open(fn_out + SUFF_INF + str(ch) + '.tmp', 'wt')
        print >>ff[ch], n_featdim 
        print >>fi[ch], n_featdim 
    fmt = '%e\t' * n_featdim    # format string

    if full: gen = get_full_spk(fn_nev, reject=reject, apply_new_thr=apply_new_thr)
    else: gen = get_img_spk(fn_mwk, fn_nev, reject=reject, apply_new_thr=apply_new_thr)

    # -- do the job for all spikes
    for wav_info in gen:
        ch = wav_info['id']
        wav = np.array(wav_info['waveform'])
        if metd == 'wav':                      # if wavelet decomposition
            ndim = len(wav)
            dec = wt.wavedec(wav, 'haar', level=wavedec_lev)
            wav = np.concatenate(dec)[:ndim]
        pc = tuple(np.dot(wav, T[ch]))

        print >>ff[ch], fmt % pc
        print >>fi[ch], wav_info['file_pos'], wav_info['timestamp'] 

    # -- clean up
    for ch in ff:
        ff[ch].close()
        fi[ch].close()

    # -- remove potential duplicates
    for ch in T:
        fn_inf = fn_out + SUFF_INF + str(ch)
        fn_fet = fn_out + SUFF_FET + str(ch)
        inf = np.loadtxt(fn_inf + '.tmp', skiprows=1)
        fet = np.loadtxt(fn_fet + '.tmp', skiprows=1)

        ninf0 = len(inf)
        infs, fet = sort_uniq(inf[:,0].astype('int'), inf, fet)
        ninf1 = len(infs)

        # if rewriting is unnecessary:
        if ninf0 == ninf1 and np.all(infs == inf):
            sh.move(fn_inf + '.tmp', fn_inf)
            sh.move(fn_fet + '.tmp', fn_fet)
            return
        
        # otherwise, rewriting is needed due to the change
        inf = infs
        print 'Ch %d: rewriting (size change: %d -> %d)' % (ch, ninf0, ninf1)
        assert len(inf) == len(fet)

        ff = open(fn_fet, 'wt')
        fi = open(fn_inf, 'wt')
        print >>ff, n_featdim
        print >>fi, n_featdim

        for i in xrange(ninf1):
            print >>ff, fmt % tuple(fet[i])
            print >>fi, int(inf[i,0]), inf[i,1]

        ff.close()
        fi.close()
        os.unlink(fn_fet + '.tmp')
        os.unlink(fn_inf + '.tmp')


# --------------------------------------------------------------------------------
def get_full_spk(fn_nev, apply_new_thr=APPLY_NEW_THR, reject=M_REJECT):   
    """get all spikes regardless of images"""
    br = load_spike_data(fn_nev)
    br.open()
    valid_elec_ids = br.spike_id
    while True:
        wav_info = br.read_once(proc_wav=True)
        if wav_info == None: break        # reached EOF
        ch = wav_info['id']
        if ch not in valid_elec_ids: continue

        if apply_new_thr:
            lthr = br.chn_info[ch]['low_thr']; hthr = br.chn_info[ch]['high_thr']
            if lthr == 0: thr0 = hthr
            else: thr0 = lthr

            thr = thr0 * reject
            wf0 = wav_info['waveform']
            wf = set_new_threshold(wf0, thr)
            
            if wf == None: continue   # `wf0` is smaller than `thr`
            wav_info['waveform'] = wf

        yield wav_info


def get_img_spk(fn_mwk, fn_nev, reject=M_REJECT, apply_new_thr=APPLY_NEW_THR):
    """get spikes around the image presentation"""
    prev_iimg = -1
    buf0 = []

    if apply_new_thr: new_thr = {'mult':reject}
    else: new_thr = None

    for arg, iimg in getspk(fn_mwk, fn_nev, full=False, new_thr=new_thr):
        # new image is presented. yield the contents in the buffer
        if prev_iimg != iimg:
            prev_iimg = iimg
            print 'At:', iimg, '  \r',
            sys.stdout.flush()
            if len(buf0) != 0:
                buf = invalidate_artifacts(buf0)
                for b in buf: yield b
                del buf, buf0
            if arg['id'] > 0: buf0 = []

        # delayed processing: store'em first
        buf0.append(arg)

    # still need to yield the contents
    if len(buf0) != 0:
        buf = invalidate_artifacts(buf0)
        for b in buf: yield b


# --------------------------------------------------------------------------------
def main(fn_mwks, fn_nevs, fn_outs, full=False, n_featdim=N_FEATDIM, reject=M_REJECT, \
        apply_new_thr=APPLY_NEW_THR, metd=FEAT_METD, nmax=N_SNIPPET_PER_IMG, wavedec_lev=WAVEDEC_LEV, ncpu=NCPU):
    # -- prepare
    wavs = get_few_waveforms(fn_mwks, fn_nevs, \
            reject=reject, apply_new_thr=apply_new_thr, nmax=nmax, wavedec_lev=wavedec_lev)  
    T = get_feat_transf(wavs, n_featdim=n_featdim, metd=metd)       # PCA
    del wavs
    for ch in T: print '* Shape:', ch, T[ch].shape

    # -- do the job
    print 'Writing features...'
    #for fn_mwk, fn_nev, fn_out in zip(fn_mwks, fn_nevs, fn_outs):
    #    write_features(T, fn_mwk, fn_nev, fn_out, full=full, n_featdim=n_featdim, metd=metd, \
    #            reject=reject, apply_new_thr=apply_new_thr, wavedec_lev=wavedec_lev)
    r = Parallel(n_jobs=ncpu, verbose=1)(delayed(parrun_push)((T, fn_mwk, fn_nev, fn_out, full, n_featdim, metd, reject, apply_new_thr, wavedec_lev)) for fn_mwk, fn_nev, fn_out in zip(fn_mwks, fn_nevs, fn_outs))


def parrun_push(args):
    fd, tmpf = tempfile.mkstemp()
    os.close(fd)
    pk.dump(args, open(tmpf, 'wb'))
    os.system('./prep_sorting.py %s %s' % (PARRUN, tmpf))

def parrun_pop(argfile):
    T, fn_mwk, fn_nev, fn_out, full, n_featdim, metd, reject, apply_new_thr, wavedec_lev = pk.load(open(argfile))
    write_features(T, fn_mwk, fn_nev, fn_out, full=full, n_featdim=n_featdim, metd=metd, reject=reject, apply_new_thr=apply_new_thr, wavedec_lev=wavedec_lev)
    os.unlink(argfile)


def prep_files(fn_mwks, fn_nevs, fn_outs):
    def load_set(flist, extchk=True):
        if flist[0][0] == '+':
            flist = [f.strip() for f in open(flist[0][1:]).readlines()]
        if extchk: assert all([os.path.exists(f) for f in flist])
        return flist

    fn_mwks = load_set(fn_mwks.split(SEP))
    fn_nevs = load_set(fn_nevs.split(SEP))
    fn_outs = load_set(fn_outs.split(SEP), extchk=False)

    assert len(fn_mwks) == len(fn_nevs)
    assert len(fn_mwks) == len(fn_outs)
    return fn_mwks, fn_nevs, fn_outs


if __name__ == '__main__':
    if len(sys.argv) == 3 and sys.argv[1] == PARRUN:
        parrun_pop(sys.argv[2])
        sys.exit(0)

    if len(sys.argv) < 4:
        print 'prep_sorting.py <mwk> <nev> <output file prefix> [options 1] [options 2] [...]'
        print 'Prepares input files for spike sorting, or simply rejects invalid spikes.'
        print 'Prefixes of "mwk", "nev", and "output" files can be comma separated'
        print 'multiple entries.'
        print
        print 'Options:'
        print '   featdim=#    - sets the number of features'
        print '   reject=#     - re-sets the threshold in terms of multiples of the' 
        print '                  previous threshold (no sign flip)'
        print '   metd=<metd>  - sets the waveform feature calculation method (wav or pca)'
        print '   nmax=#       - maximum number of spikes/image during the sampling phase'
        print '   full         - sort on all spikes, rather than during img presentations'
        print '                  (not recommended)'
        sys.exit(1)

    fn_mwks, fn_nevs, fn_outs = sys.argv[1:4]
    fn_mwks, fn_nevs, fn_outs = prep_files(fn_mwks, fn_nevs, fn_outs)

    opts = parse_opts(sys.argv[4:])
    print '* Processing:', fn_mwks, fn_nevs, fn_outs

    # -----
    if 'full' in opts: full = True
    else: full = False

    if 'featdim' in opts: n_featdim = int(opts['featdim'])
    elif 'npca' in opts: n_featdim = int(opts['npca'])   # for compatibility
    else: n_featdim = N_FEATDIM
    
    apply_new_thr = APPLY_NEW_THR
    if 'reject' in opts:
        if opts['reject'] == 'noadj':
            apply_new_thr = False
            reject = 1
        else:
            reject = float(opts['reject'])
    else: reject = M_REJECT

    if 'metd' in opts: metd = opts['metd']
    else: metd = FEAT_METD
    assert metd in ['pca', 'wav']

    if 'nmax' in opts: nmax = int(opts['nmax'])
    else: nmax = N_SNIPPET_PER_IMG

    if 'ncpu' in opts: ncpu = int(opts['ncpu'])
    else: ncpu = NCPU

    if 'wavedec_lev' in opts: wavedec_lev = int(opts['wavedec_lev'])
    else: wavedec_lev = WAVEDEC_LEV
    # -----
    print '* Variables: (full, n_featdim, reject, apply_new_thr, metd, nmax, wavedec_lev, ncpu) =', \
            (full, n_featdim, reject, apply_new_thr, metd, nmax, wavedec_lev, ncpu)

    main(fn_mwks, fn_nevs, fn_outs, full=full, n_featdim=n_featdim, reject=reject, \
            apply_new_thr=apply_new_thr, metd=metd, nmax=nmax, wavedec_lev=wavedec_lev, ncpu=ncpu)

