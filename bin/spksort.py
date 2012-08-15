#!/usr/bin/env python

import numpy as np
import cPickle as pk
import tables as tbl
import pywt as wt
import sys
from scipy import interpolate as ipl
from scipy import stats as st
from joblib import Parallel, delayed
from sklearn.cluster import AffinityPropagation
from common_fn import parse_opts2, detect_cpus

# -- defaults for feature computation
RETHRESHOLD_MULT = -2.5
ALIGN_SUBSMP = 256
ALIGN_MAXDT = 2.5
ALIGN_PEAKLOC = 12
ALIGN_FINDBWD = 3
ALIGN_FINDFWD = 6
ALIGN_CUTAT = 3
ALIGN_OUTDIM = 32
ALIGN_PEAKFUNC = 'argmin'

SKIMSPK_TB = 100000   # beginning relative time for collecting examplar spikes
SKIMSPK_TE = 180000
EXTRACT_NPERIMG = 2

FEAT_METHOD = 'wavelet'
FEAT_WAVL_LEV = 4
FEAT_KSSORT = False
FEAT_OUTDIM = 10

# -- defaults for clustering
CLUSTERING_ALG = 'affinity_prop'
AFFINITYPRP_COMMONP = 'min'
QC_MINSNR = 1.        # minimum SNR to be qualified as a cluster
QC_KS_PLEVEL = .05    # a cluster should have "KS-test p" > QC_KS_PLEVEL

# -- other defaults
NCPU = detect_cpus()  # all CPUs
NCPU_HIGHMEM = 4      # for memory intensive tasks

USAGE = \
"""spksort.py: a spike sorting swiss army knife

Feature Computing Mode
======================
Computes feautures of spikes for later clustering after the
re-thresholding (Quiroga et al., 2004) and the spike alignment.

spksort.py [options] feature <input.psf.h5> <output.c1.feat.h5>

Options:
   --with=<reference.h5> Use the parameters used in the reference file
                         (another feature output.c1.feat.h5 file).
   --rethreshold_mult=#  The multiplier for Quiroga's thresholding method
   --align_subsmp=#
   --align_maxdt=#
   --align_peakloc=#
   --align_findbwd=#
   --align_findfwd=#
   --align_cutat=#
   --align_outdim=#
   --feat_metd=<str>     The method used to extract features.  Available:
                         wavelet
   --feat_kssort=<bool>
   --feat_outdim=#
   --feat_wavelet_lev=#
   --skimspk_tb=#        Beginning relative time for collecting examplar spikes
   --skimspk_te=#        End time for examplar spikes
   --extract_nperimg=#   The max number of spikes collected for each stimulus


Clustering Mode
===============
Cluster spikes using pre-computed features.

spksort.py [options] cluster <input.c1.feat.h5> <output.c2.clu.h5>

Options:
   --with=<reference.h5>
   --cluster_alg=<str>   The clustering algorithm to use.  Avaliable:
                         affinity_prop
   --qc_minsnr=#         Minimum SNR to be qualified as a cluster
   --qc_ks_plevel=#      Desired significance level for the KS-test


Collation Mode
==============
(Optional) Find out common clusters in multiple .c2.clu.h5 files.

spksort.py [options] collate <input1.c2.clu.h5> [input2.c2.clu.h5] ...
"""


# -- House-keeping stuffs -----------------------------------------------------
def par_comp(func, spks, n_jobs=NCPU, **kwargs):
    n = spks.shape[0]
    n0 = max(n / n_jobs, 1)
    ibs = range(0, n, n0)
    ies = range(n0, n + n0, n0)
    r = Parallel(n_jobs=n_jobs, verbose=0)(delayed(func)(spks[ib: \
            ie], **kwargs) for ib, ie in zip(ibs, ies))
    return np.concatenate(r)


def KS_all(spks_wt, full=False):
    dev = []
    ps = []
    spks_wt = spks_wt - spks_wt.mean(axis=0)
    spks_wt /= spks_wt.std(ddof=1, axis=0)

    for i in xrange(spks_wt.shape[1]):
        w = spks_wt[:, i]
        #w = (w - w.mean()) / w.std(ddof=1)
        D, p = st.kstest(w, 'norm')
        dev.append(D)
        ps.append(p)

    if full:
        return np.array(dev), np.array(ps), spks_wt
    return np.array(dev), spks_wt


# -- feature extraction related -----------------------------------------------
def rethreshold_by_multiplier_core(Msnp, Msnp_ch, ch, mult=RETHRESHOLD_MULT):
    selected = np.zeros(Msnp_ch.shape, dtype='bool')

    i_ch = np.nonzero(Msnp_ch == ch)[0]
    M = Msnp[i_ch]
    thr = mult * np.median(np.abs(M) / 0.6745)

    if thr < 0:
        i_pass = np.min(M, axis=1) < thr
    else:
        i_pass = np.max(M, axis=1) > thr

    selected[i_ch[i_pass]] = True

    return selected, thr


def rethreshold_by_multiplier_par(Msnp, Msnp_ch, target_chs, \
        mult=RETHRESHOLD_MULT, ncpu=NCPU_HIGHMEM):
    """Experimental parralelized version of rethreshold_by_multiplier()"""
    core = rethreshold_by_multiplier_core
    r = Parallel(n_jobs=ncpu, verbose=0)(delayed(core)(Msnp, \
            Msnp_ch, ch, mult=mult) for ch in target_chs)
    thrs = [e[1] for e in r]
    selected = np.array([e[0] for e in r]).any(axis=0)

    return selected, thrs


def rethreshold_by_multiplier(Msnp, Msnp_ch, target_chs, \
        mult=RETHRESHOLD_MULT):
    thrs = []
    selected = np.zeros(Msnp_ch.shape, dtype='bool')
    for ch in target_chs:
        i_ch = np.nonzero(Msnp_ch == ch)[0]
        M = Msnp[i_ch]
        thr = mult * np.median(np.abs(M) / 0.6745)
        thrs.append(thr)

        if thr < 0:
            i_pass = np.min(M, axis=1) < thr
        else:
            i_pass = np.max(M, axis=1) > thr

        selected[i_ch[i_pass]] = True

    return selected, thrs


def align_core(spks, subsmp=ALIGN_SUBSMP, maxdt=ALIGN_MAXDT, \
        peakloc=ALIGN_PEAKLOC, findbwd=ALIGN_FINDBWD, \
        findfwd=ALIGN_FINDFWD, cutat=ALIGN_CUTAT, outdim=ALIGN_OUTDIM, \
        peakfunc=ALIGN_PEAKFUNC):
    """Alignment algorithm based on (Quiroga et al., 2004)"""

    if peakfunc == 'argmin':
        peakfunc = np.argmin
    else:
        raise ValueError('Not recognized "peakfunc"')

    R = np.empty((spks.shape[0], outdim), dtype='int16')
    n = spks.shape[1]
    x0 = np.arange(n)

    for i_spk, spk in enumerate(spks):
        tck = ipl.splrep(x0, spk, s=0)

        xn = np.arange(peakloc - findbwd, peakloc + findfwd, 1. / subsmp)
        yn = ipl.splev(xn, tck)

        dt = xn[peakfunc(yn)] - peakloc
        if np.abs(dt) > maxdt:
            dt = 0

        x = x0 + dt
        y = ipl.splev(x, tck)

        R[i_spk] = np.round(y).astype('int16')[cutat: cutat + outdim]
        #dts.append(dt)
    return R


def wavelet_core(spks, level=FEAT_WAVL_LEV):
    """Wavelet transform feature extraction"""
    feat = np.array([np.concatenate(wt.wavedec(spks[i], 'haar', \
            level=level)) for i in xrange(spks.shape[0])])
    return feat


def skim_imgs(Mimg, Mimg_tabs, Msnp_tabs, t_adjust=0, tb0=SKIMSPK_TB, \
        te0=SKIMSPK_TE, n_blk=20000):
    idx_eachimg = [np.nonzero(Mimg == i_img)[0][0] for i_img \
            in np.unique(Mimg)]
    t_eachimg = Mimg_tabs[idx_eachimg]
    i_eachimg = Mimg[idx_eachimg]

    ibie = []
    ib = 0
    ie = 0
    for t0 in t_eachimg:
        tb = t0 + tb0 - t_adjust
        te = t0 + te0 - t_adjust

        xb = np.searchsorted(Msnp_tabs[ib: ib + n_blk], tb)
        if xb >= n_blk:
            xb = np.searchsorted(Msnp_tabs[ib:], tb)
        ib += xb

        xe = np.searchsorted(Msnp_tabs[ie: ie + n_blk], te)
        if xe >= n_blk:
            xe = np.searchsorted(Msnp_tabs[ie:], te)
        ie += xe
        ibie.append((ib, ie))
    return ibie, i_eachimg


def get_example_spikes(Msnp, Msnp_ch, ibie, target_chs, \
        extract_nperimg=EXTRACT_NPERIMG):
    """Extract all spikes specified in ibie"""

    # not using defaultdict to save space
    res = {}
    for ch in target_chs:
        res[ch] = []

    # sweep over each image
    for ib, ie in ibie:
        idx = range(ib, ie)
        Msnp_img = Msnp[idx]
        Msnp_ch_img = Msnp_ch[idx]

        for ch in target_chs:
            res[ch].extend(Msnp_img[Msnp_ch_img == ch][:extract_nperimg])

    return res


# ----------------------------------------------------------------------------
def get_features(fn_inp, fn_out, opts):
    config = {}
    config['rethreshold_mult'] = RETHRESHOLD_MULT
    config['align'] = {}
    config['align']['subsmp'] = ALIGN_SUBSMP
    config['align']['maxdt'] = ALIGN_MAXDT
    config['align']['peakloc'] = ALIGN_PEAKLOC
    config['align']['findbwd'] = ALIGN_FINDBWD
    config['align']['findfwd'] = ALIGN_FINDFWD
    config['align']['cutat'] = ALIGN_CUTAT
    config['align']['outdim'] = ALIGN_OUTDIM
    config['align']['peakfunc'] = ALIGN_PEAKFUNC

    config['skimspk'] = {}
    config['skimspk']['tb0'] = SKIMSPK_TB
    config['skimspk']['te0'] = SKIMSPK_TE
    config['extract_nperimg'] = EXTRACT_NPERIMG

    config['feat'] = {}
    config['feat']['metd'] = FEAT_METHOD
    config['feat']['kwargs'] = {'level': FEAT_WAVL_LEV}
    config['feat']['kssort'] = FEAT_KSSORT
    config['feat']['outdim'] = FEAT_OUTDIM

    # -- process opts
    # TODO: implement!!!

    # -- preps
    print '* Began feature computations...'
    h5 = tbl.openFile(fn_inp)
    Msnp = h5.root.Msnp.read()
    Msnp_ch = h5.root.Msnp_ch.read()
    Msnp_pos = h5.root.Msnp_pos.read()
    Msnp_tabs = h5.root.Msnp_tabs.read()
    Mimg = h5.root.Mimg.read()
    Mimg_tabs = h5.root.Mimg_tabs.read()

    t_adjust = h5.root.meta.t_adjust.read()
    t_start0 = h5.root.meta.t_start0.read()
    t_stop0 = h5.root.meta.t_stop0.read()

    idx2iid = h5.root.meta.idx2iid.read()
    iid2idx_pk = h5.root.meta.iid2idx_pk.read()
    idx2ch = h5.root.meta.idx2ch.read()
    ch2idx_pk = h5.root.meta.ch2idx_pk.read()
    all_chs = range(len(idx2ch))

    # -- re-threshold
    print '* Re-thresholding...'
    if type(config['rethreshold_mult']) is float or int:
        print '* At: re-thresholding...'
        thr_sel, thrs = rethreshold_by_multiplier(Msnp, Msnp_ch, \
                all_chs, config['rethreshold_mult'])

        Msnp = Msnp[thr_sel]
        Msnp_ch = Msnp_ch[thr_sel]
        Msnp_pos = Msnp_pos[thr_sel]
        Msnp_tabs = Msnp_tabs[thr_sel]
        Msnp_selected = np.nonzero(thr_sel)[0]

    else:
        thr_sel = None
        thrs = None
        Msnp_selected = None

    # -- align
    print '* Aligning...'
    Msnp = par_comp(align_core, Msnp, **config['align'])

    # -- feature extraction
    print '* Extracting features...'
    if config['feat']['metd'] == 'wavelet':
        Msnp_feat = par_comp(wavelet_core, Msnp, \
                **config['feat']['kwargs'])

    elif config['feat']['metd'] != 'pca':
        config['feat'].pop('kwargs')
        raise NotImplementedError('PCA not implemented yet')

    else:
        raise ValueError('Not recognized "feat_metd"')

    # -- get training examples...
    print '* Collecting snippet examples...'
    ibie, iuimg = skim_imgs(Mimg, Mimg_tabs, Msnp_tabs, \
            t_adjust, **config['skimspk'])

    Msnp_feat_train = get_example_spikes(Msnp_feat, Msnp_ch, ibie, all_chs, \
            extract_nperimg=config['extract_nperimg'])
    """DBG
    Msnp_feat_train = {}
    for ch in all_chs:
        Msnp_feat_train[ch] = Msnp_feat
    """

    # -- get final features
    outdim = config['feat']['outdim']
    if config['feat']['kssort']:
        Msnp_feat_use = []

        for ch in all_chs:
            X = np.array(Msnp_feat_train[ch])
            # get deviations from Gaussian
            devs, _ = KS_all(X)
            # got top-n deviations
            devs = np.argsort(-devs)[:outdim]
            Msnp_feat_use.append(devs)
    else:
        Msnp_feat_use = [range(outdim)] * len(all_chs)
    Msnp_feat_use = np.array(Msnp_feat_use)

    # -- done! write everything...
    print '* Writing results...'
    filters = tbl.Filters(complevel=4, complib='blosc')
    t_int16 = tbl.Int16Atom()
    t_uint32 = tbl.UInt32Atom()
    t_uint64 = tbl.UInt64Atom()
    t_float32 = tbl.Float32Atom()

    h5o = tbl.openFile(fn_out, 'w')
    CMsnp = h5o.createCArray(h5o.root, 'Msnp', t_int16, \
            Msnp.shape, filters=filters)
    CMsnp_tabs = h5o.createCArray(h5o.root, 'Msnp_tabs', t_uint64, \
            Msnp_tabs.shape, filters=filters)
    CMsnp_ch = h5o.createCArray(h5o.root, 'Msnp_ch', t_uint32, \
            Msnp_ch.shape, filters=filters)
    CMsnp_pos = h5o.createCArray(h5o.root, 'Msnp_pos', t_uint64, \
            Msnp_pos.shape, filters=filters)
    CMsnp_feat = h5o.createCArray(h5o.root, 'Msnp_feat', t_float32, \
            Msnp_feat.shape, filters=filters)
    CMsnp_feat_use = h5o.createCArray(h5o.root, 'Msnp_feat_use', t_uint64, \
            Msnp_feat_use.shape, filters=filters)
    CMsnp_selected = h5o.createCArray(h5o.root, 'Msnp_selected', \
            t_uint64, Msnp_selected.shape, filters=filters)

    CMsnp[...] = Msnp
    CMsnp_tabs[...] = Msnp_tabs
    CMsnp_ch[...] = Msnp_ch
    CMsnp_pos[...] = Msnp_pos
    CMsnp_feat[...] = Msnp_feat
    CMsnp_feat_use[...] = Msnp_feat_use
    CMsnp_selected[...] = Msnp_selected
    h5o.createArray(h5o.root, 'Msnp_feat_train_pk', pk.dumps(Msnp_feat_train))

    h5o.createArray(h5o.root, 'Mimg', Mimg)
    h5o.createArray(h5o.root, 'Mimg_tabs', Mimg_tabs)

    meta = h5o.createGroup('/', 'meta', 'Metadata')
    h5o.createArray(meta, 't_start0', t_start0)
    h5o.createArray(meta, 't_stop0', t_stop0)
    h5o.createArray(meta, 't_adjust', t_adjust)
    h5o.createArray(meta, 'config_feat_pk', pk.dumps(config))

    h5o.createArray(meta, 'idx2iid', idx2iid)
    h5o.createArray(meta, 'iid2idx_pk', iid2idx_pk)
    h5o.createArray(meta, 'idx2ch', idx2ch)
    h5o.createArray(meta, 'ch2idx_pk', ch2idx_pk)
    h5o.createArray(meta, 'thrs', thrs)

    h5o.close()


# ----------------------------------------------------------------------------


def main():
    if len(sys.argv) < 3:
        print USAGE
        return

    args, opts = parse_opts2(sys.argv[1:])
    mode = args[0]

    # -- parsing extra arguments (mainly for backward compatibility)
    if mode == 'feature':
        get_features(args[1], args[2], opts)
    elif mode == 'cluster':
        pass
    else:
        raise ValueError('Invalid mode')

    print 'Done.                                '


if __name__ == '__main__':
    main()
