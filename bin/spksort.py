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

FEAT_METHOD = 'wavelet'
FEAT_WAVL_LEV = 4

# -- defaults for clustering
CLUSTERING_ALG = 'affinity_prop'
AFFINITYPRP_COMMONP = 'min'

SKIMSPK_TB = 100000   # beginning relative time for collecting examplar spikes
SKIMSPK_TE = 250000
EXTRACT_NPERIMG = 3
EXTRACT_NMAX = 2500

FEAT_KSSORT = True
FEAT_OUTDIM = 10

QC = True
QC_MINSNR = 1.2       # minimum SNR to be qualified as a cluster
QC_KS_PLEVEL = .05    # a cluster should have "KS-test p" > QC_KS_PLEVEL
QC_MINSIZE = 60

# -- other defaults
NCPU = detect_cpus()  # all CPUs
NCPU_HIGHMEM = 4      # for memory intensive tasks
UNSORTED = 0
BADIDX = -1
ATOL = 1e-4

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


Clustering Mode
===============
Cluster spikes using pre-computed features.

spksort.py [options] cluster <input.c1.feat.h5> <output.c2.clu.h5>

Options:
   --with=<reference.h5>
   --cluster_alg=<str>   The clustering algorithm to use.  Avaliable:
                         affinity_prop
   --feat_kssort=<bool>
   --feat_outdim=#
   --feat_wavelet_lev=#
   --skimspk_tb=#        Beginning relative time for collecting examplar spikes
   --skimspk_te=#        End time for examplar spikes
   --extract_nperimg=#   The max number of spikes collected for each stimulus
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
def rethreshold_by_multiplier_core(Msnp, Msnp_ch, ch, \
        mult=RETHRESHOLD_MULT, selected=None):

    return_full = False
    if selected is None:
        selected = np.zeros(Msnp_ch.shape, dtype='bool')
        return_full = True

    i_ch = np.nonzero(Msnp_ch == ch)[0]
    M = Msnp[i_ch]
    thr = mult * np.median(np.abs(M) / 0.6745)

    if thr < 0:
        i_pass = np.min(M, axis=1) < thr
    else:
        i_pass = np.max(M, axis=1) > thr

    selected[i_ch[i_pass]] = True

    if return_full:
        return selected, thr
    return thr


def rethreshold_by_multiplier_par(Msnp, Msnp_ch, target_chs, \
        mult=RETHRESHOLD_MULT, n_jobs=NCPU_HIGHMEM):
    """Experimental parralelized version of rethreshold_by_multiplier()
    TODO: The performance is terrible.  Needs optimization."""
    core = rethreshold_by_multiplier_core
    r = Parallel(n_jobs=n_jobs, verbose=0)(delayed(core)(Msnp, \
            Msnp_ch, ch, mult=mult) for ch in target_chs)
    thrs = [e[1] for e in r]
    selected = np.array([e[0] for e in r]).any(axis=0)

    return selected, thrs


def rethreshold_by_multiplier(Msnp, Msnp_ch, target_chs, \
        mult=RETHRESHOLD_MULT):
    thrs = []
    selected = np.zeros(Msnp_ch.shape, dtype='bool')
    for ch in target_chs:
        thr = rethreshold_by_multiplier_core(Msnp, Msnp_ch, \
                ch, mult=mult, selected=selected)
        thrs.append(thr)

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


# -- clustering related ------------------------------------------------------
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
        nperimg=EXTRACT_NPERIMG, nmax=EXTRACT_NMAX):
    """Extract all spikes specified in ibie"""

    # not using defaultdict to save space
    res = []
    for ch in target_chs:
        res.append([])

    # sweep over each image
    for ib, ie in ibie:
        idx = range(ib, ie)
        Msnp_img = Msnp[idx]
        Msnp_ch_img = Msnp_ch[idx]

        for ch in target_chs:
            if len(res[ch]) >= nmax:
                continue
            res[ch].extend(Msnp_img[Msnp_ch_img == ch][:nperimg])

    for ch in target_chs:
        if len(res[ch]) < nmax:
            continue
        res[ch] = res[ch][:nmax]

    return res


def cluster_affinity_prop_core(feat, commonp=AFFINITYPRP_COMMONP):
    """Copied from the sklearn website"""
    X = np.array(feat)

    # -- Compute similarities
    X_norms = np.sum(X ** 2, axis=1)
    S = - X_norms[:, np.newaxis] - X_norms[np.newaxis, :] + 2 * np.dot(X, X.T)

    if commonp == '10med':
        p = 10 * np.median(S)
    elif commonp == 'min':
        p = np.min(S)
    else:
        raise ValueError('Not recognized commonp')

    # -- Compute Affinity Propagation
    af = AffinityPropagation().fit(S, p)
    cluster_centers_indices = af.cluster_centers_indices_
    labels = af.labels_
    n_clusters_ = len(cluster_centers_indices)

    return labels, n_clusters_, cluster_centers_indices


def find_clusters_par(Msnp_feat_train, feat_use, n_jobs=NCPU, \
        metd=CLUSTERING_ALG, **clu_cfg):
    assert len(Msnp_feat_train) == feat_use.shape[0]

    if metd == 'affinity_prop':
        func = cluster_affinity_prop_core
    else:
        raise ValueError('Not recognized clustering metd')

    r = Parallel(n_jobs=n_jobs, verbose=0)(delayed(func)(np.array(\
            Msnp_feat_train[i])[:, feat_use[i]], **clu_cfg) \
            for i in xrange(feat_use.shape[0]))

    labels_all = []
    nclu_all = []
    cluctr_all = []

    for r0 in r:
        labels, nclu, cluctr = r0

        labels_all.append(labels)
        nclu_all.append(nclu)
        cluctr_all.append(cluctr)

    return labels_all, nclu_all, cluctr_all


def quality_ctrl_core(X0, labels, ks_plevel=QC_KS_PLEVEL, \
        min_snr=QC_MINSNR, min_size=QC_MINSIZE):

    n_clu = len(np.unique(labels))
    #if n_clu == 1:
    #    return {0: UNSORTED}, True

    X0 = np.array(X0)
    lbl_conv = {}
    some_unsorted = False
    new_cids = set([UNSORTED])

    for i_cl in xrange(n_clu):
        inds = labels == i_cl
        N = inds.sum()

        if N < min_size:   # apparantly, this is not a cluster
            SNR = 0
            KS_passed = False
        else:
            m = X0[inds].mean(0)
            s = X0[inds].std(0, ddof=1)
            s[np.abs(s) < ATOL] = 1. / ATOL   # mute bad signals
            devs, ps, _ = KS_all(X0[inds], full=True)

            SNR = np.mean(np.abs(m) / s)
            KS_passed = np.all(ps > ks_plevel)

        if (N > min_size) and (SNR > min_snr) and KS_passed:
            # if passed the quality criteria...
            new_cid = len(new_cids)
            lbl_conv[i_cl] = new_cid
            new_cids.add(new_cid)
        else:
            some_unsorted = True
            lbl_conv[i_cl] = UNSORTED
    return lbl_conv, some_unsorted


def quality_ctrl_par(Msnp_train, labels_all, nclu_all, cluctr_all, \
        n_jobs=NCPU, **kwargs):

    func = quality_ctrl_core
    r = Parallel(n_jobs=n_jobs, verbose=0)(delayed(func)(\
            X, lbl, **kwargs) \
            for X, lbl in zip(Msnp_train, labels_all))

    for i_ch, (lbl_conv, some_unsorted) in enumerate(r):
        labels_all[i_ch] = np.array([lbl_conv[e] for e in labels_all[i_ch]])
        nclu_all[i_ch] = len(np.unique(labels_all[i_ch]))

        cluctr = cluctr_all[i_ch]
        old_cids = range(len(cluctr))
        new_cids = [lbl_conv[old_cid] for old_cid in old_cids]
        s = sorted([(new_cid, cluctr[old_cid]) for new_cid, old_cid in \
                zip(new_cids, old_cids) if new_cid != UNSORTED])
        new_cluctr = [BADIDX] + [e[1] for e in s]

        base = 1
        if some_unsorted:
            base = 0
        assert len(new_cluctr) == nclu_all[i_ch] + base
        cluctr_all[i_ch] = new_cluctr


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

    config['feat'] = {}
    config['feat']['metd'] = FEAT_METHOD
    config['feat']['kwargs'] = {'level': FEAT_WAVL_LEV}

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
        config['feat']['kwargs'].pop('level')
        raise NotImplementedError('PCA not implemented yet')

    else:
        raise ValueError('Not recognized "feat_metd"')

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
    # TODO: support when thr_sel is None
    CMsnp_selected = h5o.createCArray(h5o.root, 'Msnp_selected', \
            t_uint64, Msnp_selected.shape, filters=filters)

    CMsnp[...] = Msnp
    CMsnp_tabs[...] = Msnp_tabs
    CMsnp_ch[...] = Msnp_ch
    CMsnp_pos[...] = Msnp_pos
    CMsnp_feat[...] = Msnp_feat
    CMsnp_selected[...] = Msnp_selected

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
    h5o.createArray(meta, 'fn_inp', fn_inp)

    h5o.close()


# ----------------------------------------------------------------------------
def cluster(fn_inp, fn_out, opts):
    config = {}
    config['skimspk'] = {}
    config['skimspk']['tb0'] = SKIMSPK_TB
    config['skimspk']['te0'] = SKIMSPK_TE
    config['extract'] = {}
    config['extract']['nperimg'] = EXTRACT_NPERIMG
    config['extract']['nmax'] = EXTRACT_NMAX

    config['feat'] = {}
    config['feat']['kssort'] = FEAT_KSSORT
    config['feat']['outdim'] = FEAT_OUTDIM

    config['cluster'] = {}
    config['cluster']['metd'] = CLUSTERING_ALG
    config['cluster']['commonp'] = AFFINITYPRP_COMMONP

    config['qc'] = {}
    config['qc']['qc'] = QC
    config['qc']['kwargs'] = {}
    config['qc']['kwargs']['min_snr'] = QC_MINSNR
    config['qc']['kwargs']['ks_plevel'] = QC_KS_PLEVEL
    config['qc']['kwargs']['min_size'] = QC_MINSIZE

    # -- process opts
    # TODO: implement!!!

    # -- preps
    print '* Began clustering...'
    h5 = tbl.openFile(fn_inp)
    Msnp = h5.root.Msnp.read()
    Msnp_feat = h5.root.Msnp_feat.read()
    Msnp_ch = h5.root.Msnp_ch.read()
    Msnp_pos = h5.root.Msnp_pos.read()
    Msnp_tabs = h5.root.Msnp_tabs.read()
    Msnp_selected = h5.root.Msnp_selected.read()
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

    # -- get training examples...
    print '* Collecting snippet examples...'
    ibie, iuimg = skim_imgs(Mimg, Mimg_tabs, Msnp_tabs, \
            t_adjust, **config['skimspk'])

    clu_feat_train = get_example_spikes(Msnp_feat, Msnp_ch, ibie, all_chs, \
            **config['extract'])
    clu_train = get_example_spikes(Msnp, Msnp_ch, ibie, all_chs, \
            **config['extract'])

    # -- get feature indices to use...
    print '* Finding useful axes...'
    outdim = config['feat']['outdim']
    if config['feat']['kssort']:
        Msnp_feat_use = []

        for ch in all_chs:
            X = np.array(clu_feat_train[ch])
            # get deviations from Gaussian
            devs, _ = KS_all(X)
            # got top-n deviations
            devs = np.argsort(-devs)[:outdim]
            Msnp_feat_use.append(devs)
    else:
        Msnp_feat_use = [range(outdim)] * len(all_chs)
    Msnp_feat_use = np.array(Msnp_feat_use)

    # -- get clusters...
    print '* Clustering...'
    # DBG
    #clu_labels, clu_nclus, clu_centers = \
    #        find_clusters_par(clu_feat_train[:4], \
    #        Msnp_feat_use[:4], **config['cluster'])
    clu_labels, clu_nclus, clu_centers = \
            find_clusters_par(clu_feat_train, \
            Msnp_feat_use, **config['cluster'])

    # -- quality control
    if config['qc']['qc']:
        print '* Run quality-control screening...'
        quality_ctrl_par(clu_train, clu_labels, clu_nclus, clu_centers, \
                **config['qc']['kwargs'])

    # -- done! write everything...
    print '* Writing results...'
    filters = tbl.Filters(complevel=4, complib='blosc')
    t_uint32 = tbl.UInt32Atom()
    t_uint64 = tbl.UInt64Atom()

    h5o = tbl.openFile(fn_out, 'w')
    CMsnp_tabs = h5o.createCArray(h5o.root, 'Msnp_tabs', t_uint64, \
            Msnp_tabs.shape, filters=filters)
    CMsnp_ch = h5o.createCArray(h5o.root, 'Msnp_ch', t_uint32, \
            Msnp_ch.shape, filters=filters)
    CMsnp_pos = h5o.createCArray(h5o.root, 'Msnp_pos', t_uint64, \
            Msnp_pos.shape, filters=filters)
    CMsnp_feat_use = h5o.createCArray(h5o.root, 'Msnp_feat_use', t_uint64, \
            Msnp_feat_use.shape, filters=filters)
    CMsnp_selected = h5o.createCArray(h5o.root, 'Msnp_selected', \
            t_uint64, Msnp_selected.shape, filters=filters)

    CMsnp_tabs[...] = Msnp_tabs
    CMsnp_ch[...] = Msnp_ch
    CMsnp_pos[...] = Msnp_pos
    CMsnp_feat_use[...] = Msnp_feat_use
    CMsnp_selected[...] = Msnp_selected

    h5o.createArray(h5o.root, 'Mimg', Mimg)
    h5o.createArray(h5o.root, 'Mimg_tabs', Mimg_tabs)

    meta = h5o.createGroup('/', 'meta', 'Metadata')
    h5o.createArray(meta, 't_start0', t_start0)
    h5o.createArray(meta, 't_stop0', t_stop0)
    h5o.createArray(meta, 't_adjust', t_adjust)
    h5o.createArray(meta, 'config_clu_pk', pk.dumps(config))

    h5o.createArray(meta, 'idx2iid', idx2iid)
    h5o.createArray(meta, 'iid2idx_pk', iid2idx_pk)
    h5o.createArray(meta, 'idx2ch', idx2ch)
    h5o.createArray(meta, 'ch2idx_pk', ch2idx_pk)

    h5o.createArray(meta, 'fn_inp', fn_inp)

    clupk = h5o.createGroup('/', 'clu_pk', 'Pickles')
    h5o.createArray(clupk, 'clu_feat_train_pk', pk.dumps(clu_feat_train))
    h5o.createArray(clupk, 'clu_train_pk', pk.dumps(clu_train))
    h5o.createArray(clupk, 'clu_labels_pk', pk.dumps(clu_labels))
    h5o.createArray(clupk, 'clu_nclus_pk', pk.dumps(clu_nclus))
    h5o.createArray(clupk, 'clu_centers_pk', pk.dumps(clu_centers))

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
        cluster(args[1], args[2], opts)
    else:
        raise ValueError('Invalid mode')

    print 'Done.                                '


if __name__ == '__main__':
    main()
