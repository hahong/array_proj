#!/usr/bin/env python

import sys
import os
import Image
import cPickle as pk
import pylab as pl
import numpy as np
import scipy.stats as st
import scipy.optimize as opt
#from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rc, cm
from collections import defaultdict
from expr_anal import get_one_ch, get_stats, plot_one_PSTH
from joblib import Memory

NARRAYS = 3   # A, M, P
BLANKS = ['OS_53', 'Nat300_101', 'Nat300_102', 'Nat300_103', 'Nat300_104',
          'Nat300_105', 'Var03_2562', 'Var03_2563', 'Var03_2564', 'Var03_2565', 'Var03_2566', 
          'Var03_2567', 'Var03_2568', 'Var03_2569', 'Var03_2570', 'Var00_641', 'Var00_642', 
          'Var00_643', 'Var00_644', 'Var00_645', 'Var00_646', 'Var00_647', 'Var00_648', 'Var00_649', 'Var00_650', 
          'Chou_stimuli_1_0.0_0.0_0.0_2.0_2.0', 'Chou_stimuli_1_0.0_0.0_0.0_4.0_4.0', 'Chou_stimuli_1_0.0_0.0_0.0_8.0_8.0', 'Chou_stimuli_1_2.0_0.0_0.0_4.0_4.0', 
          'Chou_stimuli_1_4.0_0.0_0.0_4.0_4.0', 'Nicole10x6_61', 'Nicole10x6_62', 'Nat300_F_301', 'Nat300_F_302', 'Nat300_F_303', 'Nat300_F_304', 'Nat300_F_305']
DEF_ON_RANGE_US = [75000, 175000]
IMGDIR = {'Nat300':'../../../Desktop/ChaboSetup/NATURAL300/Nat300_%s.png', 'OS':'../../../Desktop/ChaboSetup/OnlineSwap_Group/OSImage_%s.png'}
mem = Memory(cachedir='.tmp_chg_anal', verbose=0)

# ----------------------------------------------------------------------------
def read_header(f):
    import csv
    csvr = csv.reader(open(f, 'rb'))
    left = []
    for i, row in enumerate(csvr):
        if i == 0: top = row
        left.append(row[0])
    return top, left

# get the slope with errors-in-variables model regression
def slope_err_var_regr(x, y):
    num = np.sum((x - x.mean()) * (y - y.mean())**2)
    dnm = np.sum((x - x.mean())**2 * (y - y.mean()))
    return num / dnm


def load_spk(fn_pk):
    data = pk.load(open(fn_pk))
    # all_spike[`channel`][`iid`] = [[trial1_spike, trial1_spike, ...], ...]
    all_spike = data['all_spike']
    return all_spike


def _load_all_spk(files, i_128):
    all_spk = {}
    for el in i_128: all_spk[el] = defaultdict(list)

    for f in files:
        dir = os.path.dirname(f); f = os.path.basename(f) 
        animal, date, prc, arr, trial = f.split('.')[0].split('_')

        # d' data:
        if arr == 'A': offset = 0
        elif arr == 'M': offset = 96
        elif arr == 'P': offset = 192
        else: raise ValueError, 'spurious array.'

        # f_pk: pk file name: e.g., Chabo_20100818_RSVPNicole_M_001.psf.pk
        f_pk = os.path.join(dir, '%s_%s_%s_%s_%s.psf.pk' % \
                (animal, date, prc, arr, trial))   
        spk = load_spk(f_pk)
        for el in spk:
            # 0-based unified electrode id across all A, M, P arrays
            el_unified = el + offset - 1
            if not el_unified in i_128: continue
            for iid in spk[el]:
                all_spk[el_unified][iid].extend(spk[el][iid])
    return all_spk
load_all_spk = mem.cache(_load_all_spk)

# ----------------------
def anal_one_expr(f, all_iids, skipcols=1, maxp=0.1):
    data = np.loadtxt(f, delimiter=',', skiprows=1)[:,skipcols:]
    maxn = int(np.round(len(all_iids) * maxp))
    # we want to have an array sorted by absolute value.
    sorted_data = []
    ir_all = []
    for r0 in data:
        ra = [(np.abs(e), e) for e in r0]
        ra0 = [np.abs(e) for e in r0]
        r = sorted(ra)
        ir = np.argsort(ra0)
        sorted_data.append([e for _, e in r])
        ir_all.append(ir)
    sorted_data = np.array(sorted_data)

    avg = data.mean(axis=1)
    std = data.std(axis=1, ddof=1)
    # sorted[:,-maxn:] gives last `maxn` columnes
    avg_maxn = np.abs(sorted_data[:,-maxn:]).mean(axis=1)
    std_maxn = np.abs(sorted_data[:,-maxn:]).std(axis=1, ddof=1)
    avg_maxr = sorted_data[:,-maxn:].mean(axis=1)
    std_maxr = sorted_data[:,-maxn:].std(axis=1, ddof=1)

    return avg, std, avg_maxn, std_maxn, avg_maxr, std_maxr, ir_all, data


def load_all_dps(files, ydate):
    top, left = read_header(files[0])
    all_elecids = left[1:]; n_row = len(all_elecids)
    dates = {'A':[], 'M':[], 'P':[]}
    # all_dps[0-based ch: 0 to 287] = {`iid`: [d', d', ...], ...}
    all_dps = {}
    for el in range(NARRAYS * n_row): all_dps[el] = defaultdict(list)

    # collect all data
    Amaxn = np.empty((n_row, 0))
    Mmaxn = np.empty((n_row, 0))
    Pmaxn = np.empty((n_row, 0))
    Amaxr = np.empty((n_row, 0))
    Mmaxr = np.empty((n_row, 0))
    Pmaxr = np.empty((n_row, 0))
    label = []

    for f in files:
        top, left = read_header(f)
        all_iids = top[1:]        # do not mess-up with the order
        avg, std, avg_maxn, std_maxn, avg_maxr, std_maxr, used, all_dp = anal_one_expr(f, all_iids)
        dir = os.path.dirname(f); f = os.path.basename(f) 
        animal, date, prc, arr, trial = f.split('.')[0].split('_')
        dates[arr].append(date)
        label.append(date)

        # d' data:
        if arr == 'A': 
            Amaxn = np.c_[Amaxn, avg_maxn]
            Amaxr = np.c_[Amaxr, avg_maxr]
            offset = 0
        elif arr == 'M': 
            Mmaxn = np.c_[Mmaxn, avg_maxn]
            Mmaxr = np.c_[Mmaxr, avg_maxr]
            offset = 96
        elif arr == 'P': 
            Pmaxn = np.c_[Pmaxn, avg_maxn]
            Pmaxr = np.c_[Pmaxr, avg_maxr]
            offset = 192
        else: raise ValueError, 'spurious array.'

        for r in range(all_dp.shape[0]):
            for c in range(all_dp.shape[1]):
                all_dps[r + offset][all_iids[c]].append(all_dp[r,c])

    return all_dps, Amaxn, Mmaxn, Pmaxn, Amaxr, Mmaxr, Pmaxr, dates

# ----------------------------------------------------------------------------
def _prepare_map(files, i_128, rng):
    # all_spk[0-based ch: 0 to 287] = {`iid`: [[t1_trial1, t_2trial1, ...], [...], ...], ...}
    all_spk = load_all_spk(sorted(files), i_128)   # due to the cache, it's better to be here.

    n_col = len(set(all_spk[i_128[0]].keys()) - set(BLANKS))
    n_row = len(i_128)
    data = np.zeros((n_row, n_col))

    for r, i in enumerate(i_128):
        c = -1
        el = i
        for iid in all_spk[el]:
            if iid in BLANKS: continue
            c += 1
            m, s, _ = get_stats(all_spk[el][iid], rng)
            data[r, c] = m
        M = []
        for iid in BLANKS:
            m, s, _ = get_stats(all_spk[el][iid], rng)
            M.append(m)
        data[r,:] -= np.mean(M)   # evoked
    return data
prepare_map = mem.cache(_prepare_map)

def plot_map(files, i_128, dp_abs_sorted, rng=DEF_ON_RANGE_US, save=True, sel=None, aspect='auto'):
    data0 = prepare_map(sorted(files), i_128, rng)
    if sel is None:
        data = data0
    else:
        data = np.zeros((data0.shape[0], len(sel)))
        for c, c0 in enumerate(sel):
            data[:,c] = data0[:,c0]

    d_rng = [2.5, 2., 1.5, 1., 0.5]
    hlines = []
    for i in range(len(i_128) - 1):
        for dp in d_rng:
            if dp_abs_sorted[i] < dp and dp < dp_abs_sorted[i+1]:
                hlines.append(len(i_128) - i - 0.688)
    
    # convert to Hz
    ws = (rng[1] - rng[0]) / 1000000.
    data /= ws

    # evoked
    pl.figure()
    pl.imshow(data, interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Evoked spike rate/Hz')
    if save: pl.savefig('fig2_raw.pdf')

    # evoked normalized
    pl.figure()
    pl.imshow(data/np.max(data), interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Normalized evoked spike rate')
    if save: pl.savefig('fig2_norm.pdf')

    # |evoked|
    pl.figure()
    pl.imshow(np.abs(data), interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Absolute evoked spike rate/Hz')
    if save: pl.savefig('fig2_raw_abs.pdf')

    # |evoked normalized|
    pl.figure()
    pl.imshow(np.abs(data)/np.max(data), interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Normalized absolute evoked spike rate')
    if save: pl.savefig('fig2_norm_abs.pdf')

    # evoked normalized by each site
    pl.figure()
    data_indv = np.array(data)
    data_max = np.max(data, axis=1)
    for r in range(data_indv.shape[0]): data_indv[r,:] /= data_max[r]

    im = pl.imshow(data_indv, interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), vmin=-1, aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Normalized evoked spike rate')
    #cb = pl.colorbar(im, ticks=[-4, 0, 0.25, 0.5, 0.75, 1])
    #cb.ax.set_yticklabels([0, 0.25, 0.5, 0.75, 1])
    if save: pl.savefig('fig2_norm_indv.pdf')

    # |evoked normalized| by each site
    pl.figure()
    data_indv = np.abs(data)
    data_max = np.max(data_indv, axis=1)
    for r in range(data_indv.shape[0]): data_indv[r,:] /= data_max[r]
    pl.imshow(data_indv, interpolation='nearest', origin='lower', extent=(0.5, data.shape[1] + 0.5, data.shape[0] + 0.5, 0.5), aspect=aspect)
    pl.xlabel('Stimulus ID')
    pl.ylabel(r"$\langle |d^{'}|_\mathrm{10\%} \rangle$ rank")
    [pl.axhline(y=h, color='k', lw=1.2) for h in hlines]
    cb = pl.colorbar()
    cb.set_label('Normalized absolute evoked spike rate')
    if save: pl.savefig('fig2_norm_abs_indv.pdf')

# ----------------------------------------------------------------------------
def plot_sc_once(dat, dp, rng, nsub=10, ntr=60, nticks=4, agg=10, imgdir=IMGDIR):
    # get good iids
    dps = sorted([(np.mean(np.abs(dp[iid])), iid) for iid in dp], reverse=True)
    giids = [iid for _, iid in dps[:nsub]] 
    
    # blank information
    n_tr = 0
    trials = []
    for iid in BLANKS:
        for tr in dat[iid]: 
            if len(tr) < 2: continue
            trials.extend(tr)
            n_tr += 1
    pl.figure()   # XXX: quick-and-dirty
    hinfo = plot_one_PSTH(trials, n_tr, pl.gca(), rng=rng, xrng=rng, aggregate=agg,
            nticks=nticks, visible=True, ymax=300)
    xarr = []; yarr = []
    for i, y in enumerate(hinfo[0]):
        yarr.append(y); yarr.append(y)
        xarr.append(hinfo[1][i]); xarr.append(hinfo[1][i + 1])
    pl.close()

    # actual plotting
    for i_plot, iid in enumerate(giids):
        ax = pl.subplot(3, nsub, i_plot + 1)
        y = 1
        for trial in dat[iid][:ntr]:
            tr = np.array(trial) / 1000.
            tr = tr[(tr >= rng[0]) & (tr <= rng[1])]
            pl.plot(tr, [y] * len(tr), 'k.', markersize=2, alpha=0.3)
            y += 1
        pl.axvline(x=0, color='r', lw=0.5)
        ax.xaxis.set_major_locator(pl.MaxNLocator(nticks)) 
        ax.yaxis.set_major_locator(pl.NullLocator()) 

    # PSTH
    for i_plot, iid in enumerate(giids):
        ax = pl.subplot(3, nsub, i_plot + 1 + nsub)
        trials = []
        for tr in dat[iid]: trials.extend(tr)
        plot_one_PSTH(trials, len(dat[iid]), ax, rng=rng, xrng=rng, aggregate=agg,
                nticks=nticks, visible=True, ymax=300)
        pl.plot(xarr, yarr, color='g', lw=0.6, alpha=0.4)
        if i_plot == 0:
            pl.ylabel('Spikes/s')
        #print '   - PSTH:', i_plot

    # inlet image
    for i_plot, iid in enumerate(giids):
        ax = pl.subplot(3, nsub, i_plot + 1 + nsub*2, frameon=False)
        ax.xaxis.set_major_locator(pl.NullLocator()) 
        ax.yaxis.set_major_locator(pl.NullLocator()) 
        type = iid.split('_')[0]
        fimg = imgdir[type] % iid.split('_')[1]
        for line in ax.get_xticklines():
            line.set_visible(False)
        for line in ax.get_yticklines():
            line.set_visible(False)

        #print fimg
        if os.path.isfile(fimg):
            img = Image.open(fimg)
            # if img.mode != 'L': img = img.convert('L')
            pl.imshow(img, origin='lower', vmin=0, vmax=255, cmap=cm.gray)
            pl.xlabel(iid)

def plot_sc(files, all_dps, target, rng=[-50, 250], save=True):
    data = load_all_spk(sorted(files), target)   # due to the cache, it's better to be here.
        
    for el in target:
        pl.figure(figsize=(16, 5))
        plot_sc_once(data[el], all_dps[el], rng)
        pl.subplots_adjust(left=0.04, right=0.99, top=0.98, bottom=0.06, wspace=0.33, hspace=0.4)
   
        if save:
            if el < 96: arr = 'A'; el0 = el + 1
            elif el < 192: arr = 'M'; el0 = el - 95
            else: arr = 'P'; el0 = el - 191
            fsav = 'fig3_%s_%03d.pdf' % (arr, el0)
            pl.savefig(fsav)
            print '   + %s -> %s' % (el, fsav)
            pl.close()

# ----------------------------------------------------------------------------
def plot_dvd(dates, A, M, P, i_128, ydate='20100827'):
    def get_stats_all(dates, data, ydate):
        ix = []; iy = []
        for i, date in enumerate(dates):
            if date == ydate: iy.append(i)
            else: ix.append(i)
        x0 = data[:,ix]; x = x0.mean(axis=1)
        y0 = data[:,iy]; y = y0.mean(axis=1)
        xe = st.sem(x0, axis=1)
        ye = st.sem(y0, axis=1)
        return x, y, xe, ye, x0, y0
    X = np.empty(0); XE = np.empty(0); X0 = []
    Y = np.empty(0); YE = np.empty(0); Y0 = []

    for arr, data in zip(['A', 'M', 'P'], [A, M, P]):
        x, y, xe, ye, x0, y0 = get_stats_all(dates[arr], data, ydate)
        X = np.r_[X, x]; Y = np.r_[Y, y];
        XE = np.r_[XE, xe]; YE = np.r_[YE, ye]
        # X0[arr] = x0; Y0[arr] = y0
        for r in range(x0.shape[0]):
            for x00 in x0[r]:
                for y00 in y0[r]:
                    X0.append(x00)
                    Y0.append(y00)
    X0 = np.array(X0)
    Y0 = np.array(Y0)

    i_n128 = list(set(range(len(X))) - set(i_128))
    pl.figure()
    pl.axvline(x=0, c='#cccccc')
    pl.axhline(y=0, c='#cccccc')
    pl.plot([-5,8],[-5,8], color='#eaeaea')
    pl.errorbar(X[i_128], Y[i_128], xerr=XE[i_128], yerr=YE[i_128], c='r',
                capsize=1, ls='None', lw=0.6)
    pl.errorbar(X[i_n128], Y[i_n128], xerr=XE[i_n128], yerr=YE[i_n128], c='b',
                capsize=1, ls='None', lw=0.6)
    pl.xlabel(r"$\langle d^{'}_\mathrm{10\%} \rangle_\mathrm{\cdots, D\!-\!2}$")
    pl.ylabel(r"$\langle d^{'}_\mathrm{10\%} \rangle_\mathrm{D\!-\!1}$")
    pl.xlim([-3,5])
    pl.ylim([-3,5])

    # rank order
    rho, p = st.spearmanr(X, Y)
    # slope with errors-in-variables model
    s = slope_err_var_regr(X0, Y0)

    pl.title(r"$\rho = %1.4f; \ s=%1.4f$" % (rho, s))

    # which one is the outlier?
    # XXX: not used anymore.
    """
    def transf(iter):
        r = []
        for i in iter:
            if i <= 95: r.append(('A', i+1))
            elif i <= 191: r.append(('M', i-95))
            elif i <= 287: r.append(('M', i-191))
            else: raise ValueError, '!!!'
        return r

    l_high = transf(np.nonzero((X > 2.5) & (Y > 2.5))[0])
    l_low = transf(np.nonzero((X < -0.3) & (Y < -0.3))[0])
    l_highodd = transf(np.nonzero((X < 1) & (Y > 3))[0])
    l_lowodd = transf(np.nonzero((X > 2.7) & (Y < 1.3))[0])

    print '* high:', l_high 
    print '* low:', l_low
    print '* high odd:', l_highodd
    print '* low odd:', l_lowodd"""


# ----------------------------------------------------------------------------
def plot_dvi(all_dps, i_128, dp_abs_sorted, \
             rng=[(2.5, 6.), (2., 2.5), (1.5, 2.), (1., 1.5), (0.8, 1)], capsz=2):
    n_splot = len(rng)
    pl.figure()
    for i_splot in range(n_splot):
        pl.subplot(n_splot, 1, i_splot + 1)
        # "collapsed" d's for the range. dp_coll[`iid`] = [d', d', ...]
        dp_coll = defaultdict(list)
        ch_rng = i_128[(dp_abs_sorted >= rng[i_splot][0]) \
                       & (dp_abs_sorted < rng[i_splot][1])]
        for ch in ch_rng:
            for iid in all_dps[ch]:
                dp_coll[iid].extend(all_dps[ch][iid])

        stims = list(set(dp_coll.keys()) - set(BLANKS))
        iid_sorted0 = sorted([(iid.split('_')[0], int(iid.split('_')[1]), iid) for iid in stims])
        iid_sorted = [iid[-1] for iid in iid_sorted0]

        X = np.array(range(len(iid_sorted))) + 0.65
        Y = []
        E = []
        for iid in iid_sorted:
            Y.append(np.mean(np.abs(dp_coll[iid])))
            E.append(st.sem(dp_coll[iid]))
        Y = np.array(Y)
        E = np.array(E)

        i_top10 = np.argsort(-Y)[:10]
        i_ntop10 = np.argsort(-Y)[10:]

        pl.bar(X[i_ntop10], Y[i_ntop10], yerr=E[i_ntop10], width=0.7, color='#222222', edgecolor='#222222',
               ecolor='#888888', capsize=capsz)
        pl.bar(X[i_top10], Y[i_top10], yerr=E[i_top10], width=0.7, color='r', edgecolor='r',
               ecolor='#ff7777', capsize=capsz)
        if i_splot == 0: rtxt = r"\langle |d^{'}|_\mathrm{10\%%} \rangle \,\geq\, %1.1f" % rng[i_splot][0]
        else: rtxt = r"%1.1f \, \leq \, \langle |d^{'}|_\mathrm{10\%%} \rangle \, < \, %1.1f" % rng[i_splot]
        pl.text(X[-1] - 0.2, 2.37, '$%s; \ n_\mathrm{sites} = %d$' % (rtxt, len(ch_rng)),
                size=13, va='top', ha='right')
        pl.ylabel(r"$\langle|d^{'}|\rangle$")
        pl.xlim([0.1, len(iid_sorted) + 1])
        pl.ylim([0, 2.5])

    pl.subplots_adjust(left=0.04, right=0.99, top=0.985, bottom=0.02, hspace=0.33)

# ----------------------------------------------------------------------------
def plot_dvi_indv(all_dps, sel, sel_dp_abs, capsz=2, sqz=False, vabs=False):
    n_splot = len(sel)
    pl.figure()
    for i_splot, el in enumerate(sel):
        ax = pl.subplot(n_splot, 1, i_splot + 1)

        stims = list(set(all_dps[el].keys()) - set(BLANKS))
        iid_sorted0 = sorted([(iid.split('_')[0], int(iid.split('_')[1]), iid) for iid in stims])
        iid_sorted = [iid[-1] for iid in iid_sorted0]

        X = np.array(range(len(iid_sorted))) + 0.65
        Y = []
        E = []
        for iid in iid_sorted:
            if vabs:
                Y.append(np.mean(np.abs(all_dps[el][iid])))   # XXX
            else:
                Y.append(np.mean(all_dps[el][iid]))   # XXX
            E.append(st.sem(all_dps[el][iid]))
        Y = np.array(Y)
        E = np.array(E)

        i_top10 = np.argsort(-np.abs(Y))[:10]
        i_ntop10 = np.argsort(-np.abs(Y))[10:]

        pl.bar(X[i_ntop10], Y[i_ntop10], yerr=E[i_ntop10], width=0.7, color='#222222', edgecolor='#222222',
               ecolor='#888888', capsize=capsz)
        pl.bar(X[i_top10], Y[i_top10], yerr=E[i_top10], width=0.7, color='r', edgecolor='r',
               ecolor='#ff7777', capsize=capsz)

        rtxt = r"\langle |d^{'}|_\mathrm{10\%%} \rangle \,=\, %1.2f" % sel_dp_abs[i_splot]
        if el < 96: arr = 'A'; el0 = el + 1
        elif el < 192: arr = 'M'; el0 = el - 95
        else: arr = 'P'; el0 = el - 191
        pl.text(X[-1] + 0.9, 5.95, 'Channel %s%03d: $%s$' % (arr, el0, rtxt),
                size=10, va='top', ha='right')

        if vabs:
            pl.ylabel(r"$\langle|d^{'}|\rangle$")
            ax.xaxis.set_major_locator(pl.MaxNLocator(15)) 
            ax.yaxis.set_major_locator(pl.MaxNLocator(3)) 
            
            if sqz:  # XXX: quick-and-dirty..
                if i_splot == 0:
                    ax.set_xticklabels([''] * 15)
                    ax.set_yticklabels(['', '2', '4', '6'])
                elif i_splot != n_splot - 1:
                    ax.set_xticklabels([''] * 15)
                    ax.set_yticklabels(['', '2', '4', ''])
                elif i_splot == n_splot - 1:
                    ax.set_yticklabels(['0', '2', '4', ''])

            pl.xlim([0.1, len(iid_sorted) + 1])
            pl.ylim([0, 6])
        else:
            pl.ylabel(r"$\langle d^{'} \rangle$")
            pl.axhline(y=0, color='#999999')
            ax.xaxis.set_major_locator(pl.MaxNLocator(15)) 
            ax.yaxis.set_major_locator(pl.MaxNLocator(4)) 

            if sqz:  # XXX: quick-and-dirty..
                if i_splot == 0:
                    ax.set_xticklabels([''] * 15)
                    ax.set_yticklabels(['', '0', '2', '4', '6'])
                elif i_splot != n_splot - 1:
                    ax.set_xticklabels([''] * 15)
                    ax.set_yticklabels(['', '0', '2', '4', ''])
                elif i_splot == n_splot - 1:
                    ax.set_yticklabels(['-2', '0', '2', '4', ''])
            
            pl.xlim([0.1, len(iid_sorted) + 1])
            pl.ylim([-2, 6])

    if sqz:
        pl.subplots_adjust(left=0.022, right=0.995, top=0.995, bottom=0.02, hspace=0.001)
    else:
        pl.subplots_adjust(left=0.022, right=0.995, top=0.995, bottom=0.02, hspace=0.33)


# ----------------------------------------------------------------------------
def get_selected128_abs(Amaxn, Mmaxn, Pmaxn):
    m_A = Amaxn.mean(axis=1); m_M = Mmaxn.mean(axis=1); m_P = Pmaxn.mean(axis=1)
    M = np.r_[m_A, m_M, m_P]; i_128 = np.argsort(np.abs(M))[-128:]
    s_A = i_128[i_128 <= 95]; s_M = i_128[(i_128 > 95) & (i_128 <= 191)]; s_P = i_128[i_128 > 191]
    print '* selected A (0-based; n=%d):' % len(s_A), s_A
    print '* selected M (0-based; n=%d):' % len(s_M), s_M
    print '* selected P (0-based; n=%d):' % len(s_P), s_P
    print "* selected 128 channel's <|d'| of max 10%>:", np.abs(M)[i_128]
    print "* selected 128:", i_128
    return i_128, np.abs(M)[i_128]

# ----------------------------------------------------------------------------
def init_rc(opts=None):
    # init plot fonts
    rc('font', **{
        'family'          : 'serif',
        'serif'           : ['Times New Roman']})

def main():
    if len(sys.argv) < 3:
        print 'chg_anal.py <last date pattern (e.g., 20100827)> <file 1.csv> [file 2.csv] ...'
        print 'chg_anal.py: analyze some statistics across several days.'
        return

    ydate = sys.argv[1]
    files = sys.argv[2:]

    # all_dps[0-based ch: 0 to 287] = {`iid`: [d', d', ...], ...}
    all_dps, Amaxn, Mmaxn, Pmaxn, Amaxr, Mmaxr, Pmaxr, dates = load_all_dps(sorted(files), ydate)
    # i_128: index of top-128, dp_abs_sorted: actual d' values of the 128
    i_128, dp_abs_sorted = get_selected128_abs(Amaxn, Mmaxn, Pmaxn)
    print '* begin plotting.'

    # plot them
    init_rc()
    
    # NOTE: comment and uncomment for the plots. ----------------------------------------
    # XXX: kinda kludge. should clean up..

    # Fig 1:
    #plot_dvi(all_dps, i_128, dp_abs_sorted)

    # Selection of the channel for the PSTH and individual d'-vs-image ID plot (Fig 1)
    #sel = [127, 126, 125, 123, 118, 105, 103, 98, 78, 53, 28, 3, 0]
    #sel = [127, 126, 125, 118, 105, 103, 98, 78, 68, 53, 28, 3, 0]
    #plot_dvi_indv(all_dps, i_128[sel], dp_abs_sorted[sel])

    # Fig 2:
    # Plotting all d's
    # plot_map(files, i_128, dp_abs_sorted)
    # Plotting for some Image ID's
    #plot_map(files, i_128, dp_abs_sorted, sel=[104, 113, 114, 115, 119, 131, 132, 137, 144, 150], aspect=0.2)

    # Fig 3:
    #plot_sc(files, all_dps, i_128[sel])

    # Figure 4:
    # plot_dvd(dates, Amaxr, Mmaxr, Pmaxr, i_128, ydate=ydate)
    # pl.show()


if __name__ == '__main__':
    main()
