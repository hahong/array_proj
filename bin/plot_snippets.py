#!/usr/bin/env python

import numpy as np
import cPickle as pk
import pylab as pl
import struct
import sys
sys.path.append('lib')
from mworks.data import MWKFile
from mergeutil import Merge
from common_fn import get_stim_info, xget_events, load_spike_data
from collections import defaultdict
from matplotlib import rc

C_STIM = '#announceStimulus'
I_STIM_ID = 2
T_START = -100000
T_STOP = 250000
DEFAULT_ELECS = range(1, 97)
C_SUCCESS = 'number_of_stm_shown'   # one visual stimulus should have one "success" event within
T_SUCCESS = 250000                  # 250ms time window in order to be considered as valid one.
PROC_CLUSTER = False
MAX_CLUS = 5                          # number of clusters per channel
CH_PER_FIG = 16
MAX_SPK = 300
IMG_MAX = 800
EXCLUDE_IMG = None                  # exclude image by its name

# ----------------------------------------------------------------------------
def plot_all(fn_mwks, fn_nevs, fn_ref, override_delay_us=None, override_elecs=None, \
        verbose=True, c_success=C_SUCCESS, extinfo=False, prefix='out',  t_success_lim=T_SUCCESS, \
        exclude_img=EXCLUDE_IMG):
    # ref = pk.load(open(fn_ref))
    # clus_info = ref['clus_info']
    # del ref
    rc('font', **{
        'family'          : 'serif',
        'serif'           : ['Times New Roman']})
    rc('figure', **{
        'subplot.left'    : '0.02',
        'subplot.right'   : '0.98',
        'subplot.bottom'  : '0.035',
        'subplot.top'     : '0.945',
        'subplot.wspace'  : '0.3',
        'subplot.hspace'  : '0.25',
    })

    # -----------------------
    clus_info = defaultdict(list)

    for fn_mwk, fn_nev in zip(fn_mwks, fn_nevs):
        print '*', fn_mwk, fn_nev
        get_data(fn_mwk, fn_nev, clus_info, override_delay_us=override_delay_us, \
                override_elecs=override_elecs, verbose=verbose, c_success=c_success, \
                t_success_lim=t_success_lim, extinfo=extinfo, exclude_img=exclude_img)

    print '* Collection done.              '
    figs = {}
    for k in clus_info: 
        ch, cid = k
        i_fig = ((ch - 1) / CH_PER_FIG) + 1
        i_subpl = (ch - 1) % CH_PER_FIG + cid * CH_PER_FIG + 1

        if i_fig not in figs: 
            figs[i_fig] = pl.figure(i_fig, figsize=(17.96, 10.25))
        else:
            pl.figure(i_fig)

        ax = pl.subplot(MAX_CLUS, CH_PER_FIG, i_subpl)
        wavs = np.array(clus_info[k]) 
        y = wavs.mean(axis=0)
        x = range(len(y))
        pl.plot(x, y, 'b-')

        if len(wavs) > 2:
            s = wavs.std(axis=0)
            pl.fill_between(x, y + s, y - s, color='b', alpha=0.2)
        pl.xlim([x[0], x[-1]])

        pl.title('(%d,%d)' % k, fontsize=7)
        fontsize = 3
        for tick in ax.xaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)
        for tick in ax.yaxis.get_major_ticks():
            tick.label1.set_fontsize(fontsize)

    # kludge
    for i in range(1, len(override_elecs)/CH_PER_FIG + 1):
        pl.figure(i)
        pl.savefig(prefix + '%02d.pdf' % i)
        pl.close()



def get_data(fn_mwk, fn_nev, clus_info, override_delay_us=None, override_elecs=None, \
        verbose=True, c_success=C_SUCCESS, t_success_lim=T_SUCCESS, max_clus=MAX_CLUS, \
        extinfo=False, imgmax=IMG_MAX, exclude_img=EXCLUDE_IMG):
    mf = MWKFile(fn_mwk)
    mf.open()
    nf = load_spike_data(fn_nev)
    nf.open()
    cnt = defaultdict(int)

    # read TOC info from the "merged" mwk file
    toc = xget_events(mf, codes=[Merge.C_MAGIC])[0].value
    c_spikes = toc[Merge.K_SPIKE]                  # get the code name from the toc

    # when the visual stimuli presented is valid?
    t_success = [ev.time for ev in xget_events(mf, codes=[c_success])]
    t_success = np.array(t_success)

    # get the active electrodes
    if override_elecs is None: actvelecs = toc[Merge.K_SPIKEWAV].keys() 
    else: actvelecs = override_elecs               # e.g, range(1, 97)

    img_onset, img_id = get_stim_info(mf, extinfo=extinfo, exclude_img=exclude_img)
    n_stim = len(img_onset)

    # MAC-NSP time translation
    if override_delay_us != None: 
        t_delay = toc['align_info']['delay']
        t_adjust = int(np.round(override_delay_us - t_delay))
    else: 
        t_adjust = 0
    t_start = T_START - t_adjust
    t_stop = T_STOP - t_adjust

    # actual calculation -------------------------------
    # all_spike[chn_id][img_id]: when the neurons spiked?
    for i in range(n_stim):
        if i > imgmax: break
        t0 = img_onset[i]; iid = img_id[i]
        # -- check if this presentation is successful. if it's not ignore this.
        if np.sum((t_success > t0) & (t_success < (t0 + t_success_lim))) < 1: continue

        if verbose: 
            print 'At', (i + 1), 'out of', n_stim, '         \r',
            sys.stdout.flush()
   
        spikes = xget_events(mf, codes=[c_spikes], time_range=[t0 + t_start, t0 + t_stop])
        for s in spikes:
            ch = s.value['id']
            cid = s.value['cluster_id']
            pos = s.value['foffset']
            key = (ch, cid)

            if cnt[key] >= MAX_SPK: continue
            cnt[key] += 1

            dat = nf.read_once(pos=pos, proc_wav=True)
            wav = dat['waveform']
            clus_info[key].append(np.array(wav))


# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 3:
        print 'plot_snippets.py <mwk1,mwk2,...> <nev1,nev2,...> <ref.psf.pk> [override delay in us] [number of electrodes] [opts]'
        print 'Plot representative snippets for each unit'
        print 
        print 'Options:'
        print 'c_success=<string>   - code name for "success" signal'
        print 'extinfo              - collects extra stimuli information in addition to the names'
        print 'prefix               - output pdf prefix'
        return

    fn_mwks = sys.argv[1].split(',')
    fn_nevs = sys.argv[2].split(',')
    fn_ref = sys.argv[3]
    opts0 = sys.argv[6:]
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
    extinfo = False

    # parsing extra arguments
    if len(sys.argv) >= 5:
        override_delay_us = long(sys.argv[4])
        if override_delay_us < 0: override_delay_us = None
        print 'Delay override:', override_delay_us
    else: override_delay_us = None

    if len(sys.argv) >= 6:
        override_elecs = range(1, int(sys.argv[5]) + 1)
        print 'Active electrodes override: [%d..%d]' % (min(override_elecs), max(override_elecs))
    else: override_elecs = DEFAULT_ELECS

    # handle "opts"
    if 'c_success' in opts:
        c_success = opts['c_success']
        print 'c_success:', c_success
    else:
        c_success = C_SUCCESS

    if 't_success' in opts:
        t_success = int(opts['t_success'])
        print '* t_success:', t_success
    else:
        t_success = T_SUCCESS

    exclude_img = EXCLUDE_IMG
    if 'exclude_img' in opts:
        exclude_img = opts['exclude_img'].split(',')
        print '* Exclude unwanted images:', exclude_img


    extinfo = False
    if 'extinfo' in opts:
        extinfo = True
        print 'Collecting extra information of the stimuli'

    prefix = 'out'
    if 'prefix' in opts:
        prefix = opts['prefix']
        print 'Output prefix:', prefix

    # go go go
    plot_all(fn_mwks, fn_nevs, fn_ref, override_delay_us=override_delay_us, \
            override_elecs=override_elecs, c_success=c_success, extinfo=extinfo, \
            prefix=prefix, t_success_lim=t_success, exclude_img=exclude_img)
    print 'Done.     '


if __name__ == '__main__':
    main()
