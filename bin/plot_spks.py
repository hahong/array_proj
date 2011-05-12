#!/usr/bin/env python

import numpy as np
import pylab as pl
import sys
sys.path.append('lib')
from mworks.data import MWKFile
from mergeutil import BRReader
from common_fn import set_new_threshold, xget_events

#START = -2000
#END = 2000
#START = 100000
#END = 105000
START = -100000
END = 200000
OVERRIDE_DELAY_US = 300
CH = 93        # 93, 114
MOVIE = False

if MOVIE:
    ONLY_THIS = '_frame_0'
    C_SUCCESS = 'success'   # one visual stimulus should have one "success" event within
    T_SUCCESS = 2500000     # 250ms time window in order to be considered as valid one.
else:
    ONLY_THIS = None
    # Nat300_269F (ladybug) -- 1-based 93, **115**, 24, 18: good | 89: so so | 29: bad  # Nat300_F_181
    # IID = 'Nat300_F_269'
    IID = ''
    T_SUCCESS = 250000     # 250ms time window in order to be considered as valid one.
    C_SUCCESS = 'number_of_stm_shown' 

NPERIOD = 1000. / 85.
MAX_SPK = 500
DIM = 48
GET_ALL_NOISE = False

PLOT_ONLY_DOWN = False
PLOT_ONLY_EXCHG = False
PLOT_ONLY_ABV = False #1.133

C_STIM = '#announceStimulus'
#PLOT_ONLY_DOWN = True
#PLOT_ONLY_EXCHG = True
#PLOT_ONLY_ABV = 1.133


# ---------------
def myhist(x, bins=10, norm=None, cutoff=-100):
    x = x[x > cutoff]
    y0, x0 = np.histogram(x, bins=bins)
    x = (x0[1:] + x0[:-1]) / 2.

    if norm is None: y = y0
    elif norm == 'peak': y = y0 / float(y0.max())

    return x, y

# ----------------------------------------------------------------------------
def get_stim_info(mf, c_stim=C_STIM, extinfo=False, only_this=ONLY_THIS):
    stims = mf.get_events(codes=[c_stim])
    if only_this is None:
        stims = [x for x in stims if type(x.value) != int and x.value['type'] == 'image' and IID in x.value['name'] and 'blank' not in x.value['name'] and 'mask' not in x.value['name']]
    else:
        stims = [x for x in stims if type(x.value) != int and x.value['type'] == 'image_directory_movie' and \
                type(x.value['current_stimulus']) != int and only_this in x.value['current_stimulus']['name']]
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


def main(c_success=C_SUCCESS):
    #mf = MWKFile('../analysis/data_merged/Chabo_20110426_MovieGallant110413_S110204_RefAA_001.mwk')
    mf = MWKFile('../analysis/data_merged/Chabo_20110331_RSVPNicole305_S110204_RefAA_001.mwk')
    mf.open()
    #br = BRReader('../analysis/data_nev/Chabo_20110426_MovieGallant110413_S110204_RefAA_001.nev')
    br = BRReader('../analysis/data_nev/Chabo_20110331_RSVPNicole305_S110204_RefAA_001.nev')
    br.open()

    if PLOT_ONLY_ABV:
        adj_reject = 5./3.
        new_thr = {}
        for ch in br.chn_info:
            lthr = br.chn_info[ch]['low_thr']; hthr = br.chn_info[ch]['high_thr']
            if lthr == 0: thr0 = hthr
            else: thr0 = lthr
            new_thr[ch] = thr0 * adj_reject * 0.249
            print '*', ch, new_thr[ch]

    toc = xget_events(mf, codes=['#merged_data_toc'])[0].value
    # MAC-NSP time translation
    if OVERRIDE_DELAY_US != None: 
        t_delay = toc['align_info']['delay']
        t_adjust = int(np.round(OVERRIDE_DELAY_US - t_delay))
    else: 
        t_adjust = 0

    # when the visual stimuli presented is valid?
    t_success = [ev.time for ev in mf.get_events(codes=[c_success])]
    t_success = np.array(t_success)

    img_onset, img_id = get_stim_info(mf)
    print 'Len =', len(img_onset)

    i_plot_ac = 0
    i_spk = 0
    i_spk2 = 0
    i_spk_nvis = 0
    i_spk_vis = 0
    M = np.zeros((MAX_SPK, DIM))
    M2 = np.zeros((MAX_SPK, DIM))
    Mnvis = np.zeros((MAX_SPK, DIM))
    Mvis  = np.zeros((MAX_SPK, DIM))
    t = (np.arange(DIM) - 11) * 1000. / 30000.   # in ms
    #t = np.arange(len(DIM))    # in index

    for t0 in img_onset:
        if i_spk >= MAX_SPK: break
        if i_spk_nvis >= MAX_SPK: break
        if i_spk_vis >= MAX_SPK: break
        if np.sum((t_success > t0) & (t_success < (t0 + T_SUCCESS))) < 1: continue

        spks = mf.get_events(codes=['merged_spikes'], time_range=[t0 +START, t0 + END])

        for spk in spks:
            ch = spk.value['id']
            ts = spk.time
            if ch != CH: continue

            offset = spk.value['foffset']
            wav = br.read_once(pos=offset, proc_wav=True)['waveform']
            y = np.array(wav) * 0.249  # in uV

            if PLOT_ONLY_DOWN:
                if y[12] > y[11]: continue
            if PLOT_ONLY_ABV:
                wav = set_new_threshold(wav, new_thr[ch], rng=(11, 13), i_chg=32) 
                if wav == None: continue
            if PLOT_ONLY_EXCHG:
                if np.max(y[:32]) < 0: continue

            t_rel = ts + t_adjust - t0
            # print t_rel
            # -- monitor noise?
            if (t_rel/1000.) % NPERIOD < 1. or (t_rel/1000.) % NPERIOD > (NPERIOD - 1):
                if t_rel > -50000 and t_rel < 50000:
                    M[i_spk] = y
                    i_spk += 1
                    pl.figure(ch)
                    pl.plot(t, y, 'k-')
                    pl.title('Noise responses (OFF region)')
                    pl.xlabel('Time/ms')
                    pl.ylabel(r'Response/$\mu$V')
                elif t_rel > 70000 and t_rel < 170000:
                    M2[i_spk2] = y
                    i_spk2 += 1
                    pl.figure(1000 + ch)
                    pl.plot(t, y, 'k-')
                    pl.title('Noise responses (ON region)')
                    pl.xlabel('Time/ms')
                    pl.ylabel(r'Response/$\mu$V')

            elif (t_rel/1000.) % NPERIOD > 2. and (t_rel/1000.) % NPERIOD < (NPERIOD - 2):
                if t_rel > -50000 and t_rel < 50000:
                    Mnvis[i_spk_nvis] = y
                    i_spk_nvis += 1
                    pl.figure(2000 + ch)
                    pl.plot(t, y, 'k-')
                    pl.title('Non-noise blank responses')
                    pl.xlabel('Time/ms')
                    pl.ylabel(r'Response/$\mu$V')
                elif t_rel > 70000 and t_rel < 170000:
                    Mvis[i_spk_vis] = y
                    i_spk_vis += 1
                    pl.figure(3000 + ch)
                    pl.plot(t, y, 'k-')
                    pl.title('Non-noise visual responses')
                    pl.xlabel('Time/ms')
                    pl.ylabel(r'Response/$\mu$V')

            i_plot_ac += 1

    M = M[:i_spk]
    # pl.figure()
    # pl.hist(np.ravel(M[:,12] - M[:,11]), bins=20)
    print 'i_spk =', i_spk
    print 'i_spk2 =', i_spk2
    print 'i_spk_nvis =', i_spk_nvis
    print 'i_spk_vis =', i_spk_vis
    print 'i_plot_ac =', i_plot_ac
    M = M[:i_spk,:]
    M2 = M2[:i_spk2,:]
    Mnvis = Mnvis[:i_spk_nvis,:]
    Mvis  = Mvis[:i_spk_vis,:] 

    pl.figure()
    xb, yb = myhist(np.min(M,axis=1), norm='peak')
    xb2, yb2 = myhist(np.min(M2,axis=1), norm='peak')
    xg, yg = myhist(np.min(Mnvis,axis=1), norm='peak')
    xv, yv = myhist(np.min(Mvis,axis=1), norm='peak')
    pl.plot(xb, yb, 'r-', label='Noise responses (OFF)')
    pl.plot(xb2, yb2, 'm-', label='Noise responses (ON)')
    pl.plot(xg, yg, 'g-', label='Non-noise responses (OFF)')
    pl.plot(xv, yv, 'b-', label='Non-noise responses (ON)')
    #pl.axvline(ptbad.mean(), color='r',alpha=0.3)
    #pl.axvline(ptgood.mean(), color='b',alpha=0.3)
    #pl.axvline(ptvis.mean(), color='g',alpha=0.3)
    pl.xlabel(r'Peak response/$\mu$V')
    pl.ylabel('Normalized probability')
    pl.legend(loc='upper left')

    pl.figure()
    # --
    m = M.mean(axis=0); s = M.std(axis=0, ddof=1)
    pl.plot(t, m, 'r-', label='Noise responses (OFF)')
    pl.fill_between(t, m-s, m+s, color='r', alpha=0.2)
    # --
    m = Mnvis.mean(axis=0); s = Mnvis.std(axis=0, ddof=1)
    pl.plot(t, m, 'g-', label='Non-noise responses (OFF)')
    pl.fill_between(t, m-s, m+s, color='g', alpha=0.2)
    # -- 
    pl.xlabel('Time/ms')
    pl.ylabel(r'Response/$\mu$V')
    pl.legend(loc='lower right')

    pl.figure()
    # --
    m = M2.mean(axis=0); s = M2.std(axis=0, ddof=1)
    pl.plot(t, m, 'm-', label='Noise responses (ON)')
    pl.fill_between(t, m-s, m+s, color='m', alpha=0.2)
    # --
    m = Mvis.mean(axis=0); s = Mvis.std(axis=0, ddof=1)
    pl.plot(t, m, 'b-', label='Non-noise responses (ON)')
    pl.fill_between(t, m-s, m+s, color='b', alpha=0.2)
    # --
    pl.xlabel('Time/ms')
    pl.ylabel(r'Response/$\mu$V')
    pl.legend(loc='lower right')

    pl.show()

main()

