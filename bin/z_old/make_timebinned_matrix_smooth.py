#! /usr/bin/env python

import cPickle as pk
import scipy as sp
from scipy import io
import numpy as np
import sys

STIMONDUR = 100000   # stimuls on duration

def get_movie_mat(fsorted, fout, tw=100000, iid_blank=101, stminfo=None, stimondur=STIMONDUR):
    # fsorted: sorted pk file. the "temp_data.pk" thing
    # fout: output file name
    # tw: timewindow size in us
    # ii_blank: image id assigned for the off
    sdb = pk.load(open(fsorted))
    sactvch = sorted([i for i in sdb if len(sdb[i]['good_sp']) != 0])   # 1-based

    t_abs_min = min([sdb[ch]['b_sp_t_abs'][0] for ch in sdb])
    t_abs_max = max([sdb[ch]['b_sp_t_abs'][-1] for ch in sdb])
    t_abs_dur = t_abs_max - t_abs_min   # duration in us

    #nw = int(np.ceil(float(t_abs_dur) / tw))   # number of windows
    nw = range(t_abs_min,t_abs_max,10000)
    nc = len(sdb)                  # number of channels
    
    frate = np.zeros((nc, len(nw)))         # frate[0-base channel, time window index] = firing rate
    lstm = np.zeros(len(nw)).astype('int')  # lstim[time window index] = 0-based image id (e.g., Nat300_54 => 54)

    if stminfo != None:
        lst_tabs = np.array([s[0] for s in stminfo])
        lst_iid = [s[1] for s in stminfo]

    for iw in range(len(nw)):
        if iw % 10 == 0: print '* At:', iw, '/', nw
        #tb = t_abs_min + tw * iw         # begin of the window
        #te = t_abs_min + tw * (iw + 1)   # end of the window
        tb = nw[iw]
        te = tb + tw
        # -- get the image id for the window ----------------------
        if stminfo != None:
            tm = (tb + te) / 2.
            dt = lst_tabs - tb
            ic = np.argsort(np.abs(dt))[0]   # closest stimuls index
            dtc = dt[ic]                     # time difference to the closest stim
            ton = lst_tabs[ic]               # onset of the stim
            toff = ton + stimondur           # offset of the stim
            if (dtc < 0 and tm < toff) or \
               (dtc > 0 and ton < tm):
                    iid = int(lst_iid[ic].split('_')[-1])
            else:
                iid = iid_blank
        else:
            ch = 1   # now, just search for the first channel.
            # get the bad spikes during the time
            t_abs_bad = np.array(sdb[ch]['b_sp_t_abs'])
            ib = np.nonzero((t_abs_bad > tb) & (t_abs_bad < te))[0]
            # get the representative 0-based image id for the window
            # i.e., get the most frequency image id for the window
            if len(ib) == 0: iid = iid_blank
            else:
                iids = []
                for i0 in ib:
                    trel = sdb[ch]['b_sp_t_rel'][i0]
                    if trel < -75000 or trel > 175000: iid = iid_blank
                    else: iid = int(sdb[ch]['b_sp_im_id'][i0].split('_')[-1])
                    iids.append(iid)
                uiids = list(set(iids))
                p = [(iids.count(uiid), uiid) for uiid in uiids]
                iid = sorted(p, reverse=True)[0][1]

        # DBG: print iid
        lstm[iw] = iid

        # -- get the firing rate for the window ----------------------
        for ch in sactvch:
            ch0 = ch - 1   # 0-based channel
            # get the good spikes during the time
            t_abs_good = np.array(sdb[ch]['g_sp_t_abs'])
            ig = np.nonzero((t_abs_good > tb) & (t_abs_good < te))[0]
            cnt = len(ig)
            fr = float(cnt) / float(tw) * 1000000.   # convert to firing rate
            frate[ch0, iw] = fr

    pk.dump({'frate':frate, 'stm_index':lstm}, open(fout, 'wb'))


def main(args):
    if len(args) != 2:
        print 'make_timebinned_matrix.py <window_width in ms>'
        return

    tw = int(args[1])
    stminfo = pk.load(open('Chabo_20100818_RSVPNicole_A_001.stiminfo.pk'))
    print '* time window width (ms) = ', tw
    get_movie_mat('temp_data.pk',
                  'Chabo_20100818_RSVPNicole_A_001.movie.%dms.smooth.pk' % tw, tw=1000*tw, stminfo=stminfo)

if __name__ == '__main__':
    main(sys.argv)
