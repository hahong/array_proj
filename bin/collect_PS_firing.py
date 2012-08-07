#!/usr/bin/env python

import numpy as np
import cPickle as pk
import struct
import sys
sys.path.append('lib')
from mworks.data import MWKFile
from mergeutil import Merge
#from common_fn import get_stim_info, xget_events, xget_events_readahead, C_MSG, C_STIM, ERR_UTIME_MSG, ERR_UTIME_TYPE
from common_fn import *

DEFAULT_ELECS = range(1, 97)
PROC_CLUSTER = False
MAX_CLUS = 5                        # number of clusters per channel
REJECT_SLOPPY = False               # by default, do not reject sloppy (time to present > 2 frames) stimuli
EXCLUDE_IMG = None                  # exclude image by its name
DIFF_CUTOFF = 3000000               # if the animal is not working longer than 3s, then don't do readahead
MAX_READAHEAD = 6000000             # maximum readahead: 6s

CH_SHIFT = {}                       # shifting channels based on rules: CH_SHIFT[rule_name] = {src_1_based_ch:new_1_based_ch}
CH_SHIFT[None] = None
# for 1-to-1 cards
CH_SHIFT['1to1'] = {}
for ch1 in xrange(1, 49): CH_SHIFT['1to1'][ch1] = ch1 
for ch1 in xrange(81, 129): CH_SHIFT['1to1'][ch1] = ch1 - 32

# for 20110720A: assign all 40 A channels to 1-40 and all 70 M channels to 41-110
CH_SHIFT['20110720A'] = {1: 41, 2: 42, 3: 43, 4: 44, 5: 45, 6: 46, 7: 47, 8: 48, 9: 49, 10: 50, \
        11: 51, 12: 52, 13: 53, 14: 54, 15: 55, 16: 56, 17: 57, 18: 58, 19: 59, 20: 60, 21: 61, \
        22: 62, 23: 63, 24: 64, 25: 65, 26: 66, 27: 67, 28: 68, 29: 69, 30: 70, 31: 71, 32: 72, \
        33: 73, 34: 74, 35: 75, 44: 1, 45: 2, 46: 3, 47: 4, 48: 5, 49: 6, 50: 7, 51: 8, 52: 9, \
        53: 10, 54: 11, 55: 12, 56: 13, 57: 14, 58: 15, 59: 16, 60: 17, 61: 18, 62: 19, 63: 20, \
        64: 21, 65: 22, 66: 23, 67: 24, 68: 25, 69: 26, 70: 27, 71: 28, 72: 29, 73: 30, 74: 31, \
        75: 32, 76: 33, 77: 34, 78: 35, 79: 36, 80: 37, 81: 38, 82: 39, 83: 40, 94: 76, 95: 77, \
        96: 78, 97: 79, 98: 80, 99: 81, 100: 82, 101: 83, 102: 84, 103: 85, 104: 86, 105: 87, \
        106: 88, 107: 89, 108: 90, 109: 91, 110: 92, 111: 93, 112: 94, 113: 95, 114: 96, 115: 97, \
        116: 98, 117: 99, 118: 100, 119: 101, 120: 102, 121: 103, 122: 104, 123: 105, 124: 106, \
        125: 107, 126: 108, 127: 109, 128: 110}

# ----------------------------------------------------------------------------
def firrate(fn_mwk, fn_out, override_delay_us=None, override_elecs=None, verbose=2, \
        extinfo=False, c_success=C_SUCCESS, t_success_lim=T_SUCCESS, proc_cluster=PROC_CLUSTER, max_clus=MAX_CLUS, \
        t_start0=T_START, t_stop0=T_STOP, c_msg=C_MSG, c_stim=C_STIM, exclude_img=EXCLUDE_IMG, \
        reject_sloppy=REJECT_SLOPPY, err_utime_msg=ERR_UTIME_MSG, err_utime_type=ERR_UTIME_TYPE, \
        movie_begin_fname=None, ign_unregistered=False, ch_shift=None):
    """TODO: merge with get_spk() in common_fn.py"""

    mf = MWKFile(fn_mwk)
    mf.open()

    # read TOC info from the "merged" mwk file
    toc = xget_events(mf, codes=[Merge.C_MAGIC])[0].value
    c_spikes = toc[Merge.K_SPIKE]                  # get the code name from the toc

    # when the visual stimuli presented is valid?
    t_success = [ev.time for ev in xget_events(mf, codes=[c_success])]
    t_success = np.array(t_success)

    # get the active electrodes
    if override_elecs is None: actvelecs = toc[Merge.K_SPIKEWAV].keys() 
    else: actvelecs = override_elecs               # e.g, range(1, 97)
    n_actvelec = len(actvelecs)                    # number of active spike electrodes

    img_onset, img_id = get_stim_info(mf, extinfo=extinfo, exclude_img=exclude_img)

    # if requested, remove all sloppy (time spent during update main window > 2 frames)
    if reject_sloppy:
        # since get_stim_info ignores fixation point, all stimuli info must be retrived.
        all_stims = xget_events(mf, codes=[c_stim])
        all_times = np.array([s.time for s in all_stims])
        msgs = xget_events(mf, codes=[c_msg])
        errs = [m for m in msgs if m.value['type'] == err_utime_type and \
                err_utime_msg in m.value['message']]

        for e in errs:
            t0 = e.time
            rel_t = all_times - t0
            # index to the closest prior stimulus
            ci = int(np.argsort(rel_t[rel_t < 0])[-1])
            # ...and its presented MWK time
            tc = all_stims[ci].time
            # get all affected sloppy stimuli
            ss = list(np.nonzero(np.array(img_onset) == tc)[0])

            new_img_onset = []
            new_img_id = []

            # I know this is kinda O(n^2), but since ss is short, it's essentially O(n)
            for i, (io, ii) in enumerate(zip(img_onset, img_id)):
                if i in ss:
                    if verbose > 1:
                        print '** Removing sloppy:', img_id[i]
                    continue      # if i is sloppy stimuli, remove it.
                new_img_onset.append(io)
                new_img_id.append(ii)

            # trimmed the bad guys..
            img_onset = new_img_onset
            img_id = new_img_id
        assert len(img_onset) == len(img_id)

    n_stim = len(img_onset)

    # MAC-NSP time translation
    if override_delay_us != None: 
        t_delay = toc['align_info']['delay']
        t_adjust = int(np.round(override_delay_us - t_delay))
    else: 
        t_adjust = 0
    t_start = t_start0 - t_adjust
    t_stop = t_stop0 - t_adjust

    # actual calculation -------------------------------
    # all_spike[chn_id][img_id]: when the neurons spiked?
    all_spike = {}
    all_foffset = {}
    clus_info = {}

    frame_onset = {}
    movie_iid = None
    movie_onsets = []
    movie_onset0 = 0


    t0_valid = []
    iid_valid = []
    for i in xrange(n_stim):
        t0 = img_onset[i]; iid = img_id[i]
        # -- check if this presentation is successful. if it's not ignore this.
        if np.sum((t_success > t0) & (t_success < (t0 + t_success_lim))) < 1: continue
        t0_valid.append(t0)
        iid_valid.append(iid)
    n_stim_valid = len(t0_valid)

    t_slack = t_stop - t_start
    readaheads = np.zeros(n_stim_valid, 'int')
    i_cnkbegin = 0            # beginning of the chunk
    for i in xrange(1, n_stim_valid):
        t0 = t0_valid[i]
        t0p = t0_valid[i - 1]
        t0b = t0_valid[i_cnkbegin]

        if (t0 - t0p > DIFF_CUTOFF) or (t0 - t0b > MAX_READAHEAD):
            readaheads[i_cnkbegin:i] = t0p - t0b + t_slack
            i_cnkbegin = i
            continue
    readaheads[i_cnkbegin:] = t0 - t0b + t_slack

    for i in xrange(n_stim_valid):
        t0 = t0_valid[i]; iid = iid_valid[i]; readahead=int(readaheads[i])
        # -- process movie?
        if movie_begin_fname != None:
            # begin new clip?
            if movie_begin_fname in iid:
                # was there previous clip?
                if movie_iid != None:
                    if movie_iid not in frame_onset: frame_onset[movie_iid] = []
                    frame_onset[movie_iid].append(movie_onsets)
                # init for new clip
                movie_onsets = []
                iid = movie_iid = iid.replace(movie_begin_fname, '')
                movie_onset0 = t0
                movie_onsets.append(0)
            elif movie_iid != None:
                movie_onsets.append(t0 - movie_onset0)
                continue

        if verbose > 0: 
            print 'At', (i + 1), 'out of', n_stim_valid, '         \r',
            sys.stdout.flush()
   
        #spikes = xget_events(mf, codes=[c_spikes], time_range=[t0 + t_start, t0 + t_stop])
        spikes = xget_events_readahead(mf, c_spikes, (t0 + t_start, t0 + t_stop), readahead=readahead)
        actvunits = {}
        t_rel = {}
        foffset = {}
        # -- prepare the t_rel & foffset
        for ch in actvelecs:
            # if no clustering info is used...
            if not proc_cluster:
                t_rel[ch] = []
                foffset[ch] = []
                continue
            # if clustering info is used...
            cids = range(max_clus)
            actvunits[ch] = cids
            for cid in cids:
                t_rel[(ch,cid)] = []
                foffset[(ch,cid)] = []

        # -- put actual spiking info
        for s in spikes:
            ch = s.value['id']
            if ch_shift != None:   # if mapping is requested
                if ch not in ch_shift: continue
                ch = ch_shift[ch]

            if proc_cluster: 
                cid = s.value['cluster_id']
                key = (ch, cid)
            else:
                key = ch
          
            # put the relative time
            if ign_unregistered and key not in t_rel: continue
            t_rel[key].append(int(s.time + t_adjust - t0))
            foffset[key].append(int(s.value['foffset']))
            # update the clus_info and n_cluster
            if proc_cluster and key not in clus_info:
                clus_info[key] = s.value
                if s.value['nclusters'] > max_clus:
                    raise ValueError, '*** Unit %s: max_cluster(=%d) is smaller than the actual number of clusters(=%d)!' \
                            % (str(key), max_clus, s.value['nclusters'])

        # -- combine all
        for el in actvelecs:
            if proc_cluster: cids = actvunits[el]
            else: cids = [-1]
            for cid in cids:
                if proc_cluster: key = (el, cid)
                else: key = el
                if key not in all_spike: 
                    # not using defaultdict here:
                    # all_spike[key] = defaultdict(list)
                    all_spike[key] = {}
                    all_foffset[key] = {}
                if iid not in all_spike[key]:
                    all_spike[key][iid] = []
                    all_foffset[key][iid] = []
                all_spike[key][iid].append(t_rel[key])
                all_foffset[key][iid].append(foffset[key])

    # flush movie data
    if movie_iid != None:
        if movie_iid not in frame_onset: frame_onset[movie_iid] = []
        frame_onset[movie_iid].append(movie_onsets)

    # finished calculation....
    f = open(fn_out, 'w')
    out =  {'all_spike': all_spike, 
            't_start': t_start0,
            't_stop': t_stop0,
            't_adjust': t_adjust,
            'actvelecs': actvelecs}
    if proc_cluster:
        out['clus_info'] = clus_info
        out['max_clus'] = max_clus
    if movie_begin_fname != None:
        out['frame_onset'] = frame_onset
    pk.dump(out, f)

    # put all_foffset into the 2nd half to speed up reading
    out2 = {'all_foffset': all_foffset,}
    pk.dump(out2, f)

    f.close()


# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 3:
        print 'collect_PS_firing.py <mwk> <output.pk> [override delay in us] [number of electrodes] [opts]'
        print 'Collects spike timing around visual stimuli'
        print 
        print 'Options:'
        print 'extinfo              - collects extra stimuli information in addition to the names'
        print 'c_success=string     - code name for "success" signal'
        print 'proc_cluster         - process extra spike sorting information'
        print 'max_cluster=#        - maximum number of clusters per channel'
        return

    fn_mwk = sys.argv[1]
    fn_out = sys.argv[2]
    opts0 = sys.argv[5:]
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
    if len(sys.argv) >= 4:
        override_delay_us = long(sys.argv[3])
        if override_delay_us < 0: override_delay_us = None
        print '* Delay override:', override_delay_us
    else: override_delay_us = None

    if len(sys.argv) >= 5:
        override_elecs = range(1, int(sys.argv[4]) + 1)
        print '* Active electrodes override: [%d..%d]' % (min(override_elecs), max(override_elecs))
    else: override_elecs = DEFAULT_ELECS

    # handle "opts"
    if 'extinfo' in opts:
        extinfo = True
        print '* Collecting extra information of the stimuli'

    if 'c_success' in opts:
        c_success = opts['c_success']
        print '* c_success:', c_success
    else:
        c_success = C_SUCCESS

    if 't_success' in opts:
        t_success = int(opts['t_success'])
        print '* t_success:', t_success
    else:
        t_success = T_SUCCESS

    proc_cluster = PROC_CLUSTER
    if 'proc_cluster' in opts or 'proc_spksorting' in opts:
        print '* Collecting spike sorting information'
        proc_cluster = True

    max_clus = MAX_CLUS
    if 'max_cluster' in opts:
        max_clus = int(opts['max_cluster'])
        print '* Maximum number of clusters per channel:', max_clus

    t_start0 = T_START
    if 't_start' in opts:
        t_start0 = int(opts['t_start'])
        print '* t_start =', t_start0

    t_stop0 = T_STOP
    if 't_stop' in opts:
        t_stop0 = int(opts['t_stop'])
        print '* t_stop =', t_stop0

    reject_sloppy = REJECT_SLOPPY
    if 'reject_sloppy' in opts:
        reject_sloppy = True
        print '* Rejecting sloppy stimuli'

    exclude_img = EXCLUDE_IMG
    if 'exclude_img' in opts:
        exclude_img = opts['exclude_img'].split(',')
        print '* Exclude unwanted images:', exclude_img

    movie_begin_fname = None
    if 'movie_begin_fname' in opts:
        movie_begin_fname = opts['movie_begin_fname']
        print '* movie_begin_fname:', movie_begin_fname

    ign_unregistered = False
    if 'ign_unregistered' in opts:
        ign_unregistered = True
        print '* Ignore unregistered keys'

    ch_shift = None
    if 'ch_shift' in opts:
        ch_shift = opts['ch_shift']
        print '* Shifting based on this rule:', ch_shift

    # go go go
    firrate(fn_mwk, fn_out, override_delay_us=override_delay_us, override_elecs=override_elecs, \
            extinfo=extinfo, c_success=c_success, t_success_lim=t_success, proc_cluster=proc_cluster, max_clus=max_clus, \
            t_start0=t_start0, t_stop0=t_stop0, reject_sloppy=reject_sloppy, exclude_img=exclude_img, movie_begin_fname=movie_begin_fname, \
            ign_unregistered=ign_unregistered, ch_shift=CH_SHIFT[ch_shift]) 
    print 'Done.                                '


if __name__ == '__main__':
    main()
