#!/usr/bin/env python

import numpy as np
import cPickle as pk
import struct
import sys
sys.path.append('lib')
from mworks.data import MWKFile
from mergeutil import Merge
from common_fn import get_stim_info, xget_events

C_STIM = '#announceStimulus'
I_STIM_ID = 2
T_START = -100000
T_STOP = 250000
DEFAULT_ELECS = range(1, 97)
C_SUCCESS = 'number_of_stm_shown'   # one visual stimulus should have one "success" event within
T_SUCCESS = 250000                  # 250ms time window in order to be considered as valid one.
PROC_CLUSTER = False
MAX_CLUS = 5                          # number of clusters per channel

# ----------------------------------------------------------------------------
def firrate(fn_mwk, fn_out, override_delay_us=None, override_elecs=None, verbose=True, \
        extinfo=False, c_success=C_SUCCESS, proc_cluster=PROC_CLUSTER, max_clus=MAX_CLUS):
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

    img_onset, img_id = get_stim_info(mf, extinfo=extinfo)
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
    all_spike = {}
    clus_info = {}
    for i in range(n_stim):
        t0 = img_onset[i]; iid = img_id[i]
        # -- check if this presentation is successful. if it's not ignore this.
        if np.sum((t_success > t0) & (t_success < (t0 + T_SUCCESS))) < 1: continue

        if verbose: 
            print 'At', (i + 1), 'out of', n_stim, '         \r',
            sys.stdout.flush()
   
        spikes = xget_events(mf, codes=[c_spikes], time_range=[t0 + t_start, t0 + t_stop])
        actvunits = {}
        t_rel = {}
        # -- prepare the t_rel
        for ch in actvelecs:
            # if no clustering info is used...
            if not proc_cluster:
                t_rel[ch] = []
                continue
            # if clustering info is used...
            cids = range(max_clus)
            actvunits[ch] = cids
            for cid in cids:
                t_rel[(ch,cid)] = []

        # -- put actual spiking info
        for s in spikes:
            if proc_cluster: 
                ch = s.value['id']
                cid = s.value['cluster_id']
                key = (ch, cid)
            else:
                key = s.value['id']
          
            # put the relative time
            t_rel[key].append(int(s.time + t_adjust - t0))
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
                if iid not in all_spike[key]:
                    all_spike[key][iid] = []
                all_spike[key][iid].append(t_rel[key])

    # finished calculation....
    f = open(fn_out, 'w')
    out =  {'all_spike': all_spike, 
            't_start': T_START,
            't_stop': T_STOP,
            't_adjust': t_adjust,
            'actvelecs': actvelecs}
    if proc_cluster: out['clus_info'] = clus_info
    pk.dump(out, f)
    f.close()



# ----------------------------------------------------------------------------
def main():
    if len(sys.argv) < 3:
        print 'collect_PS_firing.py <mwk> <output.pk> [override delay in us] [number of electrodes] [opts]'
        print 'Collects spike timing around visual stimuli'
        print 
        print 'Options:'
        print 'extinfo              - collects extra stimuli information in addition to the names'
        print 'c_success=<string>   - code name for "success" signal'
        print 'proc_cluster         - process extra spike sorting information'
        print 'max_cluster=<#>      - maximum number of clusters per channel'
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
        print 'Delay override:', override_delay_us
    else: override_delay_us = None

    if len(sys.argv) >= 5:
        override_elecs = range(1, int(sys.argv[4]) + 1)
        print 'Active electrodes override: [%d..%d]' % (min(override_elecs), max(override_elecs))
    else: override_elecs = DEFAULT_ELECS

    # handle "opts"
    if 'extinfo' in opts:
        extinfo = True
        print 'Collecting extra information of the stimuli'

    if 'c_success' in opts:
        c_success = opts['c_success']
        print 'c_success:', c_success
    else:
        c_success = C_SUCCESS

    proc_cluster = PROC_CLUSTER
    if 'proc_cluster' in opts or 'proc_spksorting' in opts:
        print 'Collecting spike sorting information'
        proc_cluster = True

    max_clus = MAX_CLUS
    if 'max_cluster' in opts:
        max_clus = int(opts['max_cluster'])
        print 'Maximum number of clusters per channel:', max_clus

    # go go go
    firrate(fn_mwk, fn_out, override_delay_us=override_delay_us, override_elecs=override_elecs, \
            extinfo=extinfo, c_success=c_success, proc_cluster=proc_cluster, max_clus=max_clus) 
    print 'Done.                                '


if __name__ == '__main__':
    main()
