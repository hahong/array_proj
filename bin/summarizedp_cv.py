#!/usr/bin/env python

import sys
import cPickle as pk
import numpy as np
import time
from joblib import Parallel, delayed

BLANKS = 'blanks'
RMAX = 0.1     # calculate d'-CV-max10%
NITER = 16     # from on 16 iterations
RSPLIT = 0.5   # with 0.5/0.5 splitting
RSEED = 0
ABS = True     # compare absolute values
NJOBS = -1     # maximum worker threads

# ---------------------------------------------------------------------------------------
def dprime(m_y, s_y, m_n, s_n):
    return (m_y - m_n)/np.sqrt(((s_y * s_y) + (s_n * s_n)) / 2.)

def dp(img, blank):
    m_y = np.mean(img)
    s_y = np.std(img, ddof=1)
    m_n = np.mean(blank)
    s_n = np.std(blank, ddof=1)
    return dprime(m_y, s_y, m_n, s_n)

def dp_thread(dat_ch, nblank, babs=ABS, rsplit=RSPLIT, rmax=RMAX):
    blanks = list(dat_ch[nblank])
    np.random.shuffle(blanks)
    lb = len(blanks)
    blank1 = blanks[:int(lb*rsplit)]
    blank2 = blanks[int(lb*rsplit):]
    dps = []

    for n in dat_ch:
        imgs = list(dat_ch[n])
        np.random.shuffle(imgs)
        li = len(imgs)
        img1 = imgs[:int(li*rsplit)]
        img2 = imgs[int(li*rsplit):]

        # compute cross-validated d's
        dp1 = dp(img1, blank1)
        dp2 = dp(img2, blank2)

        # take absolute if needed
        if babs:
            dp1 = abs(dp1)
            dp2 = abs(dp2)

        dps.append((dp1, dp2, n)) 

    # take rmax portion only
    sdps = sorted(dps, reverse=True)
    ls = len(sdps)
    return sdps[:int(ls*rmax)]


def dp_cv(db, niter=NITER, verbose=1):
    dat = db['dat']
    iid2num = db['iid2num']
    num2iid = db['num2iid']
    nblank = iid2num[BLANKS]
    res = {}   # res[ch] = {sig_iids:sig_iids, m_dp:m_cv_rmax_dp, s_dp:s_cv_rmax_dp}

    for ch in dat:
        t0 = time.time()
        max_dps = Parallel(n_jobs=NJOBS, verbose=0)(delayed(dp_thread)(dat[ch], nblank) for _ in range(niter))

        cv_rmax_dp = [d[1] for sublst in max_dps for d in sublst]
        sig_iids = list(set([num2iid[d[2]] for sublst in max_dps for d in sublst]))
        m_cv_rmax_dp = np.mean(cv_rmax_dp)
        s_cv_rmax_dp = np.std(cv_rmax_dp, ddof=1)
        
        # save
        res[ch] = {'sig_iids':sig_iids, 'm_dp':m_cv_rmax_dp, 's_dp':s_cv_rmax_dp, 'blanks':dat[ch][nblank]}
        t1 = time.time()

        if verbose:
            print '  Ch %d: m_cv_rmax_dp = %1.4f   (std = %1.4f; dt=%f)' % (ch, m_cv_rmax_dp, s_cv_rmax_dp, t1-t0)

    return res


# ---------------------------------------------------------------------------------------
def main(rseed=RSEED):
    if len(sys.argv) < 2:
        print 'summarizedp_cv.py prep <d\'-CV output.pk> <combined frate.pk>'
        print 'NOTE: the pickle files are made by combine_firrate_data.py'
        return

    np.random.seed(rseed)    
    mode = sys.argv[1]

    # prepare d'-cv computation
    if mode == 'prep':
        fo = sys.argv[2]
        fi = sys.argv[3]
        print '* Output d\'-CV result:', fo

        db = pk.load(open(fi))
        res = dp_cv(db)
        pk.dump(res, open(fo, 'wb'))
    else:
        print 'Unrecognized mode. Aborting.'


if __name__ == '__main__':
    main()
