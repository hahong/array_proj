#!/usr/bin/env python

import sys
import cPickle as pk
import numpy as np
from joblib import Parallel, delayed

BLANKS = 'blanks'
RMAX = 0.1     # calculate d'-CV-max10%
NITER = 10     # from on 10 iterations
RSPLIT = 0.5   # with 0.5/0.5 splitting
RSEED = 0
ABS = True     # compare absolute values

# ---------------------------------------------------------------------------------------
def dprime(m_y, s_y, m_n, s_n):
    return (m_y - m_n)/np.sqrt(((s_y * s_y) + (s_n * s_n)) / 2.)

def dp(img, blank):
    m_y = np.mean(img)
    s_y = np.std(img, ddof=1)
    m_n = np.mean(blank)
    s_n = np.std(blank, ddof=1)
    return dprime(m_y, s_y, m_n, s_n)



def dp_cv(db, rsplit=RSPLIT, rmax=RMAX, niter=NITER, babs=ABS, verbose=1):
    dat = db['dat']
    iid2num = db['iid2num']
    num2iid = db['num2iid']
    nblank = iid2num[BLANKS]
    res = {}   # res[ch] = {sig_iids:sig_iids, m_dp:m_cv_rmax_dp, s_dp:s_cv_rmax_dp}

    for ch in dat:
        blanks = dat[ch][nblank]
        lb = len(blanks)
        max_dps = []

        for _ in range(niter):
            np.random.shuffle(blanks)
            blank1 = blanks[:int(lb*rsplit)]
            blank2 = blanks[int(lb*rsplit):]
            dps = []

            for n in dat[ch]:
                imgs = dat[ch][n]
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
            max_dps.extend(sdps[:int(ls*rmax)])

        cv_rmax_dp = [d[1] for d in max_dps]
        sig_iids = list(set([num2iid[d[2]] for d in max_dps]))
        m_cv_rmax_dp = np.mean(cv_rmax_dp)
        s_cv_rmax_dp = np.std(cv_rmax_dp, ddof=1)
        
        # save
        res[ch] = {'sig_iids':sig_iids, 'm_dp':m_cv_rmax_dp, 's_dp':s_cv_rmax_dp}

        if verbose:
            print '  Ch %d: m_cv_rmax_dp = %1.4f   (std = %1.4f)' % (ch, m_cv_rmax_dp, s_cv_rmax_dp)

    return res


# ---------------------------------------------------------------------------------------
def main(rseed=RSEED):
    if len(sys.argv) < 2:
        print 'summarizedp_cv.py <combined frate 1.pk> [combined frate 2.pk] ...'
        print 'NOTE: the pickle files are made by combine_firrate_data.py'
        return

    np.random.seed(rseed)

    fis = sys.argv[1:]
    db = pk.load(open(fis[0]))
    dp_cv(db)


if __name__ == '__main__':
    main()
