#!/usr/bin/env python

import sys
import numpy as np
from collections import defaultdict as dd


def getdup(x):
    cnt = dd(int)

    for e in x:
        cnt[e] += 1

    r = []
    for k in cnt:
        if cnt[k] > 1: r.append(k)

    return sorted(r)




def main(fn):
    L = [e.strip().split() for e in open(fn).readlines() if not e.startswith('#')]
    assert len(L) == 96

    L2 = []
    for i, e in enumerate(L):
        assert len(e) == 5, i
        L2.append((int(e[0]), int(e[1]), e[2], int(e[3]), e[4]))

    assert all([e[2] in ['A', 'B', 'C'] for e in L])

    loc = [(e[0], e[1]) for e in L2]
    arr = [(e[2], e[3]) for e in L2]
    lbl = [e[4] for e in L2]

    assert len(loc) == len(set(loc)), getdup(loc)
    assert len(arr) == len(set(arr)), getdup(arr)
    assert len(lbl) == len(set(lbl)), getdup(lbl)

    diff0 = np.diff([e[0] for e in L2])
    diff1 = np.diff([e[1] for e in L2])

    assert diff0.min() == -1 or -2
    assert len(diff0) == np.sum(diff0 == -1) + np.sum(diff0 == -2) + np.sum(diff0 > 5)
    assert len(diff1) == np.sum(diff1 == -1) + np.sum(diff1 == 0) #+ np.sum(diff0 > 5)

    print 'Done.'


if __name__ == '__main__':
    main(sys.argv[1])
