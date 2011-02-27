#!/usr/bin/env python

import sys
import numpy as np

ELECS = range(1, 97)


def load_impedance(fi, elecs=ELECS):
    im = np.zeros(len(elecs))    # XXX: kludge...

    for line in open(fi, 'rt').readlines():
        tokens = line.strip().split()
        if len(tokens) != 3 or tokens[2] != 'kOhm': continue
        
        ch = int(tokens[0].split('elec')[-1])
        if ch not in elecs: continue

        im[ch - 1] = int(tokens[1])
        
    return im


def cmp_impedance(fi1, fi2):
    im1 = load_impedance(fi1)
    im2 = load_impedance(fi2)

    denom = np.max([im1, im2], axis=0)
    nom = np.abs(im1 - im2)
    ratio = nom / denom

    return ratio


def cmp_impedance_all(files):
    fi1s = files[0::2]
    fi2s = files[1::2]
    res = []
    
    for fi1, fi2 in zip(fi1s, fi2s):
        ratio = cmp_impedance(fi1, fi2)
        res.append(sorted(ratio))
    n_ch = len(res[0])

    for ch in range(n_ch):
        row = []
        for col in range(len(res)):
            row.append(res[col][ch])

        fmt = '%3d' + '\t%1.3f' * len(row)
        print fmt % ((ch + 1,) + tuple(row))


def main():
    if len(sys.argv) < 3:
        print 'cmp_impedance.py <impedance log 1a.txt> <impedance log 1b.txt> [<log 2a.txt> <log 2b.txt>] ...'
        print 'cmp_impedance.py: compare two impedance log files.'
        return

    files = sys.argv[1:]
    cmp_impedance_all(files)


if __name__ == '__main__':
    main()
