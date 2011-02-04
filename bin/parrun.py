#!/usr/bin/env python

import os
import sys
from joblib import Parallel, delayed

def parrun(jobs):
    r = Parallel(n_jobs=-1, verbose=1)(delayed(os.system)(j) for j in jobs)

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'parrun: no joblist file given. reading from stdin.'
        f = sys.stdin
    else: f = open(sys.argv[1])

    jobs = [line.strip() for line in f.readlines()]
    parrun(jobs)
