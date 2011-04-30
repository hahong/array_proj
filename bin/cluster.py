#!/usr/bin/env python
# encoding: utf-8

import sys
import os
import glob
import shutil as sh
import collections as cl
import numpy as np
from common_fn import *
from joblib import Parallel, delayed
from scipy.spatial import KDTree

SUFF_FET = '.fet.'
SUFF_INF = '.inf.'
SUFF_CLU = '.clu.'
METD_CLUSTER = 'spc'
MAX_PTS = 10000
BIN = os.getcwd() + os.sep + 'cluster_maci.exe'
ERR_CLUSTER = -1
MIN_CLUS = 60
SEP = ','
NELEC = 128
NCPU = -1


# --------------------------------------------------------------------------------
# direct translation from "find_temp.m" of Wave_clus
def spc_findtemp(tree, min_clus=MIN_CLUS):
    num_temp = len(tree)
            
    aux  = np.diff(tree[:,4])   # Changes in the first cluster size                                       
    aux1 = np.diff(tree[:,5])   # Changes in the second cluster size
    aux2 = np.diff(tree[:,6])   # Changes in the third cluster size                                       
    aux3 = np.diff(tree[:,7])   # Changes in the fourth cluster size                                       
                
    temp = 0         # Initial value
                    
    for t in range(num_temp - 1):
        # Looks for changes in the cluster size of any cluster larger than min_clus.                  
        if aux[t] > min_clus or aux1[t] > min_clus or aux2[t] > min_clus or aux3[t] > min_clus:
            temp = t + 1
            
    # In case the second cluster is too small, then raise the temperature a little bit                 
    if temp == 0 and tree[temp, 5] < min_clus:
        temp = 1

    return temp


def spc_template_match(fn_fet, fn_base, osuff='', verbose=True):
    if verbose: print 'At:', fn_fet
    # load parameters
    try:
        fet = np.loadtxt(fn_fet, skiprows=1)
        clu = np.loadtxt(fn_base + '.out.dg_01.lab')
        tree = np.loadtxt(fn_base + '.out.dg_01')
        pts = np.loadtxt(fn_base + '.inp')
    except Exception, e:
        print '** Not clustered properly:', fn_fet
        # --
        fn_cluout = fn_fet.replace(SUFF_FET, osuff + SUFF_CLU)
        fp = open(fn_cluout, 'wt')
        print >>fp, ERR_CLUSTER
        """for pt in fet:
            dist, idx = K.query(pt)
            print >>fp, lbl[idx]"""
        fp.close()
        return

    # index of temperature to use
    iT = spc_findtemp(tree)
    lbl = clu[iT,2:].astype('int')
    K = KDTree(pts)
    #dist, idx = K.query(pt)
    # due to the pickling issue, doesn't work!!: 
    # XXX: r = Parallel(n_jobs=ncpu, verbose=0)(delayed(K.query)(pt) for pt in fet)

    # finished template matching: writing cluster info
    nclusters = int(tree[iT,3])
    fn_cluout = fn_fet.replace(SUFF_FET, osuff + SUFF_CLU)
    fp = open(fn_cluout, 'wt')
    print >>fp, nclusters
    for pt in fet:
        dist, idx = K.query(pt)
        print >>fp, lbl[idx]
    fp.close()


def spc_write_run(fn_run, fn_inp, fn_out, npts, ndim, rseed=None):
    # -- writing parameter files for SPC clustering (copied from wave_clus)
    frun = open(fn_run, 'wt')
    print >>frun, 'NumberOfPoints:', npts
    print >>frun, 'DataFile:',       fn_inp
    print >>frun, 'OutFile:',        fn_out
    print >>frun, 'Dimensions:',     ndim 
    # the followings are default constants. see Wave_clus
    print >>frun, 'MinTemp: 0'
    print >>frun, 'MaxTemp: 0.21'
    print >>frun, 'TempStep: 0.0025'
    print >>frun, 'SWCycles: 100'
    print >>frun, 'KNearestNeighbours: 11'
    print >>frun, 'MSTree|'
    print >>frun, 'DirectedGrowth|'
    print >>frun, 'SaveSuscept|'
    print >>frun, 'WriteLables|'
    print >>frun, 'WriteCorFile~'
    if rseed != None: 
        print >>frun, 'ForceRandomSeed:', rseed
    frun.close()


def spc_cluster(fn_prefix, rseed=None, nmax=MAX_PTS, ncpu=NCPU, bin=BIN, rm_tmp_level=1, osuff=''):
    wd = os.getcwd()
    pid = os.getpid()
    uniq = 0

    bname0 = os.path.basename(fn_prefix[0])
    basedir0 = os.path.dirname(fn_prefix[0])
    if basedir0 == '': basedir0 = '.'
    os.chdir(basedir0)

    fn_fet0s = glob.glob(bname0 + SUFF_FET + '*')
    job_spc = []       # list for the SPC runs
    job_match = []     # list for the template matching
    job_rm = []        # list for the temp file removal

    # get unique `uniq` to prevent pid collision (in a cluster)
    while os.path.exists('tmp.%d.%05d.000.inp' % (uniq, pid)):
        uniq += 1

    # -- make the input file and joblist
    for ifn, fn_fet0 in enumerate(fn_fet0s):
        # -- write SPC compatible input files
        tmp_base = 'tmp.%d.%05d.%03d' % (uniq, pid, ifn)
        fn_inp = tmp_base + '.inp'
        fn_run = tmp_base + '.run'
        fn_out = tmp_base + '.out'

        finp = open(fn_inp, 'wt')
        npts = 0    # total number of points to cluster
        # iterate over the prefixes
        for fpx in fn_prefix:
            bname = os.path.basename(fpx)
            fn_fet = fn_fet0.replace(bname0, bname)
            src = open(fn_fet).readlines()

            ndim = int(src[0].strip())
            src.pop(0)                  # remove the header
            np.random.shuffle(src)
            out = src[:nmax] 
            npts += len(out)
            finp.writelines(out)

            job_match.append((fn_fet, tmp_base))
            del src, out
        finp.close()

        spc_write_run(fn_run, fn_inp, fn_out, npts, ndim, rseed)
        cmd = bin + ' ' + fn_run
        job_spc.append(cmd)
        job_rm.append((tmp_base, fn_fet0))

    # -- clustering/template matching. (n.b.: still within `basedir`)
    # run the SPC clustering parallely
    r = Parallel(n_jobs=ncpu, verbose=1)(delayed(os.system)(cmd) for cmd in job_spc)
    # template matching with NN
    r = Parallel(n_jobs=ncpu, verbose=1)(delayed(spc_template_match)(fn_fet, tmp_base, osuff=osuff) for fn_fet, tmp_base in job_match)
        
    # -- remove temporary files
    if rm_tmp_level > 0:
        for tmp_base, fn_fet0 in job_rm:
            if rm_tmp_level == 1: 
                for ext in ['.inp', '.out.dg_01', '.out.dg_01.lab', '.run']:
                    try:
                        sh.move(tmp_base + ext, fn_fet0.replace(SUFF_FET, osuff + '.spc.') + ext)
                    except OSError, e: pass
            for ext in ['.out.mag', '.out.mst11.edges', '.out.param', \
                    '.inp', '.out.dg_01', '.out.dg_01.lab', '.run']:
                try:
                    os.unlink(tmp_base + ext)
                except OSError, e: pass


# --------------------------------------------------------------------------------
def main(fn_prefix, metd, bin, osuffix, nsample):
    if metd == 'spc': spc_cluster(fn_prefix, bin=bin, osuff=osuffix, nmax=nsample)
    else:
        print 'cluster.py: invalid metd'

# --------------------------------------------------------------------------------
def prep_files(flist, nelec=NELEC):
    flist = flist.split(SEP)
    if flist[0][0] == '+':
        flist = [f.strip() for f in open(flist[0][1:]).readlines()]
    assert all([os.path.exists(f + SUFF_FET + str(ch)) for f in flist for ch in range(1, nelec+1)])
    return flist

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'cluster.py <.fet prefixes separated by ","> [options 1] [options 2] [...]'
        print 'Do unsupervised clustering on wave features.  Compound prefixes must be'
        print 'in the same directory.'
        print
        print 'Options:'
        print '   metd        - clustering method (e.g, spc)'
        print '   bin         - clustering binary (must be abs path)'
        print '   nsample     - number of sampled points/prefix'
        print '   osuffix     - suffix of the output files'
        print '   nelec       - number of electrodes'
        sys.exit(1)

    opts = parse_opts(sys.argv[2:])
    fn_prefix = sys.argv[1]

    print '* Processing:', fn_prefix
    
    metd = METD_CLUSTER
    if 'metd' in opts: metd = opts['metd']

    bin = BIN
    if 'bin' in opts: bin = opts['bin']

    nelec = NELEC
    if 'nelec' in opts: nelec = int(opts['nelec'])

    ncpu = NCPU
    if 'ncpu' in opts: ncpu = int(opts['ncpu'])

    osuffix = ''
    if 'osuffix' in opts: osuffix = opts['osuffix']

    nsample = MAX_PTS
    if 'nsample' in opts: nsample = int(opts['nsample'])

    fn_prefix = prep_files(fn_prefix, nelec)
    assert os.path.exists(bin)

    print '* Variables: (metd, bin, osuffix, nsample, nelec, ncpu) =', (metd, bin, osuffix, nsample, nelec, ncpu)

    main(fn_prefix, metd, bin, osuffix, nsample)

