#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/home/array/array/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
REMOTE_LOCK=array/analysis/scheduled/log/02_analyze.sh.lock
JOBNAME=joblist/`date +%Y%m%d`_PSTH.mh17.sh

###################################################################
# -- Merge and collect
cd $PROJROOT/analysis/
# ./04_par_merge+collect_PS_firing.py > $JOBNAME
# parrun.py $JOBNAME 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_analysis.log
