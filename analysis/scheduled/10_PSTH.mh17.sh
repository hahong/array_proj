#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/home/array/array/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
REMOTEFILER=dicarlo2
REMOTEUSER=array
REMOTELOCK=array/analysis/scheduled/log/02_analyze.sh.lock
JOBNAME=joblist/`date +%Y%m%d`_PSTH.sh

###################################################################
# if still processing the data at dicarlo2, then quits.
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOCK" && exit   
#cd $PROJROOT/analysis/



# ./04_par_merge+collect_PS_firing.py > $JOBNAME
# parrun.py $JOBNAME 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_analysis.log
