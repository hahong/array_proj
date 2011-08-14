#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/home/array/array/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
LOCK=$LOGDIR/02_analyze.sh.lock
JOBNAME=joblist/`date +%Y%m%d`_04_merge+collect.sh

###################################################################
# -- Merge and collect

if [ -f $LOCK ]; then
	# -- if locked, terminates immediately
	echo "Locked:" $LOCK
	exit
fi

cd $PROJROOT/analysis/
touch $LOCK   # create a lock file so that no data transfer occurs from mh17
./04_par_merge+collect_PS_firing.py > $JOBNAME
parrun.py $JOBNAME 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_analysis.log
rm -f $LOCK
