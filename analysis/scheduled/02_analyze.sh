#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/home/array/array/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
LOCK=$LOGDIR/02_analyze.sh.lock
JOBNAME=joblist/`date +%Y%m%d`_04_merge+collect_`date +%H%M`.sh
JOBNAMEPSTH=joblist/`date +%Y%m%d`_10_plot_PSTH_`date +%H%M`.sh
LOGNAME=$LOGDIR/`date +%Y%m%d_%H%M%S`_analysis.log
HTMLPSTH=/home/array/public_html/psth/

###################################################################
# -- Merge and collect

if [ -f $LOCK ]; then
	# -- if locked, terminates immediately
	echo "Locked:" $LOCK
	exit
fi

cd $PROJROOT/analysis/
touch $LOCK   # create a lock file so that no data transfer occurs from mh17

# -- 1. Tito
# Merge and collect
./04_par_merge+collect_PS_firing.py > $JOBNAME
dirnev=data/d003_Tito/neudat_NSP2/ ./04_par_merge+collect_PS_firing.py --merge_opts='multinsp NSP2' >> $JOBNAME
NJOBS=5 parrun.py $JOBNAME 2>&1 | tee -a $LOGNAME

# Plot PSTHs
./10_plot_PSTH.sh > $JOBNAMEPSTH
NJOBS=5 parrun.py $JOBNAMEPSTH 2>&1 | tee -a $LOGNAME 
rsync -a $PROJROOT/data/d003_Tito/subproj/100_PSTH/*.pdf $HTMLPSTH

# -- END
rm -f $LOCK
