#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/mindhive/dicarlolab/proj/array_data/array/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
LCLBIN=$PROJROOT/analysis/scheduled/
REMOTEFILER=dicarlo2
REMOTEUSER=array
REMOTEBIN=array/analysis/scheduled/
REMOTELOG=array/analysis/scheduled/log/
REMOTEDATA=array/data/
LOCK=$LOGDIR/11_getdata.mh17.sh.lock

###################################################################
# if still processing the data at dicarlo2, then quits.
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOG/01_getdata.sh.lock" && exit 1  
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOG/02_analyze.sh.lock" && exit 1  

# -- Check lock
if [ -f $LOCK ]; then
        # -- if locked, terminates immediately
        echo "Locked:" $LOCK
        exit 1
fi
touch $LOCK

##################################################################
# Start fetch
LCLGET=$LCLBIN/12_getbadlist.mh17.sh
LCLSET=$LCLBIN/13_setbadlist.mh17.sh
RMTGET=$REMOTEBIN/12_getbadlist.mh17.sh
RMTSET=$REMOTEBIN/13_setbadlist.mh17.sh

LCLBADLST=$LOGDIR/12_getbadlist.mh17.sh.badlist
RMTBADLST=$REMOTELOG/12_getbadlist.mh17.sh.badlist

function syncall {
	rmtdir=`dirname $1`/`basename $1`/
	lcldir=`dirname $2`/`basename $2`/

	# syncbad: local -> remote
	$LCLGET $lcldir > $LCLBADLST
	scp $LCLBADLST $REMOTEUSER@$REMOTEFILER:$RMTBADLST 2>&1
	ssh $REMOTEUSER@$REMOTEFILER "$RMTSET $RMTBADLST $rmtdir; $RMTGET $rmtdir > $RMTBADLST" 
	# syncbad: remote -> local
	scp $REMOTEUSER@$REMOTEFILER:$RMTBADLST $LCLBADLST 2>&1
	$LCLSET $LCLBADLST $lcldir
	# file sync
	rsync -avzuH --exclude='*.ns5' --exclude='*.ns5.*' --exclude='.*ns5*' --exclude='.*nev*' --exclude='.*ccf*' --exclude='.*mwk*' --exclude='*cluster_wd*' $REMOTEUSER@$REMOTEFILER:$rmtdir $lcldir 2>&1
}

# -- 1. Tito
syncall $REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/ | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2_d002_Tito.log
syncall $REMOTEDATA/d003_Tito/ $PROJROOT/data/d003_Tito/ | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2_d003_Tito.log
# rsync -avzuH --exclude='*.ns5' --exclude='*.ns5.*' --exclude='*cluster_wd*' $REMOTEUSER@$REMOTEFILER:$REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/ 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2.log &
# wait
 
# -- Done
rm -f $LOCK
