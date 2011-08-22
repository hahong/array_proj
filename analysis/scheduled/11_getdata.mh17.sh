#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/mindhive/dicarlolab/proj/array_data/
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
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOG/01_getdata.sh.lock" && exit   
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOG/02_analyze.sh.lock" && exit   

# -- Check lock
if [ -f $LOCK ]; then
        # -- if locked, terminates immediately
        echo "Locked:" $LOCK
        exit
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
	rmtdir=$1
	lcldir=$2

	# syncbad: local -> remote
	$LCLGET $lcldir > $LCLBADLST
	scp $LCLBADLST $REMOTEUSER@$REMOTEFILER:$RMTBADLST
	ssh $REMOTEUSER@$REMOTEFILER "$RMTSET $RMTBADLST $rmtdir; echo $RMTGET $rmtdir " #> $RMTBADLST" 
	scp $REMOTEUSER@$REMOTEFILER:$RMTBADLST $LCLBADLST
	$LCLSET $LCLBADLST $lcldir

	# rsync -avzuH --exclude='*.ns5' --exclude='*.ns5.*' --exclude='*cluster_wd*' $REMOTEUSER@$REMOTEFILER:$REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/ 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2.log &
}

# -- 1. Tito
syncall $REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/
# rsync -avzuH --exclude='*.ns5' --exclude='*.ns5.*' --exclude='*cluster_wd*' $REMOTEUSER@$REMOTEFILER:$REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/ 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2.log &
wait

# -- Done
rm -f $LOCK
