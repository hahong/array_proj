#!/bin/bash
# This script is meant to be called by cronjobmaster.sh
test "$PROJROOT" || PROJROOT=/mindhive/dicarlolab/proj/array_data/
test "$LOGDIR" || LOGDIR=$PROJROOT/analysis/scheduled/log/
REMOTEFILER=dicarlo2
REMOTEUSER=array
REMOTELOCK=array/analysis/scheduled/log/
REMOTEDATA=array/data/
LOCK=$LOGDIR/11_getdata.mh17.sh.lock

###################################################################
# if still processing the data at dicarlo2, then quits.
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOCK/01_getdata.sh.lock" && exit   
ssh $REMOTEUSER@$REMOTEFILER "test -f $REMOTELOCK/02_analyze.sh.lock" && exit   

# -- Check lock
if [ -f $LOCK ]; then
        # -- if locked, terminates immediately
        echo "Locked:" $LOCK
        exit
fi
touch $LOCK

# -- Start fetch
# fetch Tito
rsync -avzuH --exclude='*.ns5' --exclude='*.ns5.*' --exclude='*cluster_wd*' $REMOTEUSER@$REMOTEFILER:$REMOTEDATA/d002_Tito/ $PROJROOT/data/d002_Tito/ 2>&1 | tee -a $LOGDIR/`date +%Y%m%d_%H%M%S`_Tito_dicarlo2.log &
wait

# -- Done
rm -f $LOCK
