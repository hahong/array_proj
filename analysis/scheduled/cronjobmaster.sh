#!/bin/bash
test "$HOME" || HOME=/home/array/
export PROJROOT=$HOME/array/
export CRONDIR=$PROJROOT/analysis/scheduled/
export LOGDIR=$CRONDIR/log/

. $HOME/.profile 
$CRONDIR/01_getdata.sh && $CRONDIR/02_analyze.sh
