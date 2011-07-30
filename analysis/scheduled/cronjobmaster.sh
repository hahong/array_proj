#!/bin/bash
test "$HOME" || HOME=/home/array/
PROJROOT=$HOME/array/
CRONDIR=$PROJROOT/analysis/scheduled/
LOGDIR=$CRONDIR/log/

. $HOME/.profile 
. $CRONDIR/01_getdata.sh
. $CRONDIR/02_analyze.sh
