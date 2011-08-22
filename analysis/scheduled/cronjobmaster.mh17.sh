#!/bin/bash
test "$HOME" || HOME=/mindhive/dicarlolab/u/hahong/
export PROJROOT=$HOME/teleport/array/
export CRONDIR=$PROJROOT/analysis/scheduled/
export LOGDIR=$CRONDIR/log/

. $HOME/.bash_profile
$CRONDIR/11_getdata.mh17.sh
