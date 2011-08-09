#!/bin/bash
test "$HOME" || HOME=/mindhive/dicarlolab/u/hahong/
PROJROOT=$HOME/teleport/array/
CRONDIR=$PROJROOT/analysis/scheduled/
LOGDIR=$CRONDIR/log/

. $HOME/.bash_profile
. $CRONDIR/11_getdata.mh17.sh
