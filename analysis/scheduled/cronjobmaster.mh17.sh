#!/bin/bash
test "$HOME" || HOME=/mindhive/dicarlolab/u/hahong/
PROJROOT=$HOME/teleport/array/
CRONDIR=$PROJROOT/analysis/scheduled/
LOGDIR=$CRONDIR/log/

. $HOME/.bash_profile
. $CRONDIR/10_PSTH.mh17.sh
