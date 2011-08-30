#!/bin/bash

test "$pridat" || pridat="d003_Tito"     # set the primary data path (in data/)
test "$patt" || patt="*.mwk"             # pattern of interested mworks data

test "$dirdat" || dirdat="data"
test $dirmwk || dirmwk="$dirdat/$pridat/mwk"
test $dirnev || dirnev="$dirdat/$pridat/neudat"
test $dirmg || dirmg="$dirdat/$pridat/mwk_merged"
test $dirpp || dirpp="$dirdat/$pridat/data_postproc"
sep=", "
bin="bin"

defdelay=300   # default delay
defnelec=96    # default number of electrodes active
allnelec=128   # all electrodes
