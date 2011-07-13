#!/bin/bash

#patt="Chabo_????????_RSVP*.mwk"
test "$patt" || patt="Tito*.mwk"
dirmwk="data_mwk"
test $dirnev || dirnev="data_nev_tmp"
test $dirmg || dirmg="data_merged_tmp"
test $dirpp || dirpp="data_postproc_tmp"
sep=", "
bin="bin"

defdelay=300   # default delay
defnelec=96    # default number of electrodes active
allnelec=128   # all electrodes
