#!/bin/bash

#patt="Chabo_????????_RSVP*.mwk"
patt="Chabo*.mwk"
dirmwk="data_mwk"
dirnev="data_nev"
test $dirmg || dirmg="data_merged"
test $dirpp || dirpp="data_postproc"
sep=", "
bin="bin"

defdelay=300   # default delay
defnelec=96    # default number of electrodes active
allnelec=128   # all electrodes
