#!/bin/bash
# This script is meant to be called by 11_getdata.mh17.sh
lcldir=`dirname $1`/`basename $1`/

for f0 in `find $lcldir -name '*.bad'`; do
	f=`echo $f0 | sed -e "s:${lcldir}::"`
	echo $f 
done
