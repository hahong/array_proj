#!/bin/bash
# This script is meant to be called by 11_getdata.mh17.sh
badlist=$1
wd=`dirname $2`/`basename $2`

for fbadshort in `cat $badlist`; do
	forgshort=`dirname $fbadshort`/`basename $fbadshort .bad`
	forg="$wd/$forgshort"
	fdst="$wd/$fbadshort"
	
	if [ ! -f $forg ]; then
		continue
	fi
	if [ -f $fdst ]; then
		echo Already exists: $fdst
		continue
	fi

	# process only if forg exists and fdest doesn't exist.
	echo Apply: $forg
done
