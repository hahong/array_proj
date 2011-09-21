#!/bin/bash
# load all common variables
. ./common_var.sh

test "$dirout" || dirout="subproj/100_PSTH"
wd=$dirdat/$pridat/$dirout

for pkf in `'ls' -1 $dirpp/*.psf.pk`; do
	bname=`basename $pkf .psf.pk`
	if [ -f $wd/$bname.pdf ]; then 
		continue 
	fi

	if [[ $bname == *MovieGallant110413* ]]; then
		echo bin/plot_PSTH.py $wd/$bname $pkf --n_bins=2600
	else
		echo bin/plot_PSTH.py $wd/$bname $pkf
	fi
done
