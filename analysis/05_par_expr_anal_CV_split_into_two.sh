#!/bin/bash

# load all common variables
. ./common_var.sh

for fn in `'ls' -1r ${dirpp}/*.psf.pk`; do
	bname=`basename ${fn} .psf.pk`

	# if there's already .PSTH.pdf in ${dirpp}, ignore it.
	if [ -f ${dirpp}/${bname}.dp_all_sp1.csv ]; then
		continue
	fi

	# if it's nicole's data
	if [[ $bname == *Nicole* ]]; then
		# first half
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=100,106 prange=0,0.5 suffix=_sp0"
		# second half
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=100,106 prange=0.5,1 suffix=_sp1"
	# if it's nuo's data
	elif [[ $bname == *Nuo* ]]; then
		# first half
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} prange=0,0.5 suffix=_sp0"
		# second half
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} prange=0.5,1 suffix=_sp1"
	fi
done
