#!/bin/bash

# load all common variables
. ./common_var.sh

myname=`basename $0`
tmpdp=`mktemp .${myname}.XXXXXX` || { echo "Failed to create temp file"; exit 1; }
tmpps=`mktemp .${myname}.XXXXXX` || { echo "Failed to create temp file"; exit 1; }

for fn in `'ls' -1r ${dirpp}/*.psf.pk`; do
	bname=`basename ${fn} .psf.pk`

	# if there's already .PSTH.pdf in ${dirpp}, ignore it.
	if [ -f ${dirpp}/${bname}.PSTH.pdf ]; then
		continue
	fi

	# if it's nicole's data
	if [[ $bname == *Nicole* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=100,106 corr" >> $tmpdp
		echo "${bin}/expr_anal.py ps ${fn} ${dirpp}/${bname} blanks=100,106" >> $tmpps
	# if it's nuo's data
	elif [[ $bname == *Nuo* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} corr" >> $tmpdp
		echo "${bin}/expr_anal.py ps ${fn} ${dirpp}/${bname}" >> $tmpps
	fi
done

cat $tmpdp $tmpps

rm -f $tmpdp
rm -f $tmpps
