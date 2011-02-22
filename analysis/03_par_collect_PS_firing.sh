#!/bin/bash

if [ $# -ne 1 ]; then
	echo 'Usgage: 03_par_collect_PS_firing.sh <joblist output file name | stdout | exec>'
	echo 'Special argument:'
	echo '   - print:  print joblist to stdout'
	echo '   - exec:   execute the joblist directly'
	exit 1
fi

# load all common variables
. ./common_var.sh
# get tmp file
myname=`basename $0`
fntmp=`mktemp .${myname}.XXXXXX` || { echo "Failed to create temp file"; exit 1; }

for fn in `'ls' -1d ${dirmg}/*.mwk`; do
	bname=`basename ${fn}`
	fnpsf=`basename ${fn} .mwk`.psf.pk

	# if there's already .psf in ${dirpp}, ignore it.
	# it is likely to be processed already.
	if [ -f ${dirpp}/${fnpsf} ]; then
		continue
	fi

	if [[ $bname == *pos* || $bname == *POS* ]]; then
		echo "${bin}/collect_PS_firing.py ${fn} ${dirpp}/${fnpsf} ${defdelay} ${defnelec} extinfo c_success=images_shown" >> $fntmp
	elif [[ $bname == *Chou* ]]; then
		echo "${bin}/collect_PS_firing.py ${fn} ${dirpp}/${fnpsf} ${defdelay} ${defnelec} extinfo" >> $fntmp
	elif [[ $bname == *RF* ]]; then
		echo "${bin}/collect_PS_firing.py ${fn} ${dirpp}/${fnpsf} ${defdelay} ${defnelec} extinfo" >> $fntmp
	elif [[ $bname == *on*off* ]]; then
		echo "${bin}/collect_PS_firing.py ${fn} ${dirpp}/${fnpsf} ${defdelay} ${defnelec} reject_sloppy exclude_img=circ_mask t_stop=450000" >> $fntmp
	else
		echo "${bin}/collect_PS_firing.py ${fn} ${dirpp}/${fnpsf} ${defdelay} ${defnelec}" >> $fntmp
	fi
done

if [ $1 == 'print' ]; then
	cat $fntmp
elif [ $1 == 'exec' ]; then
	${bin}/parrun.py $fntmp
else
	mv $fntmp $1
fi

# finished all jobs. remove tmp file
rm -f $fntmp
