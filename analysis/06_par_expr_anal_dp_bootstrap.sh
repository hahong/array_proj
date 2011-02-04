#!/bin/bash

# load all common variables
. ./common_var.sh

for fn in `'ls' -1r ${dirpp}/*.psf.pk`; do
	bname=`basename ${fn} .psf.pk`

	if [ -f ${dirpp}/${bname}.dp_all*.csv ]; then
		continue
	fi

	# bin/expr_anal.py dp prefix.psf.pk prefix bootstrap=200 suffix=_bstrp extbsinfo=data_postproc/Chabo_20100822_RSVPNuo_P_001.dp_bstrp.pk
	# if it's nicole's invariance data
	if [[ $bname == *Nicole10x6* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=60,62 bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# chou's invariance data
	elif [[ $bname == *Chou* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=0,4 bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# RF mapping data
	elif [[ $bname == *RF* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} extfirinfo"
	# if it's nicole's full data
	elif [[ $bname == *Nicole305* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=300,306 bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# if it's nicole's data
	elif [[ $bname == *Nicole* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=100,106 bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# if it's var00'
	elif [[ $bname == *Var00* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=640,650 bootstrap=200 suffix=_bstrp extfirinfo"
	# if it's var03'
	elif [[ $bname == *Var03* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=2560,2570 bootstrap=200 suffix=_bstrp extfirinfo"
	# if it's nuo's data
	elif [[ $bname == *Nuo* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# position leraning task based-on Nuo's images
	elif [[ $bname == *pos* || $bname == *POS* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=15,18 bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# Swap task
	elif [[ $bname == *SWAP* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} bootstrap=200 suffix=_bstrp extbsinfo=${dirpp}/${bname}.dp_extbsinfo.csv extfirinfo"
	# Mojo and Pogo old data
	elif [[ $bname == *Pogo* || $bname == *Mojo* ]]; then
		echo "${bin}/expr_anal.py dp ${fn} ${dirpp}/${bname} blanks=100,106 extfirinfo proc_cluster"
	fi
done
