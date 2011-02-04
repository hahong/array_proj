#!/bin/bash

if [ $# -ne 1 ]; then
	echo 'Usgage: 02_par_merge.sh <joblist output file name | stdout | exec>'
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

for fn in `'ls' -1 ${dirmwk}/${patt}`; do
	bname=`basename ${fn}`
	# if there's already mwk directory in ${dirmg}, ignore it.
	# it is likely to be merged already.
	if [ -d ${dirmg}/${bname} ]; then
		continue
	fi

	if [ -d ${fn} ]; then
		fnmwk=${fn}/${bname}
		if [ ! -f ${fnmwk} ]; then
			continue
		fi
	else
		fnmwk=${fn}
	fi

	fnnev="${dirnev}/`basename ${bname} .mwk`.nev"
	if [ ! -f ${fnnev} ]; then
		continue
	fi

	# basenames
	bnmwk=`basename ${fnmwk}`
	bnnev=`basename ${fnnev}`

	echo "cp ${fnmwk} ${dirmg}/${bnmwk} && ${bin}/merge.py ${dirmg}/${bnmwk} ${fnnev} nowav && rm -f ${dirmg}/${bnmwk}/*.bak" >> $fntmp
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
