#!/bin/bash

if [ $# -lt 1 ]; then
	echo 'Usgage: 02_par_merge.sh <joblist output file name | print | exec> [options]'
	echo 'Special argument:'
	echo '   - print:  print joblist to stdout'
	echo '   - exec:   execute the joblist directly'
	echo
	echo 'Options:'
	echo '   - cluster <mwk file list> <cluster prefix file list> <suffix> [options to be passed]: merge with cluster info'
	echo '   - multinsp <NSP name>'
	exit 1
fi

# load all common variables
. ./common_var.sh
# get tmp file
myname=`basename $0`
fntmp=`mktemp .${myname}.XXXXXX` || { echo "Failed to create temp file"; exit 1; }
files=(`'ls' -1 ${dirmwk}/${patt}`)
opts="none"

# process options
if [ $# -ge 5 -a "$2" = "cluster" ]; then
	files=(`cat $3`)
	cfiles=(`cat $4`)
	cfileslst="$4"
	csuffix="$5"
	cextopts="$6"
	opts="cluster"

	if [ "${#files[@]}" != "${#cfiles[@]}" ]; then
		echo 'Number of files mismatch' >&2
		exit 1
	fi
elif [ $# -ge 3 -a "$2" = "multinsp" ]; then
	nsp="$3"
	opts="multinsp"
	files0=(`'ls' -1 ${dirnev}/*${nsp}*.nev`)
	files=()
	for index in "${!files0[@]}"; do
		files[$index]="${dirmwk}/`basename ${files0[$index]} .nev`.mwk"
	done
fi

for index in "${!files[@]}"; do
        mwksrc=${files[$index]}
	mwksrcb=`basename ${mwksrc}`
	mwkdstb=$mwksrcb
	nevsrc="${dirnev}/`basename ${mwksrcb} .mwk`.nev"

	if [ "$opts" == "cluster" ]; then
		# -- clustering
		mwksrc="${dirmwk}/`basename ${mwksrc}`"
		mwkdstb="`basename ${mwksrc} .mwk`${csuffix}.mwk"
		nevsrc="${dirnev}/`basename ${mwksrcb} .mwk`.nev"
	elif [ "$opts" == "multinsp" ]; then
		# -- multiple NSPs
		if [[ $mwksrc != *$nsp* ]]; then
			continue    # discard non-multi-nsp data
		fi
		# removes $nsp and the first preceeding suffix
		# e.g., Tito_20110809_RSVPVar03a_S110720A_NSP2_P_001 -> Tito_20110809_RSVPVar03a_S110720A_001
		mwksrc=`python -c "s0 = '${mwksrc}'; s1 = s0.split('_${nsp}_'); s1a = s1[0]; s1b = '_'.join(s1[1].split('_')[1:]); print s1a + '_' + s1b"`
		mwksrcb=`basename ${mwksrc}`
		mwkdstbhost=$mwksrcb
	fi


	# if there's already mwk directory in ${dirmg}, ignore it.
	# it is likely to be merged already.
	if [ -d ${dirmg}/${mwkdstb} ]; then
		continue
	fi

	# already indexed?
	if [ -d ${mwksrc} ]; then
		mwksrc=${mwksrc}/${mwksrcb}
		if [ ! -f ${mwksrc} ]; then
			continue
		fi
	fi

	# check nev files
	if [ ! -f ${nevsrc} ]; then
		continue
	fi
	
	oargs=""
	# other options
	if [ "$opts" = "cluster" ]; then
		cfile0=${cfiles[$index]}
		oargs="cluster=${cfile0} cluster_all=+${cfileslst}"
	#elif [ "$opts" = "multinsp" ]; then
		#oargs="timetransf=${dirmg}/${mwkdstbhost}"
	fi

	echo "cp ${mwksrc} ${dirmg}/${mwkdstb} && ${bin}/merge.py ${dirmg}/${mwkdstb} ${nevsrc} nowav ${oargs} ${cextopts} && rm -f ${dirmg}/${mwkdstb}/*.bak" >> $fntmp
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
