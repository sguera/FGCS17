#!/bin/bash

#Environment variable $CURRENT_NT contains the number of threads to use
#output file name passed as first argument

module_home=`pwd`

export OMP_NUM_THREADS=$CURRENT_NT

read -ra PARAMS_NAME <<<"$(gawk -F '[""]' '/list/{print $4}' ../.project)"

#keeps track of tuning
file_bench_time=0

if `cmp --silent src/bench bin/bench` ; then
	#echo "Exec files not differ. Use old bin/bench!"
	file_bench_time=`stat -c "%Z" ./bin/bench`
else
	#echo "Exec files differ. Updated bin/bench!"
	cp -R src/bench bin/bench
	file_bench_time=`stat -c "%Z" ./src/bench`
fi

#create name of tune file
tune_file=$module_home/src/tune-${PARAM_VALUES//" "/"-"}-$CURRENT_NT

if [ ! -f $tune_file ] || [ `stat -c "%Z" $tune_file` -le $file_bench_time ]; then
	#AUTOTUNING PHASE
	cd src
	echo '------- Autotuning STARTED! -------'
	arg=""
	count=0

	for param in $PARAM_VALUES
	do
		arg="$arg ${PARAMS_NAME[$count]}=$param"
		count=$((count + 1))
	done

#echo "ARGS = "$arg
	#remove old file
	rm -f tune*.txt

	#run autotuning
	make tune$arg

	params=$(gawk '/Optimal parameter configuration found/{getline; print $0}' tune*.txt)

	#update tune_file
	if [ "$params" == "" ]; then
		echo ""
		echo "(ERROR) Tuning error with $t thread(s), execution aborted!"
		echo ""
#		continue
	fi
	echo "$params" > $tune_file

	echo '------- Autotuning FINISHED! -------'
fi

tuned_params=`less $tune_file`
echo "tuned params: $tuned_params"

if [ "$1" == "" ]; then
	echo "Output file not specified. Output on screen."
	$module_home/bin/bench $tuned_params
else
	echo "tuned params: $tuned_params" >> $1
	$module_home/bin/bench $tuned_params >> $1
fi

