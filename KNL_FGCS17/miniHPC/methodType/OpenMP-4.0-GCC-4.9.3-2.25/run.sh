#!/bin/bash

#Environment variable $CURRENT_NT contains the number of threads to use
#output file name passed as first argument

module_home=`pwd`

export OMP_NUM_THREADS=$CURRENT_NT

#LIKWID_CMD is exported in the environment by the main run script

if [ "$1" == "" ]; then
	echo "Output file not specified. Output on screen."
	$LIKWID_CMD $module_home/bin/project
else
	$LIKWID_CMD $module_home/bin/project >> $1
fi

