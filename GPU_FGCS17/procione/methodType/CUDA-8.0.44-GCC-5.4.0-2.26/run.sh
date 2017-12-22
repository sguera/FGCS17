#!/bin/bash

#Environment variable $CURRENT_NT contains the number of threads to use
#output file name passed as first argument

module_home=`pwd`

#LIKWID_CMD is exported in the environment by the main run script

if [ "$1" == "" ]; then
	echo "Output file not specified. Output on screen."
	CUDA_VISIBLE_DEVICES=0 $LIKWID_CMD $module_home/bin/project
else
	CUDA_VISIBLE_DEVICES=0 $LIKWID_CMD $module_home/bin/project >> $1
fi

