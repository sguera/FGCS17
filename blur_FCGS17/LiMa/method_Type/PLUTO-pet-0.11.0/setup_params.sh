#!/bin/bash

#the environment PARAM_VALUES contains the values of the parameters
#project descriptor passed as first argument

#project descriptor
if [ "$1" == "" ]; then
	PD=../.project
else
	PD=$1
fi

#read parameters names and put into an array
read -ra PARAMS_NAME <<<"$(gawk -F '[""]' '/list/{print $4}' $PD)"

count=0

rm -rf src/dim_input.h

touch src/dim_input.h

echo "Compiling with parameters: $PARAM_VALUES"

for param in $PARAM_VALUES
do
	echo "#define ${PARAMS_NAME[$count]} $param" >> src/dim_input.h
	count=$((count + 1))
done
