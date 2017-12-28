#!/bin/bash

module_home=`pwd`

#create folder where to put the binary files and the outputs
if [ ! -d bin ]; then
	mkdir bin
fi
if [ ! -d out ]; then
	mkdir out
fi

#to compile an openmp implementation you need only to call the makefile
cd src
#Compile source
make clean
make all

cd $module_home
