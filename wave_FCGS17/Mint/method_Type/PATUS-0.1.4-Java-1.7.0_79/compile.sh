#!/bin/bash

module_home=`pwd`

stencil_name='example.stc'

#to compile a patus implementation you need to generate the code starting from the stencil definition and then call the makefile
cd $EBROOTPATUS
source util/patusvars.sh
#go back
cd $module_home

#create folder where to put the binary files and the outputs
if [ ! -d bin ]; then
	mkdir bin
fi
if [ ! -d out ]; then
	mkdir out
fi

#Generate Code in src directory. You can customize this command
patus ./src/$stencil_name --outdir=./src --architecture="x86_64 AVX"
#patus ./src/$stencil_name --outdir=./src --architecture="x86_64 SSE"
cd src
#Compile source
make

cd $module_home
