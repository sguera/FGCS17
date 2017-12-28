#CC=icc
CC = gcc

PLCFLAGS += --parallel --pet
#general
CFLAGS = -DTIME -DVERIFY 
LDFLAGS = -lm

#Pluto
PLC = polycc

#NTHREADS=4

ifeq ($(CC), icc)
	OPT_FLAGS     := -O3 -xHost -ansi-alias -ipo -fp-model precise
	PAR_FLAGS     := -parallel
	OMP_FLAGS     := -openmp
else
	# for gcc
	OPT_FLAGS     := -O3 -march=native -mtune=native -ftree-vectorize -frecord-gcc-switches
	PAR_FLAGS     := -ftree-parallelize-loops=$(NTHREADS)
	OMP_FLAGS     := -fopenmp
endif

TILEFLAGS +=
LIKWID_FLAGS +=


SRC = example
PROJECT_NAME = project

all: module

$(SRC).par.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel $(TILEFLAGS) $(PLCFLAGS)  -o $@
#diamond tiling
$(SRC).lbpar.c:  $(SRC).c
	$(PLC) $(SRC).c --tile --parallel --partlbtile $(TILEFLAGS) $(PLCFLAGS) -o $@

par: $(SRC).par.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).par.c -o ../bin/$(PROJECT_NAME) $(LDFLAGS)
module: $(SRC).lbpar.c
	$(CC) $(OPT_FLAGS) $(CFLAGS) $(OMP_FLAGS) $(SRC).lbpar.c -o ../bin/$(PROJECT_NAME) $(LDFLAGS)	

.PHONY: clean

clean:
	rm -rf ../bin/$(PROJECT_NAME) *~ *.lbpar.c
