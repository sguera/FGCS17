GCC = gcc
LDFLAGS= -O3 -mavx512f -mavx512er -mavx512cd -mavx512pf -lm 
# gnu
OMP_FLAG = -fopenmp 

SRC = example.c 
PROJECT_NAME = project

all: module

module: $(SRC)
	$(GCC) $(OMP_FLAG) $+ ${LDFLAGS} -o ../bin/$(PROJECT_NAME)

.PHONY: clean

clean:
	rm -rf ../bin/$(PROJECT_NAME) *~
