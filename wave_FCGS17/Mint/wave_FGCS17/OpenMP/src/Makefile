GCC = gcc
LDFLAGS= -lm
# gnu
OMP_FLAG = -O3 -march=native -mtune=native -ffast-math -funroll-loops -fopenmp 

SRC = example.c 
PROJECT_NAME = project

all: module

module: $(SRC)
	$(GCC) $(OMP_FLAG) $+ ${LDFLAGS} -o ../bin/$(PROJECT_NAME)

.PHONY: clean

clean:
	rm -rf ../bin/$(PROJECT_NAME) *~
