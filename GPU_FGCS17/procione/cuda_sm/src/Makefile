CXX=g++
CUDA=nvcc
CXXFLAGS=-O3 -fopenmp
CUDAFLAGS=-O3 -std=c++11 -arch=sm_30

SOURCES = kernel.cu main.cu
HEADERS = definitions.cuh
OBJ     = kernel.o
PROJECT_NAME = project

all: main

main.o: main.cu $(HEADERS)
	$(CUDA) $(CUDAFLAGS) -c main.cu

kernel.o: kernel.cu $(HEADERS)
	$(CUDA) $(CUDAFLAGS) -c kernel.cu

main: $(OBJ) main.o
	$(CXX) $(CXXFLAGS) $(OBJ) main.o -o ../bin/$(PROJECT_NAME) -lcublas -lcudart

.PHONY: clean

clean:
	rm -f *.o ../bin/$(PROJECT_NAME) *~
