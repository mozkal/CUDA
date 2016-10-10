CC = g++
NVCC = nvcc

all: wn_gpu wn_cpu

wn_cpu: wn_cpu.cc declarations.h wunsch_needleman.o wn_cpu.o
	$(CC) -std=c++11 -g -Wall -pedantic -o wn_cpu wn_cpu.o wunsch_needleman.o

wn_gpu: declarations.h wunsch_needleman.o wn_gpu.o
	$(NVCC) -o wn_gpu wn_gpu.o wunsch_needleman.o

wunsch_needleman.o: declarations.h wunsch_needleman.cc
	$(CC) -c -std=c++11 -Wall -pedantic wunsch_needleman.cc

wn_cpu.o: wn_cpu.cc declarations.h
	$(CC) -c -std=c++11 -Wall -pedantic wn_cpu.cc

wn_gpu.o: wn_gpu.cu
	$(NVCC) -o wn_gpu.o -c wn_gpu.cu
