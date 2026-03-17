#!/bin/bash


# compila o codigo main_embarcado.c
gcc -O3 -march=native -ffast-math -fexpensive-optimizations -funroll-loops -fopenmp -o main_embarcado.x main_embarcado.c -lm

# executa o programa com 4 threads
export OMP_NUM_THREADS=4
./main_embarcado.x
