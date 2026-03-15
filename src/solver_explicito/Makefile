# ==========================================
# MAKEFILE - BENCHMARK DE ALTA PERFORMANCE
# ==========================================

CC = gcc
LIBS = -lm

# As flags extremas de otimização
FLAGS = -O3 -march=native -ffast-math -funroll-loops -fexpensive-optimizations -foptimize-register-move

# Alvos
all: serial.x serial_opt.x serial_omp.x

# Regra do código ingênuo
serial.x: serial_pca.c
	$(CC) $(FLAGS) serial_pca.c -o serial.x $(LIBS)

# Regra do código otimizado
serial_opt.x: pca_serial_opt.c
	$(CC) $(FLAGS) pca_serial_opt.c -o serial_opt.x $(LIBS)

# Regra do código paralelo
serial_omp.x: pca_omp.c
	$(CC) $(FLAGS) -fopenmp pca_omp.c -o serial_omp.x $(LIBS)

# Limpeza
clean:
	rm -f serial.x serial_opt.x serial_omp.x	
