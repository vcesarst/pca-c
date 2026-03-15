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
serial.x: main_prof.c
	$(CC) $(FLAGS) main_prof.c -o serial.x $(LIBS)

# Regra do código otimizado
serial_opt.x: pca_opt_serial.c
	$(CC) $(FLAGS) pca_opt_serial.c -o serial_opt.x $(LIBS)

# Regra do código paralelo
serial_omp.x: pca_omp.c
	$(CC) $(FLAGS) -fopenmp pca_omp.c -o serial_omp.x $(LIBS)

# Limpeza
clean:
	rm -f serial.x serial_opt.x serial_omp.x	
