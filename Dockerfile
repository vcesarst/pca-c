FROM alpine:latest

RUN apk add --no-cache build-base gcc libgomp

WORKDIR /app
COPY . .

# 1. Compila-se o arquivo e salva-se o .x lá dentro da pasta solver_explicito
RUN gcc -O3 -ffast-math -funroll-loops -fexpensive-optimizations -march=native -fopenmp src/solver_explicito/main_embarcado.c -o src/solver_explicito/pca_embarcado.x -lm

ENV OMP_NUM_THREADS=4

# 2. Entra-se na pasta do solver_explicito ANTES de rodar
WORKDIR /app/src/solver_explicito

# 3. Executa-se o arquivo .x
CMD ["./pca_embarcado.x"]
