#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "funcoes/funcoes.h"

int main() {
    int m = 100000; // linhas (amostras)
    int n = 500;  // colunas (variáveis)

    printf("--- BENCHMARK SERIAL DO PCA ---\n");
    printf("Matriz: %d x %d\n\n", m, n);

    // Variáveis para cronometrar
    clock_t inicio, fim;
    double tempo_gasto;

    // --- 1. ALOCAÇÃO E GERAÇÃO DE DADOS ---
    inicio = clock();
    double **X = aloca_matriz(m, n);
    srand(42); // Semente fixa para os testes serem sempre idênticos
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = ((double)rand() / RAND_MAX) * 100.0;
        }
    }
    fim = clock();
    tempo_gasto = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    printf("[1] Geracao dos Dados:    %f segundos\n", tempo_gasto);

    // --- 2. NORMALIZAÇÃO ---
    inicio = clock();
    double **X_norm = normaliza_matriz(X, m, n);
    fim = clock();
    tempo_gasto = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    printf("[2] Normalizacao:         %f segundos\n", tempo_gasto);

    // Preparação para o SVD
    double **V = aloca_matriz(n, n);
    double *autovalores = inicializa_vetor(n);

    // --- 3. SVD DE JACOBI UNILATERAL (O Motor Pesado) ---
    inicio = clock();
    svd_jacobi_unilateral(X_norm, m, n, V, autovalores);
    fim = clock();
    tempo_gasto = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    printf("[3] SVD de Jacobi:        %f segundos\n", tempo_gasto);

    // --- 4. ORDENAÇÃO ---
    inicio = clock();
    ordena_componentes(autovalores, V, n);
    fim = clock();
    tempo_gasto = ((double)(fim - inicio)) / CLOCKS_PER_SEC;
    printf("[4] Ordenacao:            %f segundos\n", tempo_gasto);

    printf("\nProcessamento concluido!\n");

    // Limpeza
    libera_matriz(X, m);
    libera_matriz(X_norm, m);
    libera_matriz(V, n);
    libera_vetor(autovalores);

    return 0;
}
