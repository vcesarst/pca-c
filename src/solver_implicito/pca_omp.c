#include <stdio.h>
#include <stdlib.h>
#include <omp.h> // Substitui o <time.h> para benchmarks paralelos
#include "funcoes/funcoes_omp.h"

// Função para ter a matriz de output sem travar o terminal
void print_amostra_matriz(double *matriz, int m, int n, char *label) {
    int viz_linhas = 5; 
    int viz_cols = 3;   

    printf("\n--- Amostra da Matriz: %s ---\n", label);
    for (int i = 0; i < m; i++) {
        if (i < viz_linhas || i >= m - viz_linhas) {
            for (int j = 0; j < n; j++) {
                if (j < viz_cols || j >= n - viz_cols) {
                    printf("%8.4f ", matriz[j * m + i]);
                } else if (j == viz_cols) {
                    printf("  ...   ");
                }
            }
            printf("\n");
        } else if (i == viz_linhas) {
            printf("   ... (omitindo %d linhas) ...\n", m - (viz_linhas * 2));
        }
    }
    printf("------------------------------------------\n");
}

int main() {
    int m = 100000; 
    int n = 500;    
    
    printf("==========================================\n");
    printf("   PCA BENCHMARK: C PARALELO (OPENMP)     \n");
    printf("==========================================\n");
    printf("Dimensoes: %d linhas x %d colunas\n", m, n);
    
    // Painel de Controle HPC: Confirma se o OpenMP acordou os núcleos
    printf("Threads Disponiveis: %d\n", omp_get_max_threads()); 
    printf("==========================================\n");

    // Variaveis de tempo do OpenMP (já retornam em segundos)
    double start, end; 
    
    // --- 1. Alocação ---
    double *X = aloca_matriz(m, n);
    double *V = aloca_matriz(n, n);
    double *S = (double*)malloc(n * sizeof(double));

    // --- 2. Geração de Dados (Determinística) ---
    start = omp_get_wtime();
    for (int j = 0; j < n; j++) {
        int offset = j * m;
        for (int i = 0; i < m; i++) {
            X[offset + i] = (double)((i + j) % 100); 
        }
    }
    end = omp_get_wtime();
    printf("[1 & 2] Alocacao e Geracao: %.4f s\n", end - start);

    // --- 3. Normalização ---
    start = omp_get_wtime();
    normaliza_matriz(X, m, n);
    end = omp_get_wtime();
    printf("[3] Normalizacao:           %.4f s\n", end - start);

    // --- 4. SVD de Jacobi ---
    printf("[4] Executando SVD Paralelo... (Aguarde)\n");
    start = omp_get_wtime();
    svd_jacobi_unilateral(X, m, n, V, S);
    end = omp_get_wtime();
    double tempo_svd = end - start;
    printf("[4] SVD Finalizado em:      %.4f s (%.2f min)\n", tempo_svd, tempo_svd / 60.0);

    // --- 5. Ordenação ---
    ordena_componentes(S, V, n);

    // --- RESULTADOS PARA VALIDAÇÃO ---
    printf("\n=== VALIDACAO DOS RESULTADOS ===\n");
    printf("Primeiros 5 Autovalores (Variancia Explicada):\n");
    for(int i = 0; i < 5; i++) {
        printf("PC%d: %.4f\n", i+1, S[i]);
    }
    printf("==========================================\n");

    libera_matriz(X);
    libera_matriz(V);
    free(S);

    return 0;
}
