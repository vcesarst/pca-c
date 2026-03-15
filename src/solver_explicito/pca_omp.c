#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "funcoes/funcoes_omp.h"

int main() {
    int m = 100000; 
    int n = 500;    
    
    printf("==========================================\n");
    printf(" PCA EXPLÍCITO EXTREMO: OPENMP + COVARIÂNCIA\n");
    printf("==========================================\n");
    printf("Dimensoes X: %d x %d\n", m, n);
    printf("Threads Ativas: %d\n", omp_get_max_threads());
    printf("==========================================\n");

    double start, end;
    
    // --- 1. Alocação ---
    double *X = aloca_matriz(m, n);
    double *C = aloca_matriz(n, n); 
    double *V = aloca_matriz(n, n);
    double *S = (double*)malloc(n * sizeof(double));

    // --- 2. Geração de Dados ---
    start = omp_get_wtime();
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < n; j++) {
        int offset = j * m;
        for (int i = 0; i < m; i++) {
            X[offset + i] = (double)((i + j) % 100); 
        }
    }
    end = omp_get_wtime();
    printf("[1] Geracao de Dados:       %.4f s\n", end - start);

    // --- 3. Normalização OMP ---
    start = omp_get_wtime();
    normaliza_matriz_omp(X, m, n);
    end = omp_get_wtime();
    printf("[2] Normalizacao OMP:       %.4f s\n", end - start);

    // --- 4. Construção da Matriz de Covariância OMP ---
    printf("[3] Calculando X^T * X...   (Voando!)\n");
    start = omp_get_wtime();
    calcula_covariancia_omp(X, m, n, C);
    end = omp_get_wtime();
    double tempo_cov = end - start;
    printf("    -> Covariancia pronta em: %.4f s\n", tempo_cov);

    // --- 5. SVD de Jacobi (Cache L3) ---
    printf("[4] Executando Autovetores... \n");
    start = omp_get_wtime();
    jacobi_covariancia_serial(C, n, V, S);
    end = omp_get_wtime();
    double tempo_jacobi = end - start;
    printf("    -> Jacobi finalizado em:  %.4f s\n", tempo_jacobi);

    // --- 6. Ordenação ---
    ordena_componentes(S, V, n);

    // --- RESULTADOS TOTAIS ---
    printf("\n=== RESUMO DE PERFORMANCE PARALELA ===\n");
    printf("Tempo Matriz Covariancia: %.4f s\n", tempo_cov);
    printf("Tempo Calculo Autovetores: %.4f s\n", tempo_jacobi);
    printf("TEMPO MATEMATICO TOTAL:   %.4f s\n", tempo_cov + tempo_jacobi);

    // --- VALIDACAO ---
    printf("\n=== VALIDACAO NUMERICA ===\n");
    for(int i = 0; i < 5; i++) {
        printf("PC%d: %.4f\n", i+1, S[i]);
    }
    printf("==========================================\n");

    libera_matriz(X);
    libera_matriz(C);
    libera_matriz(V);
    free(S);

    return 0;
}
