#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "funcoes/funcoes.h"

int main() {
    int m = 100000; 
    int n = 500;    
    
    printf("==========================================\n");
    printf(" PCA EXPLÍCITO: MATRIZ DE COVARIÂNCIA (INGENUO)\n");
    printf("==========================================\n");

    clock_t start, end;
    
    // --- 1. Alocação Ingênua ---
    double **X = aloca_matriz_ingenua(m, n);
    double **C = aloca_matriz_ingenua(n, n); 
    double **V = aloca_matriz_ingenua(n, n);
    double *S = (double*)malloc(n * sizeof(double));

    // --- 2. Geração de Dados ---
    start = clock();
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < m; i++) {
            X[i][j] = (double)((i + j) % 100); 
        }
    }
    end = clock();
    printf("[1] Geracao de Dados:       %.4f s\n", (double)(end - start) / CLOCKS_PER_SEC);

    // --- 3. Normalização ---
    start = clock();
    normaliza_matriz_ingenua(X, m, n);
    end = clock();
    printf("[2] Normalizacao:           %.4f s\n", (double)(end - start) / CLOCKS_PER_SEC);

    // --- 4. Construção da Matriz de Covariância ---
    printf("[3] Calculando X^T * X...   (Vai demorar!)\n");
    start = clock();
    calcula_covariancia_ingenua(X, m, n, C);
    end = clock();
    double tempo_cov = (double)(end - start) / CLOCKS_PER_SEC;
    printf("    -> Covariancia pronta em: %.4f s\n", tempo_cov);

    // --- 5. SVD de Jacobi ---
    printf("[4] Executando Autovetores... (Aguarde)\n");
    start = clock();
    jacobi_covariancia_ingenua(C, n, V, S);
    end = clock();
    double tempo_jacobi = (double)(end - start) / CLOCKS_PER_SEC;
    printf("    -> Jacobi finalizado em:  %.4f s\n", tempo_jacobi);

    // --- 6. Ordenação ---
    ordena_componentes_ingenua(S, V, n);

    printf("\n=== RESUMO DE PERFORMANCE SERIAL ===\n");
    printf("Tempo Matriz Covariancia: %.4f s\n", tempo_cov);
    printf("Tempo Calculo Autovetores: %.4f s\n", tempo_jacobi);
    printf("TEMPO MATEMATICO TOTAL:   %.4f s\n", tempo_cov + tempo_jacobi);
    printf("==========================================\n");

    libera_matriz_ingenua(X, m);
    libera_matriz_ingenua(C, n);
    libera_matriz_ingenua(V, n);
    free(S);

    return 0;
}
