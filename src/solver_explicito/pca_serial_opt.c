#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "funcoes/funcoes_opt.h"

int main() {
    int m = 100000; 
    int n = 500;    
    
    printf("==========================================\n");
    printf(" PCA EXPLÍCITO: MATRIZ DE COVARIÂNCIA (SERIAL)\n");
    printf("==========================================\n");
    printf("Dimensoes X: %d x %d\n", m, n);
    printf("Dimensoes C: %d x %d\n", n, n);
    printf("==========================================\n");

    clock_t start, end;
    
    // --- 1. Alocação ---
    double *X = aloca_matriz(m, n);
    double *C = aloca_matriz(n, n); // A matriz de Covariância
    double *V = aloca_matriz(n, n);
    double *S = (double*)malloc(n * sizeof(double));

    // --- 2. Geração de Dados ---
    start = clock();
    for (int j = 0; j < n; j++) {
        int offset = j * m;
        for (int i = 0; i < m; i++) {
            X[offset + i] = (double)((i + j) % 100); 
        }
    }
    end = clock();
    printf("[1] Geracao de Dados:       %.4f s\n", (double)(end - start) / CLOCKS_PER_SEC);

    // --- 3. Normalização ---
    start = clock();
    normaliza_matriz(X, m, n);
    end = clock();
    printf("[2] Normalizacao:           %.4f s\n", (double)(end - start) / CLOCKS_PER_SEC);

    // --- 4. Construção da Matriz de Covariância ---
    printf("[3] Calculando X^T * X...   (Aguarde)\n");
    start = clock();
    calcula_covariancia(X, m, n, C);
    end = clock();
    double tempo_cov = (double)(end - start) / CLOCKS_PER_SEC;
    printf("    -> Covariancia pronta em: %.4f s\n", tempo_cov);

    // --- 5. SVD de Jacobi no Cache L3 ---
    printf("[4] Executando Autovetores... (Aguarde)\n");
    start = clock();
    jacobi_covariancia(C, n, V, S);
    end = clock();
    double tempo_jacobi = (double)(end - start) / CLOCKS_PER_SEC;
    printf("    -> Jacobi finalizado em:  %.4f s\n", tempo_jacobi);

    // --- 6. Ordenação ---
    ordena_componentes(S, V, n);

    // --- RESULTADOS TOTAIS ---
    printf("\n=== RESUMO DE PERFORMANCE SERIAL ===\n");
    printf("Tempo Matriz Covariancia: %.4f s\n", tempo_cov);
    printf("Tempo Calculo Autovetores: %.4f s\n", tempo_jacobi);
    printf("TEMPO MATEMATICO TOTAL:   %.4f s\n", tempo_cov + tempo_jacobi);

    // --- RESULTADOS PARA VALIDAÇÃO ---
    printf("\n=== VALIDACAO NUMERICA ===\n");
    printf("Primeiros 5 Autovalores (Variancia Explicada):\n");
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
