#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "funcoes/funcoes_omp.h" 

int main() {
    printf("==========================================\n");
    printf(" PCA EMBARCADO: DADOS REAIS (Breast Cancer)\n");
    printf("==========================================\n");

    // --- 1. Abrir o arquivo gerado pelo Python ---
    FILE *file = fopen("../../examples/breast_cancer_data.txt", "r");
    if (file == NULL) {
        printf("[ERRO] Nao foi possivel abrir o arquivo de dados em ../../examples/\n");
        return 1;
    }

    // Leitura das dimensões dinamicamente
    int m, n;
    fscanf(file, "%d %d", &m, &n);
    printf("Dimensoes detectadas: %d linhas x %d colunas\n", m, n);

    // --- 2. Alocação (Usando a sua função!) ---
    double *X = aloca_matriz(m, n);
    double *C = aloca_matriz(n, n); 
    double *V = aloca_matriz(n, n);
    double *S = (double*)malloc(n * sizeof(double));

    // --- 3. Leitura e Conversão para Column-Major ---
    // O txt está em linhas (Row-Major), mas o C exige Colunas (Column-Major).
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            double valor;
            fscanf(file, "%lf", &valor);
            X[j * m + i] = valor; // Inserindo no formato otimizado do seu código!
        }
    }
    fclose(file);
    printf("[1] Dados carregados na RAM (Convertidos para Column-Major).\n");

    // ====================================================================
    // --- PCA ---
    // ====================================================================
    
    // Normalização
    double start = omp_get_wtime();
    normaliza_matriz_omp(X, m, n);
    double end = omp_get_wtime();
    printf("[2] Normalizacao OMP:       %.6f s\n", end - start);

    // Matriz de Covariância
    start = omp_get_wtime();
    calcula_covariancia_omp(X, m, n, C);
    end = omp_get_wtime();
    printf("[3] Covariancia OMP:        %.6f s\n", end - start);

    // SVD de Jacobi
    start = omp_get_wtime();
    jacobi_covariancia_serial(C, n, V, S);
    end = omp_get_wtime();
    printf("[4] Jacobi Finalizado:      %.6f s\n", end - start);

    // Ordenação
    ordena_componentes(S, V, n);

    // ====================================================================

    // --- 4. Exportar os Autovetores para o Gabarito ---
    FILE *out = fopen("../../examples/pca_c_autovetores.txt", "w");
    if (out != NULL) {
        // Salvando em formato Row-Major para o Python ler facilmente
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                // Acessando V em Column-Major (j * n + i)
                fprintf(out, "%.8f ", V[j * n + i]); 
            }
            fprintf(out, "\n");
        }
        fclose(out);
        printf("\n[SUCESSO] Autovetores salvos em 'examples/pca_c_autovetores.txt'.\n");
    }

    // --- 5. Validação Rápida no Terminal ---
    printf("\n=== TOP 5 AUTOVALORES (Variancia Explicada) ===\n");
    for(int i = 0; i < 5; i++) {
        printf("PC%d: %.4f\n", i+1, S[i]);
    }
    printf("==========================================\n");

    // --- 6. Limpeza ---
    libera_matriz(X);
    libera_matriz(C);
    libera_matriz(V);
    free(S);

    return 0;
}