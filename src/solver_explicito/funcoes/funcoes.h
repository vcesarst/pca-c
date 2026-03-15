#ifndef FUNCOES_H
#define FUNCOES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// CRIME 1: Alocação Fragmentada (Ponteiro de Ponteiros)
double** aloca_matriz_ingenua(int m, int n) {
    double **matriz = (double**)malloc(m * sizeof(double*));
    for(int i = 0; i < m; i++) {
        matriz[i] = (double*)calloc(n, sizeof(double));
    }
    return matriz;
}

void libera_matriz_ingenua(double **matriz, int m) {
    for(int i = 0; i < m; i++) {
        free(matriz[i]);
    }
    free(matriz);
}

// Acesso linha a linha para normalizar colunas (Cache Misses)
void normaliza_matriz_ingenua(double **X, int m, int n) {
    double inv_m = 1.0 / (double)m; 
    for (int j = 0; j < n; j++) {
        double soma = 0.0;
        for (int i = 0; i < m; i++) {
            soma += X[i][j];
        }
        double media = soma * inv_m; 
        for (int i = 0; i < m; i++) {
            X[i][j] -= media;
        }
    }
}

// CRIME 2 e 3: O Pesadelo do Prefetcher
void calcula_covariancia_ingenua(double **X, int m, int n, double **C) {
    double inv_m_minus_1 = 1.0 / (double)(m - 1);
    
    // Ignorando a simetria: calculando n*n posições em vez de (n*n)/2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double soma = 0.0;
            
            // Loop Interno: Onde o hardware chora.
            // X[k] obriga a CPU a consultar uma tabela de ponteiros na RAM
            // para achar onde a linha 'k' está alocada, e depois pegar as colunas 'i' e 'j'.
            for (int k = 0; k < m; k++) {
                soma += X[k][i] * X[k][j];
            }
            
            C[i][j] = soma * inv_m_minus_1;
        }
    }
}

// Jacobi rodando numa matriz fragmentada (Mesmo sendo 500x500, o acesso é ineficiente)
void jacobi_covariancia_ingenua(double **C, int n, double **V, double *S) {
    int max_iter = 100;
    double tol = 1e-9;
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) V[i][j] = 0.0;
        V[i][i] = 1.0;
    }

    for (int iter = 0; iter < max_iter; iter++) {
        double max_off = 0.0;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double abs_val = fabs(C[i][j]);
                if (abs_val > max_off) max_off = abs_val;
            }
        }
        if (max_off < tol) break;

        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double Cpq = C[p][q];
                if (fabs(Cpq) > 1e-12) { 
                    double Cpp = C[p][p];
                    double Cqq = C[q][q];
                    
                    double theta = (Cqq - Cpp) / (2.0 * Cpq);
                    double t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                    if (theta < 0.0) t = -t;
                    
                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = c * t;

                    for (int k = 0; k < n; k++) {
                        if (k != p && k != q) {
                            double Cpk = C[p][k];
                            double Cqk = C[q][k];
                            C[p][k] = c * Cpk - s * Cqk;
                            C[k][p] = C[p][k]; 
                            
                            C[q][k] = s * Cpk + c * Cqk;
                            C[k][q] = C[q][k]; 
                        }
                    }
                    
                    C[p][p] = c * c * Cpp - 2.0 * s * c * Cpq + s * s * Cqq;
                    C[q][q] = s * s * Cpp + 2.0 * s * c * Cpq + c * c * Cqq;
                    C[p][q] = 0.0;
                    C[q][p] = 0.0;

                    for (int k = 0; k < n; k++) {
                        double Vkp = V[k][p];
                        double Vkq = V[k][q];
                        V[k][p] = c * Vkp - s * Vkq;
                        V[k][q] = s * Vkp + c * Vkq;
                    }
                }
            }
        }
    }
    for (int i = 0; i < n; i++) S[i] = C[i][i];
}

void ordena_componentes_ingenua(double *S, double **V, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (S[i] < S[j]) {
                double tempS = S[i];
                S[i] = S[j];
                S[j] = tempS;
                for (int k = 0; k < n; k++) {
                    double tempV = V[k][i];
                    V[k][i] = V[k][j];
                    V[k][j] = tempV;
                }
            }
        }
    }
}
#endif // FUNCOES_H
