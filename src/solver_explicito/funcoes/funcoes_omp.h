#ifndef FUNCOES_OMP_H
#define FUNCOES_OMP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

// Alocação em bloco único contíguo (Memory Friendly)
double* aloca_matriz(int m, int n) {
    return (double*)calloc(m * n, sizeof(double));
}

void libera_matriz(double* matriz) {
    free(matriz);
}

// Normalização Paralela (Totalmente independente por coluna)
void normaliza_matriz_omp(double* X, int m, int n) {
    double inv_m = 1.0 / (double)m; 
    
    #pragma omp parallel for schedule(static)
    for (int j = 0; j < n; j++) {
        int offset = j * m;
        double soma = 0.0;
        
        // Loop de soma
        for (int i = 0; i < m; i++) {
            soma += X[offset + i];
        }
        
        double media = soma * inv_m; 
        
        // Loop de subtração
        for (int i = 0; i < m; i++) {
            X[offset + i] -= media;
        }
    }
}

// O MOTOR V8: Multiplicação X^T * X com OpenMP
void calcula_covariancia_omp(double *X, int m, int n, double *C) {
    double inv_m_minus_1 = 1.0 / (double)(m - 1);
    
    // O schedule(dynamic) é vital aqui porque o loop interno diminui a cada passo (triângulo)
    #pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < n; j++) {
        int off_j = j * m;
        
        for (int i = j; i < n; i++) {
            int off_i = i * m;
            double soma = 0.0;
            
            // O paraíso do Prefetcher e do AVX2: Dois arrays perfeitamente contíguos
            for (int k = 0; k < m; k++) {
                soma += X[off_i + k] * X[off_j + k];
            }
            
            double valor = soma * inv_m_minus_1;
            
            // Escrita thread-safe na matriz C
            C[i * n + j] = valor;
            if (i != j) {
                C[j * n + i] = valor; // Espelha a simetria
            }
        }
    }
}

// Jacobi mantido serial porque roda inteiramente dentro do Cache L2/L3 (Sem gargalo de RAM)
void jacobi_covariancia_serial(double *C, int n, double *V, double *S) {
    int max_iter = 100;
    double tol = 1e-9;
    
    for (int i = 0; i < n * n; i++) V[i] = 0.0;
    for (int i = 0; i < n; i++) V[i * n + i] = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double max_off = 0.0;
        
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                double abs_val = fabs(C[i * n + j]);
                if (abs_val > max_off) max_off = abs_val;
            }
        }
        
        if (max_off < tol) break;

        for (int p = 0; p < n - 1; p++) {
            for (int q = p + 1; q < n; q++) {
                double Cpq = C[p * n + q];
                
                if (fabs(Cpq) > 1e-12) { 
                    double Cpp = C[p * n + p];
                    double Cqq = C[q * n + q];
                    
                    double theta = (Cqq - Cpp) / (2.0 * Cpq);
                    double t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
                    if (theta < 0.0) t = -t;
                    
                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = c * t;

                    for (int k = 0; k < n; k++) {
                        if (k != p && k != q) {
                            double Cpk = C[p * n + k];
                            double Cqk = C[q * n + k];
                            C[p * n + k] = c * Cpk - s * Cqk;
                            C[k * n + p] = C[p * n + k]; 
                            
                            C[q * n + k] = s * Cpk + c * Cqk;
                            C[k * n + q] = C[q * n + k]; 
                        }
                    }
                    
                    C[p * n + p] = c * c * Cpp - 2.0 * s * c * Cpq + s * s * Cqq;
                    C[q * n + q] = s * s * Cpp + 2.0 * s * c * Cpq + c * c * Cqq;
                    C[p * n + q] = 0.0;
                    C[q * n + p] = 0.0;

                    for (int k = 0; k < n; k++) {
                        double Vkp = V[k * n + p];
                        double Vkq = V[k * n + q];
                        V[k * n + p] = c * Vkp - s * Vkq;
                        V[k * n + q] = s * Vkp + c * Vkq;
                    }
                }
            }
        }
    }

    for (int i = 0; i < n; i++) {
        S[i] = C[i * n + i];
    }
}

// Ordenação rápida na memória
void ordena_componentes(double *S, double *V, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (S[i] < S[j]) {
                double tempS = S[i];
                S[i] = S[j];
                S[j] = tempS;

                for (int k = 0; k < n; k++) {
                    double tempV = V[k * n + i];
                    V[k * n + i] = V[k * n + j];
                    V[k * n + j] = tempV;
                }
            }
        }
    }
}

#endif // FUNCOES_OMP_H
