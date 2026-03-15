#ifndef FUNCOES_OPT_H
#define FUNCOES_OPT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double* aloca_matriz(int m, int n) {
    return (double*)calloc(m * n, sizeof(double));
}

void libera_matriz(double* matriz) {
    free(matriz);
}

void normaliza_matriz(double* X, int m, int n) {
    double inv_m = 1.0 / (double)m; 
    for (int j = 0; j < n; j++) {
        int offset = j * m;
        double soma = 0.0;
        for (int i = 0; i < m; i++) {
            soma += X[offset + i];
        }
        double media = soma * inv_m; 
        for (int i = 0; i < m; i++) {
            X[offset + i] -= media;
        }
    }
}

void svd_jacobi_unilateral(double *X, int m, int n, double *V, double *S) {
    int max_iter = 50;
    double tol = 1e-9;
    double tol_sq = tol * tol; 

    for (int i = 0; i < n * n; i++) V[i] = 0.0;
    for (int i = 0; i < n; i++) V[i * n + i] = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        double maior_convergencia = 0.0;

        for (int i = 0; i < n - 1; i++) {
            int off_i = i * m; 
            for (int j = i + 1; j < n; j++) {
                int off_j = j * m;
                
                double alpha = 0.0, beta = 0.0, gamma = 0.0;

                int k = 0;
                for (; k <= m - 4; k += 4) {
                    double xi0 = X[off_i + k],   xj0 = X[off_j + k];
                    double xi1 = X[off_i + k+1], xj1 = X[off_j + k+1];
                    double xi2 = X[off_i + k+2], xj2 = X[off_j + k+2];
                    double xi3 = X[off_i + k+3], xj3 = X[off_j + k+3];

                    alpha += xi0*xi0 + xi1*xi1 + xi2*xi2 + xi3*xi3;
                    beta  += xj0*xj0 + xj1*xj1 + xj2*xj2 + xj3*xj3;
                    gamma += xi0*xj0 + xi1*xj1 + xi2*xj2 + xi3*xj3;
                }
                for (; k < m; k++) {
                    double xi = X[off_i + k], xj = X[off_j + k];
                    alpha += xi * xi;
                    beta  += xj * xj;
                    gamma += xi * xj;
                }

                double convergencia_atual = (gamma * gamma) / (alpha * beta);
                if (convergencia_atual > maior_convergencia) maior_convergencia = convergencia_atual;

                if (convergencia_atual > tol_sq) {
                    double zeta = (beta - alpha) / (2.0 * gamma);
                    double t = (zeta > 0 ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = c * t;

                    int l = 0;
                    for (; l <= m - 4; l += 4) {
                        double xi0 = X[off_i + l], xj0 = X[off_j + l];
                        double xi1 = X[off_i + l+1], xj1 = X[off_j + l+1];
                        double xi2 = X[off_i + l+2], xj2 = X[off_j + l+2];
                        double xi3 = X[off_i + l+3], xj3 = X[off_j + l+3];

                        X[off_i+l]   = c*xi0 - s*xj0; X[off_j+l]   = s*xi0 + c*xj0;
                        X[off_i+l+1] = c*xi1 - s*xj1; X[off_j+l+1] = s*xi1 + c*xj1;
                        X[off_i+l+2] = c*xi2 - s*xj2; X[off_j+l+2] = s*xi2 + c*xj2;
                        X[off_i+l+3] = c*xi3 - s*xj3; X[off_j+l+3] = s*xi3 + c*xj3;
                    }
                    for (; l < m; l++) {
                        double xi = X[off_i + l], xj = X[off_j + l];
                        X[off_i + l] = c * xi - s * xj;
                        X[off_j + l] = s * xi + c * xj;
                    }

                    int v_off_i = i * n;
                    int v_off_j = j * n;
                    for (int l = 0; l < n; l++) {
                        double vi = V[v_off_i + l], vj = V[v_off_j + l];
                        V[v_off_i + l] = c * vi - s * vj;
                        V[v_off_j + l] = s * vi + c * vj;
                    }
                }
            }
        }
        if (maior_convergencia < tol_sq) break;
    }

    double inv_m_minus_1 = 1.0 / (double)(m - 1);
    for (int i = 0; i < n; i++) {
        double soma = 0.0;
        int offset = i * m;
        for (int k = 0; k < m; k++) soma += X[offset + k] * X[offset + k];
        S[i] = soma * inv_m_minus_1;
    }
}


// Função de ordenação para layout 1D
void ordena_componentes(double *S, double *V, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (S[i] < S[j]) {
                double tempS = S[i];
                S[i] = S[j];
                S[j] = tempS;

                int off_i = i * n;
                int off_j = j * n;
                for (int k = 0; k < n; k++) {
                    double tempV = V[off_i + k];
                    V[off_i + k] = V[off_j + k];
                    V[off_j + k] = tempV;
                }
            }
        }
    }
}

#endif // FUNCOES_OPT_H
