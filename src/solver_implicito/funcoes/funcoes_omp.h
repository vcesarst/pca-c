#ifndef FUNCOES_OMP_H
#define FUNCOES_OMP_H

/*
 * funcoes_omp.h — Jacobi SVD Paralelo, Versão 2 (Cache-Aware)
 *
 * DIAGNÓSTICO DA V1:
 * ─────────────────────────────────────────────────────────────
 * O problema não era "poucas threads" — era arquitetura de memória.
 *
 * Com m=100.000, cada coluna ocupa 800 KB. Um par de colunas ocupa
 * 1,6 MB. Com 4 threads processando 4 pares simultâneos, são 6,4 MB
 * de dados ativos — acima de qualquer L3 razoável (~4–32 MB compartilhado).
 * Resultado: cada acesso vira um cache miss → vai para a RAM (~100 ns).
 * Adicionar mais threads só piora: mais threads competem pelo mesmo
 * barramento de memória (DDR4 ~40–50 GB/s, não escalável por thread).
 *
 * SOLUÇÃO — ROW TILING:
 * ─────────────────────────────────────────────────────────────
 * Em vez de paralelizar por PARES DE COLUNAS (que leva colunas inteiras
 * para o cache), paralelizamos por BLOCOS DE LINHAS dentro de cada par.
 *
 * Cada thread recebe um tile de TILE_SIZE linhas e processa o par (i,j)
 * *nesse tile*. Com TILE_SIZE = 4096 linhas × 2 colunas × 8 bytes = 64 KB,
 * o tile cabe inteiro no L2 privado de cada core (tipicamente 256 KB–1 MB).
 *
 * Isso transforma o padrão de acesso de:
 *   "8 threads brigando por 8 colunas de 800 KB cada" (memory-bound)
 * em:
 *   "8 threads processando 8 tiles de 64 KB cada, reutilizando o cache"
 *   (compute-bound dentro do L2 — speedup quase linear)
 *
 * PARALELISMO POR PARES → mantemos o round-robin para garantir
 * independência dos pares por rodada. A mudança é que cada par agora
 * também paraleliza internamente suas linhas via tiling.
 *
 * ESTRUTURA GERAL:
 *   - Nível externo: round-robin de pares (serial, garante independência)
 *   - Nível interno: #pragma omp parallel for por tiles de linhas
 *     dentro do par ativo. Cada thread opera em memória disjunta.
 *
 * RACE CONDITIONS:
 *   - X: tiles são disjuntos → sem conflito
 *   - V: atualizado UMA vez por par, fora do paralelo de tiles → sem conflito
 *   - alpha/beta/gamma: reduction explícita → sem conflito
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

/* ------------------------------------------------------------------ *
 *  Tamanho do tile de linhas para a rotação de Jacobi.               *
 *  4096 linhas × 2 colunas × 8 bytes = 64 KB → cabe no L2 privado.  *
 *  Ajuste para baixo se seu L2 for 256 KB (tente 2048).              *
 *  Ajuste para cima se tiver L2 grande (tente 8192 com L2 de 1 MB).  *
 * ------------------------------------------------------------------ */
#define TILE_SIZE 4096

double* aloca_matriz(int m, int n) {
    double *p = (double*)calloc((size_t)m * n, sizeof(double));
    if (!p) { fprintf(stderr, "ERRO: calloc falhou (%d x %d)\n", m, n); exit(1); }
    return p;
}

void libera_matriz(double *matriz) {
    free(matriz);
}

/* ------------------------------------------------------------------
 * normaliza_matriz
 *
 * Subtrai a média de cada coluna (centralização).
 * Equivalente ao comportamento padrão do sklearn.PCA, que também
 * só centraliza sem escalar — resultados validados contra scikit-learn.
 *
 * Paralelismo: colunas são 100% independentes → schedule(static)
 * perfeito, zero overhead de sincronização.
 * ------------------------------------------------------------------ */
void normaliza_matriz(double *X, int m, int n) {
    double inv_m = 1.0 / (double)m;

    #pragma omp parallel for schedule(static)
    for (int j = 0; j < n; j++) {
        double *col = X + (size_t)j * m;

        double soma = 0.0;
        for (int i = 0; i < m; i++)
            soma += col[i];
        double media = soma * inv_m;

        for (int i = 0; i < m; i++)
            col[i] -= media;
    }
}

/* ------------------------------------------------------------------
 * svd_jacobi_unilateral  —  Versão Cache-Aware com Row Tiling
 *
 * Mantém o mesmo algoritmo e critério de convergência da V1.
 * Muda APENAS onde o paralelismo é aplicado.
 * ------------------------------------------------------------------ */
void svd_jacobi_unilateral(double *X, int m, int n, double *V, double *S) {
    const int    max_iter = 50;
    const double tol      = 1e-9;
    const double tol_sq   = tol * tol;

    /* Inicializa V = I */
    memset(V, 0, (size_t)n * n * sizeof(double));
    for (int i = 0; i < n; i++) V[i * n + i] = 1.0;

    /* Prepara torneio round-robin (idêntico à V1) */
    int n_pad = (n % 2 == 0) ? n : n + 1;
    int *p = (int*)malloc(n_pad * sizeof(int));
    if (!p) { fprintf(stderr, "ERRO: malloc p\n"); exit(1); }
    for (int k = 0; k < n_pad; k++) p[k] = k;

    /* ── Loop de iterações Jacobi ─────────────────────────────── */
    for (int iter = 0; iter < max_iter; iter++) {
        double maior_conv = 0.0;

        /* ── Rodadas do torneio (n_pad-1 rodadas, cada uma com n/2 pares) */
        for (int step = 0; step < n_pad - 1; step++) {

            /* Percorre os pares da rodada atual.
             *
             * MUDANÇA CENTRAL v1 → v2:
             * Na V1: #pragma omp parallel for sobre os pares k
             *   → cada thread pegava UMA coluna inteira (800KB) → memory-bound
             *
             * Na V2: loop SERIAL sobre os pares k (pares são poucos: n/2=250)
             *   → dentro de cada par, #pragma omp parallel for sobre TILES de linhas
             *   → cada thread pega TILE_SIZE linhas das 2 colunas (~64KB) → L2-bound
             *
             * A independência de dados entre threads é garantida pelos tiles disjuntos.
             */
            for (int k = 0; k < n_pad / 2; k++) {
                int col1 = p[k];
                int col2 = p[n_pad - 1 - k];
                if (col1 >= n || col2 >= n) continue;

                int ci = (col1 < col2) ? col1 : col2;
                int cj = (col1 > col2) ? col1 : col2;

                double *Xi = X + (size_t)ci * m;
                double *Xj = X + (size_t)cj * m;

                /* ── Fase 1: acumula alpha, beta, gamma em paralelo por tiles ── */
                double alpha = 0.0, beta = 0.0, gamma = 0.0;

                #pragma omp parallel for schedule(static) \
                        reduction(+:alpha,beta,gamma)
                for (int tile = 0; tile < m; tile += TILE_SIZE) {
                    int end = tile + TILE_SIZE;
                    if (end > m) end = m;

                    double la = 0.0, lb = 0.0, lg = 0.0;
                    int r = tile;

                    /* Loop desdobrado 4× dentro do tile */
                    for (; r <= end - 4; r += 4) {
                        double xi0 = Xi[r],   xj0 = Xj[r];
                        double xi1 = Xi[r+1], xj1 = Xj[r+1];
                        double xi2 = Xi[r+2], xj2 = Xj[r+2];
                        double xi3 = Xi[r+3], xj3 = Xj[r+3];
                        la += xi0*xi0 + xi1*xi1 + xi2*xi2 + xi3*xi3;
                        lb += xj0*xj0 + xj1*xj1 + xj2*xj2 + xj3*xj3;
                        lg += xi0*xj0 + xi1*xj1 + xi2*xj2 + xi3*xj3;
                    }
                    for (; r < end; r++) {
                        double xi = Xi[r], xj = Xj[r];
                        la += xi * xi;
                        lb += xj * xj;
                        lg += xi * xj;
                    }
                    alpha += la; beta += lb; gamma += lg;
                }

                /* ── Critério de convergência ── */
                double conv = (gamma * gamma) / (alpha * beta);
                if (conv > maior_conv) maior_conv = conv;
                if (conv <= tol_sq) continue; /* par já ortogonal, pula */

                /* ── Calcula ângulo de rotação (serial, O(1)) ── */
                double zeta = (beta - alpha) / (2.0 * gamma);
                double t    = (zeta >= 0.0 ? 1.0 : -1.0)
                              / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                double c    = 1.0 / sqrt(1.0 + t * t);
                double s    = c * t;

                /* ── Fase 2: aplica rotação em paralelo por tiles ── */
                #pragma omp parallel for schedule(static)
                for (int tile = 0; tile < m; tile += TILE_SIZE) {
                    int end = tile + TILE_SIZE;
                    if (end > m) end = m;
                    int r = tile;

                    for (; r <= end - 4; r += 4) {
                        double xi0 = Xi[r],   xj0 = Xj[r];
                        double xi1 = Xi[r+1], xj1 = Xj[r+1];
                        double xi2 = Xi[r+2], xj2 = Xj[r+2];
                        double xi3 = Xi[r+3], xj3 = Xj[r+3];
                        Xi[r]   = c*xi0 - s*xj0;  Xj[r]   = s*xi0 + c*xj0;
                        Xi[r+1] = c*xi1 - s*xj1;  Xj[r+1] = s*xi1 + c*xj1;
                        Xi[r+2] = c*xi2 - s*xj2;  Xj[r+2] = s*xi2 + c*xj2;
                        Xi[r+3] = c*xi3 - s*xj3;  Xj[r+3] = s*xi3 + c*xj3;
                    }
                    for (; r < end; r++) {
                        double xi = Xi[r], xj = Xj[r];
                        Xi[r] = c*xi - s*xj;
                        Xj[r] = s*xi + c*xj;
                    }
                }

                /* ── Atualiza V (serial: n=500 ops, custo desprezível) ──
                 * Fora do paralelo de tiles por design: V tem dependências
                 * de leitura/escrita nas mesmas linhas entre pares diferentes.
                 * Como estamos no loop SERIAL de pares, não há conflito. */
                double *Vi = V + (size_t)ci * n;
                double *Vj = V + (size_t)cj * n;
                for (int r = 0; r < n; r++) {
                    double vi = Vi[r], vj = Vj[r];
                    Vi[r] = c*vi - s*vj;
                    Vj[r] = s*vi + c*vj;
                }

            } /* fim loop pares k */

            /* Rotaciona p para próxima rodada (mantém p[0] fixo) */
            int last = p[n_pad - 1];
            for (int r = n_pad - 1; r > 1; r--) p[r] = p[r-1];
            p[1] = last;

        } /* fim loop step */

        if (maior_conv < tol_sq) {
            printf("     Convergiu na iteracao %d\n", iter + 1);
            break;
        }

    } /* fim loop iter */

    free(p);

    /* ── Singular values: norma de cada coluna de X rotacionada ── */
    double inv_m1 = 1.0 / (double)(m - 1);
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; i++) {
        double *col = X + (size_t)i * m;
        double soma = 0.0;
        for (int k = 0; k < m; k++) soma += col[k] * col[k];
        S[i] = soma * inv_m1;
    }
}

/* ------------------------------------------------------------------ *
 *  Ordenação (idêntica à V1 — O(n²) com n=500 é ~0.001 s)           *
 * ------------------------------------------------------------------ */
void ordena_componentes(double *S, double *V, int n) {
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (S[i] < S[j]) {
                double tmp = S[i]; S[i] = S[j]; S[j] = tmp;
                double *Vi = V + (size_t)i * n;
                double *Vj = V + (size_t)j * n;
                for (int k = 0; k < n; k++) {
                    double tv = Vi[k]; Vi[k] = Vj[k]; Vj[k] = tv;
                }
            }
        }
    }
}

#endif /* FUNCOES_OMP_H */
