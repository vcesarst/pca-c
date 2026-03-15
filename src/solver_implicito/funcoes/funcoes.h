#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double media_vetor(double * v, int n){
    double media = 0.0;
    for (int i = 0; i < n; i++)
    {
        media += v[i];
    }
    return media / n;
}

double* inicializa_vetor(int n){
    double *v = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
    {
        v[i] = 0.0;
    }
    return v;
}

void libera_vetor(double *v){
    free(v);
}

void imprime_vetor(double *v, int n){
    for (int i = 0; i < n; i++)
    {
        printf("%f ", v[i]);
    }
    printf("\n");
}

double ** aloca_matriz(int m, int n)
{
    double **matriz = (double **)malloc(m * sizeof(double *));
    for (int i = 0; i < m; i++)
    {
        matriz[i] = (double *)malloc(n * sizeof(double));
    }
    return matriz;
}

void libera_matriz(double **matriz, int m)
{
    for (int i = 0; i < m; i++)
    {
        free(matriz[i]);
    }
    free(matriz);
}

void imprime_matriz(double **matriz, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            printf("%f ", matriz[i][j]);
        }
        printf("\n");
    }
}

double ** normaliza_matriz(double **matriz, int m, int n)
{
    // Aloca a matriz m x n (linhas x colunas)
    double **matriz_normalizada = aloca_matriz(m, n);

    // Navega coluna por coluna (j vai de 0 até n)
    for (int j = 0; j < n; j++)
    {
        // 1. Calcula a média da coluna 'j' atual
        double soma = 0.0;
        for (int i = 0; i < m; i++)
        {
            soma += matriz[i][j]; // Somando as linhas daquela coluna
        }
        double media = soma / m; // Divide pelo número de linhas (amostras)

        // 2. Subtrai a média de cada elemento da coluna 'j'
        for (int i = 0; i < m; i++)
        {
            matriz_normalizada[i][j] = matriz[i][j] - media;
        }
    }
    
    return matriz_normalizada;
}

double ** transposta_matriz(double **matriz, int m, int n)
{
    double **matriz_transposta = aloca_matriz(n, m);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            matriz_transposta[j][i] = matriz[i][j];
        }
    }
    return matriz_transposta;
}

double ** multiplicacao_matrizes(double **A, double **B, int m, int n, int p)
{
    double **C = aloca_matriz(m, p);
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < p; j++)
        {
            C[i][j] = 0.0;
            for (int k = 0; k < n; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    return C;
}

double sign(double x) {
    return (x >= 0.0) ? 1.0 : -1.0;
}

void svd_jacobi_unilateral(double **X, int m, int n, double **V, double *autovalores) {
    double tol = 1e-9;   // Tolerância para considerar as colunas ortogonais
    int max_sweeps = 100; // Limite de segurança de iterações
    int convergiu;

    // 1. Inicializa a matriz V como uma Matriz Identidade
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) V[i][j] = 1.0;
            else V[i][j] = 0.0;
        }
    }

    // 2. O Loop Principal de Varreduras (Sweeps)
    for (int sweep = 0; sweep < max_sweeps; sweep++) {
        convergiu = 1; // Assume que convergiu, a menos que precisemos rotacionar algo

        // Combina todos os pares de colunas (i, j)
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                
                // a) Calcula os produtos escalares
                double alpha = 0.0; // Norma ao quadrado da coluna i
                double beta  = 0.0; // Norma ao quadrado da coluna j
                double gamma = 0.0; // Produto escalar entre coluna i e j
                
                for (int k = 0; k < m; k++) {
                    alpha += X[k][i] * X[k][i];
                    beta  += X[k][j] * X[k][j];
                    gamma += X[k][i] * X[k][j];
                }

                // b) Verifica se as colunas já são ortogonais
                // Multiplicamos a tolerância pelas normas para evitar problemas de escala
                if (fabs(gamma) > tol * sqrt(alpha * beta)) {
                    convergiu = 0; // Se entrou aqui, o algoritmo ainda precisa rodar

                    // c) Calcula o ângulo de rotação de Jacobi (cosseno e seno)
                    double zeta = (beta - alpha) / (2.0 * gamma);
                    double t = sign(zeta) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
                    double c = 1.0 / sqrt(1.0 + t * t);
                    double s = c * t;

                    // d) Aplica a rotação simultaneamente nas colunas i e j das matrizes X e V
                    for (int k = 0; k < m; k++) {
                        double x_ki = X[k][i];
                        double x_kj = X[k][j];
                        X[k][i] = c * x_ki - s * x_kj;
                        X[k][j] = s * x_ki + c * x_kj;
                    }
                    
                    for (int k = 0; k < n; k++) {
                        double v_ki = V[k][i];
                        double v_kj = V[k][j];
                        V[k][i] = c * v_ki - s * v_kj;
                        V[k][j] = s * v_ki + c * v_kj;
                    }
                }
            }
        }
        
        // Se percorreu todas as colunas e o `convergiu` continuou 1, terminamos!
        if (convergiu) {
            break; 
        }
    }

    // 3. Extração final dos Autovalores
    // A variância de cada componente é a norma ao quadrado da coluna de X resultante, dividida por (m - 1)
    for (int i = 0; i < n; i++) {
        double norma_quadrada = 0.0;
        for (int k = 0; k < m; k++) {
            norma_quadrada += X[k][i] * X[k][i];
        }
        autovalores[i] = norma_quadrada / (m - 1.0);
    }
}

/*
  Função para ordenar os Autovalores (decrescente) e seus Autovetores correspondentes
  Parâmetros:
  - autovalores: Vetor (tamanho n)
  - V: Matriz de autovetores (n x n)
  - n: Número de variáveis (colunas)
*/
void ordena_componentes(double *autovalores, double **V, int n) {
    for (int i = 0; i < n - 1; i++) {
        // Encontra o índice do maior autovalor restante
        int max_idx = i;
        for (int j = i + 1; j < n; j++) {
            if (autovalores[j] > autovalores[max_idx]) {
                max_idx = j;
            }
        }

        // Se o maior não estiver na posição atual (i), fazemos a troca (Swap)
        if (max_idx != i) {
            // 1. Troca os autovalores de lugar
            double temp_val = autovalores[i];
            autovalores[i] = autovalores[max_idx];
            autovalores[max_idx] = temp_val;

            // 2. Troca as colunas inteiras na matriz de Autovetores V
            for (int k = 0; k < n; k++) {
                double temp_vec = V[k][i];
                V[k][i] = V[k][max_idx];
                V[k][max_idx] = temp_vec;
            }
        }
    }
}
