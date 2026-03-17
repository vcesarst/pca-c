# 📊 Benchmarks Reais e Testes para Sistemas Embarcados (IoT)

Este diretório contém a suíte de testes focada em duas missões críticas: 
1. **Validação Matemática:** Comprovar que a implementação em C gera resultados idênticos à principal biblioteca de mercado (Scikit-Learn).
2. **Stress-Test de Hardware:** Provar que a biblioteca é leve o suficiente para rodar algoritmos de Machine Learning em dispositivos de borda (*Edge Computing* / IoT) com severas restrições de memória RAM e CPU.

## 1. Validação Numérica contra o Scikit-Learn (Python)

Para abandonar os números aleatórios e testar com dados reais, utilizamos o clássico **Breast Cancer Wisconsin Dataset** (Matriz de 569 amostras x 30 características).

* O script `breast_cancer_dataset.py`, na pasta examplesx faz o download do dataset real, roda a implementação padrão do `sklearn.decomposition.PCA` e exporta o "gabarito" dos autovetores (`sklearn_autovetores.txt`).
* **Resultado da Validação:** Ao cruzar os autovetores gerados pelo código em C com os do Python, o desvio máximo encontrado foi na casa de $8.0 \times 10^{-8}$. Isso garante **exatidão matemática absoluta** no limite da precisão de ponto flutuante do processador, provando que a lógica de Covariância Explícita + SVD de Jacobi está perfeita.

## 2. Execução Local Rápida (`run_embarcado.sh`)

Para compilar e testar o código C nativamente no seu sistema operacional, está disponibilizado  o script `run_embarcado.sh`. 

Ele automatiza a compilação do arquivo `main_embarcado.c` aplicando as flags de otimização extrema do GCC (`-O3`, `-ffast-math`, `-funroll-loops`, `-march=native` e `-fopenmp`) e executa o binário na sequência, lendo a matriz de testes e cuspindo o resultado em `pca_c_autovetores.txt`.

**Como rodar:**
```bash
./run_embarcado.sh

## 3. Teste de Performance em Dispositivos de Borda (IoT)

O script 'run_docker_benchmark.sh' é a chave para validar a performance da biblioteca em um ambiente controlado de IoT. Ele utiliza uma imagem Docker leve baseada em Alpine Linux, que simula as restrições de hardware típicas de dispositivos de borda.

