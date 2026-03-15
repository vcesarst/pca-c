#!/bin/bash

# ==========================================
# CONFIGURAÇÕES INICIAIS
# ==========================================
# Cria a pasta para organizar as saídas (se não existir, o -p resolve)
PASTA_SAIDA="saida_threads"
mkdir -p $PASTA_SAIDA

# Arquivo que vai guardar o resumo limpo dos tempos
ARQ_RESUMO="tempos_escalabilidade_omp.txt"

echo "=========================================" > $ARQ_RESUMO
echo "   BENCHMARK OPENMP - ESCALABILIDADE     " >> $ARQ_RESUMO
echo "=========================================" >> $ARQ_RESUMO


# ==========================================
# EXECUÇÃO DO BENCHMARK (1 a 8 THREADS)
# ==========================================
for i in {1..8}
do
    echo "-----------------------------------------" | tee -a $ARQ_RESUMO
    echo ">> Testando com $i thread(s)..." | tee -a $ARQ_RESUMO
    
    # Esta é a mágica: dita quantas threads o binário C vai poder acordar
    export OMP_NUM_THREADS=$i
    
    ARQ_SAIDA="$PASTA_SAIDA/saida_thread_${i}.txt"
    
    # Roda o programa e redireciona os prints da matriz para a pasta
    ./pca_omp.x > $ARQ_SAIDA
    
    # Usa o grep para "pescar" apenas a linha do tempo do SVD de dentro do arquivo de saída
    TEMPO_SVD=$(grep "SVD Finalizado em:" $ARQ_SAIDA)
    
    # Salva o tempo encontrado no arquivo de resumo e mostra na tela
    echo "   -> $TEMPO_SVD" | tee -a $ARQ_RESUMO
done

echo "=========================================" | tee -a $ARQ_RESUMO
echo " Benchmark finalizado com sucesso!       " | tee -a $ARQ_RESUMO
echo " Os prints completos estao em: $PASTA_SAIDA/ " | tee -a $ARQ_RESUMO
echo "=========================================" | tee -a $ARQ_RESUMO
