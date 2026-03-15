#!/bin/bash

# Configurações
LIMITE_TEMPO="1.5h"
ARQ_TEMPOS="tempos_execucao.txt"

# Limpa o arquivo de tempos de execuções anteriores
echo "=========================================" > $ARQ_TEMPOS
echo "         RELATÓRIO DE TEMPOS (CPU)       " >> $ARQ_TEMPOS
echo "=========================================" >> $ARQ_TEMPOS

# Compilação
echo "Limpando e compilando executáveis..."
make clean > /dev/null
make all > /dev/null

if [ $? -ne 0 ]; then
    echo "Erro na compilação! Verifique o Makefile e as flags."
    exit 1
fi

# Função de Teste
rodar_teste() {
    EXEC=$1
    ARQ_SAIDA="saida_${EXEC}.txt"
    
    echo "-----------------------------------------" >> $ARQ_TEMPOS
    echo "Iniciando teste: $EXEC"
    
    # O comando time pega o %U (Tempo User/CPU) e anexa (-a -o) no arquivo de tempos.
    # O timeout corta a execução em 1h30m. 
    # O > $ARQ_SAIDA joga os prints da matriz no arquivo de saída.
    /usr/bin/time -f "-> Tempo exato de CPU: %U segundos" -a -o $ARQ_TEMPOS timeout $LIMITE_TEMPO ./$EXEC > $ARQ_SAIDA 2>&1
    
    STATUS=$?
    
    if [ $STATUS -eq 124 ]; then
        echo "-> RESULTADO: TIMEOUT (Atingiu o teto de 1h30m)" >> $ARQ_TEMPOS
        echo "Processo $EXEC cortado por timeout."
    elif [ $STATUS -eq 0 ]; then
        echo "-> RESULTADO: SUCESSO" >> $ARQ_TEMPOS
        echo "Processo $EXEC finalizado com sucesso."
    else
        echo "-> RESULTADO: ERRO (Código: $STATUS)" >> $ARQ_TEMPOS
        echo "Processo $EXEC falhou."
    fi
}

# Executa os dois programas
rodar_teste "serial.x"
rodar_teste "serial_opt.x"

echo "=========================================" >> $ARQ_TEMPOS
echo "Benchmark concluído! Verifique os arquivos gerados."
