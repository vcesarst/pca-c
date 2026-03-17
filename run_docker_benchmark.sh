#!/bin/bash

# Para o script se algum comando falhar
set -e

echo "=========================================================="
echo "  PREPARANDO AMBIENTE EMBARCADO (DOCKER)                "
echo "=========================================================="
echo "Construindo a imagem Alpine com o código C Otimizado..."
docker build -t pca-embarcado .
echo "[OK] Imagem 'pca-embarcado' construída com sucesso!"
echo ""

echo "=========================================================="
echo "  CENÁRIO 1: PLACA PODEROSA (Sem limites de hardware)   "
echo "=========================================================="
# O --rm faz o container se autodestruir após o uso para não sujar o PC
docker run --rm pca-embarcado
echo ""

echo "=========================================================="
echo "  CENÁRIO 2: IOT PADRÃO (Restrito a 50MB de RAM)        "
echo "=========================================================="
docker run --rm --memory="50m" pca-embarcado
echo ""

echo "=========================================================="
echo "  CENÁRIO 3: IOT EXTREMO (Apenas 15MB RAM e 0.5 CPU)    "
echo "=========================================================="
# Restringe a meia CPU lógica e quase nada de RAM
docker run --rm --memory="15m" --cpus="0.5" pca-embarcado
echo ""

echo "=========================================================="
echo " BENCHMARKS FINALIZADOS                                "
echo " O PCA sobreviveu a ambientes de altíssima restrição!     "
echo "=========================================================="
