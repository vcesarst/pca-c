import numpy as np
import time
from sklearn.decomposition import PCA

m = 100000
n = 500

print("==========================================")
print("   PCA BENCHMARK: PYTHON / SCIKIT-LEARN   ")
print("==========================================")
print(f"Dimensoes: {m} linhas x {n} colunas")

# --- 1 & 2. ALOCAÇÃO E GERAÇÃO DE DADOS ---
start = time.time()
# Cria vetores de índices e usa broadcasting para simular (i + j) % 100
i_idx = np.arange(m).reshape(-1, 1)
j_idx = np.arange(n).reshape(1, -1)
X = ((i_idx + j_idx) % 100).astype(np.float64)
end = time.time()
print(f"[1 & 2] Alocacao e Geracao: {end - start:.4f} s")

# --- 4 & 5. PCA (Scikit-Learn faz a normalização internamente) ---
print("[4] Executando PCA (LAPACK)... (Aguarde)")
start = time.time()

# svd_solver='full' chama a biblioteca Fortran otimizada (Divide & Conquer)
pca = PCA(svd_solver='full')
pca.fit(X)

end = time.time()
tempo_pca = end - start
print(f"[4 & 5] SVD + Ordenacao:    {tempo_pca:.4f} s")

# --- RESULTADOS PARA VALIDAÇÃO ---
print("\n=== VALIDACAO DOS RESULTADOS ===")
print("Primeiros 5 Autovalores (Variancia Explicada):")
for idx, val in enumerate(pca.explained_variance_[:5]):
    print(f"PC{idx+1}: {val:.4f}")
print("==========================================")
