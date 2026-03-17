import numpy as np
from sklearn.datasets import load_breast_cancer
from sklearn.decomposition import PCA

def main():
    print("=== Extração de Dados Reais para Teste HPC ===")
    
    # 1. Carregar o dataset real
    dataset = load_breast_cancer()
    X = dataset.data
    m, n = X.shape
    
    print(f"Dataset carregado: {m} amostras (linhas) x {n} features (colunas).")

    # 2. Exportar a Matriz X para o seu código C
    # Colocando as dimensões na primeira linha para facilitar o malloc no C
    caminho_dados = 'breast_cancer_data.txt'
    with open(caminho_dados, 'w') as f:
        f.write(f"{m} {n}\n") # Cabeçalho
        np.savetxt(f, X, fmt='%.8f', delimiter=' ')
    print(f"[OK] Dados exportados para '{caminho_dados}'.")

    # 3. Rodar o "Gabarito" (Scikit-Learn)
    # O Scikit-Learn subtrai a média internamente
    pca_sklearn = PCA()
    pca_sklearn.fit(X)

    # 4. Exportar os Autovetores do Scikit-Learn
    # pca.components_ no Python retorna as linhas como componentes.
    # Transpondo para que as colunas sejam os autovetores
    caminho_gabarito = 'sklearn_autovetores.txt'
    with open(caminho_gabarito, 'w') as f:
        np.savetxt(f, pca_sklearn.components_.T, fmt='%.8f', delimiter=' ')
    print(f"[OK] Gabarito do Scikit-Learn exportado para '{caminho_gabarito}'.")
    

if __name__ == "__main__":
    main()
