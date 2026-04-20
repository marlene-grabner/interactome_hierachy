import networkx as nx
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.preprocessing import normalize


def compute_eigenvalues_and_eigenvectors(
    G: nx.Graph, k: int = 50
) -> (np.ndarray, np.ndarray):
    L = nx.normalized_laplacian_matrix(G)
    # Calculate the k-smalles eigenvalues using SciPy's sparse solver
    # 'which=SM' means the smallest magnitude eigenvalues
    eigenvalues, eigenvectors = eigsh(L, k=k, which="SM")
    # Sort them
    eigenvalues = np.sort(eigenvalues)

    # Plot them
    plt.figure(figsize=(10, 6))
    plt.plot(range(1, k + 1), eigenvalues, marker="o", linestyle="-", color="b")
    plt.title(f"Spectral Gap Analysis of the first {k} Eigenvalues")
    plt.xlabel("Eigenvalue Index (k)")
    plt.ylabel("Eigenvalue Magnitude")
    plt.axhline(0, color="black", linewidth=0.5)
    plt.grid(True, linestyle="--", alpha=0.7)
    plt.show()

    return eigenvalues, eigenvectors


def calculate_k_cluster_from_spectrum(
    eigenvectors: np.ndarray,  # The matrix of eigenvectors (n x k)
    k: int,  # The number of clusters to form and dimensions to use
    node_list: list,  # The list of node names corresponding to the rows of the eigenvector matrix
    output_folder: str = None,  #
) -> dict:
    """
    Calculates k cluster assignments for a given k using the first k eigenvectors
    from the spectral decomposition of the graph Laplacian.
    """

    # Extract the first k eigenvectors
    U_k = eigenvectors[:, :k]

    # Row normalize the rows of U_k
    # This projects the coordinates onto a unit sphere, making K-means much more effective
    U_k_normalized = normalize(U_k, axis=1, norm="l2")

    # Run k-means clustering
    # n_init=10 runs K-means 10 times with different seeds to find the best fit
    kmeans = KMeans(n_clusters=k, n_init=10, random_state=42)
    labels = kmeans.fit_predict(U_k_normalized)

    # Map the labels back to the original node names
    clusters = {i: [] for i in range(k)}
    for node_idx, label in enumerate(labels):
        protein_name = node_list[node_idx]
        clusters[label].append(protein_name)

    # Save the cluster assignments to a file
    if output_folder is not None:
        filename = os.path.join(output_folder, f"clusters_k{k}.txt")
        with open(filename, "w") as f:
            for cluster_id, proteins in clusters.items():
                f.write(f">Cluster_{cluster_id} | Size: {len(proteins)}\n")
                f.write("\n".join(proteins) + "\n\n")

    return clusters
