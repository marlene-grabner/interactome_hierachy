# %%
import interactome_hierachy as ih
import networkx as nx
import numpy as np
import gseapy as gp
import time
import os
import re

# %%
# =========================================================================
# Load Interactome
# =========================================================================

# Load interactome
ppi_path = "/Users/marlene/Documents/data/Projects/23_noise_in_networks_systematic/Code/NoiseInNetworks/paper_notebooks/comparison_network_creation/networks_2026_03_05/originals/chloe_ppi_ncbi/chloe_ppi_ncbi.tsv"
G = nx.read_edgelist(ppi_path, delimiter="\t")

# Extract only the LCC to avoid 0 eigenvalues
largest_cc = max(nx.connected_components(G), key=len)
G_core = G.subgraph(largest_cc)

# %%
# =========================================================================
# Check out the eigenvalue spectrum
# =========================================================================
# First 10 eigenvalues and eigenvectors
ih.spectral_clustering.compute_eigenvalues_and_eigenvectors(G_core, k=10)
# First 50 eigenvalues and eigenvectors
ih.spectral_clustering.compute_eigenvalues_and_eigenvectors(G_core, k=50)
# First 100 eigenvalues and eigenvectors
ih.spectral_clustering.compute_eigenvalues_and_eigenvectors(G_core, k=100)
# First 500 eigenvalues and eigenvectors
eigenvalues_500, eigenvectors_500 = (
    ih.spectral_clustering.compute_eigenvalues_and_eigenvectors(G_core, k=500)
)
# %%

# =========================================================================
# Check for the largest gaps in the eigenvalue spectrum
# =========================================================================

# Calculate the differences between the consecutive eigenvalues to identify the largest gaps
differences = np.diff(eigenvalues_500)
top_k_candidates = np.argsort(differences)[-30:] + 1
print(
    f"The top 30 candidates for k (based on largest gaps) are: {top_k_candidates[::-1]}"
)

# =========================================================================
# Cluster the graph using different values of k
# =========================================================================

# Sort eigenvectors according to order of the eigenvalues
sort_indices = np.argsort(eigenvalues_500)
sorted_eigenvectors = eigenvectors_500[:, sort_indices]

# Order of nodes for the Laplacian
node_order = list(G_core.nodes())

# %%
# Location to save the cluster assignments
output_folder = "/Users/marlene/Documents/data/Projects/29_Interactome_Hierachy/data/processed/hierachical_clustering/cluster_assignments"
# %%
# Hierachical clustering pipleine
hierachical_clusters = {}
for k in range(1, 25):
    clusters = ih.spectral_clustering.calculate_k_cluster_from_spectrum(
        eigenvectors=sorted_eigenvectors,
        k=k,
        node_list=node_order,
        output_folder=output_folder,
    )
    hierachical_clusters[k] = clusters
# %%
# =========================================================================
# Analyse the enrichment of these clusters
# =========================================================================

# Get a dictionary mapping NCBI IDs to HGNC symbols for all nodes in the graph
# This can then be used for enrichment analysis of the clusters
all_ncbi_ids = list(G_core.nodes())
ncbi_to_hgnc = ih.utils.dict_ncbi_to_hgnc(all_ncbi_ids)


# Fetch the cluster assignemnts to a dictionary
cluster_assignemnt_files = os.listdir(output_folder)
clusters_per_level = {}
for file in cluster_assignemnt_files:
    pattern = r"clusters_k(\d+)\.txt"
    match = re.match(pattern, file)
    cluster = ih.utils.clusters_to_dict(os.path.join(output_folder, file))
    clusters_per_level[int(match.group(1))] = cluster

# Use enrichr to calculate the enrichment for each cluster and each level of the hierachy
gene_set_libraries = ["GO_Biological_Process_2023", "KEGG_2021_Human"]

# Calculate the enrichment for each cluster and each level of the hierachy
cluster_per_level_enrichment_dict = {}
for hierachical_level, clusters in clusters_per_level.items():
    level_enrichment_dict = {}
    for cluster_id, nodes in clusters.items():
        print(
            f"Calculating enrichment for cluster {cluster_id} in level {hierachical_level} with {len(nodes)} nodes..."
        )
        symbol_list = [ncbi_to_hgnc[ncbi] for ncbi in nodes if ncbi in ncbi_to_hgnc]
        enr = gp.enrichr(
            gene_list=symbol_list,
            gene_sets=gene_set_libraries,
            organism="human",
            outdir=None,
        )
        time.sleep(3)
        level_enrichment_dict[cluster_id] = enr.results
    cluster_per_level_enrichment_dict[hierachical_level] = level_enrichment_dict

# Save a summary of the enrichment for each level of the hierachy
output_folder = "/Users/marlene/Documents/data/Projects/29_Interactome_Hierachy/data/processed/hierachical_clustering/cluster_enrichment"

for level, enrichment_dict in cluster_per_level_enrichment_dict.items():
    ih.spectral_clustering.summarize_level_enrichment(
        cluster_enrichment_dict=enrichment_dict,
        level_name=f"k_{level}",
        output_folder=output_folder,
    )

# %%
