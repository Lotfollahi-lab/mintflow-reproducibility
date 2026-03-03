# Addressing R3C1 and R3C7 - Kindey (graphST V2) - TLS core 1

import os
import random
import tempfile

import numpy as np  # type: ignore
import scanpy as sc  # type: ignore
import scvi  # type: ignore
import pandas as pd

import squidpy as sq
import anndata as ad

import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import entropy

initDir = '/lustre/scratch126/cellgen/lotfollahi/dv8/mint_flow_bench/GraphST/results/xenium_old_kidney/'
save_dir = '/nfs/team361/ms83/data/kidney/xenium_graphst/analysis/'

adata = sc.read_h5ad(initDir + 'xenium_old_kidney_graphst.h5ad')
adata

adata = adata[adata.obs['level_2_cell_type'] == 'T lymphocyte'].copy()
adata

# Load the annotations
spatial_niche_df = pd.read_csv("/nfs/users/nfs_m/ms83/mintflow-revision/benchmarking/spatial_kidney_labels_annotations.csv")

# Create mapping from cell_id to INIGEN spatial labels
spatial_niche_map = dict(zip(spatial_niche_df['cell_id'], spatial_niche_df['INIGEN spatial labels']))

# Apply to other AnnData using cell_id column
adata.obs['INIGEN spatial labels'] = adata.obs['cell_id'].map(spatial_niche_map)

# Check how many matched
matched = adata.obs['INIGEN spatial labels'].notna().sum()
print(f"Successfully mapped {matched} out of {adata.n_obs} cells")

COLORS = {
 'CD4+ kidney resident': '#3B76AF',
 'kidney-tumour shared': '#EF8636',
 'CD8+ tumour infiltrating': '#529E3F',   
 'CD8+ TLS core 1': '#C53A32',
 'CD8+ TLS core 2': '#9169B8',
 'CD8+ TLS border': '#85594E',
 }
colors=COLORS

adata.obs['INIGEN spatial labels']=adata.obs['INIGEN spatial labels'].astype('category')
adata.uns['spatial_labels_colors'] = [COLORS[cat] for cat in adata.obs['INIGEN spatial labels'].cat.categories]

def calculate_metrics_for_niche(adata, niche_name, cluster_column):
    """
    Calculate Recall, Precision, F1, and Split Entropy for a specific niche
    """
    # Get true niche labels
    true_niche = (adata.obs['INIGEN spatial labels'] == niche_name).values
    n_niche_cells = true_niche.sum()
    
    # Get cluster assignments
    clusters = adata.obs[cluster_column].values
    unique_clusters = np.unique(clusters)
    
    # Find the cluster with maximum overlap with the niche
    best_cluster = None
    max_overlap = 0
    
    cluster_overlaps = {}
    for cluster in unique_clusters:
        cluster_mask = (clusters == cluster)
        overlap = (true_niche & cluster_mask).sum()
        cluster_overlaps[cluster] = overlap
        if overlap > max_overlap:
            max_overlap = overlap
            best_cluster = cluster
    
    # Calculate metrics using the best matching cluster
    pred_cluster = (clusters == best_cluster)
    
    # True Positives, False Positives, False Negatives
    TP = (true_niche & pred_cluster).sum()
    FP = (~true_niche & pred_cluster).sum()
    FN = (true_niche & ~pred_cluster).sum()
    
    # Recall: what proportion of true niche cells are captured
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    
    # Precision: what proportion of the cluster is actually the niche
    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    
    # F1 score
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0
    
    # Split Entropy: how fragmented is the niche across clusters?
    # Calculate the distribution of niche cells across all clusters
    niche_distribution = []
    for cluster in unique_clusters:
        cluster_mask = (clusters == cluster)
        n_niche_in_cluster = (true_niche & cluster_mask).sum()
        if n_niche_in_cluster > 0:
            niche_distribution.append(n_niche_in_cluster)
    
    if len(niche_distribution) > 0:
        niche_distribution = np.array(niche_distribution) / n_niche_cells
        split_entropy = entropy(niche_distribution, base=2)  # bits
    else:
        split_entropy = 0
    
    return {
        'recall': recall,
        'precision': precision,
        'f1': f1,
        'split_entropy': split_entropy,
        'best_cluster': best_cluster,
        'n_clusters_with_niche': len(niche_distribution)
    }

# Run analysis
res = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
knn = [3, 4, 5, 6, 7, 8, 9, 10]

results = []

for i in res: 
    for k in knn:
        neighbors_key = f"neighbors_{k}"
        cluster_key = f"{i}_{k}_cluster_graphST"
        
        sc.pp.neighbors(adata, use_rep="emb", n_neighbors=k, key_added=neighbors_key)
        sc.tl.leiden(adata, resolution=i, key_added=cluster_key, neighbors_key=neighbors_key)
        
        # Calculate metrics for TLS core 1 niche
        metrics = calculate_metrics_for_niche(adata, 'CD8+ TLS core 1', cluster_key)
        metrics['resolution'] = i
        metrics['k_neighbors'] = k
        results.append(metrics)

# Convert to DataFrame
results_df = pd.DataFrame(results)
results_df
results_df.to_csv('/nfs/users/nfs_m/ms83/mintflow-revision/benchmarking/graphST/analysis/kidney_TLSscore1_metrics.csv')

# Display best parameters
print("Top 5 parameter combinations by F1 score:")
print(results_df.nlargest(5, 'f1')[['resolution', 'k_neighbors', 'recall', 'precision', 'f1', 'split_entropy']])

# Create summary visualizations
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Pivot data for heatmaps
for idx, metric in enumerate(['recall', 'precision', 'f1', 'split_entropy']):
    ax = axes[idx // 2, idx % 2]
    pivot_data = results_df.pivot(index='resolution', columns='k_neighbors', values=metric)
    
    if metric == 'split_entropy':
        sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlGn_r', ax=ax, vmin=0)
    else:
        sns.heatmap(pivot_data, annot=True, fmt='.2f', cmap='RdYlGn', ax=ax, vmin=0, vmax=1)
    
    ax.set_title(f'{metric.replace("_", " ").title()} for TLS core 1 Recovery')
    ax.set_xlabel('k-neighbors')
    ax.set_ylabel('Resolution')

plt.tight_layout()
plt.savefig('TLScore1_recovery_metrics.pdf', dpi=300, bbox_inches='tight')
plt.show()

# Create figure with two subplots - main metrics and entropy
fig, axes = plt.subplots(1, 2, figsize=(6, 4), gridspec_kw={'width_ratios': [3, 1]})

# Left panel: Recall, Precision, F1
metrics_01 = ['recall', 'precision', 'f1']
plot_data_01 = []

for metric in metrics_01:
    for val in results_df[metric]:
        plot_data_01.append({
            'Metric': metric.replace('_', ' ').title(),
            'Value': val
        })

plot_df_01 = pd.DataFrame(plot_data_01)

sns.boxplot(data=plot_df_01, x='Metric', y='Value', ax=axes[0], palette='Set2')
sns.stripplot(data=plot_df_01, x='Metric', y='Value', ax=axes[0], 
             color='black', alpha=0.3, size=3)

axes[0].set_ylabel('Value', fontsize=12, fontweight='bold')
axes[0].set_xlabel('', fontsize=12)
axes[0].set_ylim(0, 1)
axes[0].grid(axis='y', alpha=0.3, linestyle='--')

# Right panel: Split Entropy
plot_data_entropy = []
for val in results_df['split_entropy']:
    plot_data_entropy.append({
        'Metric': 'Split Entropy',
        'Value': val
    })

plot_df_entropy = pd.DataFrame(plot_data_entropy)

sns.boxplot(data=plot_df_entropy, x='Metric', y='Value', ax=axes[1], color='pink')
sns.stripplot(data=plot_df_entropy, x='Metric', y='Value', ax=axes[1], 
             color='black', alpha=0.3, size=3)

axes[1].set_ylabel('Entropy', fontsize=12, fontweight='bold')
axes[1].set_xlabel('', fontsize=12)
axes[1].grid(axis='y', alpha=0.3, linestyle='--')

plt.suptitle('TLS core 1 Recovery Metrics Across All Parameter Combinations', 
             fontsize=12, fontweight='bold', y=1.00)
plt.tight_layout()
plt.savefig('TLScore1_recovery_boxplots.pdf', dpi=300, bbox_inches='tight')
plt.show()

adata.write_h5ad(save_dir + 'adata_graphST_spl_clusters_TLScore1_V2.h5ad')

# Get top 3 parameter combinations by F1 score
top_params = results_df.nlargest(3, 'f1')[['resolution', 'k_neighbors']]

# Create custom color palette
niches_to_highlight = ["CD8+ TLS core 1"]
custom_colors = []
for cat in adata.obs['INIGEN spatial labels'].cat.categories:
    if cat in niches_to_highlight:
        custom_colors.append(COLORS[cat])
    else:
        custom_colors.append('lightgrey')

# Generate UMAPs only for top parameter combinations
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

for idx, (_, row) in enumerate(top_params.iterrows()):
    i = row['resolution']
    k = int(row['k_neighbors'])
    
    neighbors_key = f"neighbors_{k}"
    cluster_key = f"{i}_{k}_cluster_graphST"
    
    # Compute UMAP for this specific neighbors_key
    sc.tl.umap(adata, neighbors_key=neighbors_key)
    
    # Get F1 score for title
    f1_score = results_df[(results_df['resolution']==i) & (results_df['k_neighbors']==k)]['f1'].values[0]
    
    # Plot
    sc.pl.umap(adata, color="INIGEN spatial labels", 
               legend_loc='on data', frameon=False,
               title=f"Res={i}, k={k}\nF1={f1_score:.3f}",
               palette=custom_colors,
               ax=axes[idx], show=False, size=20, alpha=0.8)

plt.tight_layout()
plt.savefig('TLScore1_top_recoveries.pdf', dpi=300, bbox_inches='tight')
plt.show()

adata.write_h5ad(save_dir + 'adata_graphST_spl_clusters_TLScore1_V2.h5ad')