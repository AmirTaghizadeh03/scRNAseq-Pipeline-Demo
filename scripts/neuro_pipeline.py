# scripts/neuro_pipeline.py

import os
import scanpy as sc
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

# -------------------------------
# Use safe backend for saving plots
# -------------------------------
matplotlib.use('Agg')

# -------------------------------
# Project paths
# -------------------------------
project_root = os.path.dirname(os.path.dirname(__file__))
raw_dir = os.path.join(project_root, "data", "raw")
processed_dir = os.path.join(project_root, "data", "processed")
figures_dir = os.path.join(project_root, "figures")

os.makedirs(raw_dir, exist_ok=True)
os.makedirs(processed_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)

# -------------------------------
# Load PBMC3k sample dataset
# -------------------------------
print("Loading sample PBMC3k dataset...")
adata = sc.datasets.pbmc3k()

# -------------------------------
# Preprocess data
# -------------------------------
print("Preprocessing data...")
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)

# -------------------------------
# PCA and UMAP
# -------------------------------
print("Computing PCA and UMAP...")
sc.tl.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# -------------------------------
# Leiden clustering
# -------------------------------
print("Running Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5)

# -------------------------------
# Plot UMAP clusters
# -------------------------------
umap_clusters_path = os.path.join(figures_dir, "umap_clusters.png")
sc.pl.umap(adata, color=["leiden"], show=False)
plt.savefig(umap_clusters_path, dpi=150)
plt.close()
print(f"UMAP clusters plot saved → {umap_clusters_path}")

# -------------------------------
# Find marker genes
# -------------------------------
print("Finding marker genes...")
sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')

marker_plot_path = os.path.join(figures_dir, "marker_genes.png")
sc.pl.rank_genes_groups(adata, n_genes=5, sharey=False, show=False)
plt.savefig(marker_plot_path, dpi=150)
plt.close()
print(f"Marker genes plot saved → {marker_plot_path}")

# Save marker genes to CSV
result = adata.uns['rank_genes_groups']
groups = result['names'].dtype.names
marker_df = pd.DataFrame({group: result['names'][group] for group in groups})
marker_csv_path = os.path.join(processed_dir, "marker_genes.csv")
marker_df.to_csv(marker_csv_path, index=False)
print(f"Marker genes CSV saved → {marker_csv_path}")

# -------------------------------
# Plot multiple marker genes on UMAP
# -------------------------------
marker_genes = ['IL7R', 'CD14', 'MS4A1']
existing_genes = [gene for gene in marker_genes if gene in adata.var_names]

if existing_genes:
    multi_marker_plot_path = os.path.join(figures_dir, "umap_multi_marker.png")
    sc.pl.umap(adata, color=existing_genes, wspace=0.4, show=False)
    plt.savefig(multi_marker_plot_path, dpi=150)
    plt.close()
    print(f"UMAP multi-marker plot saved → {multi_marker_plot_path}")

    # Save expression of marker genes to CSV
    marker_expr = pd.DataFrame({gene: adata[:, gene].X.toarray().flatten() if hasattr(adata[:, gene].X, "toarray") else adata[:, gene].X.flatten() for gene in existing_genes})
    multi_marker_csv_path = os.path.join(processed_dir, "multi_marker_expression.csv")
    marker_expr.to_csv(multi_marker_csv_path, index=False)
    print(f"Multi-marker expression CSV saved → {multi_marker_csv_path}")
else:
    print("No marker genes found in dataset!")

# -------------------------------
# Violin and Dot plots
# -------------------------------
if existing_genes:
    violin_path = os.path.join(figures_dir, "violin_plot.png")
    dotplot_path = os.path.join(figures_dir, "dotplot.png")
    sc.pl.violin(adata, keys=existing_genes, groupby='leiden', rotation=45, show=False)
    plt.savefig(violin_path, dpi=150)
    plt.close()
    sc.pl.dotplot(adata, var_names=existing_genes, groupby='leiden', show=False)
    plt.savefig(dotplot_path, dpi=150)
    plt.close()
    print(f"Violin and Dot plots saved → {violin_path}, {dotplot_path}")

# -------------------------------
# Save processed data
# -------------------------------
processed_file = os.path.join(processed_dir, "pbmc3k_processed.h5ad")
adata.write(processed_file)
print(f"Processed data saved → {processed_file}")

print("Pipeline complete ✅")
