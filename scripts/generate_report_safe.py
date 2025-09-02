# scripts/generate_report_safe.py

import os
import pandas as pd
import scanpy as sc

# -------------------------------
# Paths
# -------------------------------
project_root = os.path.dirname(os.path.dirname(__file__))
processed_dir = os.path.join(project_root, "data", "processed")
figures_dir = os.path.join(project_root, "figures")
report_file = os.path.join(project_root, "NeuroProject_Report.md")

processed_file = os.path.join(processed_dir, "pbmc3k_processed.h5ad")

# -------------------------------
# Check if processed file exists
# -------------------------------
if not os.path.exists(processed_file):
    print(f"[ERROR] Processed dataset not found: {processed_file}")
    print("Run neuro_pipeline.py first!")
    exit(1)

# -------------------------------
# Load processed dataset
# -------------------------------
adata = sc.read_h5ad(processed_file)

# -------------------------------
# Clustering info
# -------------------------------
clusters = adata.obs['leiden'].unique() if 'leiden' in adata.obs else []
n_clusters = len(clusters)

# -------------------------------
# Top marker genes (pandas DataFrame)
# -------------------------------
top_markers_df = pd.DataFrame()
if 'rank_genes_groups' in adata.uns:
    result = adata.uns['rank_genes_groups']
    groups = result['names'].dtype.names
    top_markers_df = pd.DataFrame({group: result['names'][group][:5] for group in groups})

# -------------------------------
# Generate Markdown
# -------------------------------
with open(report_file, "w") as f:
    f.write("# NeuroProject Analysis Report\n\n")
    
    f.write("## Dataset\n")
    f.write(f"- PBMC3k (processed): {adata.n_obs} cells × {adata.n_vars} genes\n\n")
    
    f.write("## Clustering\n")
    f.write(f"- Leiden clusters found: {n_clusters}\n\n")
    
    f.write("## Top Marker Genes\n")
    if not top_markers_df.empty:
        f.write(top_markers_df.to_markdown(index=True))
    else:
        f.write("No marker genes found.\n")
    f.write("\n\n")
    
    f.write("## Visualizations\n")
    for plot_name in ["umap_clusters.png", "umap_multi_marker.png", "violin_plot.png", "dotplot.png"]:
        path = os.path.join('figures', plot_name)
        f.write(f"- {plot_name}: `{path}`\n")
    
    f.write("\n## Notes\n")
    f.write(f"- Processed dataset: `{processed_file}`\n")
    marker_csv = os.path.join(processed_dir, "marker_genes.csv")
    f.write(f"- Marker genes CSV: `{marker_csv}`\n")

print(f"[INFO] Report generated → {report_file}")
