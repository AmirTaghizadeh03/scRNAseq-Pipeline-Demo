# NeuroProject Analysis Report

## Dataset
- PBMC3k (processed): 2700 cells × 2000 genes

## Clustering
- Leiden clusters found: 8

## Top Marker Genes
|    | 0      | 1      | 2    | 3        | 4    | 5      | 6        | 7     |
|---:|:-------|:-------|:-----|:---------|:-----|:-------|:---------|:------|
|  0 | LTB    | LYZ    | CCL5 | CD74     | NKG7 | LST1   | CD74     | PF4   |
|  1 | RPS27A | CST3   | NKG7 | HLA-DRA  | GNLY | AIF1   | HLA-DRB1 | PPBP  |
|  2 | RPL13  | TYROBP | IL32 | HLA-DPB1 | CTSW | FCER1G | ARPC1B   | GPX1  |
|  3 | IL32   | FTL    | B2M  | CD79A    | GZMB | FTL    | HLA-DRA  | GNG11 |
|  4 | RPS3A  | S100A9 | GZMA | HLA-DRB1 | PRF1 | COTL1  | HLA-DPB1 | CCL5  |

## Visualizations
- umap_clusters.png: `figures\umap_clusters.png`
- umap_multi_marker.png: `figures\umap_multi_marker.png`
- violin_plot.png: `figures\violin_plot.png`
- dotplot.png: `figures\dotplot.png`

## Notes
- Processed dataset: `c:\users\asus\desktop\neuroproject\data\processed\pbmc3k_processed.h5ad`
- Marker genes CSV: `c:\users\asus\desktop\neuroproject\data\processed\marker_genes.csv`
