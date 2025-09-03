# NeuroProject - Single-Cell RNA-Seq Analysis Pipeline

## Overview
This project demonstrates a complete single-cell RNA-seq analysis pipeline using the PBMC3k sample dataset from Scanpy. 
It includes preprocessing, dimensionality reduction, clustering, marker gene identification, and visualization.

## Pipeline Steps
1. **Preprocessing**: filter cells/genes, normalization, log transformation, highly variable gene selection  
2. **PCA and UMAP**: dimensionality reduction  
3. **Leiden clustering**: unsupervised clustering  
4. **Marker gene identification**: top genes for each cluster  
5. **Visualization**:
   - UMAP clusters
   - Multi-marker UMAP
   - Violin plots
   - Dot plots  
6. **Save processed data**: H5AD files and CSV of marker genes  

## Project Structure

NeuroProject/
│── data/ # Raw and processed data
│── figures/ # Plots (UMAPs, violin plots, dot plots, etc.)
│── notebooks/ # Jupyter notebooks for analysis
│── scripts/ # Python scripts
│── README.md # Project documentation


## Example Plots
### UMAP Clusters
![UMAP clusters](figures/umap_clusters.png)

### Multi-marker UMAP
![Multi-marker UMAP](figures/multi_marker_umap.png)

### Violin Plots
![Violin plots](figures/violin_plots.png)

### Dot Plots
![Dot plots](figures/dot_plot.png)

## Requirements
- Python 3.9+
- Scanpy
- Matplotlib
- Seaborn
- Pandas
- Numpy

Install with:
```bash
pip install scanpy matplotlib seaborn pandas numpy
```
How to Run

    Clone this repository:

git clone https://github.com/AmirTaghizadeh03/scRNAseq-Pipeline-Demo.git
cd scRNAseq-Pipeline-Demo

Open Jupyter Notebook or run the pipeline script:

    jupyter notebook notebooks/pipeline.ipynb

Notes

    Plots are stored in figures/.

    Processed datasets are stored in data/processed/.

    Raw PBMC3k data is stored in data/raw/.

