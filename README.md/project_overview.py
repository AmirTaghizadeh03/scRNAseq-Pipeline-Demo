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
