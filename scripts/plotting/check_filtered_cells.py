#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# check_filtered_cells.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Highlight where do potentially to be filtered cells end up in umap, how many
#   would that be and which cell types do they belong to
#
#
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/python/combine_adatas.py \
        -i {input.complete_Xenium} \
	    --input_dir {params.in_dir} \
        -o {output.combined_adatas} \
        --output_plot {params.out_plot_dir}
"""


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
import argparse
import os

# Read data
input = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/combined_adatas.h5ad'
adata = sc.read_h5ad(input)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Specify filtering criteria
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
min_transcripts_per_cell = 20
max_transcripts_per_cell = 98  # percentile
min_cells_per_gene = 100

thres = np.quantile(adata.obs['total_counts'], 0.98)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Highlight which cells would be filtered out
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a column in adata.obs that specifies whether a cell would be filtered out
adata.obs['filtered'] = (
    (adata.obs['transcript_counts'] < min_transcripts_per_cell))
#|   (adata.obs['transcript_counts'] > thres))

# Print the number of cells that would be filtered out
print(f"Number of cells to be filtered out: {adata.obs['filtered'].sum()} out of {adata.n_obs}")

# Print the number of cells that would be filtered out by cell type
filtered_by_cell_type = adata.obs[adata.obs['filtered']].groupby('celltype').size()
print("Number of cells to be filtered out by cell type:")
print(filtered_by_cell_type)

# Plot UMAP with filtered cells highlighted
sc.pl.umap(adata, color='filtered', title='Filtered Cells Highlighted')
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/to_be_filtered/umap_filtered_cells.png')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Print genes that would be filtered out
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate the number of cells each gene is expressed in
gene_cell_counts = np.array((adata.X > 0).sum(axis=0)).flatten()
genes_to_keep = gene_cell_counts >= min_cells_per_gene
filtered_genes = adata.var_names[~genes_to_keep]
print(f"Number of genes to be filtered out: {len(filtered_genes)} out of {adata.n_vars}")
print("Genes to be filtered out:")
print(filtered_genes)
