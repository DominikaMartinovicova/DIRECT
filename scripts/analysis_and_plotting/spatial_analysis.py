#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# spatial_analysis.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze spatial proximity of cell types.
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Define analysis functions
#   3 
#   4 Choose analyses to perform:
#       a. Compute centrality scores
#           i. closeness centrality - measure of how close the group is to other nodes
#           ii. clustering coefficient - measure of the degree to which nodes cluster together
#           iii.degree centrality - fraction of non-group members connected to group members
#   5 Save
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
output_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/spatial/'
exclude_v17 = False
# Set aesthetics
sns.set_style("whitegrid")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define analysis functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build a spatial neighborhood graph
print('Building spatial neighbors graph...')
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

# Compute centrality scores
#--------------------------------------------------------------------------------
print('Computing centrality scores...')
sq.gr.centrality_scores(adata, cluster_key=celltype_key)
sq.pl.centrality_scores(adata, cluster_key=celltype_key, 
                        scores=["closeness_centrality", "clustering_coefficient", "degree_centrality"],
                        figsize=(16, 5), show=False)
plt.tight_layout()
plt.savefig(output_dir + 'centrality_scores.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

# Compute co-occurrence probability
#--------------------------------------------------------------------------------
print('Computing co-occurrence probabilities...')
sq.gr.co_occurrence(adata, cluster_key=celltype_key)
sq.pl.co_occurrence(adata, cluster_key=celltype_key, figsize=(8, 6), show=False)
plt.tight_layout()
plt.savefig(output_dir + 'co_occurrence_probabilities.svg',format='svg', dpi=300, bbox_inchces='tight')
plt.close()

# Compute neighbors enrichment
#--------------------------------------------------------------------------------
print('Computing neighbors enrichment...')
sq.gr.nhood_enrichment(adata, cluster_key=celltype_key)
sq.pl.nhood_enrichment(adata, cluster_key=celltype_key)
plt.tight_layout()
plt.savefig(output_dir + 'neighbors_enrichment.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

# Compute Moran's I spatial autocorrelation
#--------------------------------------------------------------------------------
print("Computing Moran's I spatial autocorrelation...")
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)
adata.uns['moranI'].head(10)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
