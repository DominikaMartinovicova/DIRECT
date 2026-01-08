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
#   3 Choose analyses to perform:
#       a. Compute centrality scores
#           i. closeness centrality - measure of how close the group is to other nodes
#           ii. degree centrality - fraction of non-group members connected to group members
#           iii. betweenness centrality - measure of how often a node appears on shortest paths between other nodes
#           iv. clustering coefficient - measure of the degree to which nodes cluster together
#       b. Compute neighbors enrichment - measure of how often cell types are neighbors compared to random expectation
#       c. Compute Moran's I spatial autocorrelation - measure of how gene expression is correlated with spatial location
#       d. Compute interaction matrix - measure of interactions between cell types
#       e. Compute co-occurrence probabilities - measure of how often cell types co-occur in the same neighborhood
#   4 Save
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import pickle
import argparse

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 spatial_analysis.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Spatially analyze each sample.  ')
    parser.add_argument('-i', help='path to adata sample subset',
                        dest='input',
                        type=str)
    parser.add_argument('-o_results', help='path to output dir for results',
                        dest='output_dir_results',
                        type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    parser.add_argument('--celltype_key', help='phenotyping level',
                        dest='celltype_key',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

celltype_key = args.celltype_key
output_dir_results=args.output_dir_results
output_dir_plots=args.output_dir_plots
os.makedirs(output_dir_results, exist_ok=True)
os.makedirs(output_dir_plots, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad(args.input)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Perform spatial analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build a spatial neighborhood graph
print('Building spatial neighbors graph...')
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

# Compute centrality scores
#--------------------------------------------------------------------------------
print('Computing centrality scores...')
sq.gr.centrality_scores(adata, cluster_key=celltype_key)
scores = ["average_clustering", "closeness_centrality", "degree_centrality", "betweenness_centrality"]
centrality_scores_result = adata.uns[celltype_key + '_centrality_scores']
centrality_scores_result.to_csv(output_dir_results + '/centrality_scores.csv')
sq.pl.centrality_scores(adata, cluster_key=celltype_key, 
                        score=scores)
#plt.savefig(output_dir_plots + '/centrality_scores.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

# Compute neighbors enrichment
#--------------------------------------------------------------------------------
print('Computing neighbors enrichment...')
sq.gr.nhood_enrichment(adata, cluster_key=celltype_key)
nhood_enrichment_result = adata.uns[celltype_key + '_nhood_enrichment']
with open(output_dir_results + "/neighbors_enrichment.pkl", "wb") as f:
    pickle.dump(nhood_enrichment_result, f)
sq.pl.nhood_enrichment(adata, cluster_key=celltype_key, figsize=(8, 8), cmap='bwr')
plt.rcParams["font.size"] = 45
#plt.savefig(output_dir_plots + '/neighbors_enrichment.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

# Compute Moran's I spatial autocorrelation
#--------------------------------------------------------------------------------
print("Computing Moran's I spatial autocorrelation...")
sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)
#print(adata.uns['moranI'].head(10))

# Compute interaction matrix
#--------------------------------------------------------------------------------
print('Computing interaction matrix...')
sq.gr.interaction_matrix(adata, cluster_key=celltype_key)
interaction_matrix_result = adata.uns[celltype_key + '_interactions']
with open(output_dir_results + "/interaction_matrix.pkl", "wb") as f:
    pickle.dump(interaction_matrix_result, f)

# Compute co-occurrence probability
#--------------------------------------------------------------------------------
print('Computing co-occurrence probabilities...')
sq.gr.co_occurrence(adata, cluster_key=celltype_key)
co_occurrence_result = adata.uns[celltype_key + '_co_occurrence']
with open(output_dir_results + "/co_occurrence_probabilities.pkl", "wb") as f:
    pickle.dump(co_occurrence_result, f)
for celltype in adata.obs[celltype_key].unique().tolist():
#    print(f'Plotting cell type: {celltype}')
    sq.pl.co_occurrence(adata, cluster_key=celltype_key, clusters=celltype, figsize=(12, 9))
#    plt.savefig(output_dir_plots + f'/co_occurrence_probabilities_{celltype}.svg',format='svg', dpi=300, bbox_inches='tight')
    plt.close()
#sq.gr.co_occurrence(adata, cluster_key=celltype_key)
#sq.pl.co_occurrence(adata, cluster_key=celltype_key) #, figsize=(8, 6))
#plt.savefig(output_dir + 'co_occurrence_probabilities.svg',format='svg', dpi=300, bbox_inches='tight')
#plt.close()

