#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# spatial_analysis_old.py
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
#           ii. clustering coefficient - measure of the degree to which nodes cluster together
#           iii. degree centrality - fraction of non-group members connected to group members
#       b. Compute neighbors enrichment - measure of how often cell types are neighbors compared to random expectation
#       c. Compute Moran's I spatial autocorrelation - measure of how gene expression is correlated with spatial location
#       d. Compute co-occurrence probabilities - measure of how often cell types co-occur in the same neighborhood
#   4 Save
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
import pickle

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
celltype_key = 'Neutro_Epi_extImm_pooled_A_EM_N'
adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{celltype_key}_combined_adatas.h5ad')
output_dir=f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{celltype_key}/spatial/'
output_dir_results=f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/{celltype_key}/spatial/'
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_dir_results, exist_ok=True)
exclude_v17 = True
if exclude_v17:
    print('Excluding v1.7 samples...')
    adata = adata[adata.obs['treatment_scheme'] != 'v1.7', :].copy()
# Set aesthetics
sns.set_style("whitegrid")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define analysis functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Build a spatial neighborhood graph
def spatial_analysis(adata, celltype_key=celltype_key, output_dir=output_dir, output_dir_results=output_dir_results):
    print('Building spatial neighbors graph...')
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)

    # Compute centrality scores
    #--------------------------------------------------------------------------------
    print('Computing centrality scores...')
    sq.gr.centrality_scores(adata, cluster_key=celltype_key)
    scores = ["average_clustering", "closeness_centrality", "degree_centrality"]
    centrality_scores_result = adata.uns[celltype_key + '_centrality_scores']
    centrality_scores_result.to_csv(output_dir_results + '/centrality_scores.csv')
    sq.pl.centrality_scores(adata, cluster_key=celltype_key, 
                            score=scores)
    plt.savefig(output_dir + '/centrality_scores.svg',format='svg', dpi=300, bbox_inches='tight')
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
    plt.savefig(output_dir + '/neighbors_enrichment.svg',format='svg', dpi=300, bbox_inches='tight')
    plt.close()

    # Compute Moran's I spatial autocorrelation
    #--------------------------------------------------------------------------------
    print("Computing Moran's I spatial autocorrelation...")
    sq.gr.spatial_autocorr(adata, mode='moran', n_perms=100, n_jobs=1)
    print(adata.uns['moranI'].head(10))

    # Compute co-occurrence probability
    #--------------------------------------------------------------------------------
    print('Computing co-occurrence probabilities...')
    sq.gr.co_occurrence(adata, cluster_key=celltype_key)
    co_occurrence_result = adata.uns[celltype_key + '_co_occurrence']
    with open(output_dir_results + "/co_occurrence_probabilities.pkl", "wb") as f:
        pickle.dump(co_occurrence_result, f)
    for celltype in adata.obs[celltype_key].unique().tolist():
        print(f'Plotting cell type: {celltype}')
        sq.pl.co_occurrence(adata, cluster_key=celltype_key, clusters=celltype, figsize=(12, 9))
        plt.savefig(output_dir + f'/co_occurrence_probabilities_{celltype}.svg',format='svg', dpi=300, bbox_inches='tight')
        plt.close()
    #sq.gr.co_occurrence(adata, cluster_key=celltype_key)
    #sq.pl.co_occurrence(adata, cluster_key=celltype_key) #, figsize=(8, 6))
    #plt.savefig(output_dir + 'co_occurrence_probabilities.svg',format='svg', dpi=300, bbox_inches='tight')
    #plt.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
categories = [None, 'MPR'] # 'treatment']
for category in categories:
    biopsy_adata = adata[adata.obs['sample_type'] == 'Biopsy', :].copy()
    resection_adata = adata[adata.obs['sample_type'] == 'Resection', :].copy()
    if category == None:
        print('====================== Analyzing all data together only split on sample type (Biopsy vs Resection). ======================')
        for adata_subset, sample_type in zip([biopsy_adata, resection_adata], ['Biopsy', 'Resection']):
            print(f'-------------------- Analyzing sample type: {sample_type} ----------------------')
            output_dir_sample = os.path.join(output_dir, 'None', sample_type)
            output_dir_results_sample = os.path.join(output_dir_results, 'None', sample_type)
            os.makedirs(output_dir_sample, exist_ok=True)
            os.makedirs(output_dir_results_sample, exist_ok=True)
            spatial_analysis(adata_subset, celltype_key=celltype_key, output_dir=output_dir_sample, output_dir_results=output_dir_results_sample)
    else:
        print(f'====================== Analyzing category: {category} ======================')
        for adata_subset, sample_type in zip([biopsy_adata, resection_adata], ['Biopsy', 'Resection']):
            for category_value in adata_subset.obs[category].unique().tolist():
                print(f'-------------------- Analyzing sample type: {sample_type}, category value: {category_value} ----------------------')
                adata_category = adata_subset[adata_subset.obs[category] == category_value, :].copy()
                output_dir_sample = os.path.join(output_dir, category, category_value, sample_type)
                output_dir_results_sample = os.path.join(output_dir_results, category, category_value, sample_type)
                os.makedirs(output_dir_sample, exist_ok=True)
                os.makedirs(output_dir_results_sample, exist_ok=True)
                spatial_analysis(adata_category, celltype_key=celltype_key, output_dir=output_dir_sample, output_dir_results=output_dir_results_sample)

#adata.write_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/analyzed/{celltype_key}_analyzed.h5ad')
print('Analysis complete. Results saved.')