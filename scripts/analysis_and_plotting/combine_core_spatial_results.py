#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# combine_spatial_results.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Combine the results of spatial analysis from all cores into one file 
#   per analysis type.
#
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
import squidpy as sq

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 combine_spatial_results.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Combine the results of spatial analysis from all samples into one file per analysis type. ') 
    parser.add_argument('-i', help='path to input directory with spatial analysis results per sample',
                        dest='input',
                        type=str)
    parser.add_argument('--adata', help='path to adata object with all samples combined',
                        dest='adata',
                        type=str)
    parser.add_argument('-o_results', help='path to output dir for results',
                        dest='output_dir_results',
                        type=str)
    # parser.add_argument('--samples_list', help='list of samples to combine',
    #                     dest='samples_list',
    #                     type=str)
    args = parser.parse_args()
    return args

args = parse_args()

output_dir_results=args.output_dir_results
os.makedirs(output_dir_results, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Sort samples to groups for comparison 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a dataframe with sample information
#--------------------------------------------------------------------------------
adata = sc.read_h5ad(args.adata)
print(adata)
core_info = adata.obs[['sample', 'sample_type', 'regression', "MPR", "treatment_scheme", 'T_number','pt_id']].drop_duplicates().set_index('sample')
print(core_info)
cores_list = core_info.index.tolist()  #list of cores to combine

# Create a dictionary for each analysis
#--------------------------------------------------------------------------------
list_dir = []
for core in cores_list:
    list_dir.append(args.input + core + "/")  #list of directories with analyses results per sample

# Create dictionaries of results for each analysis
centrality_results={}
nhood_results={}
interaction_results={}
cooccurrence_results={}
for i, core_dir in enumerate(list_dir):   # Loop over all folders and files in the directory
    core = core_dir.split("/")[-2]
    #print(f"Reading results from {core}")
    for file in os.listdir(core_dir):
        if file.endswith("centrality_scores.csv"):
            file = pd.read_csv(os.path.join(core_dir, file), index_col=0)
            centrality_results[core] = file
        elif file.endswith("neighbors_enrichment.pkl"):
            with open(os.path.join(core_dir, file), 'rb') as f:
                file = pickle.load(f)
            nhood_results[core] = file
        elif file.endswith("interaction_matrix.pkl"):
            with open(os.path.join(core_dir, file), 'rb') as f:
                file = pickle.load(f)
            interaction_results[core] = file  
        elif file.endswith("co_occurrence_probabilities.pkl"):
            with open(os.path.join(core_dir, file), 'rb') as f:
                file = pickle.load(f)
            cooccurrence_results[core] = file
        #else:
            #print(f"Skipping over {file}")

print("Centrality scores - " + str(len(centrality_results)))
print("Neighborhood enrichment - " + str(len(nhood_results)))
print("Interaction matrix - " + str(len(interaction_results)))
print("Cooccurrence probabilities - " + str(len(cooccurrence_results)))



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Combine results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combine centrality scores
#--------------------------------------------------------------------------------
print('Combining centrality scores...')
degree_centrality_scores = pd.DataFrame()
average_clustering_scores = pd.DataFrame()
closeness_centrality_scores = pd.DataFrame()

for core in centrality_results.keys():
    for analysis in centrality_results[core].columns:
        if analysis == 'degree_centrality':
            degree_centrality_scores[core] = centrality_results[core][analysis]
        elif analysis == 'average_clustering':
            average_clustering_scores[core] = centrality_results[core][analysis]
        elif analysis == 'closeness_centrality':
            closeness_centrality_scores[core] = centrality_results[core][analysis]
        else:
            print(f'Unknown analysis type: {analysis}')

# Add sample info
degree_centrality_scores = degree_centrality_scores.transpose()
degree_centrality_scores = degree_centrality_scores.join(core_info[['sample_type', 'MPR', 'pt_id']], how='left')

average_clustering_scores = average_clustering_scores.transpose()
average_clustering_scores = average_clustering_scores.join(core_info[['sample_type', 'MPR', 'pt_id']], how='left')

closeness_centrality_scores = closeness_centrality_scores.transpose()
closeness_centrality_scores = closeness_centrality_scores.join(core_info[['sample_type', 'MPR', 'pt_id']], how='left')

combined_centrality_scores = {'degree_centrality': degree_centrality_scores,'average_clustering': average_clustering_scores,'closeness_centrality': closeness_centrality_scores}
with open(os.path.join(output_dir_results, 'combined_centrality_scores.pkl'), 'wb') as f:
    pickle.dump(combined_centrality_scores, f)

# Combine neighborhood enrichment
#--------------------------------------------------------------------------------
print('Combining neighborhood enrichments...')
combined_nhood_enrichment = {}

for core in nhood_results.keys():
    combined_nhood_enrichment[core] = nhood_results[core]
    combined_nhood_enrichment[core]['sample_type'] = core_info.loc[core, 'sample_type']
    combined_nhood_enrichment[core]['MPR'] = core_info.loc[core, 'MPR']
    combined_nhood_enrichment[core]['pt_id'] = core_info.loc[core, 'pt_id']
#print(combined_nhood_enrichment)
with open(os.path.join(output_dir_results, 'combined_neighbors_enrichment.pkl'), 'wb') as f:
    pickle.dump(combined_nhood_enrichment, f)

# Combine interaction matrix 
#--------------------------------------------------------------------------------
print('Combining interaction matrices...')
combined_interaction_matrices = {}
for core in interaction_results.keys():
    combined_interaction_matrices[core] = {'matrix': interaction_results[core],
                                             'sample_type': core_info.loc[core, 'sample_type'],
                                             'MPR': core_info.loc[core, 'MPR']}
#print(combined_interaction_matrices)
with open(os.path.join(output_dir_results, 'combined_interaction_matrix.pkl'), 'wb') as f:
    pickle.dump(combined_interaction_matrices, f)

# Combine co-occurrence probabilities
#--------------------------------------------------------------------------------
print('Combining co-occurrence probabilities...')
combined_cooccurrence_probs={}
for core in cooccurrence_results.keys():
    combined_cooccurrence_probs[core] = cooccurrence_results[core]
    combined_cooccurrence_probs[core]['sample_type'] = core_info.loc[core, 'sample_type']
    combined_cooccurrence_probs[core]['MPR'] = core_info.loc[core, 'MPR']
    combined_cooccurrence_probs[core]['pt_id'] = core_info.loc[core, 'pt_id']
#print(combined_cooccurrence_probs)
with open(os.path.join(output_dir_results, 'combined_cooccurrence_probabilities.pkl'), 'wb') as f:
    pickle.dump(combined_cooccurrence_probs, f)

print('Done combining results!')

