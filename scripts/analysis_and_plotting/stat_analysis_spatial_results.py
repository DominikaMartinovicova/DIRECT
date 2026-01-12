#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_spatial_results.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze and plot combined spatial analysis results across samples.
#   
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
from matplotlib import category
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
    parser = argparse.ArgumentParser(prog = 'python3 stat_analysis_spatial_results.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Perform statistical analysis and plotting between groups of samples. ') 
    parser.add_argument('-i', help='path to input directory with spatial analysis results per sample',
                        dest='input',
                        type=str)
    parser.add_argument('-o_report', help='report with stat analysis results',
                        dest='output_dir_report',
                        type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
input_dir=args.input
output_dir_report=args.output_dir_report
output_dir_plots=args.output_dir_plots

# Set aesthetics
sns.set_style("whitegrid")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def centrality_shifts(centrality_scores, category):
    for key in centrality_scores.keys():
        scores_df = centrality_scores[key]
        if category == None:
            scores_df_melted = scores_df.melt(id_vars=['sample_type'], var_name='cell_type', value_name=key).drop(columns=['MPR'])


def stat_analysis_centrality_scores(input_file, output_dir_report, output_dir_plots, category):
    # Centrality scores analysis and plotting
    #------------------------------------------------------------------------------
    centrality_scores = input_file
    print(centrality_scores.keys())
    for key in centrality_scores.keys():
        if category == None:
            scores_df = centrality_scores[key].drop(columns=['MPR'])
            scores_df_melted = scores_df.melt(id_vars=['sample_type'], var_name='cell_type', value_name=key)
            plt.figure(figsize=(20,15))
            g=sns.catplot(scores_df_melted, x='cell_type', y=key, hue='sample_type', kind='bar', palette='tab20')
            plt.title(f'Centrality scores: {key}')
            plt.xticks(rotation=90)
            plt.xlabel('Cell type')
            plt.ylabel(f'{key} score')
            g.legend.set_title('Sample type')
            g.legend.set_loc('upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir_plots, f'centrality_scores/{key}.svg'), format='svg', bbox_inches='tight')
            plt.close()
        elif category != None:
            scores_df = centrality_scores[key]
            scores_df_melted = scores_df.melt(id_vars=['sample_type', category], var_name='cell_type', value_name=key)
            plt.figure(figsize=(20,15))
            g=sns.catplot(scores_df_melted, x='cell_type', y=key, hue='sample_type', col=category, kind='box', palette='tab20')
            g.set_xticklabels(rotation=90)
            g.set_xlabels('Cell type')
            g.set_ylabels(f'{key} score')
            plt.suptitle(f'Centrality scores: {key} by {category}', y=1.03)
            g.legend.set_title('Sample type')
            g.legend.set_loc('upper right')
            plt.tight_layout()
            plt.savefig(os.path.join(output_dir_plots, f'centrality_scores/{key}_by_{category}.svg'), format='svg', bbox_inches='tight')
            plt.close()

# def stat_analysis_nhood_enrichment(input_file, output_dir_report, output_dir_plots, category, mode, cell_types):
#     # Neighborhood enrichment analysis and plotting
#     #------------------------------------------------------------------------------
#     nhood_enrichment = input_file
#     if category == None:
#         biopsy_samples = {}
#         resection_samples = {}
#         for sample in nhood_enrichment.keys():
#             if nhood_enrichment[sample]['sample_type']=='Biopsy':
#                 biopsy_samples[sample] = nhood_enrichment[sample]
#             elif nhood_enrichment[sample]['sample_type']=='Resection':
#                 resection_samples[sample] = nhood_enrichment[sample]
        
#         # Calculate average in each group
#         zscore_biopsy = [biopsy_samples[sample][mode] for sample in biopsy_samples.keys()]
#         for i, arr in enumerate(zscore_biopsy):
#             print(f"Array {i}: shape = {arr.shape}")

#         mean_zscore_biopsy = np.mean(np.stack(zscore_biopsy, axis=0), axis=0)

#         # Calculate average in each group
#         zscore_resection = [resection_samples[sample][mode] for sample in resection_samples.keys()]
#         mean_zscore_resection = np.mean(np.stack(zscore_resection, axis=0), axis=0)

#         # Plot heatmaps
#         plt.figure(figsize=(20,15))
#         sns.heatmap(mean_zscore_biopsy, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Average neighborhood enrichment {mode} - Biopsy samples')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'neighborhood_enrichment/average_{mode}_biopsy.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#         plt.figure(figsize=(20,15))
#         sns.heatmap(mean_zscore_resection, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
#         plt.title(f'Average neighborhood enrichment {mode} - Resection samples')
#         plt.tight_layout()
#         plt.savefig(os.path.join(output_dir_plots, f'neighborhood_enrichment/average_{mode}_resection.svg'), format='svg', bbox_inches='tight')
#         plt.close()

#     elif category != None:
#         print('Statistical analysis by category not yet implemented for neighborhood enrichment.')

def stat_analysis_interaction_matrices(input_file, output_dir_report, output_dir_plots, category, cell_types):
    # Interaction matrix analysis and plotting
    #------------------------------------------------------------------------------
    interaction_matrices = input_file
    print(interaction_matrices.keys())
    print(interaction_matrices[list(interaction_matrices.keys())[0]])
    if category == None:
        biopsy_samples = {}
        resection_samples = {}
        for sample in interaction_matrices.keys():
            if interaction_matrices[sample]['sample_type']=='Biopsy':
                biopsy_samples[sample] = interaction_matrices[sample]
            elif interaction_matrices[sample]['sample_type']=='Resection':
                resection_samples[sample] = interaction_matrices[sample]
        
        # Calculate average in each group
        matrix_biopsy = [biopsy_samples[sample]['matrix'] for sample in biopsy_samples.keys()]
        mean_matrix_biopsy = np.mean(np.stack(matrix_biopsy, axis=0), axis=0)

        # Calculate average in each group
        matrix_resection = [resection_samples[sample]['matrix'] for sample in resection_samples.keys()]
        mean_matrix_resection = np.mean(np.stack(matrix_resection, axis=0), axis=0)

        # Plot heatmaps
        plt.figure(figsize=(20,15))
        sns.heatmap(mean_matrix_biopsy, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
        plt.title(f'Average interaction matrix - Biopsy samples')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/average_interaction_matrix_biopsy.svg'), format='svg', bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(20,15))
        sns.heatmap(mean_matrix_resection, cmap='vlag', center=0, xticklabels=cell_types, yticklabels=cell_types)
        plt.title(f'Average interaction matrix - Resection samples')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir_plots, f'interaction_matrices/average_interaction_matrix_resection.svg'), format='svg', bbox_inches='tight')
        plt.close()

    elif category != None:
        print('Statistical analysis by category not yet implemented for interaction matrices.')



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Prepare data and run analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
centrality_scores_path = os.path.join(input_dir, 'combined_centrality_scores.pkl')
with open(centrality_scores_path, 'rb') as f:
    centrality_scores = pickle.load(f)
cell_types = centrality_scores['degree_centrality'].columns.tolist()
#stat_analysis_centrality_scores(input_file=centrality_scores, output_dir_report=output_dir_report, output_dir_plots=output_dir_plots, category=None)

nhood_enrichment_path = os.path.join(input_dir, 'combined_neighbors_enrichment.pkl')
with open(nhood_enrichment_path, 'rb') as f:
    nhood_enrichment = pickle.load(f)
# stat_analysis_nhood_enrichment(input_file=nhood_enrichment, output_dir_report=output_dir_report, output_dir_plots=output_dir_plots, 
                            #    category=None, mode='zscore', cell_types=cell_types)

interaction_matrix_path = os.path.join(input_dir, 'combined_interaction_matrix.pkl')
with open(interaction_matrix_path, 'rb') as f:
    interaction_matrices = pickle.load(f)
stat_analysis_interaction_matrices(input_file=interaction_matrices, output_dir_report=output_dir_report, output_dir_plots=output_dir_plots, category=None, cell_types=cell_types)




