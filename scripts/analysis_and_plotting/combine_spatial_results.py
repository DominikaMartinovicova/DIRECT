#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# combine_spatial_results.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Combine the results of spatial analysis from all samples into one file 
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
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    parser.add_argument('--celltype_key', help='phenotyping level',
                        dest='celltype_key',
                        type=str)
    parser.add_argument('--samples_list', help='list of samples to combine',
                        dest='samples_list',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

celltype_key = args.celltype_key
samples_list = args.samples_list.split(",")
print(samples_list)
output_dir_results=args.output_dir_results
output_dir_plots=args.output_dir_plots
os.makedirs(output_dir_results, exist_ok=True)
os.makedirs(output_dir_plots, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Sort samples to groups for comparison 
#   - biopsy vs. resection
#   - biopsy + <90 vs. biopsy + >=90
#   - resection + <90 vs. resection + >=90
#   - biopsy + <90 vs. resection + <90
#   - biopsy + >=90 vs. resection + >=90
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a dataframe with sample information
#--------------------------------------------------------------------------------
adata = sc.read_h5ad(args.adata)
sample_info = adata.obs[['sample', 'sample_type', 'regression', "MPR", "treatment_scheme"]].drop_duplicates().set_index('sample')
print(sample_info)

# Create a dictionary for each analysis
#--------------------------------------------------------------------------------
list_dir = []
for sample in samples_list:
    list_dir.append(args.input + sample + "/")  #list of directories with analyses results per sample

# Create dictionaries of results for each analysis
centrality_results={}
nhood_results={}
interaction_results={}
cooccurrence_results={}
for i, sample_dir in enumerate(list_dir):   # Loop over all folders and files in the directory
    print(sample_dir)
    sample = sample_dir.split("/")[-2]
    print("Reading results from " + sample)
    for file in os.listdir(sample_dir):
        if file.endswith("centrality_scores.csv"):
            centrality_results[sample] = file
        elif file.endswith("neighbors_enrichment.pkl"):
            nhood_results[sample] = file
        elif file.endswith("interaction_matrix.pkl"):
            interaction_results[sample] = file  
        elif file.endswith("co_occurrence_probabilities.pkl"):
            cooccurrence_results[sample] = file
        else:
            print(f"Skipping over {file}")

print("Centrality scores - " + str(len(centrality_results)))
print("Neighborhood enrichment - " + str(len(nhood_results)))
print("Interaction matrix - " + str(len(interaction_results)))
print("Cooccurrence probabilities - " + str(len(cooccurrence_results)))

# Create groups to compare
#--------------------------------------------------------------------------------





#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Combine results
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Combine centrality scores
#--------------------------------------------------------------------------------
print('Combining centrality scores...')



# Combine neighborhood enrichment
#--------------------------------------------------------------------------------
print('Combining neighborhood enrichments...')



# Combine interaction matrix 
#--------------------------------------------------------------------------------
print('Combining interaction matrices...')


# Combine co-occurrence probabilities
#--------------------------------------------------------------------------------
print('Combining co-occurrence probabilities...')




