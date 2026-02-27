#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cell_fraction_analysis.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze shifts in cell fractions before and after treatment. Possibly split  
#   patients into groups based on chosen category.
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Define analysis functions
#   3 Create fractions dataframe
#   4 Choose analyses to perform:
#       a. Cell type fraction shifts (lineplot)
#       b. Cell type fraction shifts (boxplots)
#       c. Cell type fraction composition (boxplots)
#       d. Cell type fraction composition within sample type (boxplots)
#       * statistical testing functions - paired and independent samples
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
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import argparse
import os
import numpy as np
import pandas as pd
import math
import anndata as ad
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator
import argparse

import warnings
warnings.filterwarnings("ignore")

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 cell_fraction_analysis.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Analyze shifts in cell fractions before and after treatment.  ')
    parser.add_argument('-i', help='path to combined adata file',
                        dest='input',
                        type=str)
    parser.add_argument('--phen_level', help='key for cell type annotation in adata.obs',
                        dest='phen_level',
                        type=str)
    parser.add_argument('--exclude_v17', action='store_true',
                        help='Exclude v1.7 samples')
    parser.add_argument('-o', '--output_results', help='path to output dir with patches per sample',
                        dest='output_dir_results',
                        type=str)
    parser.add_argument('--output_dir_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

celltype_key = args.phen_level
exclude_v17 = args.exclude_v17
print(f'Excluding v1.7 samples: {exclude_v17}')
output_dir_results = args.output_dir_results
output_dir = args.output_dir_plots

# Make sure output directories exist
os.makedirs(output_dir, exist_ok=True)
os.makedirs(output_dir_results, exist_ok=True)

sns.set_style('whitegrid')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad(args.input)
print(adata)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define functions for analyses
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def celltype_fraction_composition_by_structure(df, output_dir, output_dir_results, exclude_v17, category=None, stat_test=mannwhitneyu, perform_stat_test=False, immune=False):
    print(df)
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    structures = df['structure'].unique()
    
    if category is None:
        # Plot all structures together
        plt.figure(figsize=(14, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'structure'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction', '')
        ax = sns.boxplot(data=df_melted, x="variable", y="value", hue="structure", palette='tab20')
        
        title_prefix = "Immune Cell Type" if immune else "Cell Type"
        title_suffix = " (excluding v1.7 treatment scheme)" if exclude_v17 else ""
        plt.title(f"{title_prefix} Fractions by Structure{title_suffix}")
        file_name = f'{output_dir}{"immune_" if immune else ""}celltype_fraction_by_structure{"_wo_v1.7" if exclude_v17 else "_w_v1.7"}.svg'
        
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        plt.legend(title='Structure')
        plt.tight_layout()
        plt.savefig(file_name, format='svg', bbox_inches='tight')
        plt.close()
    
    else:
        # Plot split by both structure and category
        plt.figure(figsize=(16, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'structure', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction', '')
        g = sns.catplot(data=df_melted, x="variable", y="value", hue="structure", col=category, 
                        kind='box', palette='tab20', height=6, aspect=1.5)
        
        title_prefix = "Immune Cell Type" if immune else "Cell Type"
        title_suffix = " (excluding v1.7 treatment scheme)" if exclude_v17 else ""
        g.fig.suptitle(f"{title_prefix} Fractions by Structure and {category}{title_suffix}")
        
        for ax in g.axes.flat:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_xlabel("Cell Type")
            ax.set_ylabel("Fraction")
        
        file_name = f'{output_dir}{category}_{"immune_" if immune else ""}celltype_fraction_by_structure{"_wo_v1.7" if exclude_v17 else "_w_v1.7"}.svg'
        plt.tight_layout()
        plt.savefig(file_name, format='svg', bbox_inches='tight')
        plt.close()


# Statistical testing functions
#--------------------------------------------------------------------------------
# Perform statistical testing
def paired_stat_testing(df, cell_fraction_cols, stat_test):
    # Perform statistical test for each cell type
    stat_results = []
    for celltype in cell_fraction_cols:
        biopsy_values = df[df['sample_type']=='Biopsy'][celltype]
        resection_values = df[df['sample_type']=='Resection'][celltype]
        # Ensure paired samples
        stat, p_value = stat_test(biopsy_values, resection_values)
        stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})

    stat_df = pd.DataFrame(stat_results)
    #print(stat_df)

    # Prepare stat_df for Annotator (expected format to be able to draw asterisks on plot)
    stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
    stat_df_annot["group1"] = "Biopsy"
    stat_df_annot["group2"] = "Resection"
    stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected format
    return stat_df, stat_df_annot

# Perform statistical testing for independent samples
def ind_stat_testing(df, cell_fraction_cols, stat_test, category=None):
    # Perform statistical test for each cell type
    stat_results = []
    for celltype in cell_fraction_cols:
        if category == None:
            biopsy_values = df[df['sample_type']=='Biopsy'][celltype]
            resection_values = df[df['sample_type']=='Resection'][celltype]
            # Ensure paired samples
            stat, p_value = stat_test(biopsy_values, resection_values)
            stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})
            stat_df = pd.DataFrame(stat_results)
            # Prepare stat_df for Annotator (expected format to be able to draw asterisks on plot)
            stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
            stat_df_annot["group1"] = "Biopsy"
            stat_df_annot["group2"] = "Resection"
            stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected forma
        else:
            categories = df[category].unique()
            subset_1 = df[df[category] == categories[0]][celltype]
            subset_2 = df[df[category] == categories[1]][celltype]
            stat, p_value = stat_test(subset_1, subset_2)
            stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})
            stat_df = pd.DataFrame(stat_results)   
            stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
            stat_df_annot["group1"] = categories[0]
            stat_df_annot["group2"] = categories[1]
            stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected forma
    #print(stat_df)
    return stat_df, stat_df_annot


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create fractions dataframe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate fractions per T_number
fractions_df = pd.DataFrame(dtype=object)
for i, element in enumerate(adata.obs['T_number'].unique().dropna()):
    adata_temp = adata[adata.obs['T_number'] == element, :] # Subset adata for element in T_number
    total_cells_temp = adata_temp.shape[0] # Total number of cells for this T_number
    temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions
    fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

    # Add metadata to the fractions_df
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme', 'MPR', 'treatment', 'structure']
    for meta in meta_list: 
        fractions_df.loc[meta, element] = adata_temp.obs[meta].unique()[0]

# Adjust dataframe for plotting
fractions_df = fractions_df.T.fillna(0) # Transpose for easier plotting and fill NaNs with 0
fractions_df.columns = [f'{col} fraction' if col not in meta_list else col for col in fractions_df.columns] # Add suffix to fraction columns

# Keep only patients with matched biopsy and resection samples
resection_pts = fractions_df[fractions_df['sample_type']=='Resection']['pt_id'].tolist()
biopsy_pts = fractions_df[fractions_df['sample_type']=='Biopsy']['pt_id'].tolist()
paired_pts = list(set(resection_pts) & set(biopsy_pts))
paired_fractions_df = fractions_df[fractions_df['pt_id'].isin(paired_pts)]
print(f'Number of paired patients: {len(paired_fractions_df["pt_id"].unique())}')


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if exclude_v17 == True:
    categories = [None, 'MPR'] 
else:
    categories = [None, 'MPR', 'treatment']
sample_types = ['Biopsy', 'Resection']

# All cell types
#--------------------------------------------------------------------------------

for category in categories:
    print(f'Analyzing category: {category}')
    for sample_type in sample_types:
        sample_type_df = fractions_df[fractions_df['sample_type']==sample_type]
        if sample_type=='Resection':
                celltype_fraction_composition_by_structure(sample_type_df, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=False)            

# Focus on immune cell types only
#--------------------------------------------------------------------------------
non_immune = ['Epithelial cell fraction', 'Fibroblast fraction', 'Endothelial cell fraction', 'Pericyte fraction', 'Stromal fraction', 'Tumor cells fraction']
to_exclude = set(non_immune).intersection(paired_fractions_df.columns)
print(f'Excluding non-immune cell types: {to_exclude}')
df_only_immune = paired_fractions_df.drop(labels=to_exclude, axis=1)

df_only_immune = df_only_immune[df_only_immune.columns[df_only_immune.columns.str.contains(' fraction')]]  # only keep columns that have ' fraction' in their name
df_immune = df_only_immune.div(df_only_immune.sum(axis=1), axis=0)  # Re-normalize to sum to 1
df_immune[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure']] = paired_fractions_df[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure']].values # Add metadata back
print(df_immune.head())

for category in categories:    
    for sample_type in sample_types:
        sample_type_df_immune = df_immune[df_immune['sample_type']==sample_type]
        if sample_type=='Resection':
            celltype_fraction_composition_by_structure(sample_type_df_immune, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=True)