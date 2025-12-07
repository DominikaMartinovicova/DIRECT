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
#   2 Calculate fraction
#   3 Save
#
#
#
#
# Author: Mischa Steketee (m.f.b.steketee@amsterdamumc.nl)
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
category = 'structure' # e.g., structure, treatment, response

# Create missing directories
output_plot_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/celltype_fraction/'
#os.makedirs(output_plot_dir + 'boxplots/', exist_ok=True)
#os.makedirs(output_plot_dir + 'swarmplots/', exist_ok=True)
#os.makedirs(output_plot_dir + 'lineplots/', exist_ok=True)

# Set aesthetics
sns.set_style("whitegrid")
sns.color_palette("tab20")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Create fractions dataframe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Calculate fractions per T_number
fractions_df = pd.DataFrame(dtype=object)
for i, element in enumerate(adata.obs['T_number'].unique().dropna()):
    print(f'Processing {i}. T_number: {element}')
    adata_temp = adata[adata.obs['T_number'] == element, :] # Subset adata for element in T_number
    total_cells_temp = adata_temp.shape[0] # Total number of cells for this T_number
    temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions
    fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

    # Add metadata to the fractions_df
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme'] #'structure',
    for meta in meta_list: 
        fractions_df.loc[meta, element] = adata_temp.obs[meta].unique()[0]

# Adjust dataframe for plotting
fractions_df = fractions_df.T.fillna(0) # Transpose for easier plotting and fill NaNs with 0
fractions_df.columns = [f'{col} fraction' if col not in meta_list else col for col in fractions_df.columns] # Add suffix to fraction columns
fractions_df['MPR'] = np.where(fractions_df['regression']>=90, '>=90', '<90') # Create MPR column
fractions_df['treatment'] = np.where(fractions_df['treatment_scheme'] == 'v1.7', 'aggressive', 'milder') # Create treatment column
print(fractions_df.head())

# Keep only patients with matched biopsy and resection samples
paired_pts = fractions_df['pt_id'][fractions_df['sample_type']=='biopsy'].isin(fractions_df['pt_id'][fractions_df['sample_type']=='resection'])
paired_fractions_df = fractions_df[fractions_df['pt_id'].isin(paired_pts.index.unique())]
print(f'Number of paired patients: {len(paired_fractions_df["pt_id"].unique())}')

# Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyse shifts in cell fractions before and after treatment
def celltype_fraction_shifts(df, category, output_dir):
    # Split data into pre- and post-treatment
    biopsy_df = df[df['sample_type']=='biopsy']
    resection_df = df[df['sample_type']=='resection']
    cell_fraction_keys = [col for col in df.columns if col.endswith('fraction')]
    if category == None:
           # Do not split into groups, compare biopsy vs resection for all patients





celltype_fraction_shifts(paired_fractions_df, )
# 2. Swarmplots of fractions per chosen category (e.g., structure) with paired pts connected
#plot_swarmplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)
# 3. Lineplots of fractions per chosen category (e.g., structure) with paired pts connected
#plot_lineplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)