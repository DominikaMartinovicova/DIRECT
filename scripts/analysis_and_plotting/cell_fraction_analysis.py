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
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
category = 'structure' # e.g., structure, treatment, response

# Create missing directories
output_plot_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/celltype_fraction/'
#os.makedirs(output_plot_dir + 'boxplots/', exist_ok=True)
#os.makedirs(output_plot_dir + 'swarmplots/', exist_ok=True)
#os.makedirs(output_plot_dir + 'lineplots/', exist_ok=True)

# Set aesthetics
sns.set_style("whitegrid")
sns.color_palette("Tab20")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Calculate fractions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fractions_df = pd.DataFrame()
for element in adata.obs['T_number'].unique().dropna():
    print('Processing T_number: ' + element)
    adata_temp = adata[adata.obs['T_number'] == element, :] # Subset adata for element in T_number
    total_cells_temp = adata_temp.shape[0] # Total number of cells for this T_number
    temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions


    for celltype in temp_fractions.index.to_list():
        adata_temp.obs[celltype + '_fraction'] = temp_fractions[celltype]

    # save fractions to df
    fractions_df = pd.concat([fractions_df, temp_fractions], axis=1)
    
    # add to adata_patient
    if 'adata_sample' in locals():
        adata_sample = ad.concat([adata_sample, adata_temp], join = 'outer')
    else:
        adata_sample = adata_temp


    # add metadata
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme'] #'structure', 


    
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyse shifts in cell fractions before and after treatment
def celltype_fraction_shifts(adata, fraction_columns, category, output_dir):





# 2. Swarmplots of fractions per chosen category (e.g., structure) with paired pts connected
plot_swarmplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)
# 3. Lineplots of fractions per chosen category (e.g., structure) with paired pts connected
plot_lineplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)