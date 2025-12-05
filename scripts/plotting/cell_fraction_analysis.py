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
celltype = 'Neutro_Epi_extImm'

# Create missing directories
output_plot_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/celltype_fractions/'
os.makedirs(output_plot_dir + 'boxplots/', exist_ok=True)
os.makedirs(output_plot_dir + 'swarmplots/', exist_ok=True)
os.makedirs(output_plot_dir + 'lineplots/', exist_ok=True)

# Set aesthetics
sns.set_style("whitegrid")
sns.color_palette("Tab20")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Calculate fractions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for element in adata.obs['T_number'].unique().dropna():
    print('Processing T_number: ' + element)
    adata_temp = adata[adata.obs['T_number'] == element, :] # Subset adata for element in T_number
    # Calculate fractions per samples
    temp_fractions = adata_temp.obs[celltype].value_counts()/np.sum(adata_temp.obs[celltype].value_counts())
    total_cells_temp = adata_temp.shape[0]
    
    ## save total # of cells
    adata_temp.obs['ncells'] = total_cells_temp
    ## Keep 1 row of sample and add fractions obs
    adata_temp = adata_temp[0,:]
    for celltype in temp_fractions.index.to_list():
        adata_temp.obs[celltype + '_fraction'] = temp_fractions[celltype]
    
    ## add to adata_patient
    if 'adata_sample' in locals():
        adata_sample = ad.concat([adata_sample, adata_temp], join = 'outer')
    else:
        adata_sample = adata_temp


#-------------------------------------------------------------------------------
# 3.2 Do analysis
#-------------------------------------------------------------------------------
## count biopsies/resections/paired
adata_sample.obs.sample_type.value_counts()

## Look at fractions in general
fraction_columns = [s + '_fraction' for s in adata.obs.celltype.unique().to_list()]

## Check fractions in control sample
pd.set_option('display.max_columns', None)
adata_sample.obs.loc[adata_sample.obs.structure == 'tumor_bed',fraction_columns]


adata_sample.obs[fraction_columns].boxplot()
plt.xticks(rotation=90)
plt.savefig(output_plot_dir + '/boxplots/boxplot_fractions_columns.png', dpi=300, bbox_inches='tight')
plt.close()

## look at paired fractions
paired_pts = adata_sample[adata_sample.obs.sample_type == 'tumor_bed'].obs.pt_id.to_list()

adata_sample_paired = adata_sample[adata_sample.obs.pt_id.isin(paired_pts)]

# for celltype in fraction_columns:
#     file_name = celltype.replace(' ', '_').replace('/', '_')
#     print('Plotting celltype: ' + celltype)
#     sns.boxplot(data = adata_sample_paired.obs, x= 'sample_type', y = celltype)
#     plt.savefig(output_plot_dir + 'boxplot_' + file_name + '_biopsy_resection.png', dpi=300, bbox_inches='tight')
#     plt.close()

#     sns.swarmplot(data = adata_sample_paired.obs, x= 'sample_type', y = celltype, hue = 'pt_id')
#     plt.savefig(output_plot_dir + 'swarmplot_' + file_name + '_biopsy_resection.png', dpi=300, bbox_inches='tight')
#     plt.close()

#     sns.lineplot(data = adata_sample_paired.obs, x= 'sample_type', y = celltype, units = 'pt_id', estimator = None, hue = 'pt_id')
#     plt.savefig(output_plot_dir + 'lineplot_' + file_name + '_biopsy_resection.png', dpi=300, bbox_inches='tight')
#     plt.close()

## combined boxplot
df_melted = adata_sample.obs.melt(id_vars=['structure', 'pt_id'], value_vars=fraction_columns,
                                  var_name='celltype', value_name='fraction')

print(df_melted)

plt.figure(figsize=(12, 6))
sns.boxplot(data=df_melted, x='celltype', y='fraction', hue='structure', palette='Paired')
plt.xticks(rotation=90)
plt.xlabel('Cell type')
plt.ylabel('Fraction')
plt.title('Cell Type Fractions Across Structures')
plt.legend(title='Structure', loc='upper left')
plt.tight_layout()

plt.savefig(output_plot_dir + '/boxplots/combined_boxplot_fractions_structure.png', dpi=300, bbox_inches='tight')
plt.close()
    