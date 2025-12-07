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
os.makedirs(output_plot_dir, exist_ok=True)
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
resection_pts = fractions_df[fractions_df['sample_type']=='Resection']['pt_id'].tolist()
biopsy_pts = fractions_df[fractions_df['sample_type']=='Biopsy']['pt_id'].tolist()
paired_pts = list(set(resection_pts) & set(biopsy_pts))
paired_fractions_df = fractions_df[fractions_df['pt_id'].isin(paired_pts)]
print(f'Number of paired patients: {len(paired_fractions_df["pt_id"].unique())}')
print(paired_fractions_df)

# Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyse shifts in cell fractions before and after treatment
def celltype_fraction_shifts(df, category, output_dir):
    # Split data into pre- and post-treatment
    #biopsy_df = df[df['sample_type']=='biopsy']
    #resection_df = df[df['sample_type']=='resection']
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])

    if category == None:       # Do not split into groups, compare biopsy vs resection for all patients
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        
        # Plot stripplot with lines connecting paired samples
        plt.figure(figsize=(12, 6))
        ax = sns.stripplot(data = df_melted, x = 'variable', y = 'value', hue='sample_type', dodge=True, jitter=False, size=7, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'})

        # Prepare the data for line plotting
        wide = df_melted.pivot_table(index='pt_id', columns=['variable', 'sample_type'], values='value')
        print(wide)

        # x positions of categorical axis
        categories = df_melted['variable'].unique()
        xticks = ax.get_xticks()
        x_map = dict(zip(categories, xticks))

        # offset for biopsy vs resection points
        offset = 0.18

        # Draw lines connecting paired samples
        for celltype in categories:
            sub = wide[celltype].dropna()
            for _, row in sub.iterrows():
                x_left = x_map[celltype] - offset   # biopsy x-position
                x_right = x_map[celltype] + offset  # resection x-position
                y_bio = row['Biopsy']
                y_res = row['Resection']

                color = 'blue' if y_res > y_bio else 'red'
                ax.plot([x_left, x_right],[y_bio, y_res],color=color,linewidth=1,alpha=0.8)


        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        plt.title("Cell Type Fractions in Biopsy vs Resection (per sample)")
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/celltype_fraction_shifts.svg', format='svg')



#def celltype_fraction_shifts_immune(df, category, output_dir):
    # Focus on immune cell types only
    immune_cell_types = [col for col in df.columns if col.endswith('fraction') and col not in ['Tuft cells fraction', 'Pericytes fraction', 'Epithelial cells fraction']]
    if category == None:       # Do not split into groups, compare biopsy vs resection for all patients
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=immune_cell_types)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        print(df_melted)
        plt.figure(figsize=(12, 6))
        sns.stripplot(data = df_melted, x = 'variable', y = 'value', hue='sample_type', dodge=True, jitter=False, size=5, alpha=0.7, palette='tab20')
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Immune Cell Type")
        plt.ylabel("Fraction")
        plt.title("Immune Cell Type Fractions in Biopsy vs Resection (per sample)")
        plt.legend()
        plt.tight_layout()
        plt.savefig(f'{output_dir}/immune_celltype_fraction_shifts.svg', format='svg')


celltype_fraction_shifts(paired_fractions_df, None, output_plot_dir)
#celltype_fraction_shifts_immune(paired_fractions_df, None, output_plot_dir)
# 2. Swarmplots of fractions per chosen category (e.g., structure) with paired pts connected
#plot_swarmplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)
# 3. Lineplots of fractions per chosen category (e.g., structure) with paired pts connected
#plot_lineplot(adata_sample_paired, fraction_columns, 'structure', output_plot_dir)