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
from scipy.stats import wilcoxon, ttest_rel, ttest_ind
from statannotations.Annotator import Annotator

import warnings
warnings.filterwarnings("ignore")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
output_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/celltype_fraction/'
exclude_v17 = False
# Set aesthetics
sns.set_style("whitegrid")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define functions for analyses
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Analyse shifts in cell fractions before and after treatment
def celltype_fraction_shifts_lineplot(df, output_dir, category, exclude_v17, stat_test = None, perform_stat_test = False, immune = False, ):
    # Split data into pre- and post-treatment
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])

    if category == None:       # Do not split into groups, compare biopsy vs resection for all patients
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        
        # Plot stripplot with lines connecting paired samples
        ax = sns.stripplot(data = df_melted, x = 'variable', y = 'value', hue='sample_type', dodge=True, jitter=False, size=4, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'})

        # Prepare the data for line plotting
        wide = df_melted.pivot_table(index='pt_id', columns=['variable', 'sample_type'], values='value')
        categories = df_melted['variable'].unique()
        xticks = ax.get_xticks()
        x_map = dict(zip(categories, xticks))
        offset = 0.18 # offset for biopsy vs resection points

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
        
        # Perform statistical testing if specified
        if perform_stat_test==True:
            stat_df_annot = paired_stat_testing(df, cell_fraction_cols, output_dir, stat_test)

            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            print(sig_df)
            pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
            annot = Annotator(ax,pairs,data=df_melted,x='variable', y='value',hue='sample_type')
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df['pval'])
        
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fractions in Biopsy vs Resection") 
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fractions in Biopsy vs Resection")
        else:
            plt.title("Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(f'{output_dir}celltype_fraction_shifts_lineplot.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}immune_celltype_fraction_shifts_lineplot.svg', format='svg')
        plt.close()

    elif category != None:   # Split into groups based on chosen category
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')

        # Plot stripplot with lines connecting paired samples
        g = sns.catplot(data=df_melted, x="variable", y="value", hue="sample_type", col=category, dodge=True, jitter=False, size=4, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'}, kind='strip')
        
        # Loop over each axis to add lines
        for ax, (facet_key, subdata) in zip(g.axes.flat, g.facet_data()): 
            cat_value = g.col_names[facet_key[1]]
            print(f'Processing category: {cat_value}')
            # Prepare the data for line plotting
            wide = subdata.pivot_table(index='pt_id', columns=['variable', 'sample_type'], values='value')
            # x positions of categorical axis
            categories = subdata['variable'].unique()
            xticks = ax.get_xticks()
            x_map = dict(zip(categories, xticks))
            offset = 0.18 # offset for biopsy vs resection points

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
        
            # Perform statistical testing if specified
            if perform_stat_test==True:
                subset_df = df[df[category]==cat_value]
                subset_df_melted = pd.melt(subset_df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
                subset_df_melted['variable'] = subset_df_melted['variable'].str.replace(' fraction','')
                print(subset_df_melted)
                stat_df_annot = paired_stat_testing(subset_df, cell_fraction_cols, output_dir, stat_test)

                # Generate pairs for significant comparisons only
                alpha = 0.05
                sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
                if sig_df.empty:
                    print(f"No significant results for category: {cat_value} — skipping annotation.")
                    continue

                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax, pairs, data=subset_df_melted, x='variable', y='value', hue='sample_type')
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])
        
        #sns.move_legend(g, "upper right", title='Sample Type')
        g.set_xticklabels(rotation=45, ha='right')
        g.set_xlabels("Cell Type")
        g.set_ylabels("Fraction")
        g.legend.set_title('Sample Type')
        g.legend.set_loc('upper right')
        plt.tight_layout()
        plt.savefig(f'{output_dir}{category}_celltype_fraction_shifts_lineplot.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}{category}_immune_celltype_fraction_shifts_lineplot.svg', format='svg')
        plt.close()

# Perform statistical testing
def paired_stat_testing(df, cell_fraction_cols, output_dir, exclude_v17, stat_test, immune=False):
    # Perform statistical test for each cell type
    stat_results = []
    for celltype in cell_fraction_cols:
        biopsy_values = df[df['sample_type']=='Biopsy'][celltype]
        resection_values = df[df['sample_type']=='Resection'][celltype]
        # Ensure paired samples
        stat, p_value = stat_test(biopsy_values, resection_values)
        stat_results.append({'cell_type': celltype.replace(' fraction',''), 'statistic': stat, 'p_value': p_value})

    stat_df = pd.DataFrame(stat_results)
    #stat_df.to_csv(f'{output_dir}/results/analysis/celltype_fraction/paired_celltype_fraction_statistical_results.csv', index=False) if immune==False else stat_df.to_csv(f'{output_dir}/results/analysis/celltype_fraction/paired_immune_celltype_fraction_statistical_results.csv', index=False)
    print(stat_df)

    # Prepare stat_df for Annotator (expected format to be able to draw asterisks on plot)
    stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
    stat_df_annot["group1"] = "Biopsy"
    stat_df_annot["group2"] = "Resection"
    stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]] # Reorder columns to expected format
    return stat_df_annot



def celltype_fraction_composition_box(df, output_dir, exclude_v17, category=None, stat_test=None, perform_stat_test=False, immune=False):
    print(df)
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    
    if category==None:
        plt.figure(figsize=(12, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        print(df_melted)
        sns.boxplot(data=df_melted, x="variable", y="value", hue="sample_type", palette='tab20', fill=True, gap=0.2)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fractions in Biopsy vs Resection") 
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fractions in Biopsy vs Resection")
        else:
            plt.title("Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(f'{output_dir}/celltype_fraction_composition_box.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}/immune_celltype_fraction_composition_box.svg', format='svg')
        plt.close()
        
    elif category != None:
        plt.figure(figsize=(12, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        g = sns.catplot(data=df_melted, x="variable", y="value", hue="sample_type", col=category, kind='box', palette='tab20')
        sns.move_legend(g, "upper right", title='Sample Type')
        g.set_xticklabels(rotation=45, ha='right')
        g.set_xlabels("Cell Type")
        g.set_ylabels("Fraction")
        plt.tight_layout()
        plt.savefig(f'{output_dir}{category}_celltype_fraction_composition_box.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}{category}_immune_celltype_fraction_composition_box.svg', format='svg')
        plt.close()

def celltype_fraction_shifts_box(df, output_dir, exclude_v17, category=None, stat_test=ttest_ind, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    biopsy_df = df[df['sample_type']=='Biopsy']
    biopsy_fractions = biopsy_df[cell_fraction_cols].set_index(biopsy_df['pt_id'])
    resection_df = df[df['sample_type']=='Resection']
    resection_fractions = resection_df[cell_fraction_cols].set_index(resection_df['pt_id']).reindex(biopsy_fractions.index)  # Align indices
    print(biopsy_fractions.head())
    print(resection_fractions.head())

    if category == None:
        # Calculate the difference between resection and biopsy for each pair
        diff_df = resection_fractions-biopsy_fractions #.values.set_index(resection_fractions.index)
        diff_df.columns = diff_df.columns.str.replace(' fraction','')
        plt.figure(figsize=(12, 6))
        sns.boxplot(diff_df)
        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection") 
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection")
        else:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Shift in Fraction (Resection - Biopsy)")
        plt.tight_layout()
        plt.savefig(f'{output_dir}celltype_fraction_shift_box.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}immune_celltype_fraction_shift_box.svg', format='svg')
        plt.close()

    elif category != None:
        # Calculate the difference between resection and biopsy for each pair
        diff_df = resection_fractions-biopsy_fractions 
        diff_df[category] = resection_df[category].values
        print(diff_df.head())
        diff_df_melted = pd.melt(diff_df, id_vars=[category], value_vars=cell_fraction_cols)
        diff_df_melted['variable'] = diff_df_melted['variable'].str.replace(' fraction','')

        # ensure consistent (alphabetical) order of categories
        cat_order = sorted(diff_df_melted[category].dropna().unique())
        diff_df_melted[category] = pd.Categorical(diff_df_melted[category],categories=cat_order, ordered=True)

        plt.figure(figsize=(12, 6))
        sns.boxplot(data=diff_df_melted, x="variable", y="value", hue=category, palette='tab20')
        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection") 
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection")
        else:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Shift in Fraction (Resection - Biopsy)")
        plt.tight_layout()
        plt.savefig(f'{output_dir}{category}_celltype_fraction_shift_box.svg', format='svg') if immune==False else plt.savefig(f'{output_dir}{category}_immune_celltype_fraction_shift_box.svg', format='svg')
        plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Create fractions dataframe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove cells coming from samples treated with treatment_scheme=v1.7
if exclude_v17 == True:
    print('Removing samples with treatment_scheme=v1.7...')
    adata = adata[adata.obs['treatment_scheme'] != 'v1.7', :]
    output_dir = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/celltype_fraction_wo_v1.7/'
    print(f'Number of cells after removing v1.7 treatment scheme samples: {adata.shape[0]}')

# Calculate fractions per T_number
fractions_df = pd.DataFrame(dtype=object)
for i, element in enumerate(adata.obs['T_number'].unique().dropna()):
    #print(f'Processing {i}. T_number: {element}')
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
print(paired_fractions_df.head())

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
categories = [None,'MPR', 'treatment']
for category in categories:
    print(f'Analyzing category: {category}')
    celltype_fraction_shifts_lineplot(paired_fractions_df, output_dir, category=category, stat_test=wilcoxon, perform_stat_test=False, immune=False, exclude_v17=exclude_v17)
    celltype_fraction_composition_box(fractions_df, category = category, exclude_v17=exclude_v17, output_dir=output_dir, immune=False)#, stat_test = ttest_rel, perform_stat_test=False,exclude_v17=exclude_v17)
    celltype_fraction_shifts_box(paired_fractions_df, output_dir, category=category, stat_test=wilcoxon, perform_stat_test=False, immune=False, exclude_v17=exclude_v17)


# Focus on immune cell types only
non_immune = ['Epithelial cell fraction', 'Fibroblast fraction', 'Endothelial cell fraction', 'Pericyte fraction', 'Stromal fraction', 'Tumor cells fraction']
to_exclude = set(non_immune).intersection(paired_fractions_df.columns)
print(f'Excluding non-immune cell types: {to_exclude}')
df_only_immune = paired_fractions_df.drop(labels=to_exclude, axis=1)
    
df_only_immune = df_only_immune[df_only_immune.columns[df_only_immune.columns.str.contains(' fraction')]]  # only keep columns that have ' fraction' in their name
df_immune = df_only_immune.div(df_only_immune.sum(axis=1), axis=0)  # Re-normalize to sum to 1
df_immune[['pt_id', 'sample_type', 'MPR', 'treatment']] = paired_fractions_df[['pt_id', 'sample_type', 'MPR', 'treatment']].values # Add metadata back
print(df_immune.head())

for category in categories:    
    celltype_fraction_shifts_lineplot(df_immune, output_dir, category=category, stat_test=wilcoxon, perform_stat_test=False, immune=True,exclude_v17=exclude_v17)
    celltype_fraction_composition_box(df_immune, category = category, output_dir=output_dir, immune=True, exclude_v17=exclude_v17)#, stat_test = ttest_rel, perform_stat_test=False)
    celltype_fraction_shifts_box(df_immune, output_dir, category=category, stat_test=wilcoxon, perform_stat_test=False, immune=True, exclude_v17=exclude_v17)


