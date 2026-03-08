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
    parser = argparse.ArgumentParser(prog = 'python3 cell_fraction_analysis_core_wise.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Analyze shifts in cell fractions before and after treatment.  ')
    parser.add_argument('-i', help='path to combined adata file',
                        dest='input',
                        type=str)
    parser.add_argument('--phen_level', help='key for cell type annotation in adata.obs',
                        dest='phen_level',
                        type=str)
    parser.add_argument('--celltype_list', help='list of all cell types in the data',
                        dest='celltype_list',
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
cell_type_list = pd.read_csv(args.celltype_list, header=None, sep='\t')[0].tolist()
exclude_v17 = args.exclude_v17
print(f'Excluding v1.7 samples: {exclude_v17}')
output_dir_results = args.output_dir_results
output_dir_plots = args.output_dir_plots

# Make sure output directories exist
os.makedirs(output_dir_plots, exist_ok=True)
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
# Statistical testing functions
#--------------------------------------------------------------------------------
# Perform statistical testing for two groups of samples (stat_test can be mannwhitneyu - independent samples, wilcoxon - paired samples)
def stat_testing_two_groups(df, cell_cols, stat_test, group, groups):
    stat_results = []
    for celltype in cell_cols:
        group1_values = df[df[group] == groups[0]][celltype]
        group2_values = df[df[group] == groups[1]][celltype]
        group1_name, group2_name = groups[0], groups[1]

        stat, p_value = stat_test(group1_values, group2_values)
        stat_results.append({'cell_type': celltype, 'statistic': stat, 'p_value': p_value})

    stat_df = pd.DataFrame(stat_results)
    stat_df_annot = stat_df.rename(columns={"cell_type": "variable", "p_value": "pval"})
    stat_df_annot["group1"] = group1_name
    stat_df_annot["group2"] = group2_name
    stat_df_annot = stat_df_annot[["variable", "group1", "group2", "pval"]]
    
    return stat_df, stat_df_annot



# Centrality scores analysis and plotting - boxplot
#------------------------------------------------------------------------------
def stat_analysis_cell_fraction_box(input_file, output_dir_plots, output_dir_results, group, category, exclude_v17, immune, stat_test, cell_type_list, key):
    # Common setup
    groups = sorted([g for g in input_file[group].dropna().unique() if not str(g).isdigit()])
    id_vars = ['sample', group, category] if category else ['sample', group]
    scores_df = input_file[id_vars + cell_type_list].copy()
    scores_df_melted = scores_df.melt(id_vars=id_vars, var_name='cell_type', value_name=key)
    col_order=sorted([g for g in input_file[category].dropna().unique() if not str(g).isdigit()]) if category else None
    col = category if category else None

    # Create plot
    g = sns.catplot(scores_df_melted, x='cell_type', y=key, hue=group, hue_order=groups, col=col, col_order=col_order, kind='box', palette='tab20', height=6, aspect=1.5)
    
    # File naming helper
    suffix = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    cat_suffix = f'{category}' if category else ''
    imm_suffix = 'immune' if immune else ''
    base_filename = f'{output_dir_plots}cf_box_{imm_suffix}_{group}_{cat_suffix}_{stat_test.__name__}_{suffix}.svg'
    
    # Process each facet
    axes = g.axes.flat if category else [g.ax]
    facet_data = list(g.facet_data()) if category else [(None, input_file)]
    print(facet_data)
    for ax, (_, subdata) in zip(axes, facet_data):       
        if category:
            # make the subset_df in the same format as the one without category to be able to use the same function for stat testing and annotation, from melted format to wide format
            subset_df = subdata.pivot(index=id_vars, columns='cell_type', values=key).reset_index()
            subset_df_melted = subdata
        else:
            subset_df = subdata[id_vars + cell_type_list]
            subset_df_melted = subset_df.melt(id_vars=id_vars, value_vars=cell_type_list, var_name='cell_type', value_name=key)
        _, stat_df_annot = stat_testing_two_groups(subset_df, cell_type_list, stat_test, group, groups)
        
        # Annotate significant results
        alpha = 0.05
        sig_df = stat_df_annot[stat_df_annot["pval"] < alpha].copy().reset_index(drop=True)
        if not sig_df.empty:
            pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
            annot = Annotator(ax, pairs, data=subset_df_melted, x='cell_type', y=key, hue=group)
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df['pval'])
    
    # Set labels and title
    if category:
        title = f"{key} {imm_suffix} in {groups[0]} vs {groups[1]} split on {category}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        g.set_xticklabels(rotation=45, ha='right')
        g.set_xlabels('Cell type')
        g.set_ylabels(key)
        plt.suptitle(title, y=1.03)
    else:
        title = f"{key} {imm_suffix} in {groups[0]} vs {groups[1]}"
        if exclude_v17:
            title += " (excluding v1.7 treatment scheme)"
        title += f" ({stat_test.__name__})"
        plt.title(title)
        plt.xticks(rotation=45, ha='right')
        plt.xlabel('Cell type')
        plt.ylabel(key)
    
    g.legend.set_title(group)
    g.legend.set_loc('upper right')
    plt.tight_layout()
    plt.savefig(base_filename, format='svg', bbox_inches='tight')
    plt.close()



def celltype_fraction_composition_B_vs_tumorbed(df, output_dir, output_dir_results, category=None, exclude_v17=False, stat_test=mannwhitneyu, perform_stat_test=False, immune=False):
    #df_grouped = df.groupby(['sample_type', 'pt_id'], observed=True).mean(numeric_only=True).reset_index()
    #category_map = df[[category, 'pt_id']].drop_duplicates() if category in ['MPR', 'treatment'] else df[['pt_id']].drop_duplicates()
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])

    if category == None:       # Do not split into groups, compare biopsy vs resection (tumor_bed) for all patients
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        
        # Plot stripplot with lines connecting paired samples
        plt.figure(figsize=(12,6))
        ax = sns.stripplot(data = df_melted, x = 'variable', y = 'value', hue='sample_type',hue_order=['Biopsy', 'Resection'], dodge=True, jitter=False, size=4, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'})
        
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
        stat_df, stat_df_annot = paired_stat_testing(df, cell_fraction_cols, stat_test)
        if immune==False and exclude_v17==False:
            plt.title(f"Cell Type Fractions in Biopsy vs Resection (tumor border) ({stat_test.__name__})")
            file_name =  f'{output_dir}celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
        elif immune==True and exclude_v17==False:
            plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (tumor border) ({stat_test.__name__})")
            file_name = f'{output_dir}immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            plt.title(f"Cell Type Fractions in Biopsy vs Resection (tumor border) (excluding v1.7 treatment scheme) ({stat_test.__name__})")
            file_name = f'{output_dir}celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
        elif immune==True and exclude_v17==True:
            plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (tumor border) (excluding v1.7 treatment scheme) ({stat_test.__name__})")
            file_name = f'{output_dir}immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
        
        # Generate pairs for significant comparisons only
        alpha = 0.05
        sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
        pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
        annot = Annotator(ax, pairs, data=df_melted, x='variable', y='value', hue='sample_type')
        annot.configure(text_format="star")
        annot.set_pvalues_and_annotate(sig_df['pval'])
        
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(file_name, format='svg', bbox_inches='tight')
        plt.close()

    elif category != None:   # Split into groups based on chosen category
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        # Plot stripplot with lines connecting paired samples
        g = sns.catplot(data=df_melted, x="variable", y="value", hue="sample_type",hue_order=['Biopsy', 'Resection'], col=category, dodge=True, jitter=False, size=4, alpha=0.7, palette={'Biopsy':'gray', 'Resection':'black'}, kind='strip', height=6, aspect=1.5)
        
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
            subset_df = df[df[category]==cat_value]
            subset_df_melted = pd.melt(subset_df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
            subset_df_melted['variable'] = subset_df_melted['variable'].str.replace(' fraction','')
            #print(subset_df_melted)
            stat_df, stat_df_annot = paired_stat_testing(subset_df, cell_fraction_cols, stat_test)
            if immune==False and exclude_v17==False:
                file_name = f'{output_dir}{category}_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
                title = f"Cell Type Fractions in Biopsy vs Resection (tumor border) by {category} ({stat_test.__name__})"
            elif immune==True and exclude_v17==False:
                file_name = f'{output_dir}{category}_immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
                title = f"Immune Cell Type Fractions in Biopsy vs Resection (tumor border) by {category} ({stat_test.__name__})"
            elif immune==False and exclude_v17==True:
                file_name = f'{output_dir}{category}_celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
                title = f"Cell Type Fractions in Biopsy vs Resection (tumor border) by {category} (excluding v1.7 treatment scheme) ({stat_test.__name__})"
            elif immune==True and exclude_v17==True:
                file_name = f'{output_dir}{category}_immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
                title = f"Immune Cell Type Fractions in Biopsy vs Resection (tumor border) by {category} (excluding v1.7 treatment scheme) ({stat_test.__name__})"
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
        plt.suptitle(title, y=1.02)
        plt.tight_layout()
        plt.savefig(file_name, format='svg', bbox_inches='tight')
        plt.close()









#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create fractions dataframe
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate fractions per sample
fractions_df = pd.DataFrame(dtype=object)
for i, element in enumerate(adata.obs['sample'].unique().dropna()):
    adata_temp = adata[adata.obs['sample'] == element, :] # Subset adata for element in sample
    total_cells_temp = adata_temp.shape[0] # Total number of cells for this sample
    temp_fractions = adata_temp.obs[args.phen_level].value_counts()/total_cells_temp # Calculate fractions
    fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

    # Add metadata to the fractions_df
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme', 'MPR', 'treatment', 'structure', 'structure_core']
    for meta in meta_list: 
        fractions_df.loc[meta, element] = adata_temp.obs[meta].unique()[0]

# Adjust dataframe for plotting
fractions_df = fractions_df.T.fillna(0) # Transpose for easier plotting and fill NaNs with 0

res_df = fractions_df[fractions_df['sample_type']=='Resection']
print(f'Number of resection samples: {res_df.shape[0]}')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if exclude_v17 == True:
    categories = [None, 'MPR'] 
else:
    categories = [None, 'MPR', 'treatment']

structure_columns = ['structure', 'structure_core']
cores = res_df['structure_core'].unique()
structure_types = res_df['structure'].unique()

# All cell types
#--------------------------------------------------------------------------------
for structure_col in structure_columns:    
    for category in categories:
        print(f'Analyzing category: {category}')
        stat_analysis_cell_fraction_box(input_file=res_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=structure_col, category=category, exclude_v17=exclude_v17, immune=False, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key='Cell fraction')            
        if category != None:
            stat_analysis_cell_fraction_box(input_file=res_df, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=category, category=structure_col, exclude_v17=exclude_v17, immune=False, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key='Cell fraction')            



# Focus on immune cell types only
#--------------------------------------------------------------------------------
non_immune = ['Epithelial_cell', 'Fibroblast', 'Endothelial_cell', 'Pericyte', 'Stromal', 'Tumor_cells']
to_exclude = set(non_immune).intersection(res_df.columns)
print(f'Excluding non-immune cell types: {to_exclude}')
df_only_immune = res_df.drop(labels=to_exclude, axis=1)

df_only_immune_cells = df_only_immune[[col for col in df_only_immune.columns if col in cell_type_list]]  # only keep columns that have ' fraction' in their name
df_immune_fractions = df_only_immune_cells.div(df_only_immune_cells.sum(axis=1), axis=0)  # Re-normalize to sum to 1
df_immune_fractions[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']] = df_only_immune[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']].values # Add metadata back
df_immune_fractions['sample'] = df_immune_fractions.index

res_df_immune = df_immune_fractions[df_immune_fractions['sample_type']=='Resection']
print(f'Number of resection samples (immune only): {res_df_immune.shape[0]}')

cell_type_list = [cell_type for cell_type in cell_type_list if cell_type not in to_exclude]
for structure_col in structure_columns:    
    for category in categories:
        print(f'Analyzing category: {category}')
        stat_analysis_cell_fraction_box(input_file=res_df_immune, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=structure_col, category=category, exclude_v17=exclude_v17, immune=True, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key='Cell fraction')            
        if category != None:
            stat_analysis_cell_fraction_box(input_file=res_df_immune, output_dir_plots=output_dir_plots, output_dir_results=output_dir_results, group=category, category=structure_col, exclude_v17=exclude_v17, immune=True, stat_test=mannwhitneyu, cell_type_list=cell_type_list, key='Cell fraction')            



# # Compare Biopsy vs tumor_bed (core_1 and core_2)
# #--------------------------------------------------------------------------------
# # Calculate fractions per sample excluding core_3 (to focus on tumor_bed)
# adata_no_core3 = adata[adata.obs['structure_core'] != 'core_3', :]
# fractions_df = pd.DataFrame(dtype=object)
# for i, element in enumerate(adata_no_core3.obs['T_number'].unique().dropna()):
#     adata_temp = adata_no_core3[adata_no_core3.obs['T_number'] == element, :] # Subset adata for element in sample
#     total_cells_temp = adata_temp.shape[0] # Total number of cells for this sample
#     temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions
#     fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

#     # Add metadata to the fractions_df
#     meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme', 'MPR', 'treatment', 'structure', 'structure_core']
#     for meta in meta_list: 
#         fractions_df.loc[meta, element] = adata_temp.obs[meta].unique()[0]

# # Adjust dataframe for plotting
# fractions_df = fractions_df.T.fillna(0) # Transpose for easier plotting and fill NaNs with 0
# fractions_df.columns = [f'{col} fraction' if col not in meta_list else col for col in fractions_df.columns] # Add suffix to fraction columns

# # Keep only patients with matched biopsy and resection samples
# resection_pts = fractions_df[fractions_df['sample_type']=='Resection']['pt_id'].tolist()
# biopsy_pts = fractions_df[fractions_df['sample_type']=='Biopsy']['pt_id'].tolist()
# paired_pts = list(set(resection_pts) & set(biopsy_pts))
# paired_fractions_df = fractions_df[fractions_df['pt_id'].isin(paired_pts)]
# print(f'Number of paired patients: {len(paired_fractions_df["pt_id"].unique())}')
# print(f'Number of resection samples: {paired_fractions_df[paired_fractions_df["sample_type"]=="Resection"].shape[0]}')
# print(f'Number of biopsy samples: {paired_fractions_df[paired_fractions_df["sample_type"]=="Biopsy"].shape[0]}')
# print(paired_fractions_df[paired_fractions_df["sample_type"]=="Resection"])
# print(paired_fractions_df[paired_fractions_df["sample_type"]=="Biopsy"])

# to_exclude = set(non_immune).intersection(paired_fractions_df.columns)
# paired_immune_df = paired_fractions_df.drop(labels=to_exclude, axis=1)
# df_only_immune = paired_immune_df[paired_immune_df.columns[paired_immune_df.columns.str.contains(' fraction')]]  # only keep columns that have ' fraction' in their name
# final_df_immune = df_only_immune.div(df_only_immune.sum(axis=1), axis=0)  # Re-normalize to sum to 1
# final_df_immune[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']] = paired_fractions_df[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']].values # Add metadata back



# for category in categories:
#     celltype_fraction_composition_B_vs_tumorbed(paired_fractions_df, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=wilcoxon, perform_stat_test=True, immune=False)
#     celltype_fraction_composition_B_vs_tumorbed(final_df_immune, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=wilcoxon, perform_stat_test=True, immune=True)



