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
def celltype_fraction_composition_by_structure(df, output_dir, output_dir_results, exclude_v17, category=None, structure_col=None, stat_test=mannwhitneyu, perform_stat_test=False, immune=False):
    #print(df)
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    if structure_col is 'structure_core':
        order=['core_1', 'core_2', 'core_3']
    elif structure_col is 'structure':
        order = ['tumor_bed', 'TLS']
    else:
        order = None

    if category is None:
        # Plot all structures together
        plt.figure(figsize=(14, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', structure_col], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction', '')
        ax = sns.boxplot(data=df_melted, x="variable", y="value", hue=structure_col, hue_order=order, palette='tab20')
        
        title_prefix = "Immune Cell Type" if immune else "Cell Type"
        title_suffix = " (excluding v1.7 treatment scheme)" if exclude_v17 else ""
        plt.title(f"{title_prefix} Fractions by {structure_col}{title_suffix}")
        file_name = f'{output_dir}{"immune_" if immune else ""}celltype_fraction_by_{structure_col}{"_wo_v1.7" if exclude_v17 else "_w_v1.7"}.svg'
        
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        plt.legend(title=structure_col)
        plt.tight_layout()
        plt.savefig(file_name, format='svg', bbox_inches='tight')
        plt.close()
    
    else:
        # Plot split by both structure and category
        plt.figure(figsize=(16, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', structure_col, category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction', '')
        g = sns.catplot(data=df_melted, x="variable", y="value", hue=structure_col, hue_order=order, col=category, 
                        kind='box', palette='tab20', height=6, aspect=1.5)
        
        title_prefix = "Immune Cell Type" if immune else "Cell Type"
        title_suffix = " (excluding v1.7 treatment scheme)" if exclude_v17 else ""
        g.fig.suptitle(f"{title_prefix} Fractions by Structure and {category}{title_suffix}")
        
        for ax in g.axes.flat:
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
            ax.set_xlabel("Cell Type")
            ax.set_ylabel("Fraction")
        
        file_name = f'{output_dir}{category}_{"immune_" if immune else ""}celltype_fraction_by_{structure_col}{"_wo_v1.7" if exclude_v17 else "_w_v1.7"}.svg'
        g.legend.set_loc('upper right')
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


def celltype_fraction_composition_by_structure_within_structure_type(df, output_dir, output_dir_results, category=None, structure_col=None, structure_type=None, exclude_v17=False, stat_test=mannwhitneyu, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    plt.figure(figsize=(12, 6))
    df_melted = pd.melt(df, id_vars=['pt_id', category], value_vars=cell_fraction_cols)
    df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
    ax = sns.boxplot(data=df_melted, x="variable", y="value", hue=category, palette='tab20')

    if perform_stat_test == True:
        stat_df, stat_df_annot = ind_stat_testing(df, cell_fraction_cols, stat_test, category)
        if immune==False and exclude_v17==False:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_{structure_type}_statistical_results_w_v1.7.csv', index=False)
            plt.title(f"Cell Type Fractions in {structure_type} ({stat_test.__name__})") 
            file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{structure_type}_{stat_test.__name__}_w_v1.7.svg'
        elif immune==True and exclude_v17==False:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_{structure_type}_statistical_results_w_v1.7.csv', index=False)
            plt.title(f"Immune Cell Type Fractions in {structure_type} ({stat_test.__name__})")
            file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{structure_type}_{stat_test.__name__}_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_{structure_type}_statistical_results_wo_v1.7.csv', index=False)
            plt.title(f"Cell Type Fractions in {structure_type} (excluding v1.7 treatment scheme) ({stat_test.__name__})")
            file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{structure_type}_{stat_test.__name__}_wo_v1.7.svg'
        elif immune==True and exclude_v17==True:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_{structure_type}_statistical_results_wo_v1.7.csv', index=False)
            plt.title(f"Immune Cell Type Fractions in {structure_type} (excluding v1.7 treatment scheme) ({stat_test.__name__})") 
            file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{structure_type}_{stat_test.__name__}_wo_v1.7.svg'
    
        # Generate pairs for significant comparisons only
        alpha = 0.05
        sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
        if sig_df.empty:
            print(f"No significant results for category: {category} — skipping annotation.")
        else:
            pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
            annot = Annotator(ax,pairs,data=df_melted,x='variable', y='value', hue=category)
            annot.configure(text_format="star")
            annot.set_pvalues_and_annotate(sig_df['pval'])

    plt.xticks(rotation=45, ha='right')
    plt.xlabel("Cell Type")
    plt.ylabel("Fraction")
    plt.tight_layout()
    plt.savefig(file_name, format='svg')  
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
    temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions
    fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

    # Add metadata to the fractions_df
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme', 'MPR', 'treatment', 'structure', 'structure_core']
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
cores = ['core_1', 'core_2', 'core_3']
structure_types = ['tumor_bed', 'TLS']

# All cell types
#--------------------------------------------------------------------------------
for structure_col in structure_columns:    
    for category in categories:
        print(f'Analyzing category: {category}')
        celltype_fraction_composition_by_structure(res_df, output_dir, output_dir_results, category=category, structure_col=structure_col, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=False)            

for structure_col in structure_columns:
    if structure_col == 'structure_core':
        for core in cores:
            df_core = res_df[res_df['structure_core'] == core]
            for category in categories:
                if category != None:
                    print(f'Analyzing category: {category} for {core}')
                    celltype_fraction_composition_by_structure_within_structure_type(df_core, output_dir, output_dir_results, category=category, structure_col=structure_col, structure_type=core, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=False)
    elif structure_col == 'structure':
        for structure_type in structure_types:
            df_structure = res_df[res_df['structure'] == structure_type]
            for category in categories:
                if category != None:
                    print(f'Analyzing category: {category} for {structure_type}')
                    celltype_fraction_composition_by_structure_within_structure_type(df_structure, output_dir, output_dir_results, category=category, structure_col=structure_col, structure_type=structure_type, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=False)

# Focus on immune cell types only
#--------------------------------------------------------------------------------
non_immune = ['Epithelial cell fraction', 'Fibroblast fraction', 'Endothelial cell fraction', 'Pericyte fraction', 'Stromal fraction', 'Tumor cells fraction']
to_exclude = set(non_immune).intersection(paired_fractions_df.columns)
print(f'Excluding non-immune cell types: {to_exclude}')
df_only_immune = fractions_df.drop(labels=to_exclude, axis=1)

df_only_immune = df_only_immune[df_only_immune.columns[df_only_immune.columns.str.contains(' fraction')]]  # only keep columns that have ' fraction' in their name
df_immune = df_only_immune.div(df_only_immune.sum(axis=1), axis=0)  # Re-normalize to sum to 1
df_immune[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']] = fractions_df[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']].values # Add metadata back

res_df_immune = df_immune[df_immune['sample_type']=='Resection']
print(f'Number of resection samples (immune only): {res_df_immune.shape[0]}')

for structure_col in structure_columns:
    for category in categories:    
        celltype_fraction_composition_by_structure(res_df_immune, output_dir, output_dir_results, category=category, structure_col=structure_col, exclude_v17=exclude_v17, stat_test=mannwhitneyu, perform_stat_test=True, immune=True)



# Compare Biopsy vs tumor_bed (core_1 and core_2)
#--------------------------------------------------------------------------------
# Calculate fractions per sample excluding core_3 (to focus on tumor_bed)
adata_no_core3 = adata[adata.obs['structure_core'] != 'core_3', :]
fractions_df = pd.DataFrame(dtype=object)
for i, element in enumerate(adata_no_core3.obs['T_number'].unique().dropna()):
    adata_temp = adata_no_core3[adata_no_core3.obs['T_number'] == element, :] # Subset adata for element in sample
    total_cells_temp = adata_temp.shape[0] # Total number of cells for this sample
    temp_fractions = adata_temp.obs[celltype_key].value_counts()/total_cells_temp # Calculate fractions
    fractions_df = pd.concat([fractions_df, temp_fractions.rename(element)], axis=1) # Save fractions to df

    # Add metadata to the fractions_df
    meta_list = ['sample', 'pt_id', 'sample_type', 'disease_stage', 'T_number', 'regression', 'treatment_scheme', 'MPR', 'treatment', 'structure', 'structure_core']
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
print(f'Number of resection samples: {paired_fractions_df[paired_fractions_df["sample_type"]=="Resection"].shape[0]}')
print(f'Number of biopsy samples: {paired_fractions_df[paired_fractions_df["sample_type"]=="Biopsy"].shape[0]}')
print(paired_fractions_df[paired_fractions_df["sample_type"]=="Resection"])
print(paired_fractions_df[paired_fractions_df["sample_type"]=="Biopsy"])

to_exclude = set(non_immune).intersection(paired_fractions_df.columns)
paired_immune_df = paired_fractions_df.drop(labels=to_exclude, axis=1)
df_only_immune = paired_immune_df[paired_immune_df.columns[paired_immune_df.columns.str.contains(' fraction')]]  # only keep columns that have ' fraction' in their name
final_df_immune = df_only_immune.div(df_only_immune.sum(axis=1), axis=0)  # Re-normalize to sum to 1
final_df_immune[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']] = paired_fractions_df[['pt_id', 'sample_type', 'MPR', 'treatment', 'structure','structure_core']].values # Add metadata back



for category in categories:
    celltype_fraction_composition_B_vs_tumorbed(paired_fractions_df, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=wilcoxon, perform_stat_test=True, immune=False)
    celltype_fraction_composition_B_vs_tumorbed(final_df_immune, output_dir, output_dir_results, category=category, exclude_v17=exclude_v17, stat_test=wilcoxon, perform_stat_test=True, immune=True)



