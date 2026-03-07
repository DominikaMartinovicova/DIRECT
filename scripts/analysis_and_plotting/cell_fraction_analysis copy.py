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
    parser.add_argument('--celltype_list', help='list of all cell types in the data',
                        dest='celltype_list',
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
# Analyse shifts in cell fractions before and after treatment
#--------------------------------------------------------------------------------
def celltype_fraction_shifts_lineplot(df, output_dir, output_dir_results, category, exclude_v17, stat_test, perform_stat_test = False, immune = False):
    # Split data into pre- and post-treatment
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])

    if category == None:       # Do not split into groups, compare biopsy vs resection for all patients
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
        if perform_stat_test==True:
            stat_df, stat_df_annot = paired_stat_testing(df, cell_fraction_cols, stat_test)
            if immune==False and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                file_name =  f'{output_dir}celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
            elif immune==True and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_immune_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                file_name = f'{output_dir}immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
            elif immune==False and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir}celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
            elif immune==True and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_paired_immune_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
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
            if perform_stat_test==True:
                subset_df = df[df[category]==cat_value]
                subset_df_melted = pd.melt(subset_df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
                subset_df_melted['variable'] = subset_df_melted['variable'].str.replace(' fraction','')
                #print(subset_df_melted)
                stat_df, stat_df_annot = paired_stat_testing(subset_df, cell_fraction_cols, stat_test)
                if immune==False and exclude_v17==False:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_paired_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
                    title = f"Cell Type Fractions in Biopsy vs Resection ({category}) ({stat_test.__name__})"
                elif immune==True and exclude_v17==False:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_paired_immune_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_w_v1.7.svg'
                    title = f"Immune Cell Type Fractions in Biopsy vs Resection ({category}) ({stat_test.__name__})"
                elif immune==False and exclude_v17==True:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_paired_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
                    title = f"Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({category}) ({stat_test.__name__})"
                elif immune==True and exclude_v17==True:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_paired_immune_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_immune_celltype_fraction_shifts_lineplot_{stat_test.__name__}_wo_v1.7.svg'
                    title = f"Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({category}) ({stat_test.__name__})"
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


# Analyze cell type fraction compositions in biopsy vs resection, possibly split into groups based on chosen category
#--------------------------------------------------------------------------------
def celltype_fraction_composition_box(df, output_dir, output_dir_results, exclude_v17, category=None, stat_test=None, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    
    if category==None:
        plt.figure(figsize=(12, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type'], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        ax=sns.boxplot(data=df_melted, x="variable", y="value", hue="sample_type",hue_order=['Biopsy', 'Resection'], palette='tab20', fill=True)

        if perform_stat_test == True:
            stat_df, stat_df_annot = ind_stat_testing(df, cell_fraction_cols, stat_test)
            if immune==False and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                file_name = f'{output_dir}celltype_fraction_composition_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==True and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                file_name = f'{output_dir}immune_celltype_fraction_composition_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==False and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir}celltype_fraction_composition_box_{stat_test.__name__}_wo_v1.7.svg'
            elif immune==True and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir}immune_celltype_fraction_composition_box_{stat_test.__name__}_wo_v1.7.svg'
        
            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            if sig_df.empty:
                print(f"No significant results for category: {cat_value} — skipping annotation.")
            else:
                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax,pairs,data=df_melted,x='variable', y='value',hue='sample_type')
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])
        
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(file_name, format='svg')
        
        # Save zoomed-in version
        ax.set_ylim(-0.01, 0.1)
        zoomed_file_name = file_name.replace('.svg', '_zoom_0_0.1.svg')
        plt.savefig(zoomed_file_name, format='svg')
        plt.close()

    elif category != None:
        plt.figure(figsize=(20, 10))
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        g = sns.catplot(data=df_melted, x="variable", y="value", hue="sample_type", hue_order=['Biopsy', 'Resection'], col=category, kind='box', palette='tab20', height=6, aspect=1.5)

        # Loop over each axis to add lines
        for ax, (facet_key, subdata) in zip(g.axes.flat, g.facet_data()): 
            cat_value = g.col_names[facet_key[1]]
            print(f'Processing category: {cat_value}')
        
            # Perform statistical testing if specified
            if perform_stat_test==True:
                subset_df = df[df[category]==cat_value]
                subset_df_melted = pd.melt(subset_df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
                subset_df_melted['variable'] = subset_df_melted['variable'].str.replace(' fraction','')
                stat_df, stat_df_annot = paired_stat_testing(subset_df, cell_fraction_cols, stat_test)
                
                if immune==False and exclude_v17==False:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{stat_test.__name__}_w_v1.7.svg'
                    #plt.title(f"Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                elif immune==True and exclude_v17==False:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_immune_celltype_fraction_statistical_results_w_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{stat_test.__name__}_w_v1.7.svg'
                    #plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection ({stat_test.__name__})")
                elif immune==False and exclude_v17==True:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{stat_test.__name__}_wo_v1.7.svg'
                    #plt.title(f"Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                elif immune==True and exclude_v17==True:
                    stat_df.to_csv(f'{output_dir_results}/{category}_{stat_test.__name__}_immune_celltype_fraction_statistical_results_wo_v1.7.csv', index=False)
                    file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{stat_test.__name__}_wo_v1.7.svg'
                    #plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                
                # Generate pairs for significant comparisons only
                alpha = 0.05
                sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
                if sig_df.empty:
                    print(f"No significant results for category: {cat_value} — skipping annotation.")
                else:
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
        plt.savefig(file_name, format='svg')

        # Save zoomed-in version
        ax.set_ylim(-0.01, 0.1)
        zoomed_file_name = file_name.replace('.svg', '_zoom_0_0.1.svg')
        plt.savefig(zoomed_file_name, format='svg')
        plt.close()

        # Plot compositions of cell type fractions as boxplots for all four categories into one plot
        # B-<90, B->=90, R-<90, R->=90
        print('Plotting combined boxplots for all categories...')
        plt.figure(figsize=(15, 6))
        df_melted = pd.melt(df, id_vars=['pt_id', 'sample_type', category], value_vars=cell_fraction_cols)
        df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
        df_melted[f'sample_type-{category}'] = df_melted['sample_type'] + '-' + df_melted[category].astype(str)
        ax=sns.boxplot(data=df_melted, x="variable", y="value", hue= f"sample_type-{category}", 
                        hue_order=sorted(df_melted[f"sample_type-{category}"].unique()), 
                        palette='tab20', fill=True, gap=0.3)
            
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Fraction")
        if immune==False and exclude_v17==False:
            plt.title(f"Cell Type Fractions in Biopsy vs Resection")
            file_name = f'{output_dir}{category}_celltype_fraction_combined_composition_box_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            plt.title(f"Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}{category}_celltype_fraction_combined_composition_box_wo_v1.7.svg'
        elif immune==True and exclude_v17==False:
            plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection")
            file_name = f'{output_dir}{category}_immune_celltype_fraction_combined_composition_box_w_v1.7.svg'
        else:
            plt.title(f"Immune Cell Type Fractions in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}{category}_immune_celltype_fraction_combined_composition_box_wo_v1.7.svg'
        plt.legend(title='Sample Type')
        plt.tight_layout()
        plt.savefig(file_name, format='svg')

        # Save zoomed-in version
        ax.set_ylim(-0.01, 0.1)
        zoomed_file_name = file_name.replace('.svg', '_zoom_0_0.1.svg')
        plt.savefig(zoomed_file_name, format='svg')
        plt.close()

# Analyze shifts in cell type fractions before and after treatment, possibly split into groups based on chosen category, as boxplots
#--------------------------------------------------------------------------------
def celltype_fraction_shifts_box(df, output_dir, output_dir_results, exclude_v17, category=None, stat_test=ttest_ind, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    biopsy_df = df[df['sample_type']=='Biopsy']
    biopsy_fractions = biopsy_df[cell_fraction_cols].set_index(biopsy_df['pt_id'])
    resection_df = df[df['sample_type']=='Resection']
    resection_fractions = resection_df[cell_fraction_cols].set_index(resection_df['pt_id']).reindex(biopsy_fractions.index)  # Align indices

    if category == None:
        # Calculate the difference between resection and biopsy for each pair
        diff_df = resection_fractions-biopsy_fractions #.values.set_index(resection_fractions.index)
        diff_df.columns = diff_df.columns.str.replace(' fraction','')
        plt.figure(figsize=(12, 6))
        sns.boxplot(diff_df)

        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection")
            file_name = f'{output_dir}celltype_fraction_shift_box_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}celltype_fraction_shift_box_wo_v1.7.svg'
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection")
            file_name = f'{output_dir}immune_celltype_fraction_shift_box_w_v1.7.svg'
        else:
            plt.title("Immune Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}immune_celltype_fraction_shift_box_wo_v1.7.svg'
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Shift in Fraction (Resection - Biopsy)")
        plt.tight_layout()
        plt.savefig(file_name, format='svg')
        plt.close()

    elif category != None:
        # Calculate the difference between resection and biopsy for each pair
        diff_df = resection_fractions-biopsy_fractions 
        diff_df[category] = diff_df.index.map(biopsy_df.set_index('pt_id')[category])  # Add category column based on biopsy samples, match this on pt_id
        diff_df_melted = pd.melt(diff_df, id_vars=[category], value_vars=cell_fraction_cols)
        diff_df_melted['variable'] = diff_df_melted['variable'].str.replace(' fraction','')

        # ensure consistent (alphabetical) order of categories
        cat_order = sorted(diff_df_melted[category].dropna().unique())
        diff_df_melted[category] = pd.Categorical(diff_df_melted[category],categories=cat_order, ordered=True)

        plt.figure(figsize=(12, 6))
        ax=sns.boxplot(data=diff_df_melted, x="variable", y="value", hue=category, palette='tab20')
        
        if perform_stat_test == True:
            stat_df, stat_df_annot = ind_stat_testing(diff_df, cell_fraction_cols, stat_test, category)
            if immune==False and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_shift_statistical_results_w_v1.7.csv', index=False) 
                plt.title(f"Cell Type Fraction Shift in Biopsy vs Resection, ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_celltype_fraction_shift_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==True and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_shift_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fraction Shift in Biopsy vs Resection, ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_immune_celltype_fraction_shift_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==False and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_shift_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme), ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_celltype_fraction_shift_box_{stat_test.__name__}_wo_v1.7.svg'
            elif immune==True and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_shift_statistical_results_wo_v1.7.csv', index=False) 
                plt.title(f"Immune Cell Type Fraction Shift in Biopsy vs Resection (excluding v1.7 treatment scheme), ({stat_test.__name__})") 
                file_name = f'{output_dir}{category}_immune_celltype_fraction_shift_box_{stat_test.__name__}_wo_v1.7.svg'
        
            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            if sig_df.empty:
                print(f"No significant results for category: {category} — skipping annotation.")
            else:
                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax,pairs,data=diff_df_melted,x='variable', y='value', hue=category)
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])

        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Shift in Fraction (Resection - Biopsy)")
        plt.tight_layout()
        plt.savefig(file_name, format='svg')
        plt.close()

# Calculate fold change between resection and biopsy for each pair
def celltype_fraction_shifts_foldchange(df, output_dir, output_dir_results, exclude_v17, category=None, stat_test=ttest_ind, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    biopsy_df = df[df['sample_type']=='Biopsy']
    biopsy_fractions = biopsy_df[cell_fraction_cols].set_index(biopsy_df['pt_id'])
    resection_df = df[df['sample_type']=='Resection']
    resection_fractions = resection_df[cell_fraction_cols].set_index(resection_df['pt_id']).reindex(biopsy_fractions.index)

    if category == None:
        # Calculate fold change: resection / biopsy (add pseudocount to avoid division by zero)
        fc_df = (resection_fractions + 1e-6) / (biopsy_fractions + 1e-6)
        fc_df.columns = fc_df.columns.str.replace(' fraction','')
        fc_df_log = np.log2(fc_df)  # Log2 transform for better visualization
        plt.figure(figsize=(12, 6))
        sns.boxplot(fc_df_log)

        if immune==False and exclude_v17==False:
            plt.title("Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection")
            file_name = f'{output_dir}celltype_fraction_foldchange_box_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            plt.title("Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}celltype_fraction_foldchange_box_wo_v1.7.svg'
        elif immune==True and exclude_v17==False:
            plt.title("Immune Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection")
            file_name = f'{output_dir}immune_celltype_fraction_foldchange_box_w_v1.7.svg'
        else:
            plt.title("Immune Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection (excluding v1.7 treatment scheme)")
            file_name = f'{output_dir}immune_celltype_fraction_foldchange_box_wo_v1.7.svg'
        plt.axhline(y=0, color='red', linestyle='--', linewidth=1, label='No change')
        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Log2 Fold Change (Resection / Biopsy)")
        plt.tight_layout()
        plt.savefig(file_name, format='svg')
        plt.close()

    elif category != None:
        # Calculate fold change
        fc_df = (resection_fractions + 1e-6) / (biopsy_fractions + 1e-6)
        fc_df[category] = fc_df.index.map(biopsy_df.set_index('pt_id')[category])
        fc_df_log = np.log2(fc_df[cell_fraction_cols])
        fc_df_log[category] = fc_df[category]
        fc_df_melted = pd.melt(fc_df_log, id_vars=[category], value_vars=cell_fraction_cols)
        fc_df_melted['variable'] = fc_df_melted['variable'].str.replace(' fraction','')

        cat_order = sorted(fc_df_melted[category].dropna().unique())
        fc_df_melted[category] = pd.Categorical(fc_df_melted[category], categories=cat_order, ordered=True)

        plt.figure(figsize=(12, 6))
        ax = sns.boxplot(data=fc_df_melted, x="variable", y="value", hue=category, palette='tab20')
        plt.axhline(y=0, color='red', linestyle='--', linewidth=1, label='No change')

        if perform_stat_test == True:
            stat_df, stat_df_annot = ind_stat_testing(fc_df_log, cell_fraction_cols, stat_test, category)
            if immune==False and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_foldchange_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_celltype_fraction_foldchange_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==True and exclude_v17==False:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_foldchange_statistical_results_w_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_immune_celltype_fraction_foldchange_box_{stat_test.__name__}_w_v1.7.svg'
            elif immune==False and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_foldchange_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_celltype_fraction_foldchange_box_{stat_test.__name__}_wo_v1.7.svg'
            elif immune==True and exclude_v17==True:
                stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_foldchange_statistical_results_wo_v1.7.csv', index=False)
                plt.title(f"Immune Cell Type Fraction Fold Change (Log2) in Biopsy vs Resection (excluding v1.7 treatment scheme) ({stat_test.__name__})")
                file_name = f'{output_dir}{category}_immune_celltype_fraction_foldchange_box_{stat_test.__name__}_wo_v1.7.svg'
            
            # Generate pairs for significant comparisons only
            alpha = 0.05
            sig_df = stat_df_annot[stat_df_annot["pval"] < alpha ].copy().reset_index(drop=True)
            if sig_df.empty:
                print(f"No significant results for category: {category} — skipping annotation.")
            else:
                pairs = [((row.variable, row.group1), (row.variable, row.group2)) for _, row in sig_df.iterrows()]
                annot = Annotator(ax, pairs, data=fc_df_melted, x='variable', y='value', hue=category)
                annot.configure(text_format="star")
                annot.set_pvalues_and_annotate(sig_df['pval'])

        plt.xticks(rotation=45, ha='right')
        plt.xlabel("Cell Type")
        plt.ylabel("Log2 Fold Change (Resection / Biopsy)")
        plt.tight_layout()
        plt.savefig(file_name, format='svg')
        plt.close()


# Analyze compositions of cell type fractions in biopsy/resection, split into groups based on chosen category, as boxplots
#--------------------------------------------------------------------------------
def composition_within_sampletype_box(df, output_dir, output_dir_results, exclude_v17, category, sample_type, stat_test=mannwhitneyu, perform_stat_test=False, immune=False):
    cell_fraction_cols = sorted([col for col in df.columns if col.endswith('fraction')])
    plt.figure(figsize=(12, 6))
    df_melted = pd.melt(df, id_vars=['pt_id', category], value_vars=cell_fraction_cols)
    df_melted['variable'] = df_melted['variable'].str.replace(' fraction','')
    ax = sns.boxplot(data=df_melted, x="variable", y="value", hue=category, palette='tab20')

    if perform_stat_test == True:
        stat_df, stat_df_annot = ind_stat_testing(df, cell_fraction_cols, stat_test, category)
        if immune==False and exclude_v17==False:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_{sample_type}_statistical_results_w_v1.7.csv', index=False)
            plt.title(f"Cell Type Fractions in {sample_type} ({stat_test.__name__})") 
            file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{sample_type}_{stat_test.__name__}_w_v1.7.svg'
        elif immune==True and exclude_v17==False:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_{sample_type}_statistical_results_w_v1.7.csv', index=False)
            plt.title(f"Immune Cell Type Fractions in {sample_type} ({stat_test.__name__})")
            file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{sample_type}_{stat_test.__name__}_w_v1.7.svg'
        elif immune==False and exclude_v17==True:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_celltype_fraction_{sample_type}_statistical_results_wo_v1.7.csv', index=False)
            plt.title(f"Cell Type Fractions in {sample_type} (excluding v1.7 treatment scheme) ({stat_test.__name__})")
            file_name = f'{output_dir}{category}_celltype_fraction_composition_box_{sample_type}_{stat_test.__name__}_wo_v1.7.svg'
        elif immune==True and exclude_v17==True:
            stat_df.to_csv(f'{output_dir_results}/{stat_test.__name__}_immune_celltype_fraction_{sample_type}_statistical_results_wo_v1.7.csv', index=False)
            plt.title(f"Immune Cell Type Fractions in {sample_type} (excluding v1.7 treatment scheme) ({stat_test.__name__})") 
            file_name = f'{output_dir}{category}_immune_celltype_fraction_composition_box_{sample_type}_{stat_test.__name__}_wo_v1.7.svg'
    
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

    # Save zoomed-in version
    ax.set_ylim(-0.01, 0.1)
    zoomed_file_name = file_name.replace('.svg', '_zoom_0_0.1.svg')
    plt.savefig(zoomed_file_name, format='svg')
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
    celltype_fraction_shifts_lineplot(paired_fractions_df, output_dir, output_dir_results, category=category, stat_test=wilcoxon, perform_stat_test=True, immune=False, exclude_v17=exclude_v17)
    celltype_fraction_composition_box(fractions_df, output_dir, output_dir_results, category = category, exclude_v17=exclude_v17, immune=False, stat_test = mannwhitneyu, perform_stat_test=True)
    celltype_fraction_shifts_box(paired_fractions_df, output_dir, output_dir_results, category=category, stat_test=mannwhitneyu, perform_stat_test=True, immune=False, exclude_v17=exclude_v17)
    celltype_fraction_shifts_foldchange(paired_fractions_df, output_dir, output_dir_results, category=category, stat_test=mannwhitneyu, perform_stat_test=True, immune=False, exclude_v17=exclude_v17)
    for sample_type in sample_types:
        sample_type_df = fractions_df[fractions_df['sample_type']==sample_type]
        if category != None:
                composition_within_sampletype_box(sample_type_df, output_dir, output_dir_results, category=category, sample_type=sample_type, stat_test=mannwhitneyu, perform_stat_test=True, immune=False, exclude_v17=exclude_v17)
            
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
    celltype_fraction_shifts_lineplot(df_immune, output_dir,output_dir_results, category=category, stat_test=wilcoxon, perform_stat_test=True, immune=True,exclude_v17=exclude_v17)
    celltype_fraction_composition_box(df_immune, output_dir, output_dir_results, category = category, immune=True, exclude_v17=exclude_v17, stat_test = mannwhitneyu, perform_stat_test=True)
    celltype_fraction_shifts_box(df_immune, output_dir, output_dir_results, category=category, stat_test=mannwhitneyu, perform_stat_test=True, immune=True, exclude_v17=exclude_v17)
    celltype_fraction_shifts_foldchange(df_immune, output_dir, output_dir_results, category=category, stat_test=mannwhitneyu, perform_stat_test=True, immune=True, exclude_v17=exclude_v17)
    for sample_type in sample_types:
        sample_type_df_immune = df_immune[df_immune['sample_type']==sample_type]
        if category != None:
            composition_within_sampletype_box(sample_type_df_immune, output_dir, output_dir_results, category=category, sample_type=sample_type, stat_test=mannwhitneyu, perform_stat_test=True, immune=True, exclude_v17=exclude_v17)