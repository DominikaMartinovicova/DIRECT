#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# gex_analysis.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze gene expression profiles of cell types before and after treatment. Possibly split  
#   patients into groups based on chosen category.
#
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Define analysis functions
#   3 Choose analyses to perform:
#       a. Volcano plot
#       b. Heatmap
#       c. 
#       d. 
#       * for Biopsy vs. Resection; 
#   4 Save
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
from matplotlib import category
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import pertpy as pt
import decoupler as dc
import argparse
import os
import numpy as np
import pandas as pd
import math
import anndata as ad
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu
#from statannotations.Annotator import Annotator
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
def dge_analysis_BvsR(pdata, celltype_key, exclude_v17, celltype, output_dir):
    os.makedirs(output_dir + 'BvsR/overall/', exist_ok=True)
    exv17 = '(excluding v1.7)' if exclude_v17 else ''
    excl_v17 = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    if celltype == None:
        print('Creating edgeR object and fitting model...')
        edgr = pt.tl.EdgeR(pdata, design="~sample_type")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column="sample_type", baseline="Biopsy", group_to_compare="Resection"))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE Biopsy vs. Resection samples {exv17}')
        plt.savefig(os.path.join(output_dir, f'BvsR/overall/volcano_plot_BvsR_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        edgr.plot_paired(pdata, results_df=res_df, n_top_vars=5, groupby="sample_type", pairedby="MPR")
        plt.suptitle(f'Top 5 differentially expressed genes in Biopsy vs. Resection {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'BvsR/overall/paired_top5_dge_BvsR_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        res_df = edgr.compare_groups(pdata, column="sample_type", baseline="Biopsy", groups_to_compare="Resection")
        edgr.plot_multicomparison_fc(res_df)
        plt.title(f'Log2FC comparison between Biopsy and Resection {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'BvsR/overall/log2fc_comparison_BvsR_{excl_v17}.png'), bbox_inches='tight')
        plt.close()
    else:
        print(f'Creating edgeR object and fitting model for cell type: {celltype}...')
        edgr = pt.tl.EdgeR(pdata, design="~sample_type")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column="sample_type", baseline="Biopsy", group_to_compare="Resection"))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE between Biopsy and Resection samples in {celltype} {exv17}')
        plt.savefig(os.path.join(output_dir, f'BvsR/overall/volcano_plot_BvsR_{celltype}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        edgr.plot_paired(pdata, results_df=res_df, n_top_vars=5, groupby="sample_type", pairedby="MPR")
        plt.suptitle(f'Top 5 differentially expressed genes between Biopsy and Resection in {celltype} {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'BvsR/overall/paired_top5_dge_BvsR_{celltype}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()   


def dge_analysis_within_sampletype(pdata, celltype_key, sample_type, category, exclude_v17, celltype, output_dir):
    os.makedirs(output_dir + f'{sample_type}/{category}/', exist_ok=True)
    exv17 = '(excluding v1.7)' if exclude_v17 else ''
    excl_v17 = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    
    if pdata.obs[category].nunique() != 2:
        print(f'Category {category} does not have exactly 2 groups in sample type {sample_type}. Skipping DGE analysis for this category and sample type.')
        return
    elif pdata.obs[category].nunique() == 2:
        group1, group2 = pdata.obs[category].unique()
    
    if celltype == None:
        print(f'Creating edgeR object and fitting model for sample type: {sample_type}...')
        edgr = pt.tl.EdgeR(pdata, design=f"~{category}")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column=category, baseline=group1, group_to_compare=group2))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE in {sample_type} samples in {group1} vs. {group2} {exv17}')
        plt.savefig(os.path.join(output_dir, f'{sample_type}/{category}/volcano_plot_{category}_{sample_type}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        res_df = edgr.compare_groups(pdata, column=category, baseline=group1, groups_to_compare=group2)
        edgr.plot_multicomparison_fc(res_df)
        plt.title(f'Log2FC comparison between {category} groups in {sample_type} {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{sample_type}/{category}/log2fc_comparison_{category}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

    else:
        print(f'Creating edgeR object and fitting model for sample type: {sample_type} and cell type: {celltype}...')
        edgr = pt.tl.EdgeR(pdata, design=f"~{category}")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column=category, baseline=group1, group_to_compare=group2))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE in {sample_type} samples in {group1} vs. {group2} {celltype} {exv17}')
        plt.savefig(os.path.join(output_dir, f'{sample_type}/{category}/volcano_plot_{category}_{sample_type}_{celltype}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        res_df = edgr.compare_groups(pdata, column=category, baseline=group1, groups_to_compare=group2)
        edgr.plot_multicomparison_fc(res_df)
        plt.title(f'Log2FC comparison between {category} groups in {sample_type} in {celltype} {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'{sample_type}/{category}/log2fc_comparison_{category}_{sample_type}_{celltype}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()



def dge_analysis_within_category(pdata, celltype_key, category, group, sample_type, exclude_v17, celltype, output_dir):
    os.makedirs(output_dir + f'BvsR/{category}/', exist_ok=True)
    exv17 = '(excluding v1.7)' if exclude_v17 else ''
    excl_v17 = 'wo_v1.7' if exclude_v17 else 'w_v1.7'
    
    if pdata.obs[sample_type].nunique() != 2:
        print(f'Category {category} does not have exactly 2 groups in sample type. Skipping DGE analysis for this category and sample type.')
        return
    elif pdata.obs[sample_type].nunique() == 2:
        group1, group2 = 'Biopsy', 'Resection' #pdata.obs[sample_type].unique()
    
    if celltype == None:
        print(f'Creating edgeR object and fitting model for category: {category}...')
        edgr = pt.tl.EdgeR(pdata, design=f"~{sample_type}")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column=sample_type, baseline=group1, group_to_compare=group2))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE in {category} {group} in {group1} vs. {group2} {exv17}')
        plt.savefig(os.path.join(output_dir, f'BvsR/{category}/volcano_plot_{category}_{group}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

        res_df = edgr.compare_groups(pdata, column=sample_type, baseline=group1, groups_to_compare=group2)
        edgr.plot_multicomparison_fc(res_df)
        plt.title(f'Log2FC comparison between {group1} and {group2} groups in {category} {group} {exv17}')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'BvsR/{category}/log2fc_comparison_{category}_{group}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()

    else:
        print(f'Creating edgeR object and fitting model for category: {category} and cell type: {celltype}...')
        edgr = pt.tl.EdgeR(pdata, design=f"~{sample_type}")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column=sample_type, baseline=group1, group_to_compare=group2))
        #print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'DGE in {category} {group} in {group1} vs. {group2} in  {celltype} {exv17}')
        plt.savefig(os.path.join(output_dir, f'BvsR/{category}/volcano_plot_{category}_{group}_{celltype}_{excl_v17}.png'), bbox_inches='tight')
        plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create pseudobulk adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print('Creating pseudobulk adata...')
pdata = dc.pp.pseudobulk(adata, sample_col='sample', groups_col=celltype_key, layer='raw_counts', mode='sum')
dc.pp.filter_samples(pdata, min_cells=10, min_counts=100)  # Filter samples with too low count or too few cells
print(pdata)
print(pdata.obs.head())


dge_analysis_BvsR(pdata = pdata, celltype_key=celltype_key, celltype=None, exclude_v17=exclude_v17, output_dir=output_dir)
for celltype in pdata.obs[celltype_key].unique():
    print(f'Performing DGE analysis for cell type: {celltype}...')
    pdata_ct = pdata[pdata.obs[celltype_key] == celltype, :].copy()
    dge_analysis_BvsR(pdata = pdata_ct, celltype_key=celltype_key, celltype=celltype, exclude_v17=exclude_v17, output_dir=output_dir)


categories = ['MPR', 'treatment']
sample_type = ['Biopsy', 'Resection']
for sample_t in sample_type:
    pdata_st = pdata[pdata.obs['sample_type'] == sample_t, :].copy()
    for category in categories:
        print(f'Performing DGE analysis for category: {category}...')
        dge_analysis_within_sampletype(pdata = pdata_st, celltype_key=celltype_key, sample_type=sample_t, category=category, output_dir=output_dir, exclude_v17=exclude_v17, celltype=None)
        for celltype in pdata.obs[celltype_key].unique():
            print(f'Performing DGE analysis for category: {category} in cell type: {celltype}...')
            pdata_st_ct = pdata_st[pdata_st.obs[celltype_key] == celltype, :].copy()
            dge_analysis_within_sampletype(pdata = pdata_st_ct, celltype_key=celltype_key, sample_type=sample_t, category=category, output_dir=output_dir, exclude_v17=exclude_v17, celltype=celltype)


for category in categories:
    print(f'Performing DGE analysis for category: {category}...')
    group1, group2 = pdata.obs[category].unique()
    for group in [group1, group2]:
        pdata_subset = pdata[pdata.obs[category]==group, :].copy()
        dge_analysis_within_category(pdata = pdata_subset, celltype_key=celltype_key, category=category, group=group, sample_type='sample_type', output_dir=output_dir, exclude_v17=exclude_v17, celltype=None)
        for celltype in pdata.obs[celltype_key].unique():
            print(f'Performing DGE analysis for category: {category} in cell type: {celltype}...')
            pdata_ct = pdata[pdata.obs[celltype_key] == celltype, :].copy()
            dge_analysis_within_category(pdata = pdata_ct, celltype_key=celltype_key, category=category, group=group, sample_type='sample_type', output_dir=output_dir, exclude_v17=exclude_v17, celltype=celltype)

