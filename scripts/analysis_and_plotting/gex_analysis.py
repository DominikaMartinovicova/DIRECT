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
#       * statistical testing functions 
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
def dge_analysis(pdata, celltype_key, sample_type, exclude_v17, output_dir):
    if sample_type is None:
        print('Creating edgeR object and fitting model...')
        edgr = pt.tl.EdgeR(pdata, design="~sample_type")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column="sample_type", baseline="Biopsy", group_to_compare="Resection"))
        print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title('Volcano plot of DGE between Biopsy and Resection samples')
        plt.savefig(os.path.join(output_dir, 'volcano_plot.png'), bbox_inches='tight')
        plt.close()

        edgr.plot_paired(pdata, results_df=res_df, n_top_vars=5, groupby="sample_type", pairedby="MPR")
        plt.suptitle('Top 5 differentially expressed genes between Biopsy and Resection samples')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'paired_top5_dge.png'), bbox_inches='tight')
        plt.close()

        res_df = edgr.compare_groups(pdata, column="sample_type", baseline="Biopsy", groups_to_compare="Resection")
        edgr.plot_multicomparison_fc(res_df)
        plt.title('Log2FC comparison between Biopsy and Resection samples')
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, 'log2fc_comparison.png'), bbox_inches='tight')
        plt.close()
    else:
        print(f'Performing DGE analysis for sample type: {sample_type}...')
        edgr = pt.tl.EdgeR(pdata, design="~MPR + treatment")
        edgr.fit()

        res_df = edgr.test_contrasts(edgr.contrast(column="treatment", baseline="milder", group_to_compare="aggressive"))
        print(res_df)

        edgr.plot_volcano(res_df, log2fc_thresh=0)
        plt.title(f'Volcano plot of DGE between treatment groups in {sample_type} samples')
        plt.savefig(os.path.join(output_dir, f'volcano_plot_{sample_type}.png'), bbox_inches='tight')
        plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create pseudobulk adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print('Creating pseudobulk adata...')
pdata = dc.pp.pseudobulk(adata, sample_col='sample', groups_col=celltype_key, layer='raw_counts', mode='sum')
dc.pp.filter_samples(pdata, min_cells=10, min_counts=100)  # Filter samples with too low count or too few cells
print(pdata)
print(pdata.obs.head())

categories = ['MPR', 'treatment']
sample_type = [None, 'Biopsy', 'Resection']
# for category in categories:
#     print(f'Performing DGE analysis for category: {category}...')
#     dge_analysis(pdata = pdata, celltype_key=celltype_key, category=category, output_dir=output_dir)

dge_analysis(pdata = pdata, celltype_key=celltype_key, sample_type=sample_type, output_dir=output_dir)

for sample_t in sample_type:
    if sample_t is None:
        dge_analysis(pdata = pdata, celltype_key=celltype_key, sample_type=sample_t, exclude_v17=exclude_v17, output_dir=output_dir)
    else:
        pdata = pdata[pdata.obs['sample_type'] == sample_t, :].copy()
        dge_analysis(pdata = pdata, celltype_key=celltype_key, sample_type=sample_t, exclude_v17=exclude_v17, output_dir=output_dir)



