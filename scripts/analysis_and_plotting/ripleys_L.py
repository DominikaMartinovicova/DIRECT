#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ripleys.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Calculate Ripley's statistics across samples (core-/tissue-level).
#   
#   0 Import libraries and parse arguments
#   1 Define functions for statistical testing and plotting of centrality scores
#       a. Statistical testing functions for paired and independent samples
#       b. Functions for centrality scores analysis and plotting (boxplots and lineplots)
#   2 Prepare data and run analysis
#   
#
#
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import scanpy as sc
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import pickle
import argparse
import squidpy as sq
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 ripleys.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Calculate Ripley\'s statistics across samples (core-/tissue-level). ') 
    parser.add_argument('-i', help='path to input directory with spatial analysis results per sample',
                        dest='input',
                        type=str)
    parser.add_argument('--phen_level', help='key with celltype annotations in adata.obs',
                        dest='phenotyping_level',
                        type=str)
    parser.add_argument('--exclude_v17', action='store_true',
                        help='Exclude v1.7 samples')
    parser.add_argument('-o_report', help='report with stat analysis results',
                        dest='output_dir_report',
                        type=str)
    parser.add_argument('-o_results', help='path to output dir for results',
                    dest='output_dir_results',
                    type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args


args = parse_args()
adata = sc.read_h5ad(args.input)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Calculate Ripley's statistics across samples
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
dict_ripley = sq.gr.ripley(adata, cluster_key = args.phenotyping_level, mode = 'L', copy=True)

with open(args.output_dir_results, "wb") as f:
    pickle.dump(dict_ripley, f)

print("Ripley's L statistics calculated and saved to results directory.")










