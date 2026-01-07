#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# spatial_analysis.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze spatial proximity of cell types.
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Define analysis functions
#   3 Choose analyses to perform:
#       a. Compute centrality scores
#           i. closeness centrality - measure of how close the group is to other nodes
#           ii. degree centrality - fraction of non-group members connected to group members
#           iii. betweenness centrality - measure of how often a node appears on shortest paths between other nodes
#           iv. clustering coefficient - measure of the degree to which nodes cluster together
#       b. Compute neighbors enrichment - measure of how often cell types are neighbors compared to random expectation
#       c. Compute Moran's I spatial autocorrelation - measure of how gene expression is correlated with spatial location
#       d. Compute interaction matrix - measure of interactions between cell types
#       e. Compute co-occurrence probabilities - measure of how often cell types co-occur in the same neighborhood
#   4 Save
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

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Run_tacco.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Run tacco to transfer the cell labels from reference scRNA atlas. Preprocess datasets individually for potential inspection  ')
    parser.add_argument('-i', help='path to adata sample subset',
                        dest='input',
                        type=str)
    parser.add_argument('-o_results', help='path to output dir for results',
                        dest='output_dir_results',
                        type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    parser.add_argument('--celltype_key', help='phenotyping level',
                        dest='celltype_key',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
os.makedirs(args.output_dir_results, exist_ok=True)
os.makedirs(args.output_dir_plots, exist_ok=True)

print('Running')