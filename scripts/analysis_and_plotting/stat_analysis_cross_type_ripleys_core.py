#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# stat_analysis_cross_type_ripleys.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Analyze and plot cross type Ripley's statistics for different groups and conditions.
#   
#   0 Import libraries and parse arguments
#   1 Define functions for statistical testing and plotting of centrality scores
#       a. Statistical testing functions for comparing two groups (can be mannwhitneyu - independent samples, wilcoxon - paired samples)
#       b. Functions for centrality scores analysis and plotting (boxplots and lineplots)
#           i. Boxplot - compare two groups of samples, possibly split on category (e.g. Biopsy vs. Resection, split on MPR -> separate plots for MPR B vs. R and non-MPR B vs. R)
#           ii. Lineplot - compare two groups of matched samples
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
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import os
import pickle
import argparse
from scipy.stats import wilcoxon, ttest_rel, ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 stat_analysis_centrality.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Perform statistical analysis and plotting between groups of samples for centrality scores. ') 
    parser.add_argument('-i', help='path to input ripleys analyzed adata',
                        dest='input',
                        type=str)
    parser.add_argument('--coi', nargs='+', required=True, 
                        help='list for cell types of interest')
    parser.add_argument('--exclude_v17', action='store_true',
                        help='Exclude v1.7 samples')
    parser.add_argument('-o_results', help='path to output dir for results',
                    dest='output_dir_results',
                    type=str)
    parser.add_argument('-o_plots', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

adata = args.input

exclude_v17=args.exclude_v17
print(f'Excluding v1.7 samples: {exclude_v17}')
if exclude_v17==True:
    adata = adata[adata.obs['treatment_scheme'] != 'v1.7',:].copy()

coi = pd.read_csv(args.coi)
output_dir_plots=args.output_dir_plots
output_dir_results=args.output_dir_results

#os.makedirs(os.path.join(output_dir_plots, 'centrality_scores'), exist_ok=True)

# Set aesthetics
sns.set_style("whitegrid")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Helper to create matrix of coi x coi
#--------------------------------------------------------------------------------
def interaction_to_matrix(values, interactions, celltypes):
    mat = np.full((len(celltypes), len(celltypes)), np.nan)
    for idx, name in enumerate(interactions):
        for i in celltypes:
            prefix = i + "_"
            if name.startswith(prefix):
                j = name[len(prefix):]
                if j in celltypes:
                    ii = celltypes.index(i)
                    jj = celltypes.index(j)
                    mat[ii, jj] = values[idx]
                break

    return mat


#def heatmap_zscore_per_r():


#def heatmap_integrals():


def ripleys_curve():









#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Prepare data and plot
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Extract statistics and metadata
#--------------------------------------------------------------------------------
obs_curve = adata.uns["ripley_obs_curve"]
sim_mean = adata.uns["ripley_sim_mean"]
sim_std = adata.uns["ripley_sim_std"]

interactions = adata.uns["ripley_interactions"]
n_types = len(coi)
radii = np.array(adata.uns["ripley_params"]["radii"])

print('Interactions: ', interactions)
print(n_types, 'cell types of interest: ', coi)
print('Tested radii: ', radii) 

# Choose groups and categories to compare
#--------------------------------------------------------------------------------
if exclude_v17==True:
    categories = [None, 'MPR']
else: 
    categories = [None, 'MPR', 'treatment']

groups = [None, 'sample_type']











# Z score at different radii
z_tensor = adata.uns["ripley_z"]

interaction = 'B_cell_T_cell_regulatory'
idx_inter = np.where(interactions == interaction)[0][0]

# average across cells in sample
z_mean_curve = np.nanmean(z_tensor[mask, idx_inter, :], axis=0)

plt.figure()
plt.plot(radii, z_mean_curve, marker='o')
plt.axhline(0, linestyle='--')
plt.title(f"Z-score curve: {interaction} ({core})")
plt.xlabel("Radius")
plt.ylabel("Z-score")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_zscore_{interaction}_{core}.png', bbox_inches='tight')
plt.close


# Ripley's stats curve
obs = np.nanmean(obs_curve[mask, idx_inter, :], axis=0)
sim_m = np.nanmean(sim_mean[mask, idx_inter, :], axis=0)
sim_s = np.nanmean(sim_std[mask, idx_inter, :], axis=0)
# upper = sim_m + 1.96 * sim_s
# lower = sim_m - 1.96 * sim_s

plt.figure(figsize=(6,5))

# shift by CSR
obs_shift = obs #- radii
sim_shift = sim_m #- radii
# upper_shift = upper - radii
# lower_shift = lower - radii

#plt.fill_between(radii, lower_shift, upper_shift, alpha=0.3, label="95% envelope")
plt.plot(radii, sim_shift, linestyle='--', label="Simulation mean")
plt.plot(radii, obs_shift, marker='o', label="Observed")

plt.axhline(0, linestyle=':', color='black')

plt.xlabel("Radius")
plt.ylabel("L(r) - r")
plt.title(f"{interaction} ({core})")
plt.legend()
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_curve_{interaction}_{core}.png', bbox_inches='tight')
plt.close




# Heatmaps of integrals
signed = np.nanmean(adata.obsm["ripley_signed"][mask], axis=0)
absolute = np.nanmean(adata.obsm["ripley_abs"][mask], axis=0)

signed_mat = interaction_to_matrix(signed, interactions, celltypes)
abs_mat = interaction_to_matrix(absolute, interactions, celltypes)

plt.figure(figsize=(6,5))
sns.heatmap(signed_mat, xticklabels=celltypes, yticklabels=celltypes,cmap="RdBu_r", center=0)
plt.title(f"Signed Ripley integral ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_signed_{core}.png', bbox_inches='tight')
plt.close


plt.figure(figsize=(6,5))
sns.heatmap(abs_mat, xticklabels=celltypes, yticklabels=celltypes,cmap="Reds")
plt.title(f"Absolute Ripley integral ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_abs_{core}.png', bbox_inches='tight')
plt.close


# Heatmap of mean z scores in all interactions
z_mean = np.nanmean(
    np.nanmean(z_tensor[mask], axis=0),  # cells → mean
    axis=1                               # radii → mean
)

z_mat = interaction_to_matrix(z_mean, interactions, celltypes)
plt.figure(figsize=(6,5))
sns.heatmap(z_mat, xticklabels=celltypes, yticklabels=celltypes,
            cmap="RdBu_r", center=0)
plt.title(f"Mean Z-score ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_zscore_{core}.png', bbox_inches='tight')
plt.close









