#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cross_type_ripleys.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Calculate Ripley's K function for cross-type interactions in spatial data.
#   Convert K to L.
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Define analysis functions

#   4 Choose analyses to perform:
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
import numpy as np
import scanpy as sc
import pandas as pd
from scipy.spatial import ConvexHull, cKDTree
import seaborn as sns
from numpy.random import default_rng
import argparse
import os
import pickle
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import combinations, product
from shapely.geometry import Polygon, Point

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 cross_type_ripleys.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Calculate cross-type Ripley\'s L function for spatial data.  ')
    parser.add_argument('-i', help='path to adata file',
                        dest='input',
                        type=str)
    parser.add_argument('--phen_level', help='key for cell type annotation in adata.obs',
                        dest='phen_level',
                        type=str)
    parser.add_argument('--coi', help='list for cell types of interest',
                        dest='coi', type=str, required=True, nargs='+')
    parser.add_argument('-o', '--output_results', help='path to output dir with patches per sample',
                        dest='output_dir_results',
                        type=str)
    parser.add_argument('--output_plots', help='path to output dir with patches per sample',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

adata = sc.read_h5ad(args.input)
cluster_key = args.phen_level
output_dir_results = args.output_dir_results
coi_list = args.coi
print(f"Cell types of interest: {coi_list}")

# Make sure output directories exist
os.makedirs(output_dir_results, exist_ok=True)
os.makedirs(args.output_dir_plots, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Set variables to be used throughout calculation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# minimum number of cells in source and target to continue with ripley's stat calculation 
min_i_threshold = 3
min_j_threshold = 3

coordinates = np.asarray(adata.obsm["spatial"])     # get coordinates of all the points (cells) in the sample
hull = ConvexHull(coordinates)      # identify cells that form the perimeter
window_polygon = Polygon(coordinates[hull.vertices])       # create a polygon connecting the outermost points to create a convex shape
area = window_polygon.area      # calculate area of the polygon


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define analysis functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate the distance to the hull edges
#--------------------------------------------------------------------------------
def distance_to_hull(coords, polygon):
    return np.array([polygon.boundary.distance(Point(p)) for p in coords])

# Calculate cross-type Ripley's L function with convex hull edge correction
#--------------------------------------------------------------------------------
def _cross_L(coords_i, coords_j, radii, polygon):
    # get the number of source and target cells
    n_i = coords_i.shape[0]
    n_j = coords_j.shape[0]

    area = polygon.area
    tree = cKDTree(coords_j)
    max_r = radii[-1]
    neighbors = tree.query_ball_point(coords_i, max_r)  # find neighbors within the maximum radius from source cell

    # compute distances to all target cells within the max radius for each source cell
    distances = [
        np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)
        if len(idx) > 0 else np.array([])
        for i, idx in enumerate(neighbors)]

    # distance to hull edge
    d_edge = distance_to_hull(coords_i, polygon)
    counts = np.zeros(len(radii))

    for k, r in enumerate(radii):
        total = 0   # store the number of target cells
        for i, dists in enumerate(distances):
            if len(dists) == 0:
                continue
            n_neighbors = np.sum(dists <= r)
            if n_neighbors == 0:
                continue
            # isotropic edge correction approximation
            if d_edge[i] >= r:  # if distance to edge larger than radius do not adjust the counts
                weight = 1.0
            else:   # if distance to edge smaller than radius that means the cell close to border, necessary to adjust for missing tissue beyonf the border
                weight = r / max(d_edge[i], 1e-6)   # 1e-6 to prevent division by too small numbers
                weight = min(weight, 5)     # limit the weight, if cell very close to the border than the weight might become very high adding noise
            total += n_neighbors * weight   # adjust total counts by their respective weights
        counts[k] = total   # for each radius save the weighted counts of target cells

    lambda_j = n_j / area   # calculate density for normalization
    K_ij = counts / (n_i * lambda_j)    # calculate Ripley's K
    L_ij = np.sqrt(K_ij / np.pi)    # calculate Ripley's L; possibly subtract radii to get the 0 as the baseline
    return L_ij

# Calculate integral difference between simulated L and observed L
#--------------------------------------------------------------------------------
def integral_obs_sim(radii, observed_L, sims):
    sim_mean = np.mean(sims, axis=0)    # get the mean value of all the simulations performed
    diff = observed_L - sim_mean    # get the difference between these two lines
    integral_value = np.trapz(diff,radii)   # integrate (calculate the area)
    integral_abs = np.trapz(abs(diff), radii)   # get the magnitude
    return integral_value, integral_abs


# Calculate cross-type Ripley's L function and perform permutation testing
#--------------------------------------------------------------------------------
def cross_ripley(adata, cluster_key, type_i, type_j, spatial_key="spatial", n_steps=50, max_dist=None, n_simulations=100, seed=None, copy=False):
    # Extract coordinates and labels into arrays
    coordinates = np.asarray(adata.obsm[spatial_key])
    labels = np.asarray(adata.obs[cluster_key])

    # Select coordinates for the two cell types
    coords_i = coordinates[labels == type_i]
    coords_j = coordinates[labels == type_j]

    if max_dist is None:
        max_dist = np.sqrt(area / 2)

    # Create array of radii to evaluate
    radii = np.linspace(0, max_dist, n_steps)

    # Calculate observed cross-type L statistic
    observed_L = _cross_L(coords_i, coords_j, radii, window_polygon)

    # Perform permutation testing - randomly shuffle labels and recalculate L for each permutation
    rng = default_rng(seed)
    sims = np.zeros((n_simulations, len(radii)))

    failed_density = 0
    valid_simulations = 0
    for s in range(n_simulations):
        permuted = rng.permutation(labels)
        # Get coordinates for the permuted cell types
        coords_i_perm = coordinates[permuted == type_i]
        coords_j_perm = coordinates[permuted == type_j]
        # If a permutation results in one of the types having no cells, skip this permutation
        if len(coords_i_perm) == 0 or len(coords_j_perm) == 0:  
            failed_density += 1
            continue
        # Calculate L for the permuted data
        sims[s] = _cross_L(coords_i_perm, coords_j_perm, radii, window_polygon)
        valid_simulations += 1
    
    sd_sim = np.std(sims, axis=0)

    failed_sd = False
    if np.any(sd_sim == 0) or np.any(np.isnan(sd_sim)):
        failed_sd = True
    # two-sided Monte Carlo p-values
    pvals = (np.sum(sims >= observed_L, axis=0) + 1) / (n_simulations + 1)
    pvals = np.minimum(pvals, 1 - pvals)

    integral_signed, integral_abs = integral_obs_sim(radii, observed_L, sims)

    result = {
        "r": radii,
        "L_observed": observed_L,
        "L_simulations": sims,
        "pvalues": pvals,
        "csr_expectation": radii,
        "type_i": type_i,
        "type_j": type_j,
        "n_source_cells": n_i,
        "n_target_cells": n_j,
        "integral_signed": integral_signed,
        "integral_abs": integral_abs,
        "reason": None,
        "failed_density_simulations": failed_density,
        "valid_simulations": valid_simulations,
        "failed_sd": failed_sd
        }

    if copy:
        return result

    key = f"cross_ripley_{type_i}_{type_j}"
    adata.uns[key] = result

# Create empty results dictionary to ensure consistent output for each sample and pair of celltypes
#------------------------------------------------------------------------------
def empty_ripley_result(type_i, type_j, n_steps, reason):
    radii = np.full(n_steps, np.nan)

    result = {
        "r": radii,
        "L_observed": np.full(n_steps, np.nan),
        "L_simulations": np.full((0, n_steps), np.nan),
        "pvalues": np.full(n_steps, np.nan),
        "csr_expectation": radii,
        "type_i": type_i,
        "type_j": type_j,
        "n_source_cells": n_i,
        "n_target_cells": n_j,
        "integral_signed": np.nan,
        "integral_abs": np.nan,
        "reason": reason,
        "failed_density_simulations": np.nan,
        "valid_simulations": np.nan,
        "failed_sd": np.nan}
    return result

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cell_counts = adata.obs[cluster_key].value_counts()
n_steps = 30

dict_res = defaultdict(dict)
for type_i, type_j in product(coi_list,coi_list):
    if type_i == type_j:
        continue

    n_i = cell_counts.get(type_i, 0)
    n_j = cell_counts.get(type_j, 0)

    if n_i <= min_i_threshold:
        dict_res[type_i][type_j] = empty_ripley_result(type_i, type_j, n_steps, 'limited_source_cell')
        continue

    if n_j < min_j_threshold:
        dict_res[type_i][type_j] = empty_ripley_result(type_i, type_j, n_steps, 'limited_target_cell')
        continue

    print(f"Calculating cross-type Ripley's L for {type_i} vs {type_j}...")
    
    res = cross_ripley(adata, cluster_key, type_i, type_j, spatial_key="spatial", n_steps=n_steps, max_dist=300, n_simulations=500, seed=42, copy=True)
    dict_res[type_i][type_j] = res

# Join the integrals into a matrix
#----------------------------------------------------------------------------
celltypes = coi_list

# initialize DataFrames with NaNs
signed_integral_df = pd.DataFrame(np.nan, index=celltypes, columns=celltypes)
absolute_integral_df = pd.DataFrame(np.nan, index=celltypes, columns=celltypes)

for type_i in celltypes:
    for type_j in celltypes:
        if type_i == type_j:
            continue
        res = dict_res.get(type_i, {}).get(type_j, None)
        if res is None:
            continue
        signed_integral_df.at[type_i, type_j] = res.get("integral_signed", np.nan)
        absolute_integral_df.at[type_i, type_j] = res.get("integral_abs", np.nan)

dict_res["signed_integral_matrix"] = signed_integral_df
dict_res["absolute_integral_matrix"] = absolute_integral_df
dict_res["celltypes"] = celltypes

# Save results as dictionary
with open(os.path.join(output_dir_results, f"dict_ripleys_L.pkl"), "wb") as f:
    pickle.dump(dict_res, f)

# Plot signed and absolute matrix for a check 
#----------------------------------------------------------------------------
signed_df = dict_res["signed_integral_matrix"]
absolute_df = dict_res["absolute_integral_matrix"]

# Quick heatmap of signed integral
plt.figure(figsize=(7, 6))
sns.heatmap(signed_df, cmap="RdBu_r", center=0)
plt.title("Signed Integral Heatmap (Cross-type Ripley's L)")
plt.xlabel("Target Cell Type")
plt.ylabel("Source Cell Type")
plt.tight_layout()
plt.savefig(args.output_dir_plots + '/signed_integral.svg', bbox_inches='tight')

# Quick heatmap of absolute integral
plt.figure(figsize=(7, 6))
sns.heatmap(absolute_df, cmap="Reds")
plt.title("Absolute Integral Heatmap (Cross-type Ripley's L)")
plt.xlabel("Target Cell Type")
plt.ylabel("Source Cell Type")
plt.tight_layout()
plt.savefig(args.output_dir_plots + '/absolute_integral.svg', bbox_inches='tight')





#with open(os.path.join(output_dir_results, 'cross_ripley', f"cross_ripley_{type_i}_vs_{type_j}.pkl"), "wb") as f:
#    pickle.dump(res, f)











# Calculate distance of each cell to the nearest border edge
#--------------------------------------------------------------------------------
#def distance_to_border(coords):
#    xmin, ymin = coords.min(axis=0)
#    xmax, ymax = coords.max(axis=0)
#
#    dx = np.minimum(coords[:,0] - xmin, xmax - coords[:,0])
#    dy = np.minimum(coords[:,1] - ymin, ymax - coords[:,1])
#
#    return np.minimum(dx, dy)

# Calculate cross-type Ripley's L function
#--------------------------------------------------------------------------------
#def _cross_L(coords_i, coords_j, radii, area):
#    # Calculate pairwise distances between points of type i and type j
#    n_i = coords_i.shape[0]
#    n_j = coords_j.shape[0]
#    distances = cdist(coords_i, coords_j)
#    # For each radius, count the number of pairs (i,j) with distance <= r
#    counts = np.array([np.sum(distances <= r) for r in radii])
#    # Calculate K_ij and convert to L_ij
#    lambda_j = n_j / area
#    K_ij = counts / (n_i * lambda_j)
#    L_ij = np.sqrt(K_ij / np.pi)
#
#    return L_ij