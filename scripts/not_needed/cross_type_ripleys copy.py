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
from scipy.spatial.distance import cdist
from numpy.random import default_rng
import argparse
import os
import pickle
import matplotlib.pyplot as plt
from collections import defaultdict
from itertools import combinations, product

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
    args = parser.parse_args()
    return args

args = parse_args()

adata = sc.read_h5ad(args.input)
cluster_key = args.phen_level
output_dir_results = args.output_dir_results
coi_list = args.coi
print(f"Cell types of interest: {coi_list}")

min_i_threshold = 3
min_j_threshold = 3

# Make sure output directories exist
os.makedirs(output_dir_results, exist_ok=True)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define analysis functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate cross-type Ripley's L function and perform permutation testing
#--------------------------------------------------------------------------------
def cross_ripley(adata, cluster_key, type_i, type_j, spatial_key="spatial", n_steps=50, max_dist=None, n_simulations=100, seed=None, copy=False):
    # Check if spatial coordinates and cluster labels are present
    if spatial_key not in adata.obsm:
        raise ValueError(f"{spatial_key} not found in adata.obsm")
    if cluster_key not in adata.obs:
        raise ValueError(f"{cluster_key} not found in adata.obs")
    # Extract coordinates and labels into arrays
    coordinates = np.asarray(adata.obsm[spatial_key])
    labels = np.asarray(adata.obs[cluster_key])
    # Select coordinates for the two cell types
    coords_i = coordinates[labels == type_i]
    coords_j = coordinates[labels == type_j]

    if len(coords_i) == 0 or len(coords_j) == 0:
        raise ValueError("Selected types contain no cells.")

    # Calculate area of the convex hull of all points to use for intensity estimation
    hull = ConvexHull(coordinates)
    area = hull.volume

    if max_dist is None:
        max_dist = np.sqrt(area / 2)

    # Create array of radii to evaluate
    radii = np.linspace(0, max_dist, n_steps)

    # Calculate observed cross-type L statistic
    observed_L = _cross_L(coords_i, coords_j, radii, area)

    # Perform permutation testing - randomly shuffle labels and recalculate L for each permutation
    rng = default_rng(seed)
    sims = np.zeros((n_simulations, len(radii)))

    for s in range(n_simulations):
        permuted = rng.permutation(labels)
        # Get coordinates for the permuted cell types
        coords_i_perm = coordinates[permuted == type_i]
        coords_j_perm = coordinates[permuted == type_j]
        # If a permutation results in one of the types having no cells, skip this permutation
        if len(coords_i_perm) == 0 or len(coords_j_perm) == 0:  
            continue
        # Calculate L for the permuted data
        sims[s] = _cross_L(coords_i_perm, coords_j_perm, radii, area)

    # two-sided Monte Carlo p-values
    pvals = (np.sum(sims >= observed_L, axis=0) + 1) / (n_simulations + 1)
    pvals = np.minimum(pvals, 1 - pvals)

    result = {
        "r": radii,
        "L_observed": observed_L,
        "L_simulations": sims,
        "pvalues": pvals,
        "csr_expectation": radii,
        "type_i": type_i,
        "type_j": type_j,
    }

    if copy:
        return result

    key = f"cross_ripley_{type_i}_{type_j}"
    adata.uns[key] = result

# Calculate cross-type Ripley's L function
#--------------------------------------------------------------------------------
def _cross_L(coords_i, coords_j, radii, area, window_coords):
    n_i = coords_i.shape[0]
    n_j = coords_j.shape[0]
    max_r = radii[-1]
    tree_j = cKDTree(coords_j)

    # query once at maximum radius
    neighbors = tree_j.query_ball_point(coords_i, max_r)
    # compute distances only once
    distances = [np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)
        if len(idx) > 0 else np.array([])
        for i, idx in enumerate(neighbors)]

    # bounding box for edge correction
    xmin, ymin = coords_i.min(axis=0)
    xmax, ymax = coords_i.max(axis=0)

    dx = np.minimum(coords_i[:,0] - xmin, xmax - coords_i[:,0])
    dy = np.minimum(coords_i[:,1] - ymin, ymax - coords_i[:,1])
    d_edge = np.minimum(dx, dy)

    counts = np.zeros(len(radii))

    for k, r in enumerate(radii):
        total = 0
        for i, dists in enumerate(distances):
            if len(dists) == 0:
                continue
            # neighbors within radius
            n_neighbors = np.sum(dists <= r)
            if n_neighbors == 0:
                continue
            # circumference edge correction
            if d_edge[i] >= r:
                weight = 1.0
            else:
                # isotropic correction
                weight = r / max(d_edge[i], 1e-9)
            total += n_neighbors * weight
        counts[k] = total

    lambda_j = n_j / area
    K_ij = counts / (n_i * lambda_j)
    L_ij = np.sqrt(K_ij / np.pi)

    return L_ij

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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Choose analyses to perform
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cell_counts = adata.obs[cluster_key].value_counts()


dict_res = defaultdict(dict)
for type_i, type_j in product(coi_list,coi_list):
    if type_i == type_j:
        continue

    n_i = cell_counts.get(type_i, 0)
    n_j = cell_counts.get(type_j, 0)

    if n_i <= min_i_threshold:
        dict_res[type_i][type_j] = {"reason": "limited_source_cell"}
        continue

    if n_j < min_j_threshold:
        dict_res[type_i][type_j] = {"reason": "limited_target_cell"}
        continue

    print(f"Calculating cross-type Ripley's L for {type_i} vs {type_j}...")
    
    res = cross_ripley(adata, cluster_key, type_i, type_j, spatial_key="spatial", n_steps=30, max_dist=300, n_simulations=500, seed=42, copy=True)
    dict_res[type_i][type_j] = res

#print(dict_res)


# Save results as dictionary
with open(os.path.join(output_dir_results, f"dict_ripleys_L.pkl"), "wb") as f:
    pickle.dump(dict_res, f)

#with open(os.path.join(output_dir_results, 'cross_ripley', f"cross_ripley_{type_i}_vs_{type_j}.pkl"), "wb") as f:
#    pickle.dump(res, f)
