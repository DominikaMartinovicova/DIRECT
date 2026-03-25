#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cross_type_ripleys_per_cell.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Calculate Ripley's K function for cross-type interactions in spatial data.
#   Convert K to L.
#
#   0 Import libraries and parse arguments
#   1 Set variables
#   2 Define analysis functions
#       a. Edge correction 
#       b. Cross type Ripley's L 
#       c. Process cells per sample
#   3 Main loop - call process sample function for each sample in adata
#   4 Perform statistical analysis and store results in adata
#   5 Save adata with Ripley's results
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
from scipy.spatial import ConvexHull, cKDTree
from shapely.geometry import Polygon, Point
from itertools import product
from numpy.random import default_rng
from joblib import Parallel, delayed
import argparse
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
import pickle

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(prog = 'python3 cross_type_ripleys_cell_level.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Calculate cross-type Ripley\'s L function for spatial data.  ')
    parser.add_argument('-i', dest='input', type=str, help='path to adata file')
    parser.add_argument('-o', dest='output', type=str, help = 'path to output adata that contains results of Ripley\'s stats per cell saved in adata.obsm')
    parser.add_argument('--phen_level', dest='phen_level', type=str, help='key for cell type annotation in adata.obs')
    parser.add_argument('--coi', nargs='+', required=True, help='list for cell types of interest')
    parser.add_argument('--sample_key', default="sample", help='key for pieces of tissue/cores in adata.obs')
    parser.add_argument('--max_dist', type=float, default=250, help='maximum radius to calculate Ripley\'s stats')
    parser.add_argument('--n_steps', type=int, default=10, help='number of radii to calculate the stats for (evenly distributed between 0-max_dist)')
    parser.add_argument('--n_sim', type=int, default=100, help='number of simulations to obtain null distribution')
    return parser.parse_args()

args = parse_args()
adata = sc.read_h5ad(args.input)
print(adata)
cluster_key = args.phen_level
print(adata.obs[cluster_key].value_counts())
sample_key = args.sample_key
coi_list = args.coi
print(f'{len(coi_list)} celltypes of interest  {coi_list}.')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Set variables to be used throughout calculation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# minimum number of cells in source and target to continue with ripley's stat calculation 
min_i_threshold = 3
min_j_threshold = 3

rng = default_rng(42)   # set random seed for reproducibility

coordinates = np.asarray(adata.obsm["spatial"])     # extract coordinates
print(f'Coordinates: {coordinates.shape}, \n {coordinates[:5,:5]}')
labels = np.asarray(adata.obs[cluster_key])         # extract labels 
print(f'{len(np.unique(labels))} unique labels: {np.unique(labels)}')
samples = np.asarray(adata.obs[sample_key])
print(f'{len(np.unique(samples))} unique samples: {np.unique(samples)}')


#radii = np.linspace(0, args.max_dist, args.n_steps)     # create list of radii to be tested
radii = [25, 50, 75, 100, 250]
n_r = len(radii)
print(f'{n_r} radii to be tested: {radii}')

n_cells = adata.n_obs
interactions = [f"{i}_{j}" for i, j in product(coi_list, coi_list)]# if i != j]
interaction_index = {k: i for i, k in enumerate(interactions)}
n_interactions = len(interactions)
print(f'Number of celltype pairs to be tested: {n_interactions}')
print("Coordinate range (x, y):", np.ptp(coordinates, axis=0))
print("Max radius:", radii[-1])

# Arrays to store results
#-------------------------------------------------------------------------------
obs_curve = np.full((n_cells, n_interactions, len(radii)), np.nan)
sim_sum = np.zeros_like(obs_curve)
sim_sq_sum = np.zeros_like(obs_curve)
sim_valid = np.zeros((n_cells, n_interactions))
sim_failed = np.zeros((n_cells, n_interactions))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Define analysis functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Edge correction
#-------------------------------------------------------------------------------
def compute_edge_weights(coords, polygon, radii, max_weight=2.0):
    dists = np.array([polygon.boundary.distance(Point(p)) for p in coords])
    weights = np.ones((len(coords), len(radii)))    # default weight is one, if cell close to the border weight is increased up to max_weight
    for k, r in enumerate(radii):
        mask = dists < r                            # identify cells that are closer to the border than the tested radius
        weights[mask, k] = np.minimum(r / np.maximum(dists[mask], 1e-6), max_weight)
    return weights

# Calculate cross-type Ripley's L function
#--------------------------------------------------------------------------------
def cross_L(coords_i, coords_j, radii, weights_i, lambda_j):
    tree = cKDTree(coords_j)        # identify all type j cells
    neighbors = tree.query_ball_point(coords_i, radii[-1])  # identify j neighbors of source cells within max radius
    L = np.zeros((len(coords_i), len(radii)))       # create np array with 0s to store results, default zero (if no neighbors found within max radius)
    for i, idx in enumerate(neighbors):
        if len(idx)==0:
            continue
        dists = np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)     # calculate distances to neighbors of type j
        dists.sort()      # sort distances to neighbors of type j
        counts = np.searchsorted(dists, radii)       # count how many neighbors of type j are within each radius
        L[i, :] = (counts * weights_i[i]) / lambda_j    # apply edge correction and normalize by intensity of type j to get K (variable is called L but it is actually K at this point)
    return np.sqrt(L / np.pi)   # division makes this statistic into Ripley's L (before it was K even if variable was called L)

# Process sample
#--------------------------------------------------------------------------------
def process_sample(sample_id):   
    mask = samples == sample_id     # find cells belonging to the sample being processed
    coords_core = coordinates[mask] # extract coordinates of the sample being processed
    labels_core = labels[mask]      # extract labels of the sample being processed
    global_idx = np.where(mask)[0]  # store global indices of the cells being processed to be able to put the results back in the right place in the adata.obsm after processing

    print(f"[START] Sample {sample_id} | {len(coords_core)} cells")

    hull = ConvexHull(coords_core)  # calculate convex hull of the sample to be able to perform edge correction
    polygon = Polygon(coords_core[hull.vertices])   # create polygon from convex hull vertices (cells forming the shape)

    weights_core = compute_edge_weights(coords_core, polygon, radii)   # calculate edge weights
    n_local = len(coords_core)  # number of cells in the sample being processed

    # Create arrays to store results for the sample being processed, default NaN (if not enough cells of type i or j to calculate the statistic)
    obs_local = np.full((n_local, n_interactions, n_r), np.nan)
    sim_L = np.full((args.n_sim, n_interactions, n_r), np.nan)  # store L values for each simulation for each interaction for each radius to be able to calculate mean and variance of the simulated L values later
    sim_valid_inter = np.zeros(n_interactions)   # count how many simulations were valid (enough cells of type i and j) for each interaction to be able to filter out interactions with low simulation support later
    
    # ----------------- Calculate observed Ripley's L for the sample being processed ----------------------
    for type_i, type_j in product(coi_list, coi_list):      # iterate through all pairs of cell types to be tested
        idx_inter = interaction_index[f"{type_i}_{type_j}"] # get index of the interaction in the results arrays
        mask_i = labels_core == type_i      # identify source cells of type i
        mask_j = labels_core == type_j      # identify target cells of type j

        if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
            continue    # skip if too few source or target cells

        coords_i = coords_core[mask_i]      # extract coordinates of source cells of type i
        coords_j = coords_core[mask_j]      # extract coordinates of target cells of type j
        weights_i = weights_core[mask_i]    # extract edge weights of source cells of type i

        lambda_j = len(coords_j) / polygon.area     # normalization factor for Ripley's K/L to account for different densities of target cell type j
        L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)     # calculate Ripley's L for the source cells of type i and target cells of type j

        obs_local[mask_i, idx_inter, :] = L     # store the observed L values for the source cells of type i and interaction with target type j in the results array

    # ----------------- Perform simulations to get null distribution of Ripley's L ----------------------
    for s in range(args.n_sim):
        if s % 100 == 0:
            print(f"[{sample_id}] simulation {s}/{args.n_sim}")
        perm = rng.permutation(labels_core)      # shuffle cell labels to create null distribution

        for type_i, type_j in product(coi_list, coi_list):
            idx_inter = interaction_index[f"{type_i}_{type_j}"]
            mask_i = perm == type_i
            mask_j = perm == type_j

            if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
                continue

            coords_i = coords_core[mask_i]      # extract coordinates of source cells of type i
            coords_j = coords_core[mask_j]      # extract coordinates of target cells of type j
            weights_i = weights_core[mask_i]    # extract edge weights of source cells of type i

            lambda_j = len(coords_j) / polygon.area     # normalization factor for Ripley's K/L to account for different densities of target cell type j
            L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)     # calculate Ripley's L for the source cells of type i and target cells of type j

            # average across source cells → one curve per interaction
            L_mean = np.mean(L, axis=0)   # shape: (n_r,)
            sim_L[s, idx_inter, :] = L_mean
            sim_valid_inter[idx_inter] += 1

    # create a matrix with average L for each simulation for each interaction for each radius to be able to compare with observed L later
    sim_mean_sample = np.nanmean(sim_L, axis=0)   # (interactions × radii)
    sim_std_sample  = np.nanstd(sim_L, axis=0)    # (interactions × radii)


    print(f"[DONE] Sample {sample_id}")
    return global_idx, obs_local, sim_mean_sample, sim_std_sample, sim_valid_inter

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Run in parallel
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sample_list = np.unique(samples)
results = Parallel(n_jobs=20)(delayed(process_sample)(s) for s in tqdm(sample_list, desc='Processing samples'))

# Assign and save results
#----------------------------------------------------------------------------
obs_curve = np.full((n_cells, n_interactions, n_r), np.nan)

# create dictionaries to store mean and std of simulated L values for each sample 
sim_mean_dict = {}
sim_std_dict = {}
sim_valid_dict = {}

for sample_id, (idx, o, sm, ss, sv) in zip(sample_list, results):
    obs_curve[idx] = o      # (cells × interactions × radii)
    sim_mean_dict[sample_id] = sm   # (interactions × radii)
    sim_std_dict[sample_id] = ss    # (interactions × radii)
    sim_valid_dict[sample_id] = sv  # (interactions)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Statistics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Z-score
z_scores_unfiltered = np.full_like(obs_curve, np.nan)
invalid_mask = np.zeros_like(z_scores_unfiltered, dtype=bool)

for sample_id in sample_list:
    mask = samples == sample_id
    sm = sim_mean_dict[sample_id]   # (interactions × radii)
    ss = sim_std_dict[sample_id]
    z_scores_unfiltered[mask] = (obs_curve[mask] - sm[None, :, :]) / ss[None, :, :] #!!! prevent division by very small number to prevent exploding of the z score

    invalid_mask[mask] = (np.isnan(obs_curve[mask]) | (ss[None,:,:] <= 1e-3))  # mark as invalid if observed curve is NaN or if simulated std is too low (indicating low simulation support)

z_scores = z_scores_unfiltered.copy()
z_scores[invalid_mask] = np.nan


# integrals
diff = np.full_like(obs_curve, np.nan)
for sample_id in sample_list:
    mask = samples == sample_id
    sm = sim_mean_dict[sample_id]
    diff[mask] = obs_curve[mask] - sm[None, :, :]

signed = np.trapezoid(diff, radii, axis=2)
absolute = np.trapezoid(np.abs(diff), radii, axis=2)

# Store in adata
#----------------------------------------------------------------------------
adata.obsm["ripley_z"] = z_scores.astype(np.float32)
adata.obsm["ripley_signed"] = signed.astype(np.float32)
adata.obsm["ripley_abs"] = absolute.astype(np.float32)
adata.obsm["ripley_z_max"] = np.nanmax(np.abs(z_scores), axis=2).astype(np.float32)
adata.obsm["ripley_z_mean"] = np.nanmean(z_scores, axis=2).astype(np.float32)

adata.uns["ripley_obs_curve"] = obs_curve.astype(np.float32)
adata.uns["ripley_sim_mean"] = sim_mean_dict
adata.uns["ripley_sim_std"] = sim_std_dict
adata.uns["ripley_sim_valid"] = sim_valid_dict
adata.uns["ripley_interactions"] = interactions
adata.uns["ripley_params"] = {
    "radii": radii,
    "n_sim": args.n_sim,
    "min_i_threshold": min_i_threshold,
    "min_j_threshold": min_j_threshold,
    "method": "interaction-level null, per-sample normalization"
}

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5 Save the analyzed adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("\nSaving...")
print(adata)
adata.write(args.output)
print("Done.")