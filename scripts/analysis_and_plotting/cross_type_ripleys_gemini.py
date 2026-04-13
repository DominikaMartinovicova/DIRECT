#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# cross_type_ripleys_per_cell.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Calculate Ripley's K function for cross-type interactions in spatial data.
#   Convert K to centered L (L(r) - r).
#
#   0 Import libraries and parse arguments
#   1 Set variables
#   2 Define analysis functions
#       a. Edge correction (Isotropic)
#       b. Cross type Ripley's K 
#       c. Process cells per sample
#   3 Main loop - call process sample function for each sample in adata
#   4 Perform statistical analysis and store results in adata
#   5 Save adata with Ripley's results
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
# Updated to mirror R/spatstat mathematical standards.
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

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
    parser.add_argument('--n_sim', type=int, default=500, help='number of simulations to obtain null distribution')
    return parser.parse_args()

args = parse_args()
core = 'T23_004535_110005_1'
adata = sc.read_h5ad(args.input)
adata = adata[adata.obs['sample']==core].copy()
print(adata)
cluster_key = args.phen_level
print(adata.obs[cluster_key].value_counts())
sample_key = args.sample_key
coi_list = args.coi
#coi_list=["B_cell", "Macrophage", "Macrophage_alveolar", "T_cell_regulatory"]
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

#radii = np.linspace(0, args.max_dist, args.n_steps)    # create list of radii to be tested
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
def compute_edge_weights(coords, polygon, radii, max_weight=5.0):
    dists = np.array([polygon.boundary.distance(Point(p)) for p in coords])
    radii_arr = np.array(radii)
    weights = np.ones((len(coords), len(radii_arr)))    # default weight is one
    for k, r in enumerate(radii_arr):
        mask = dists < r                            # identify cells closer to border than radius
        if np.any(mask):
            # Isotropic correction: proportion of circumference inside the boundary
            d_over_r = np.clip(dists[mask] / r, -1.0, 1.0)
            prop_inside = 1.0 - (1.0 / np.pi) * np.arccos(d_over_r)
            prop_inside = np.maximum(prop_inside, 1.0 / max_weight) # clamp to prevent extreme weights
            weights[mask, k] = 1.0 / prop_inside
    return weights

# Calculate cross-type Ripley's K function (Transformed to L in process_sample)
#--------------------------------------------------------------------------------
def cross_K(coords_i, coords_j, radii, weights_i, lambda_j):
    tree = cKDTree(coords_j)        # identify all type j cells
    radii_arr = np.array(radii)
    neighbors = tree.query_ball_point(coords_i, radii_arr[-1])  # identify j neighbors of source cells within max radius
    K = np.zeros((len(coords_i), len(radii_arr)))       # create np array with 0s to store results
    for i, idx in enumerate(neighbors):
        if len(idx)==0:
            continue
        dists = np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)     # calculate distances
        dists.sort()      # sort distances
        counts = np.searchsorted(dists, radii_arr, side='right')       # 'right' ensures <= r standard
        K[i, :] = (counts * weights_i[i]) / lambda_j    # apply edge correction and normalize
    return K   # return raw K to allow proper mathematical aggregation

# Process sample
#--------------------------------------------------------------------------------
def process_sample(sample_id):   
    mask = samples == sample_id     # find cells belonging to the sample being processed
    coords_core = coordinates[mask] # extract coordinates of the sample being processed
    labels_core = labels[mask]      # extract labels of the sample being processed
    global_idx = np.where(mask)[0]  # store global indices
    radii_arr = np.array(radii)

    print(f"[START] Sample {sample_id} | {len(coords_core)} cells")

    hull = ConvexHull(coords_core)  # calculate convex hull
    polygon = Polygon(coords_core[hull.vertices])   # create polygon

    weights_core = compute_edge_weights(coords_core, polygon, radii_arr)   # calculate edge weights
    n_local = len(coords_core)  # number of cells in the sample

    # Create arrays to store results
    obs_local = np.full((n_local, n_interactions, n_r), np.nan)
    sim_L = np.full((args.n_sim, n_interactions, n_r), np.nan)  
    sim_valid_inter = np.zeros(n_interactions)   
    
    # ----------------- Calculate observed Ripley's L for the sample being processed ----------------------
    for type_i, type_j in product(coi_list, coi_list):      # iterate through all pairs
        idx_inter = interaction_index[f"{type_i}_{type_j}"] 
        mask_i = labels_core == type_i      
        mask_j = labels_core == type_j      

        if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
            continue    # skip if too few cells

        coords_i = coords_core[mask_i]      
        coords_j = coords_core[mask_j]      
        weights_i = weights_core[mask_i]    

        lambda_j = len(coords_j) / polygon.area     # normalization factor
        
        # Calculate K, then transform to centered L(r) - r for each cell
        K_obs = cross_K(coords_i, coords_j, radii_arr, weights_i, lambda_j)     
        L_centered = np.sqrt(K_obs / np.pi) - radii_arr

        obs_local[mask_i, idx_inter, :] = L_centered     # store the observed centered L values

    print(np.nanmean(obs_local[:,idx_inter,:],axis=0))   
    
    # ----------------- Perform simulations to get null distribution of Ripley's L ----------------------
    for s in range(args.n_sim):
        if s % 100 == 0:
            print(f"[{sample_id}] simulation {s}/{args.n_sim}")
        perm = rng.permutation(labels_core)      # shuffle cell labels

        for type_i, type_j in product(coi_list, coi_list):
            idx_inter = interaction_index[f"{type_i}_{type_j}"]
            mask_i = perm == type_i
            mask_j = perm == type_j

            if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
                continue

            coords_i = coords_core[mask_i]      
            coords_j = coords_core[mask_j]      
            weights_i = weights_core[mask_i]    

            lambda_j = len(coords_j) / polygon.area     
            
            # Calculate K, aggregate to global K to fix Jensen's inequality, then transform
            K_sim = cross_K(coords_i, coords_j, radii_arr, weights_i, lambda_j)     
            K_mean = np.mean(K_sim, axis=0)   # average K across source cells -> shape: (n_r,)
            L_mean_centered = np.sqrt(K_mean / np.pi) - radii_arr
            
            sim_L[s, idx_inter, :] = L_mean_centered
            sim_valid_inter[idx_inter] += 1

    sim_mean_sample = np.nanmean(sim_L, axis=0)   # (interactions × radii)
    sim_std_sample  = np.nanstd(sim_L, axis=0)    # (interactions × radii)

    print(sim_mean_sample)
    print(sim_std_sample)

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
min_std = 0.1  # minimum std to consider the simulation valid

for sample_id in sample_list:
    mask = samples == sample_id
    sm = sim_mean_dict[sample_id]   # (interactions × radii)
    ss = sim_std_dict[sample_id]
    low_std_mask = (ss <= min_std) | (~np.isfinite(ss))  
    safe_std = np.where(low_std_mask, np.nan, ss)
    z_scores_unfiltered[mask] = (obs_curve[mask] - sm[None, :, :]) / safe_std[None, :, :] 

    invalid_mask[mask] = (np.isnan(obs_curve[mask]) | low_std_mask[None,:,:])  

z_scores = z_scores_unfiltered.copy()
z_scores[invalid_mask] = np.nan

# integrals
diff = np.full_like(obs_curve, np.nan)
for sample_id in sample_list:
    mask = samples == sample_id
    sm = sim_mean_dict[sample_id]
    diff[mask] = obs_curve[mask] - sm[None, :, :]

diff[invalid_mask] = np.nan  
signed = np.trapezoid(diff, radii, axis=2)
absolute = np.trapezoid(np.abs(diff), radii, axis=2)

# Store in adata
#----------------------------------------------------------------------------
adata.obsm["ripley_z"] = z_scores.astype(np.float32)
adata.obsm["ripley_z_unfiltered"] = z_scores_unfiltered.astype(np.float32)
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
    "method": "interaction-level null, per-sample normalization, centered L(r)-r"
}

print("NaN fraction in Z:", np.mean(np.isnan(z_scores)))
print("Min/Max Z:", np.nanmin(z_scores), np.nanmax(z_scores))

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5 Save the analyzed adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("\nSaving...")
print(adata)
adata.write(args.output)
print("Done.")