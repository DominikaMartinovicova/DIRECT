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
    parser.add_argument('--n_sim', type=int, default=500, help='number of smulations to obtain null distribution')
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
labels = np.asarray(adata.obs[cluster_key])         # extract labels 
samples = np.asarray(adata.obs[sample_key])

#radii = np.linspace(0, args.max_dist, args.n_steps)     # create list of radii to be tested
radii = [25, 50, 75, 100, 250]
n_r = len(radii)
print(f'Radii to be tested: {radii}')

n_cells = adata.n_obs
interactions = [f"{i}_{j}" for i, j in product(coi_list, coi_list) if i != j]
interaction_index = {k: i for i, k in enumerate(interactions)}
n_interactions = len(interactions)
print(f'Number of celltype pairs to be tested: {n_interactions}')

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
    weights = np.ones((len(coords), len(radii)))
    for k, r in enumerate(radii):
        mask = dists < r
        weights[mask, k] = np.minimum(r / np.maximum(dists[mask], 1e-6), max_weight)
    return weights

# Calculate cross-type Ripley's L function
#--------------------------------------------------------------------------------
def cross_L(coords_i, coords_j, radii, weights_i, lambda_j):
    tree = cKDTree(coords_j)
    neighbors = tree.query_ball_point(coords_i, radii[-1])
    L = np.zeros((len(coords_i), len(radii)))
    for i, idx in enumerate(neighbors):
        if not idx:
            continue
        dists = np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)
        dists.sort()
        counts = np.searchsorted(dists, radii)
        L[i, :] = (counts * weights_i[i]) / lambda_j
    return np.sqrt(L / np.pi)

# Process sample
#--------------------------------------------------------------------------------
def process_sample(sample_id):   
    mask = samples == sample_id
    coords_core = coordinates[mask]
    labels_core = labels[mask]
    global_idx = np.where(mask)[0]

    print(f"[START] Sample {sample_id} | {len(coords_core)} cells")

    hull = ConvexHull(coords_core)
    polygon = Polygon(coords_core[hull.vertices])

    weights_core = compute_edge_weights(coords_core, polygon, radii)

    n_local = len(coords_core)

    obs_local = np.full((n_local, n_interactions, n_r), np.nan)
    sim_sum = np.zeros_like(obs_local)
    sim_sq = np.zeros_like(obs_local)
    sim_valid = np.zeros((n_local, n_interactions))

    # ---- observed
    for type_i, type_j in product(coi_list, coi_list):
        if type_i == type_j:
            continue

        idx_inter = interaction_index[f"{type_i}_{type_j}"]

        mask_i = labels_core == type_i
        mask_j = labels_core == type_j

        if np.sum(mask_i) < 3 or np.sum(mask_j) < 3:
            continue

        coords_i = coords_core[mask_i]
        coords_j = coords_core[mask_j]
        weights_i = weights_core[mask_i]

        lambda_j = len(coords_j) / polygon.area

        L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)

        obs_local[mask_i, idx_inter, :] = L

    # ---- simulations
    for s in range(args.n_sim):
        if s % 100 == 0:
            print(f"[{sample_id}] simulation {s}/{args.n_sim}")
        perm = rng.permutation(labels_core)

        for type_i, type_j in product(coi_list, coi_list):
            if type_i == type_j:
                continue

            idx_inter = interaction_index[f"{type_i}_{type_j}"]

            mask_i = perm == type_i
            mask_j = perm == type_j

            if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
                continue

            coords_i = coords_core[mask_i]
            coords_j = coords_core[mask_j]
            weights_i = weights_core[mask_i]

            lambda_j = len(coords_j) / polygon.area

            L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)

            sim_sum[mask_i, idx_inter, :] += L
            sim_sq[mask_i, idx_inter, :] += L**2
            sim_valid[mask_i, idx_inter] += 1

    print(f"[DONE] Sample {sample_id}")
    return global_idx, obs_local, sim_sum, sim_sq, sim_valid

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Run in parallel
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sample_list = np.unique(samples)
results = Parallel(n_jobs=20)(delayed(process_sample)(s) for s in tqdm(sample_list, desc='Processing samples'))

# Merge results
#----------------------------------------------------------------------------
obs_curve = np.full((n_cells, n_interactions, n_r), np.nan)
sim_sum = np.zeros_like(obs_curve)
sim_sq = np.zeros_like(obs_curve)
sim_valid = np.zeros((n_cells, n_interactions))

for idx, o, ss, sq, sv in results:
    obs_curve[idx] = o
    sim_sum[idx] = ss
    sim_sq[idx] = sq
    sim_valid[idx] = sv


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Statistics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
valid_expanded = np.maximum(sim_valid[..., None], 1)

sim_mean = sim_sum / valid_expanded
sim_var = (sim_sq / valid_expanded) - sim_mean**2
sim_std = np.sqrt(np.maximum(sim_var, 1e-8))

# Z-score
z_scores = (obs_curve - sim_mean) / sim_std

invalid_mask = ((sim_valid == 0)[..., None] |np.isnan(obs_curve) |(sim_std < 1e-6))

z_scores[invalid_mask] = np.nan



# integrals
diff = obs_curve - sim_mean
signed = np.trapz(diff, radii, axis=2)
absolute = np.trapz(np.abs(diff), radii, axis=2)

# Store in adata
#----------------------------------------------------------------------------
adata.obsm["ripley_signed"] = signed
adata.obsm["ripley_abs"] = absolute
adata.obsm["ripley_sim_valid"] = sim_valid

adata.obsm["ripley_z_max"] = np.nanmax(np.abs(z_scores), axis=2)
adata.obsm["ripley_z_mean"] = np.nanmean(z_scores, axis=2)

adata.uns["ripley_z"] = z_scores.astype(np.float32)
adata.uns["ripley_interactions"] = interactions
adata.uns["ripley_params"] = {
    "radii": radii,
    "n_sim": args.n_sim,
    "method": "optimized cell-level Ripley L with fast isotropic correction"
}

obs_curve = obs_curve.astype(np.float32)
sim_mean = sim_mean.astype(np.float32)
sim_std = sim_std.astype(np.float32)

adata.uns["ripley_obs_curve"] = obs_curve
adata.uns["ripley_sim_mean"] = sim_mean
adata.uns["ripley_sim_std"] = sim_std

###
# adata.obsm["ripley_signed"] = signed
# adata.obsm["ripley_abs"] = absolute

# adata.obsm["ripley_sim_valid"] = sim_valid
# adata.obsm["ripley_sim_failed"] = sim_failed
# adata.obsm["ripley_failed_sd"] = failed_sd

# adata.uns["ripley_interactions"] = interactions
# adata.uns["ripley_params"] = {
#     "radii": radii,
#     "n_sim": args.n_sim,
#     "method": "cell-level cross-type Ripley L with isotropic correction",
# }

# adata.uns["ripley_z"] = z_scores
# adata.uns["ripley_z_info"] = {
#     "dimensions": ["cells", "interactions", "radii"],
#     "note": "Z-score of observed vs simulated Ripley's L"
# }
###
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 5 Save the analyzed adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("\nSaving...")
print(adata)
adata.write(args.output)
print("Done.")