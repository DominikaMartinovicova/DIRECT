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
    parser.add_argument('--n_sim', type=int, default=10, help='number of simulations to obtain null distribution')
    return parser.parse_args()

args = parse_args()
#core = 'T23_004535_110005_1'
adata = sc.read_h5ad(args.input)
#adata = adata[adata.obs['sample']==core].copy()
print(adata)
cluster_key = args.phen_level
print(adata.obs[cluster_key].value_counts())
sample_key = args.sample_key
coi_list = args.coi
#coi_list=["B_cell", "T_cell_regulatory"]
print(f'{len(coi_list)} celltypes of interest  {coi_list}.')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Set variables to be used throughout calculation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# minimum number of cells in source and target to continue with ripley's stat calculation 
min_i_threshold = 3
min_j_threshold = 3

rng = default_rng(42)   # set random seed for reproducibility

coordinates = np.asarray(adata.obsm["spatial"])     # extract coordinates
print(f'Coordinates: {coordinates.shape}, {coordinates[:10,:10]}')
labels = np.asarray(adata.obs[cluster_key])         # extract labels 
print(f'Labels: {labels.shape}, {labels[:10]}')
samples = np.asarray(adata.obs[sample_key])
print(f'Samples: {samples.shape}, {samples[:10]}')


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
    weights = np.ones((len(coords), len(radii)))
    for k, r in enumerate(radii):
        mask = dists < r
        weights[mask, k] = np.minimum(r / np.maximum(dists[mask], 1e-6), max_weight)
    return weights

# Calculate cross-type Ripley's L function
#--------------------------------------------------------------------------------
def cross_L(coords_i, coords_j, radii, weights_i, lambda_j):
    tree = cKDTree(coords_j)        # identify all type j cells
    neighbors = tree.query_ball_point(coords_i, radii[-1])  # identify j neighbors of source cells within max radius
    #print(f'Neighbors: {neighbors}, shape {neighbors.shape}')
    L = np.zeros((len(coords_i), len(radii)))       # create np array with 0s to store results
    n_skipped=0
    for i, idx in enumerate(neighbors):
        if len(idx)==0:
            n_skipped += 1
            #print(f'No neighbors found within 250 radius. Skipping source cell {i}...')
            continue
        dists = np.linalg.norm(coords_j[idx] - coords_i[i], axis=1)
        dists.sort()
        counts = np.searchsorted(dists, radii)
        L[i, :] = (counts * weights_i[i]) / lambda_j
    #print(f'Skipped {n_skipped} cells.')
    #
    # print(L)
    #print(f'Number of nans in L: {np.sum(np.isnan(L))}')
    return np.sqrt(L / np.pi)   # division makes this statistic into Ripley's L (before it was K even if variable was called L)

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
    #print(f"[WEIGHTS] Sample {sample_id} | {weights_core.shape} shape of weights_core; max weight: {np.nanmax(weights_core)}; min weight: {np.nanmin(weights_core)}; NaNs: {np.sum(np.isnan(weights_core))}")

    n_local = len(coords_core)

    obs_local = np.full((n_local, n_interactions, n_r), np.nan)
    sim_sum = np.zeros_like(obs_local)
    sim_sq = np.zeros_like(obs_local)
    sim_valid = np.zeros((n_local, n_interactions))

    # ---- observed
    for type_i, type_j in product(coi_list, coi_list):
        # if type_i == type_j:
        #     continue

        idx_inter = interaction_index[f"{type_i}_{type_j}"]

        mask_i = labels_core == type_i
        mask_j = labels_core == type_j
        #print(f"[DEBUG] Sample {sample_id} | {type_i}->{type_j} | n_source={np.sum(mask_i)}, n_target={np.sum(mask_j)}")

        if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
            #print(f"[DEBUG] Skipping {type_i}->{type_j} due to min threshold")
            continue

        coords_i = coords_core[mask_i]
        coords_j = coords_core[mask_j]
        weights_i = weights_core[mask_i]

        lambda_j = len(coords_j) / polygon.area
        #print(f'Calculating cross type L for {type_i}->{type_j}')
        L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)
        #print(f'[OBSERVED] {type_i} -> {type_j} | L shape: {L.shape}, NaNs in L: {np.sum(np.isnan(L))}')

        neighbor_fraction = np.sum(~np.isnan(L).any(axis=1)) / len(L)
        #print(f"[COVERAGE] {type_i} -> {type_j} | fraction of source cells with ≥1 neighbor: {neighbor_fraction:.2%}")

        obs_local[mask_i, idx_inter, :] = L
        #print(f'obs_local shape: {obs_local.shape}')
        #print(f'Number of nans in obs_local: {np.sum(np.isnan(obs_local))}')
    # ---- simulations
    for s in range(args.n_sim):
        if s % 100 == 0:
            print(f"[{sample_id}] simulation {s}/{args.n_sim}")
        perm = rng.permutation(labels_core)

        for type_i, type_j in product(coi_list, coi_list):
            # if type_i == type_j:
            #     continue

            idx_inter = interaction_index[f"{type_i}_{type_j}"]

            mask_i = perm == type_i
            mask_j = perm == type_j

            if np.sum(mask_i) < min_i_threshold or np.sum(mask_j) < min_j_threshold:
                #print(f"[DEBUG] Skipping simulation {type_i}->{type_j} due to min threshold")
                continue

            coords_i = coords_core[mask_i]
            coords_j = coords_core[mask_j]
            weights_i = weights_core[mask_i]

            lambda_j = len(coords_j) / polygon.area

            L = cross_L(coords_i, coords_j, radii, weights_i, lambda_j)

            sim_sum[mask_i, idx_inter, :] += L
            sim_sq[mask_i, idx_inter, :] += L**2
            sim_valid[mask_i, idx_inter] += 1
            #print(f'Simulation {s}')
            #print(f'sim_sum: \n {sim_sum[:5,:5,0]}, \n sim_sq: \n {sim_sq[:5,:5,0]}, \n sim_valid: \n {sim_valid[:5,:5]}')

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
    print(f'obs_curve shape: {obs_curve.shape}, sim_sum shape: {sim_sum.shape}, sim_sq shape: {sim_sq.shape}, sim_valid shape: {sim_valid.shape}')
    print(obs_curve[:5,:5,0])
    print(sim_sum[:5,:5,0])
    print(sim_sq[:5,:5,0])
    print(sim_valid[:5,:5])
    


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Statistics
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
valid_expanded = np.maximum(sim_valid[..., None], 1)
print(f'valid_expanded shape: {valid_expanded.shape}, sim_valid shape: {sim_valid.shape}')

sim_mean = sim_sum / valid_expanded
print(f'sim_mean shape: {sim_mean.shape}, sim_mean[:5,:5,0]: \n {sim_mean[:5,:5,0]}')
sim_var = (sim_sq / valid_expanded) - sim_mean**2
print(f'sim_var shape: {sim_var.shape}, sim_var[:5,:5,0]: \n {sim_var[:5,:5,0]}')
sim_std = np.sqrt(np.maximum(sim_var, 1e-4))
print(f'sim_std shape: {sim_std.shape}, sim_std[:5,:5,0]: \n {sim_std[:5,:5,0]}')

# Z-score
z_scores = (obs_curve - sim_mean) / sim_std

# filter low support simulations and invalid values
min_valid = int(args.n_sim * 0.3)

invalid_mask = ((sim_valid == 0)[..., None] |np.isnan(obs_curve) |(sim_std < 1e-6)|(sim_valid < min_valid)[..., None])
print(f'Invalid mask shape: {invalid_mask.shape}, sim_valid shape: {sim_valid.shape}, obs_curve shape: {obs_curve.shape}, sim_std shape: {sim_std.shape}')
print(f'Number of invalid entries: {np.sum(invalid_mask)}, total entries: {invalid_mask.size}, fraction invalid: {np.sum(invalid_mask)/invalid_mask.size:.2%}')
print(f'Invalid mask sample: \n {invalid_mask[:5,:5,0]}')
z_scores[invalid_mask] = np.nan

# plot distribution of z-scores and sim std to check if they look reasonable
n_r = sim_std.shape[-1]

fig, axes = plt.subplots(1, n_r, figsize=(4*n_r, 4), sharey=True)

for k in range(n_r):
    ax = axes[k] if n_r > 1 else axes

    std_k = sim_std[..., k].flatten()
    valid_k = sim_valid.flatten()

    mask = ~np.isnan(std_k)
    
    std_filtered = std_k[mask]
    valid_filtered = valid_k[mask]

    # 5. Proceed with log and plotting using the filtered versions
    #log_std = np.log10(std_filtered + 1e-10)
    high_support = valid_filtered >= min_valid
    low_support = valid_filtered < min_valid

    # plot
    ax.hist(std_filtered[high_support], bins=100, alpha=0.7,color='blue', label='High support')
    ax.hist(std_filtered[low_support], bins=100, alpha=0.7, color='red', label='Low support')

    ax.set_yscale('log')
    ax.set_title(f'Radius = {radii[k]} µm')
    ax.set_xlabel('Standard deviation')

    if k == 0:
        ax.set_ylabel('Frequency (log scale)')

    ax.legend()

plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N_old/spatial/patching/5000um_50um/sim_std_per_radius.png', dpi=300)
plt.close()


fig, axes = plt.subplots(1, n_r, figsize=(4*n_r, 4), sharey=True)

for k in range(n_r):
    ax = axes[k] if n_r > 1 else axes

    std_k = sim_std[..., k].flatten()
    valid_k = sim_valid.flatten()

    mask = ~np.isnan(std_k)
    
    std_filtered = std_k[mask]
    valid_filtered = valid_k[mask]

    log_std = np.log10(std_filtered + 1e-10)
    high_support = valid_filtered >= min_valid
    low_support = valid_filtered < min_valid

    ax.hist(log_std[high_support], bins=100, alpha=0.7,color='blue', label='High support')
    ax.hist(log_std[low_support], bins=100, alpha=0.7,color='red', label='Low support')

    ax.set_yscale('log')
    ax.set_title(f'Radius = {radii[k]} µm')
    ax.set_xlabel('log10(std)')

    if k == 0:
        ax.set_ylabel('Frequency (log scale)')

    ax.legend()

plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N_old/spatial/patching/5000um_50um/sim_std_log_per_radius.png', dpi=300)
plt.close()



# integrals
diff = obs_curve - sim_mean
signed = np.trapezoid(diff, radii, axis=2)
absolute = np.trapezoid(np.abs(diff), radii, axis=2)

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