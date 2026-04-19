#!/usr/bin/env python3

import numpy as np
import scanpy as sc
from scipy.spatial import cKDTree, ConvexHull
from shapely.geometry import Polygon, Point
from itertools import product
from joblib import Parallel, delayed
from tqdm import tqdm
import argparse

# =============================================================================
# 0. Arguments
# =============================================================================
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('--phen_level', required=True)
    parser.add_argument('--coi', nargs='+', required=True)
    parser.add_argument('--sample_key', default="sample")
    parser.add_argument('--radii', nargs='+', type=float, required=True)
    return parser.parse_args()

args = parse_args()

adata = sc.read_h5ad(args.i)
core = 'T23_004535_110005_1'
adata = adata[adata.obs['sample']==core].copy()

cluster_key = args.phen_level
sample_key = args.sample_key
coi_list = args.coi
radii = np.array(sorted(args.radii))

coords = np.asarray(adata.obsm["spatial"])
labels = np.asarray(adata.obs[cluster_key])
samples = np.asarray(adata.obs[sample_key])

interactions = [(i, j) for i, j in product(coi_list, coi_list)]

# =============================================================================
# 1. Window + area
# =============================================================================
def compute_window(coords):
    hull = ConvexHull(coords)
    polygon = Polygon(coords[hull.vertices])
    return polygon, polygon.area

# =============================================================================
# 2. Edge correction (pair-based)
# =============================================================================
def edge_weight(point, polygon, r):
    d = polygon.boundary.distance(Point(point))
    if d >= r:
        return 1.0
    ratio = np.clip(d / r, -1, 1)
    prop = 1 - (1 / np.pi) * np.arccos(ratio)
    return 1.0 / max(prop, 1e-6)

# =============================================================================
# 3. Compute close pairs (like spatstat::closepairs)
# =============================================================================
def compute_close_pairs(coords, labels, r_max, polygon):
    tree = cKDTree(coords)
    pairs = []

    neighbors = tree.query_ball_tree(tree, r_max)

    for i, neigh in enumerate(neighbors):
        for j in neigh:
            # include self-pairs like R (distinct=FALSE)
            d = np.linalg.norm(coords[i] - coords[j])
            if d <= r_max:
                pairs.append((i, j, d))

    # convert to structured arrays
    i_idx = np.array([p[0] for p in pairs])
    j_idx = np.array([p[1] for p in pairs])
    dists = np.array([p[2] for p in pairs])

    cell_i = labels[i_idx]
    cell_j = labels[j_idx]

    # edge weights per pair (based on i)
    edge = np.array([edge_weight(coords[i], polygon, r_max) for i in i_idx])

    return i_idx, j_idx, dists, cell_i, cell_j, edge

# =============================================================================
# 4. Ripley's L (MATCHES R)
# =============================================================================
def L_function(i_idx, j_idx, dists, cell_i, cell_j, edge,
               child1, child2, r, area):

    # filter pairs within r
    mask_r = dists <= r

    i_idx = i_idx[mask_r]
    j_idx = j_idx[mask_r]
    cell_i = cell_i[mask_r]
    cell_j = cell_j[mask_r]
    edge = edge[mask_r]

    # counts
    mask_child1 = cell_i == child1
    mask_child2 = cell_i == child2

    n_child1 = np.sum(mask_child1)
    n_child2 = np.sum(cell_i == child2)

    if n_child1 == 0 or n_child2 == 0:
        return np.nan

    # numerator: weighted pairs
    numerator = np.sum(
        edge[(cell_i == child1) & (cell_j == child2)]
    )

    # lambda
    lambda2 = n_child2 / area

    # K
    K = (numerator / lambda2) * (1 / n_child1)

    # L
    L = np.sqrt(K / np.pi) - r

    return L

# =============================================================================
# 5. Per-sample computation
# =============================================================================
def process_sample(sample_id):

    mask = samples == sample_id
    coords_s = coords[mask]
    labels_s = labels[mask]

    polygon, area = compute_window(coords_s)

    r_max = np.max(radii)

    i_idx, j_idx, dists, cell_i, cell_j, edge = compute_close_pairs(
        coords_s, labels_s, r_max, polygon
    )

    results = {}

    for (type_i, type_j) in interactions:

        L_vals = []

        for r in radii:
            L = L_function(
                i_idx, j_idx, dists,
                cell_i, cell_j, edge,
                type_i, type_j, r, area
            )
            L_vals.append(L)

        results[f"{type_i}__{type_j}"] = np.array(L_vals)

    return sample_id, results

# =============================================================================
# 6. Run
# =============================================================================
sample_list = np.unique(samples)

out = Parallel(n_jobs=-1)(
    delayed(process_sample)(s)
    for s in tqdm(sample_list)
)

# =============================================================================
# 7. Store results
# =============================================================================
adata.uns["ripley_L"] = {sid: res for sid, res in out}
adata.uns["ripley_radii"] = radii

adata.write(args.o)