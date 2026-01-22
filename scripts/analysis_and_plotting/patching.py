#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# patching.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Divide each piece of tissue into patches. 
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Divide each piece of tissue into patches
#   3 Create patch adata (patches as observations)
#       - Filter patches with fewer than 20 cells
#       a. adata features - mean gene expression
#       b. adata features - cell type fractions
#       (c. adata features - spatial stats (save adata for each patch to be able to paralellize spatial statistics calculation in another script))
#   4 Save patches individually as adata 
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
import warnings
warnings.filterwarnings("ignore")

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 patching.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Divide each piece of tissue into patches.  ')
    parser.add_argument('-i', help='path to adata sample subset',
                        dest='input',
                        type=str)
    parser.add_argument('--patch_size', help='size of the patches in microns',
                        dest='patch_size',
                        type=int)
    parser.add_argument('--overlap', help='overlap between patches in microns',
                        dest='overlap',
                        type=int)
    parser.add_argument('--celltype_key', help='key for cell type annotation in adata.obs',
                        dest='celltype_key',
                        type=str)
    parser.add_argument('-o', '--output_patches', help='path to output dir with patches per sample',
                        dest='output_dir_patches',
                        type=str)
    parser.add_argument('--plots_dir', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
overlap = args.overlap
p_size = args.patch_size
celltype_key = args.celltype_key
os.makedirs(args.output_dir_patches, exist_ok=True)
os.makedirs(args.output_dir_plots, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad(args.input)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Divide each piece of tissue into patches
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print('Dividing each piece of tissue into patches...')
sq.tl.sliding_window(adata=adata, library_key="sample", window_size=p_size, sliding_window_key='patch', overlap=overlap, copy=False)
if overlap == 0:
    adata.obs['patch'] = adata.obs['patch'].astype('category')  # Convert to categorical
elif overlap > 0:
    patch_cols = adata.obs.columns[adata.obs.columns.str.contains("window")]
    for col in patch_cols:
        adata.obs[col] = adata.obs[col].astype(bool)  # Convert to boolean

print('Writing adata with sliding window assignment...')
adata.write_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{celltype_key}_combined_adatas_w_patches_{p_size}_{overlap}.h5ad')

# Plot spatial scatter with patches
#print('Plotting spatial scatter with patches...')
#sq.pl.spatial_scatter(adata, color="pt_id", library_id="spatial", shape=None)
#plt.legend().remove()
#plt.savefig(os.path.join(args.output_dir_plots,f'overlap_{overlap}/spatial_scatter_patches.png'))
#plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create patch adata (patches as observations)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if overlap == 0:
    print('No overlap between patches. Patch size:', args.patch_size, 'Overlap:', overlap)
    patches = adata.obs['patch'].unique()
    print(f'Total number of patches: {len(patches)}')

    # Remove patches with too few cells
    #--------------------------------------------------------------------------------
    patches_to_keep = []
    for patch in patches:
        adata_patch = adata[adata.obs['patch']==patch]
        number_of_cells = adata_patch.n_obs
        if number_of_cells > 20 and adata_patch.obs[celltype_key].nunique() > 1:
            patches_to_keep.append(patch)
            adata_patch.obs.drop(columns=patch_cols, inplace=True)
            adata_patch.write_h5ad(os.path.join(args.output_dir_patches,f'{patch}.h5ad'))
        else:
            with open(os.path.join(args.output_dir_patches,"skipped_patches.txt"), "a") as f:
                f.write(f"Skipped {patch}: num cells={number_of_cells}, num cell types={adata_patch.obs[celltype_key].nunique()}\n")
    print(f'Number of patches after filtering for min number of cells (>20) and min number of cell types (>1): {len(patches_to_keep)}')

    # Mean gene expression as features
    #--------------------------------------------------------------------------------
    print('---------------------- Creating adata with mean gene expression per patch... ----------------------------')
    X_patch_list = []
    obs_patch_list = []

    for patch in patches_to_keep:
        index = (adata.obs['patch']==patch).to_numpy()   # Find cells belonging to the same patch
        X_patch_list.append(adata.layers['raw_counts'][index].mean(axis=0).A1)        # For each patch get mean expression of all genes across all cells in the patch
        obs_patch_list.append(adata.obs[index].iloc[0])         # For each patch get obs info from the first cell in the patch

    X_patch = np.vstack(X_patch_list)
    adata_patches_gex = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=adata.var.copy())
    adata_patches_gex.obs["patch"] = adata_patches_gex.obs["patch"].astype("category")


    # Cell type fractions as features
    #--------------------------------------------------------------------------------
    print('---------------------- Creating adata with cell type fractions per patch... ----------------------------')
    X_patch_list = []
    obs_patch_list = []

    for patch in patches_to_keep:
        index = (adata.obs['patch']==patch).to_numpy()   # Find cells belonging to the same patch
        ct_counts = adata.obs[index][celltype_key].value_counts(normalize=True)  # For each patch get cell type fractions
        ct_fractions = ct_counts.reindex(adata.obs[celltype_key].cat.categories, fill_value=0).values  # Ensure all cell types are represented in the same order
        X_patch_list.append(ct_fractions)
        obs_patch_list.append(adata.obs[index].iloc[0])         # For each patch get obs info from the first cell in the patch

    X_patch = np.vstack(X_patch_list)
    adata_patches_ct = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=pd.DataFrame(index=adata.obs[celltype_key].cat.categories))
    adata_patches_ct.obs["patch"] = adata_patches_ct.obs["patch"].astype("category")
    print(adata_patches_ct)

elif overlap > 0:
    print('Patches are overlapping. Patch size:', p_size, 'Overlap:', overlap)
    patches = adata.obs.columns[adata.obs.columns.str.contains("patch")]
    print(f'Total number of patches: {len(patches)}')
    
    # Remove patches with too few cells
    #--------------------------------------------------------------------------------
    patches_to_keep = []
    for patch in patches:
        adata_patch = adata[adata.obs[patch]==True].copy()
        number_of_cells = adata_patch.n_obs
        if number_of_cells > 20 and adata_patch.obs[celltype_key].nunique() > 1:
            patches_to_keep.append(patch)
            adata_patch.obs.drop(columns=patch_cols, inplace=True)
            adata_patch.write_h5ad(os.path.join(args.output_dir_patches,f'{patch}.h5ad'))
        else:
            with open(os.path.join(args.output_dir_patches,"skipped_patches.txt"), "a") as f:
                f.write(f"Skipped {patch}: num cells={number_of_cells}, num cell types={adata_patch.obs[celltype_key].nunique()}\n")
    print(f'Number of patches after filtering for min number of cells (>20) and min number of cell types (>1): {len(patches_to_keep)}')

    # Mean gene expression as features
    #--------------------------------------------------------------------------------
    print('---------------------- Creating adata with mean gene expression per patch... ----------------------------')
    X_patch_list = []
    obs_patch_list = []

    for patch in patches_to_keep:
        index = adata.obs[patch].values   # Find cells belonging to the same patch
        X_patch_list.append(adata.layers['raw_counts'][index].mean(axis=0).A1)        # For each patch get mean expression of all genes across all cells in the patch
        obs_patch_tmp = adata.obs[index].iloc[0].copy()
        obs_patch_tmp['patch'] = patch
        obs_patch_list.append(obs_patch_tmp)         # For each patch get obs info from the first cell in the patch

    X_patch = np.vstack(X_patch_list)
    adata_patches_gex = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=adata.var.copy())
    adata_patches_gex.obs["patch"] = adata_patches_gex.obs["patch"].astype("category")
    adata_patches_gex.obs.drop(columns=patch_cols, inplace=True)
    print(adata_patches_gex)

    # Cell type fractions as features
    #--------------------------------------------------------------------------------
    print('---------------------- Creating adata with cell type fractions per patch... ----------------------------')
    X_patch_list = []
    obs_patch_list = []

    for patch in patches_to_keep:
        index = adata.obs[patch].values   # Find cells belonging to the same patch
        ct_counts = adata.obs[index][celltype_key].value_counts(normalize=True)  # For each patch get cell type fractions
        ct_fractions = ct_counts.reindex(adata.obs[celltype_key].cat.categories, fill_value=0).values  # Ensure all cell types are represented in the same order
        X_patch_list.append(ct_fractions)
        obs_patch_tmp = adata.obs[index].iloc[0].copy()
        obs_patch_tmp['patch'] = patch
        obs_patch_list.append(obs_patch_tmp)         # For each patch get obs info from the first cell in the patch

    X_patch = np.vstack(X_patch_list)
    adata_patches_ct = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=pd.DataFrame(index=adata.obs[celltype_key].cat.categories))
    adata_patches_ct.obs["patch"] = adata_patches_ct.obs["patch"].astype("category")
    adata_patches_ct.obs.drop(columns=patch_cols, inplace=True)
    print(adata_patches_ct)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Plot patches adatas in umap space
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print('Plotting patches adatas in umap space...')
for adata, suffix in [(adata_patches_gex, 'gex'), (adata_patches_ct, 'ctFraction')]:
    # Normalize and log-transform
    print('Normalizing and log-transforming...')
    adata.layers['raw_counts'] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Dimension reduction
    print('Dimension reduction...')
    sc.pp.pca(adata)
    sc.pp.neighbors(adata, n_neighbors=16)
    sc.tl.umap(adata)
    sc.pl.umap(adata, color='pt_id', show=False, size=5)
    plt.legend().remove()
    plt.savefig(os.path.join(args.output_dir_plots,f'umap_patches_{suffix}.png'), dpi=300, bbox_inches='tight')
    plt.close()

    adata.write_h5ad(os.path.join(args.output_dir_patches,f'adata_patches_{suffix}.h5ad'))


print('Done.')



