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
warnings.simplefilter("ignore", FutureWarning)
warnings.simplefilter("ignore", PerformanceWarning)

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
    parser.add_argument('-o', '--output_patches', help='path to output dir with patches per sample',
                        dest='output_dir_patches',
                        type=str)
    parser.add_argument('--plots_dir', help='path to output dir for plots',
                        dest='output_dir_plots',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
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
sq.tl.sliding_window(adata=adata,library_key="sample", window_size=args.patch_size, overlap=args.overlap,copy=False)

sq.pl.spatial_scatter(adata, color="sliding_window_assignment", library_key="sample")
plt.savefig(os.path.join(args.output_dir_plots,'spatial_scatter_sliding_window_assignment.png'))
plt.close()

sample="T23_004535_110005_1"
sq.pl.spatial_scatter(adata,color="sliding_window_assignment",library_key="sample",library_id=[sample],figsize=(10, 10))
plt.savefig(os.path.join(args.output_dir_plots,f'spatial_scatter_sliding_window_assignment_{sample}.png'))
plt.close()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Create patch adata (patches as observations)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
patches = adata.obs['sliding_window_assignment'].unique()
print(f'Total number of patches: {len(patches)}')

# Remove patches with too few cells
#--------------------------------------------------------------------------------
patches_to_keep = []
for patch in patches:
    adata_patch = adata[adata.obs['sliding_window_assignment']==patch]
    number_of_cells = adata_patch.n_obs
    if number_of_cells > 20:
        patches_to_keep.append(patch)
print(f'Number of patches after filtering for min number of cells (>20): {len(patches_to_keep)}')

# Mean gene expression as features
#--------------------------------------------------------------------------------
print('---------------------- Creating adata with mean gene expression per patch... ----------------------------')
X_patch_list = []
obs_patch_list = []

for patch in patches_to_keep:
    index = adata.obs['sliding_window_assignment']==patch   # Find cells belonging to the same patch
    X_patch_list.append(adata.raw[index].mean(axis=0))        # For each patch get mean expression of all genes across all cells in the patch
    obs_patch_list.append(adata.obs[index].iloc[0])         # For each patch get obs info from the first cell in the patch

X_patch = np.vstack(X_patch_list)
adata_patches_gex = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=adata.var.copy())

# Normalize and log-transform
print('Normalizing and log-transforming...')
adata_patches_gex.layers['raw_counts'] = adata_patches_gex.X.copy()
sc.pp.normalize_total(adata_patches_gex, target_sum=1e4)
sc.pp.log1p(adata_patches_gex)

# Dimension reduction
print('Dimension reduction...')
sc.pp.pca(adata_patches_gex)
sc.pp.neighbors(adata_patches_gex, n_neighbors=16)
sc.tl.umap(adata_patches_gex)
sc.pl.umap(adata_patches_gex, color='sliding_window_assignment', show=False, size=5)
plt.savefig(os.path.join(args.output_dir_plots,'umap_patches_gex.png'), dpi=300, bbox_inches='tight')
plt.close()

# Cell type fractions as features
#--------------------------------------------------------------------------------
print('---------------------- Creating adata with cell type fractions per patch... ----------------------------')
X_patch_list = []
obs_patch_list = []

for patch in patches_to_keep:
    index = adata.obs['sliding_window_assignment']==patch   # Find cells belonging to the same patch
    ct_counts = adata.obs[index]['Neutro_Epi_extImm_pooled_A_EM_N'].value_counts(normalize=True)  # For each patch get cell type fractions
    ct_fractions = ct_counts.reindex(adata.obs['Neutro_Epi_extImm_pooled_A_EM_N'].cat.categories, fill_value=0).values  # Ensure all cell types are represented in the same order
    X_patch_list.append(ct_fractions)
    obs_patch_list.append(adata.obs[index].iloc[0])         # For each patch get obs info from the first cell in the patch

X_patch = np.vstack(X_patch_list)
adata_patches_ct = sc.AnnData(X=X_patch,obs=pd.DataFrame(obs_patch_list).reset_index(drop=True), var=pd.DataFrame(index=adata.obs['Neutro_Epi_extImm_pooled_A_EM_N'].cat.categories))

# Normalize and log-transform
print('Normalizing and log-transforming...')
adata_patches_ct.layers['raw_counts'] = adata_patches_ct.X.copy()
sc.pp.normalize_total(adata_patches_ct, target_sum=1e4)
sc.pp.log1p(adata_patches_ct)

# Dimension reduction
print('Dimension reduction...')
sc.pp.pca(adata_patches_ct)
sc.pp.neighbors(adata_patches_ct, n_neighbors=16)
sc.tl.umap(adata_patches_ct)
sc.pl.umap(adata_patches_ct, color='sliding_window_assignment', show=False, size=5)
plt.savefig(os.path.join(args.output_dir_plots,'umap_patches_ctFraction.png'), dpi=300, bbox_inches='tight')
plt.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Save adata of each patch
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
for patch in patches_to_keep:
    # Gene expression adata
    adata_patch = adata[adata.obs['sliding_window_assignment']==patch]
    adata_patch.write_h5ad(os.path.join(args.output_dir_patches,f'adata_patch_{patch}.h5ad'))

print('Done.')



