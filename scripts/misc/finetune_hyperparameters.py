#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# finetune_hyperparameters.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Find a combination of min_dist and k for kNN to see clusters in the umap
#
#
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/python/combine_adatas.py \
        -i {input.complete_Xenium} \
	    --input_dir {params.in_dir} \
        -o {output.combined_adatas} \
        --output_plot {params.out_plot_dir}
"""




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
import argparse
import os
import geosketch as gs

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Run_tacco.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Create celltype specific signature matrices  ')
    parser.add_argument('-i', help='path to preprocessed Xenium dirs metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('--input_dir', help='path to preprocessed Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('-threads', help='n threads to use',
                        dest='threads',
                        type=int)
    parser.add_argument('-o', help='path to output phenotyped xenium dirs metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_plot', help='path to output phenotyped xenium dirs metadata file',
                        dest='output_plot',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# From input directory read all zarr files and combine them into one adata
#--------------------------------------------------------------------------------
input_dir = '/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/complete/'
for folder in os.listdir(input_dir):    # Loop over all folders and files in the directory
    if folder.endswith(".zarr"):    # Only select .zarr folders 
        print("Reading " + folder)
        sdata_tmp = sd.read_zarr(input_dir + folder)    # Read the spatial data
        adata_tmp = sdata_tmp.tables["table"]   #Subset for only adata 
        adata_tmp.obs_names = adata_tmp.obs_names + '_' + folder.replace('.zarr','')    # Add identifier to each cell name to prevent duplicates
        adata_tmp.obs['slide'] = folder.replace('.zarr','')     # Add column in .obs to retain information about the data origin
        if 'adata_combined' not in locals():
            adata_combined = adata_tmp
        else:
            adata_combined = ad.concat([adata_combined, adata_tmp], join='outer')

print("Combined adata shape: " + str(adata_combined.shape))
print(adata_combined.X.toarray()[0:5,0:5])
print(adata_combined.layers['raw_counts'].toarray()[0:5,0:5])
print("Combined adata shape of raw_counts: " + str(adata_combined.layers['raw_counts'].toarray().shape))
adata_combined.X = adata_combined.layers['raw_counts'].copy()  # Set the main matrix to raw counts
print("Combined adata shape: " + str(adata_combined.shape))
print(adata_combined.X.toarray()[0:5,0:5])
adata_combined = adata_combined[~adata_combined.obs['sample'].isna(), :]
print("Combined adata shape after removing cells not in any sample: " + str(adata_combined.shape))
print(adata_combined)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Preprocessing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Normalize gene expression and save raw
#-------------------------------------------------------------------------------
print("Normalization...")
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

print("Calculating PCA...")
adata_combined.layers["lognorm"] = adata_combined.X.copy()
sc.pp.pca(adata_combined, n_comps=200)

print('Subsetting cells...')
N = 100000 # Number of cells to obtain from the data set
sketch_index = gs.gs(adata_combined.obsm['X_pca'], N, replace=False)

# sketch the anndata object using sketch_index
adata_sketch = adata_combined[sketch_index]
print("Sketched adata shape: " + str(adata_sketch.shape))
print(adata_sketch)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Hyperparameter tuning
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Test multiple hyperparameter combinations for UMAP calculation
#--------------------------------------------------------------------------------
ks = [2,5,10,15,20]
min_dists = [0,0.05,0.1,0.5,1]
colors = ['celltype', 'sample_type', 'slide', 'T_number', 'pt_id']


for k in ks:
    for min_dist in min_dists:
        os.makedirs('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/hyperparameter_tuning_200pc/UMAP_k' + str(k) + '_minDist' + str(min_dist) + '/', exist_ok=True)
        print("Calculating UMAP with k = " + str(k) + " and min_dist = " + str(min_dist))
        print("Neighbors...")
        sc.pp.neighbors(adata_sketch, n_neighbors = k)
        print("Calculating UMAP...")
        sc.tl.umap(adata_sketch, min_dist = min_dist)

        # Plot UMAP
        print("Plotting...")
        for color in colors:
            print("Plotting color: " + color)
            sc.pl.umap(adata_sketch, color = color, show = False)
            plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/hyperparameter_tuning_200pc/UMAP_k' + str(k) + '_minDist' + str(min_dist) + '/'
                         + "/UMAP_k" + str(k) + "_minDist" + str(min_dist) + "_" + color + ".png", dpi=300, bbox_inches='tight')
            plt.close()

        # Save
        adata_sketch.write_h5ad("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/sketch/combined_adatas_200pc_k" + str(k) + "_minDist" + str(min_dist) + ".h5ad")
        print("Saved combined adata with k = " + str(k) + " and min_dist = " + str(min_dist))





