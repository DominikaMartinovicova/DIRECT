#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# combine_adatas.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine raw counts adata from all the slides to one dataset and preprocess together
# to further check phenotyping and perform downstream analysis.
#
#   0 Import libraries and parse arguments
#   1 Read data and check if raw layers contain raw counts (integers)
#   2 Preprocessing
#       a. Normalization
#       b. PCA
#       c. kNN
#       d. PAGA
#       e. UMAP
#       f. Leiden clustering
#   3 Save
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/preprocessing/combine_adatas.py \
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

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 combine_adatas.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  # Combine raw counts adata from all the slides to one dataset and preprocess together to further check phenotyping and perform downstream analysis.  ')
    parser.add_argument('-i', help='path to phenotyped Xenium dirs metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('--input_dir', help='path to phenotyped Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('--phen_level', help='phenotyping level',
                        dest='phen_level')
    parser.add_argument('-o', help='path to output combined xenium dirs metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_plot', help='path to output combined xenium plots dir',
                        dest='output_plot',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
os.makedirs(args.output_plot, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From input directory read all zarr files and combine them into one adata
#--------------------------------------------------------------------------------
input_dir = args.input_dir
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

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Preprocessing
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Normalize gene expression and save raw
#-------------------------------------------------------------------------------
print("Normalization...")
sc.pp.normalize_total(adata_combined, target_sum=1e4)
sc.pp.log1p(adata_combined)

# Calculate dimensionality reduction
#-------------------------------------------------------------------------------
print("PCA...")
sc.tl.pca(adata_combined)
# Plot variance ratios of PCs
sc.pl.pca_variance_ratio(adata_combined, log=False, n_pcs = 30, show= False)
plt.savefig(args.output_plot + "/PCA_variance_explained.png")
plt.close()

print("Neighbors...")
sc.pp.neighbors(adata_combined)

print("PAGA...")
sc.tl.paga(adata_combined, groups = args.phen_level)
sc.pl.paga(adata_combined)
plt.savefig(args.output_plot + "/PAGA_combined_adatas.svg", format='svg')

print("Calculating UMAP...")
sc.tl.umap(adata_combined)

print("Leiden clustering...")
sc.tl.leiden(adata_combined)

print(adata_combined.X.toarray()[0:5,0:5])
print(adata_combined.layers['raw_counts'].toarray()[0:5,0:5])

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Save
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
adata_combined.write_h5ad(args.output)




