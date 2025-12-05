#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run_tacco.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Run tacco to transfer the cell labels from reference scRNA atlas. Preprocess
#   datasets individually for potential inspection.
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Run tacco
#   3 Dimensionality reduction
#       a. Normalization
#       b. PCA
#       c. kNN
#       d. PAGA
#       e. UMAP
#   4 Save
#
# Author: Mischa Steketee (m.f.b.steketee@amsterdamumc.nl)
#
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/preprocessing/Run_tacco.py \
        -i {input.preprocessed_Xenium} \
        --input_atlas {input.scRNAseq_atlas} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        --phen_level {params.phen_level} \
        -o {output.phenotyped_Xenium} \
        --output_plot {output.phenotyped_umap}
"""
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import tacco as tc
from scipy.sparse import csr_matrix
import argparse
import os

# Return indices of elements in l1 that are not in l2
def non_matched_index(l1, l2):
    l2 = set(l2)
    return [i for i, el in enumerate(l1) if el not in l2]


# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 Run_tacco.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Run tacco to transfer the cell labels from reference scRNA atlas. Preprocess datasets individually for potential inspection  ')
    parser.add_argument('-i', help='path to preprocessed Xenium dirs metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('--input_atlas', help='path to scRNAseq atlas',
                        dest='input_atlas',
                        type=str)
    parser.add_argument('--input_dir', help='path to preprocessed Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('--output_dir', help='path to phenotyped output dir',
                        dest='output_dir',
                        type=str)
    parser.add_argument('-threads', help='n threads to use',
                        dest='threads',
                        type=int)
    parser.add_argument('-o', help='path to output phenotyped xenium dirs metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_dir_plot', help='path to plot of UMAP',
                        dest='output_dir_plot',
                        type=str)
    parser.add_argument('--phen_level', help='phenotyping level',
                        dest='phen_level',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
os.makedirs(args.output_dir_plot, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print("Reading adata and sdata...")
sdata = sd.read_zarr(args.input_dir)
adata_atlas = sc.read_h5ad(args.input_atlas)
adata = sdata.tables["table"]

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Run tacco
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
adata_atlas = adata_atlas.raw.to_adata()
adata_atlas.var.index = adata_atlas.var['feature_name'] # Set index to gene names to match the atlas

print("Annotating data...")
tc.tl.annotate(adata, adata_atlas, annotation_key=args.phen_level, result_key='celltype_major', multi_center = 10)

adata.obs[args.phen_level] = adata.obsm['celltype_major'].idxmax(axis=1)
adata.obs[args.phen_level] = adata.obs[args.phen_level].astype('category')

# Remove cells with no assigned cell type by tacco
NAN_ct_id = adata.obs.loc[adata.obs[args.phen_level].isna(),:].cell_id.to_list()
adata = adata[non_matched_index(adata.obs.cell_id.to_list(), list(NAN_ct_id)),:]
print(f"Nan celltypes removed from {args.input_dir}: ", len(NAN_ct_id))


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Dimensionality reduction
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Normalize gene expression and save raw
#-------------------------------------------------------------------------------
print("Normalization...")
adata.layers["raw_counts"] = adata.X.copy()
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)

# Calculate dimension reductions
#-------------------------------------------------------------------------------
print("PCA...")
sc.pp.pca(adata)
sc.pp.neighbors(adata, n_neighbors=16)

print("PAGA...")
sc.tl.paga(adata, groups = args.phen_level)
sc.pl.paga(adata)
           
print("Calculating UMAP...")
sc.tl.umap(adata, init_pos = 'paga')
sc.pl.umap(adata, color=args.phen_level, show = False, size=5)
plt.savefig(args.output_dir_plot + f'/umap_{args.phen_level}.png', dpi=300, bbox_inches='tight')


# Plot scatter with cell labels
#-------------------------------------------------------------------------------
sq.pl.spatial_scatter(adata,library_id="spatial",shape=None,color=args.phen_level,size = 1)
plt.savefig(args.output_dir_plot + f'/spatial_{args.phen_level}.png', dpi=300, bbox_inches='tight')
plt.close()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Save 
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sdata.tables['table'] = adata   # Rewrite table to adata with phenotyped cells
os.rmdir(args.output_dir)
del adata.uns['celltype_major_mc10'] # del this uns as it somehow is not writeable
sdata.write(args.output_dir)
