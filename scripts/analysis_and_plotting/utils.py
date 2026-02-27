#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# utils.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Utils code for DIRECT Xenium spatial analysis snakemake
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import scanpy as sc
import numpy as np
import os
import pandas as pd

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove v1.7 samples and add MPR column
#-------------------------------------------------------------------------------
def preprocess_adata(adata_path, exclude_v17, phenotyping_level=None):
    print('Reading adata...')
    adata = sc.read_h5ad(adata_path)
    ex_v17 = 'w_v1.7'
    if exclude_v17 == True:
        print('Excluding v1.7 samples...')
        ex_v17 = 'wo_v1.7'
        adata = adata[adata.obs['treatment_scheme'] != 'v1.7', :].copy()
    
    # Add MPR column
    print('Adding MPR column...')
    reg = pd.to_numeric(adata.obs['regression'], errors="coerce")
    adata.obs['MPR'] = np.where(reg < 90.0,"<90",">=90")
    print(adata)
    #adata.write_h5ad(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{phenotyping_level}_combined_adatas_for_analysis_{ex_v17}.h5ad")

    # Add treatment column
    print('Adding treatment column...')
    adata.obs['treatment'] = np.where(adata.obs['treatment_scheme'] == 'v1.7',"aggressive","milder")
    print(adata)

    print('Adding structure_core column...')
    resections_loc = adata.obs['sample_type'] == 'Resection'
    adata.obs['structure_core'] = np.where(resections_loc, 'core_' + adata.obs['sample'].str.split('_').str[-1], np.nan)
    print(adata)
 
    adata.write_h5ad(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{phenotyping_level}_combined_adatas_for_analysis_{ex_v17}.h5ad")
    
    # Save cell type list
    print('Saving cell type list...')
    celltype_list = sorted(adata.obs[phenotyping_level].unique().tolist())
    pd.Series(celltype_list).to_csv(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{phenotyping_level}_celltype_list_{ex_v17}.txt", index=False, header=False)
    print(celltype_list)
    print('Preprocessing done.')
    samples_list = adata.obs['sample'].unique().tolist()
    return samples_list, adata, ex_v17

# Save per-sample adatas
#-------------------------------------------------------------------------------
def save_per_sample_adata(adata, output_dir):
    print('Saving per-sample adatas...')
    os.makedirs(output_dir, exist_ok=True)
    for i, sample in enumerate(adata.obs['sample'].unique()):
        print(f'Saving sample {i+1}/{len(adata.obs["sample"].unique())}: {sample}')
        adata_sample = adata[adata.obs['sample'] == sample, :].copy()
        adata_sample.write_h5ad(os.path.join(output_dir, f'{sample}_adata.h5ad'))
    print('Per-sample adatas saved.')

# Create a dictionary with samples as keys and adata subset for that sample as value
#-------------------------------------------------------------------------------
def get_sample_dict(adata):
    print('Creating sample dictionary...')
    sample_dict = {}
    for sample in adata.obs['sample'].unique():
        print(f'Processing sample: {sample}')
        sample_dict[sample] = adata[adata.obs['sample'] == sample, :].copy()
    print(sample_dict.keys())
    print('Sample dictionary created.')
    return sample_dict

