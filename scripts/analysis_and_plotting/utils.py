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


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Define functions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Remove v1.7 samples and add MPR column
#-------------------------------------------------------------------------------
def preprocess_adata(adata_path, exclude_v17=True):
    adata = sc.read_h5ad(adata_path)
    if exclude_v17 == True:
        print('Excluding v1.7 samples...')
        adata = adata[adata.obs['treatment_scheme'] != 'v1.7', :].copy()
    
    # Add MPR column
    print('Adding MPR column...')
    adata.obs['MPR'] = "<90" if adata.obs['regression'].iloc[0] < 90 else ">=90"
    print(adata)
    print('Preprocessing done.')
    return adata


# Create a dictionary with samples as keys and adata subset for that sample as value
#-------------------------------------------------------------------------------
def get_sample_dict(adata):
    sample_dict = {}
    for sample in adata.obs['sample'].unique():
        print(f'Processing sample: {sample}')
        sample_dict[sample] = adata[adata.obs['sample'] == sample, :].copy()
    print(sample_dict.keys())
    print('Sample dictionary created.')
    return sample_dict