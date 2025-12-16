import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


# Load data
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
output_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/combined/Neutro_Epi_extImm/resp_vs_nonresp'

exclude_v17 = True
if exclude_v17:
    print('Removing samples with treatment_scheme=v1.7...')
    adata = adata[adata.obs['treatment_scheme'] != 'v1.7', :]
    output_dir = output_dir + '/wo_v1.7'
    print(f'Number of cells after removing v1.7 treatment scheme samples: {adata.shape[0]}')



# Plot umap with responders vs non-responders
adata.obs['MPR'] = adata.obs['regression'].apply(lambda x: '>=90' if x >= 90 else '<90')
sc.pl.umap(adata, color='MPR', show=False)
plt.title('All samples wo v1.7' if exclude_v17 else 'All samples')
plt.savefig(f'{output_dir}umap_responders_vs_nonresponders.png', dpi=300, bbox_inches='tight')
plt.close() 

# Subset data for Biopsy and Resection samples
biopsy_adata = adata[adata.obs['sample_type'] == 'Biopsy']
sc.pl.umap(biopsy_adata, color='MPR', show=False)
plt.title('Biopsy Samples wo v1.7' if exclude_v17 else 'Biopsy Samples')
plt.savefig(f'{output_dir}umap_responders_vs_nonresponders_biopsy.png', dpi=300, bbox_inches='tight')
plt.close()

resection_adata = adata[adata.obs['sample_type'] == 'Resection']
sc.pl.umap(resection_adata, color='MPR', show=False)
plt.title('Resection Samples wo v1.7' if exclude_v17 else 'Resection Samples')
plt.savefig(f'{output_dir}umap_responders_vs_nonresponders_resection.png', dpi=300, bbox_inches='tight')
plt.close()