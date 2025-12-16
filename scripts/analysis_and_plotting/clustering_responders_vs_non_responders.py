import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt


# Load data
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/Neutro_Epi_extImm_combined_adatas.h5ad')
celltype_key = 'Neutro_Epi_extImm'
output_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/combined/Neutro_Epi_extImm/resp_vs_nonresp'

# Plot umap with responders vs non-responders
adata.obs['MPR'] = adata.obs['regression'].apply(lambda x: '>=90' if x >= 90 else '<90')
sc.pl.umap(adata, color='MPR', show=False)
plt.savefig(f'{output_dir}/umap_responders_vs_nonresponders.png', dpi=300)
plt.close() 

# Subset data for Biopsy and Resection samples
biopsy_adata = adata[adata.obs['sample_type'] == 'Biopsy']
sc.pl.umap(biopsy_adata, color='MPR', show=False)
plt.savefig(f'{output_dir}/umap_responders_vs_nonresponders_biopsy.png', dpi=300)
plt.close()

resection_adata = adata[adata.obs['sample_type'] == 'Resection']
sc.pl.umap(resection_adata, color='MPR', show=False)
plt.savefig(f'{output_dir}/umap_responders_vs_nonresponders_resection.png', dpi=300)
plt.close()