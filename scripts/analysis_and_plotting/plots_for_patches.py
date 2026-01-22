import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq


patch_size=5000
overlap=50
phenotyping_level = "Neutro_Epi_extImm_pooled_A_EM_N"
adata_type = "ctFraction"
output_dir_plots = f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/'
patch = 'T24_029745_110009_3_window_0'
# adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_patch/{phenotyping_level}/{patch_size}um_{overlap}um/adata_patches_{adata_type}.h5ad')
# print(adata)

# if adata_type == "ctFraction":
#     colors = ['sample_type', 'treatment_scheme', 'Tumor cells', 'B cell', 'NK cell']
# else:
#     colors = ['sample_type', 'treatment_scheme', 'CD3D','MS4A1']

# for color in colors:
#     sc.pl.umap(adata, color=[color], show=False)
#     plt.savefig(output_dir_plots + f'umap_patches_{adata_type}_{color}.png', dpi=300, bbox_inches='tight')
#     plt.close()


adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_patch/{phenotyping_level}/{patch_size}um_{overlap}um/patch_{patch}.h5ad')
print(adata)

sq.pl.spatial_scatter(adata, color="Neutro_Epi_extImm_pooled_A_EM_N", size=1,shape=None,figsize=(10,10))
plt.savefig(output_dir_plots + f'spatial_scatter_patch_{patch}.png', dpi=300, bbox_inches='tight')
