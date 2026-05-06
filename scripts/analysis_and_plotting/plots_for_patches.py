import scanpy as sc
import matplotlib.pyplot as plt
import squidpy as sq
import spatialdata as sd
import pandas as pd


patch_size=5000
overlap=50
phenotyping_level = "Neutro_Epi_extImm_pooled_A_EM_N"
adata_type = "ctFraction"
output_dir_plots = f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/spatial_scatter/'
patch = 'T23_004535_110005_1'
patches = [
    "T24_040744_130003_3",
    "T24_041869_130005_3",
    "T24_041865_130004_2",
    "T24_041865_130004_1",
    "T24_016798_110007_1",
    "T25_009511_110013_1",
    "T24_041865_130004_3",
    "T24_041869_130005_1",
    "T24_041869_130005_2",
    "T24_040746_130002_1",
    "T24_040744_130003_2",
    "T24_040746_130002_2",
    "T24_040744_130003_1",
    "T24_016798_110007_2",
    "T25_009511_110013_2",
    "T24_051361_110012_1",
    "T24_050721_110011_1",
    "T24_051361_110012_2",
    "T24_050721_110011_2",
    "T24_016798_110007_3",
    "T24_051361_110012_3",
    "T25_009511_110013_3",
    "T24_040746_130002_3",
    "T25_009512_110014_3",
    "T24_050721_110011_3",
    "T25_009512_110014_1",
    "T25_009512_110014_2",
    "TVU22_821_110002_2",
    "TVU21_11612_110001_1",
    "TVU22_821_110002_1",
    "TVU21_11612_110001_2",
    "TVU22_821_110002_3",
    "TVU21_11612_110001_3",
    "TVU22_1407_110003_3",
    "T23_004535_110005_3",
    "T23_004719_110006_3",
    "TVU22_1407_110003_2",
    "T24_040097_110010_1",
    "T24_040097_110010_3",
    "T24_040097_110010_2",
    "T23_004535_110005_2",
    "T24_040748_130001_3",
    "T24_022675_110008_3",
    "T23_004719_110006_2",
    "T24_022675_110008_1",
    "T24_029745_110009_1",
    "T24_040748_130001_1",
    "T24_029745_110009_3",
    "T24_029745_110009_2",
    "T24_022675_110008_2",
    "T24_040748_130001_2",
    "TVU23_3277_110004_1",
    "TVU23_3277_110004_2",
    "TVU23_3277_110004_3",
    "T23_004719_110006_1",
    "T23_004535_110005_1",
    "TVU22_1407_110003_1",
    "T24_041866_130005_1",
    "T24_040742_130003_2",
    "T24_040742_130003_3",
    "T24_040742_130003_1",
    "T24_041864_130004_1",
    "T24_040745_130002_1",
    "T25_000668_110014_2",
    "T24_044068_110012_4",
    "T24_041768_110011_2",
    "T24_008367_110007_1",
    "T25_000668_110014_1",
    "T24_044068_110012_2",
    "T24_041768_110011_1",
    "T24_044068_110012_1",
    "T24_044068_110012_3",
    "T25_002336_110005_1",
    "T24_040747_130001_1",
    "TVU23_03906_110006_2",
    "TVU23_03906_110006_1",
    "TVU21_14439_110002_2",
    "TVU21_14439_110002_1",
    "T24_012138_110008_1",
    "T24_041920_110003_1",
    "TVU21_9662_110001_4",
    "T24_012138_110008_3",
    "TVU21_9662_110001_2",
    "T24_012138_110008_2",
    "T24_032894_110010_1",
    "TVU21_9662_110001_3",
    "T24_041921_110004_1",
    "T24_022363_110009_3",
    "T24_022363_110009_1",
    "T24_041921_110004_2",
    "TVU21_9662_110001_1",
    "T24_022363_110009_2",
    "T24_032894_110010_2",
]
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


#adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_patch/{phenotyping_level}/{patch_size}um_{overlap}um/patch_{patch}.h5ad')
#print(adata)

celltype_palette = {
    "B_cell": "#1f77b4",
    "DC_mature": "#d62728",
    "Endothelial_cell": "#ff9896",
    "Epithelial_cell": "#c7c7c7",
    "Macrophage": "#ffbb78",
    "Macrophage_alveolar": "#ff7f0e",
    "Mast_cell": "#8b8b12",
    "Monocyte_classical": "#bcbd22",
    "Monocyte_non-classical": "#dbdb8d",
    "NAN": "#9edae5",
    "NK_cell": "#8c564b",
    "Plasma_cell": "#aec7e8",
    "Stromal": "#c49c94",
    "TAN": "#17becf",
    "T_cell_CD4": "#f7b6d2",
    "T_cell_CD8_functional": "#2ca02c",
    "T_cell_CD8_terminally_exhausted": "#98df8a",
    "T_cell_NK-like": "#c5b0d5",
    "T_cell_regulatory": "#e377c2",
    "Tumor_cells": "#7f7f7f",
    "cDC1": "#800020",
    "cDC2": "#CF4265",
    "pDC": "#9467bd",
}


def apply_palette_to_obs(adata, obs_key, palette, fallback="#d3d3d3"):
    # Ensure colors follow the exact category order expected by scanpy/squidpy.
    categories = adata.obs[obs_key].astype("category").cat.categories.tolist()
    adata.uns[f"{obs_key}_colors"] = [palette.get(ct, fallback) for ct in categories]

def apply_highlight_palette(adata,obs_key,palette,highlight_celltypes,fallback="#d3d3d3"):
    categories = (adata.obs[obs_key].astype("category").cat.categories.tolist())
    adata.uns[f"{obs_key}_colors"] = [palette.get(ct, fallback) if ct in highlight_celltypes else fallback for ct in categories]


#for patch in patches:
#    adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_sample/{phenotyping_level}/{patch}.h5ad')
#    apply_palette_to_obs(adata, phenotyping_level, celltype_palette)
#    sq.pl.spatial_scatter(adata, color="Neutro_Epi_extImm_pooled_A_EM_N", size=1,shape=None,figsize=(10,10))
#    plt.savefig(output_dir_plots + f'spatial_scatter_patch_{patch}.png', dpi=300, bbox_inches='tight')
#    plt.close()


# # create a plot of the whole slide
# sdata = sd.read_zarr("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/phenotyped/Neutro_Epi_extImm_pooled_A_EM_N/slide_5.zarr/")
# adata = sdata.tables['table']
# # rewrite each space in cell label to _
# adata.obs[phenotyping_level] = adata.obs[phenotyping_level].str.replace(" ", "_")
# apply_palette_to_obs(adata, phenotyping_level, celltype_palette)
# sq.pl.spatial_scatter(adata,library_id="spatial",shape=None,color=phenotyping_level,size = 2,figsize=(10,10))
# #change title to "Spatial scatter of phenotyped cells in slide 5"
# plt.title("Spatial scatter of phenotyped cells in slide 5", fontsize=18)
# # make the legend bigger and two column
# legend = plt.legend(title=phenotyping_level, title_fontsize=18, fontsize=16, loc='upper right', bbox_to_anchor=(2, 1))
# legend.get_title().set_fontsize(18)
# plt.tight_layout()
# plt.savefig(f'plots/tacco/{phenotyping_level}/slide_5/' + f'/spatial_{phenotyping_level}_new.png', dpi=300, bbox_inches='tight')
# plt.close()

# plot a scatterplot for a patch and highlight only selected celltypes in their respective colors and all the other cells in light gray
celltype_1 = 'B_cell'
celltype_2 = 'T_cell_regulatory'
highlight_celltypes = [celltype_1, celltype_2]
adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_sample/{phenotyping_level}/{patch}.h5ad')
apply_highlight_palette(adata=adata,obs_key=phenotyping_level,palette=celltype_palette,highlight_celltypes=highlight_celltypes,fallback="#d3d3d3")
# Put non-highlighted cells first, highlighted cells last
is_highlight = adata.obs[phenotyping_level].isin(highlight_celltypes)

adata_plot = adata[pd.concat([(~is_highlight).loc[~is_highlight],   # background first 
                              is_highlight.loc[is_highlight]        # highlighted last
    ]).index].copy()

sq.pl.spatial_scatter(adata_plot, color="Neutro_Epi_extImm_pooled_A_EM_N", size=1,shape=None,figsize=(10,10))
plt.title(f"{patch}", fontsize=18)
plt.savefig(output_dir_plots + f'spatial_scatter_patch_{patch}_{celltype_1}_{celltype_2}.png', dpi=300, bbox_inches='tight')
plt.close()



