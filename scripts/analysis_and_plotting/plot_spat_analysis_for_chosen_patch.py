import matplotlib.pyplot as plt
import pickle
import os
import seaborn as sns
import pandas as pd
import squidpy as sq
import scanpy as sc
import numpy as np

patch = 'T23_004535_110005_1_window_0'
core = 'T23_004535_110005_1'
overlap=50
patch_size=5000
phenotyping_level = 'Neutro_Epi_extImm_pooled_A_EM_N'
input_dir = f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/{phenotyping_level}/spatial/per_patch/{patch_size}um_{overlap}um/{patch}/'

# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# # 2 Prepare data and run analysis
# #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# centrality_scores=pd.read_csv(os.path.join(input_dir, 'centrality_scores.csv'), index_col=0)
# print(centrality_scores)
# # Plot centrality scores
# y='closeness_centrality'
# plt.figure(figsize=(8,6))
# sns.barplot(data=centrality_scores, x=centrality_scores.index, y=y)
# plt.title(f'Centrality Scores for Patch {patch}')
# plt.xticks(rotation=90)
# plt.tight_layout()
# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_{y}.png')
# plt.close()

# nhood_enrichment_path = os.path.join(input_dir, 'neighbors_enrichment.pkl')
# with open(nhood_enrichment_path, 'rb') as f:
#     nhood_enrichment = pickle.load(f)

# plt.figure(figsize=(8,8))
# sns.heatmap(nhood_enrichment['zscore'], cmap='vlag', center=0, xticklabels=centrality_scores.index, yticklabels=centrality_scores.index)
# plt.title(f'Neighborhood enrichment for Patch {patch}')
# plt.tight_layout()
# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_nhood_enrichment.png')
# plt.close()

# interaction_matrix_path = os.path.join(input_dir, 'interaction_matrix.pkl')
# with open(interaction_matrix_path, 'rb') as f:
#     interaction_matrices = pickle.load(f)
# plt.figure(figsize=(8,8))
# sns.heatmap(interaction_matrices, cmap='vlag', center=0, xticklabels=centrality_scores.index, yticklabels=centrality_scores.index)
# plt.title(f'Interaction matrix for Patch {patch}')
# plt.tight_layout()
# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_interaction_matrix.png')
# plt.close()


# coocurrence_path = os.path.join(input_dir, 'co_occurrence_probabilities.pkl')
# with open(coocurrence_path, 'rb') as f:
#     cooccurrence_matrices = pickle.load(f)

# print(cooccurrence_matrices)
# occ = cooccurrence_matrices["occ"]                 # (cell, cell, distance)
# distances = cooccurrence_matrices["interval"]
# cell_types = cooccurrence_matrices["cell_type_key"]

# b_idx = cell_types.index("B cell")

# dist_edges = cooccurrence_matrices["interval"]              # length 50
# dist_centers = 0.5 * (dist_edges[:-1] + dist_edges[1:])  # length 49

# # choose a colormap (tab20 has 20 distinct colors)
# cmap = plt.get_cmap("tab20")
# colors = [cmap(i / len(cell_types)) for i in range(len(cell_types))]

# fig, ax = plt.subplots(figsize=(8, 6))

# for j, ct in enumerate(cell_types):
#     ax.plot(
#         dist_centers,
#         occ[b_idx, j, :],
#         label=ct,
#         linewidth=1,
#         color=colors[j]
#     )

# # reference line (random expectation)
# ax.axhline(1, linestyle="--", linewidth=1)

# ax.set_xlabel("Distance (µm)")
# ax.set_ylabel("Co-occurrence score")
# ax.set_title("B cell co-occurrence vs distance (all cell types)")

# ax.legend(
#     bbox_to_anchor=(1.05, 1),
#     loc="upper left",
#     fontsize=8,
#     frameon=False
# )

# plt.tight_layout()
# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_cooccurrence_mean.png')
# plt.close()



# # Plot ripley's L statistics
# ripleys_L_path = os.path.join(input_dir, 'dict_ripleys_L.pkl')
# with open(ripleys_L_path, 'rb') as f:
#     dict_ripley = pickle.load(f)

# df = dict_ripley['L_stat']

# plt.figure(figsize=(12, 7))
# sns.lineplot(data=df,x='bins',y='stats',hue='Neutro_Epi_extImm_pooled_A_EM_N',estimator=None, palette="tab20")  # ensures raw lines are drawn

# plt.legend(title='Cell type', bbox_to_anchor=(1.02, 1), loc='upper left')
# plt.tight_layout()
# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_ripleys_L.png', bbox_inches='tight')


# # Plot cross Ripley's statistics for cells of interest
# type_i = 'B_cell'
# type_j = 'T_cell_regulatory'
# ripleys_L_path = os.path.join(input_dir, 'dict_ripleys_L.pkl')
# with open(os.path.join(ripleys_L_path, f"cross_ripley_{type_i}_vs_{type_j}.pkl"), "rb") as f:
#     res = pickle.load(f)


# r = res["r"]
# observed = res["L_observed"]
# sims = res["L_simulations"]
# csr = res["csr_expectation"]

# # simulation envelope (5–95%)
# lower = np.percentile(sims, 5, axis=0)
# upper = np.percentile(sims, 95, axis=0)

# plt.figure(figsize=(6,5))

# # simulation envelope
# plt.fill_between(r, lower, upper, alpha=0.3, label="Permutation envelope (5–95%)")

# # observed L
# plt.plot(r, observed, linewidth=2, label="Observed L(r)")

# # CSR expectation
# plt.plot(r, csr, linestyle="--", label="CSR expectation")

# plt.xlabel("Radius (r)")
# plt.ylabel("L(r)")
# plt.title(f"Cross-type Ripley's L: {res['type_i']} vs {res['type_j']}")
# plt.legend()
# plt.tight_layout()

# plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_{type_i}_vs_{type_j}.png', bbox_inches='tight')

# Plot one cell type in spatial
#-----------------------------------------------
adata = sc.read(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_patch/{phenotyping_level}/{patch_size}um_{overlap}um/patch_{patch}.h5ad')
unique_cts = adata.obs[phenotyping_level].unique()
for i, ct in enumerate(unique_cts):
    tmp_key = f"{phenotyping_level}_{ct}"
    
    # convert to object to allow assignment of 'other'
    adata.obs[tmp_key] = adata.obs[phenotyping_level].astype(str)
    
    # assign all non-target cells to 'other'
    adata.obs[tmp_key] = adata.obs[tmp_key].where(adata.obs[tmp_key] == ct, other="zother")
    
    # convert back to categorical with ordered categories ['other', ct]
    adata.obs[tmp_key] = adata.obs[tmp_key].astype("category")
    adata.obs[tmp_key].cat.reorder_categories(["zother", ct]) #, inplace=True)
    
    # assign colors
    adata.uns[f"{tmp_key}_colors"] = ["red", "lightgray"]
    fig, ax = plt.subplots()
    # plot and save
    sq.pl.spatial_scatter(
        adata,
        color=tmp_key,
        title=f"{ct}",
        size=1, shape= None,figsize=(10,10), ax=ax)
    
    # # pick a reference cell
    # cx, cy = adata.obsm["spatial"][0]

    # radii = [25, 50, 75, 100, 250]  # in microns

    # for r in radii:
    #     circle = plt.Circle((cx, cy), r, fill=False, linewidth=2)
    #     ax.add_patch(circle)   
    
    
    plt.savefig(
        f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N_old/spatial/patching/5000um_50um/spatial_scatter_patch_{patch}_{ct}.png',
        bbox_inches='tight'
    )
    plt.close()
    
    # clean up
    del adata.obs[tmp_key]
    del adata.uns[f"{tmp_key}_colors"]






# Check Ripley's calculations
#----------------------------------------------------------------------------
adata = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_adatas_ripleys_{core}.h5ad')
mask = adata.obs["sample"] == core
obs_curve = adata.uns["ripley_obs_curve"]
sim_mean = adata.uns["ripley_sim_mean"]
sim_std = adata.uns["ripley_sim_std"]

# extract metadata
interactions = adata.uns["ripley_interactions"]
celltypes = ["B_cell", "Macrophage", "Macrophage_alveolar", "NK_cell", "Stromal", "T_cell_CD4", "T_cell_CD8_functional", "T_cell_CD8_terminally_exhausted", "T_cell_regulatory", "Tumor_cells"]
print(celltypes)
n_types = len(celltypes)

radii = np.array(adata.uns["ripley_params"]["radii"])

# helper to create matrix of celltype x celltype
def interaction_to_matrix(values, interactions, celltypes):
    mat = np.full((len(celltypes), len(celltypes)), np.nan)

    for idx, name in enumerate(interactions):
        for i in celltypes:
            prefix = i + "_"
            if name.startswith(prefix):
                j = name[len(prefix):]
                if j in celltypes:
                    ii = celltypes.index(i)
                    jj = celltypes.index(j)
                    mat[ii, jj] = values[idx]
                break

    return mat

# Z score at different radii
z_tensor = adata.uns["ripley_z"]

#interactions_list = ['T_cell_regulatory_B_cell','B_cell_Stromal','Macrophage_alveolar_Macrophage','Macrophage_alveolar_Tumor_cells','Tumor_cells_Macrophage_alveolar','NK_cell_Stromal','Stromal_NK_cell']
for interaction in interactions:
    idx_inter = np.where(interactions == interaction)[0][0]

    # average across cells in sample
    z_mean_curve = np.nanmean(z_tensor[mask, idx_inter, :], axis=0)

    plt.figure()
    plt.plot(radii, z_mean_curve, marker='o')
    plt.axhline(0, linestyle='--')
    plt.title(f"Z-score curve: {interaction} ({core})")
    plt.xlabel("Radius")
    plt.ylabel("Z-score")
    plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/interactions/cross_ripley_zscore_{interaction}_{core}.png', bbox_inches='tight')
    plt.close()


    # Ripley's stats curve
    obs = np.nanmean(obs_curve[mask, idx_inter, :], axis=0)
    sim_m = np.nanmean(sim_mean[mask, idx_inter, :], axis=0)
    sim_s = np.nanmean(sim_std[mask, idx_inter, :], axis=0)
    upper = sim_m + 1.96 * sim_s
    lower = sim_m - 1.96 * sim_s

    plt.figure(figsize=(6,5))

    # shift by CSR
    obs_shift = obs #- radii
    sim_shift = sim_m #- radii
    upper_shift = upper #- radii
    lower_shift = lower #- radii

    plt.fill_between(radii, lower_shift, upper_shift, alpha=0.3, label="95% envelope")
    plt.plot(radii, sim_shift, linestyle='--', label="Simulation mean")
    plt.plot(radii, obs_shift, marker='o', label="Observed")

    plt.axhline(0, linestyle=':', color='black')

    plt.xlabel("Radius")
    plt.ylabel("L(r)") # - r)
    plt.title(f"{interaction} ({core})")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/interactions/cross_ripley_curve_{interaction}_{core}.png', bbox_inches='tight')
    plt.close()




# Heatmaps of integrals
signed = np.nanmean(adata.obsm["ripley_signed"][mask], axis=0)
absolute = np.nanmean(adata.obsm["ripley_abs"][mask], axis=0)

signed_mat = interaction_to_matrix(signed, interactions, celltypes)
abs_mat = interaction_to_matrix(absolute, interactions, celltypes)

plt.figure(figsize=(6,5))
sns.heatmap(signed_mat, xticklabels=celltypes, yticklabels=celltypes,cmap="RdBu_r", center=0,vmax=50000)
plt.title(f"Signed Ripley integral ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_signed_{core}.png', bbox_inches='tight')
plt.close()


plt.figure(figsize=(6,5))
sns.heatmap(abs_mat, xticklabels=celltypes, yticklabels=celltypes,cmap="Reds",vmax=50000)
plt.title(f"Absolute Ripley integral ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_abs_{core}.png', bbox_inches='tight')
plt.close()




# Heatmap of mean z scores in all interactions
z_mean = np.nanmean(
    np.nanmean(z_tensor[mask], axis=0),  # cells → mean
    axis=1                               # radii → mean
)

z_mat = interaction_to_matrix(z_mean, interactions, celltypes)
plt.figure(figsize=(6,5))
sns.heatmap(z_mat, xticklabels=celltypes, yticklabels=celltypes,
            cmap="RdBu_r", center=0, vmax=100, vmin=-100)
plt.title(f"Mean Z-score ({core})")
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_zscore_{core}.png', bbox_inches='tight')
plt.close()















