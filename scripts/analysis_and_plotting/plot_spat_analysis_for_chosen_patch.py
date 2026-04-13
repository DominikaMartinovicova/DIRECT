import matplotlib.pyplot as plt
import pickle
import os
import seaborn as sns
import pandas as pd
import squidpy as sq
import scanpy as sc
import numpy as np

patch = 'T23_004535_110005_2_window_0'
#patch = 'T24_041865_130004_3_window_0'
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






# # Check Ripley's calculations
# #----------------------------------------------------------------------------
# adata_s = sc.read_h5ad(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/analyzed/Neutro_Epi_extImm_pooled_A_EM_N_adatas_ripleys_{core}_gpt.h5ad')

# celltypes = ["B_cell", "Macrophage", "Macrophage_alveolar", "NK_cell", "Stromal", "T_cell_CD4", "T_cell_CD8_functional", "T_cell_CD8_terminally_exhausted", "T_cell_regulatory", "Tumor_cells"]
# interactions = adata_s.uns["ripley_interactions"]
# radii = adata_s.uns["ripley_params"]["radii"]

# obs_curve = adata_s.uns["ripley_obs_curve"]
# z_scores = adata_s.obsm["ripley_z"]
# print(z_scores.shape)

# sim_mean = adata_s.uns["ripley_sim_mean"][core]
# sim_std  = adata_s.uns["ripley_sim_std"][core]

# interactions = adata_s.uns["ripley_interactions"]
# radii = adata_s.uns["ripley_params"]["radii"]

# # aggregate
# obs_mean = np.nanmean(obs_curve, axis=0)
# z_mean   = np.nanmean(z_scores, axis=0)
# print(z_mean.shape)


# # helper to create matrix of celltype x celltype
# def interaction_to_matrix(values, interactions, celltypes):
#     mat = np.full((len(celltypes), len(celltypes)), np.nan)
#     interaction_dict = dict(zip(interactions, values))

#     for i, ci in enumerate(celltypes):
#         for j, cj in enumerate(celltypes):
#             key = f"{ci}_{cj}"
#             if key in interaction_dict:
#                 mat[i, j] = interaction_dict[key]

#     return mat

# # observed curve and simulated curve
# for i, name in enumerate(interactions):
#     plt.figure(figsize=(5,4))
    
#     plt.plot(radii, obs_mean[i], label="Observed", color='orange', marker='o')
#     plt.plot(radii, sim_mean[i], label="Simulated", color='gray', linestyle='--')
    
#     # optional confidence band
#     plt.fill_between(radii,sim_mean[i] - 1.96*sim_std[i],sim_mean[i] + 1.96*sim_std[i],alpha=0.3)
    
#     plt.xlabel("Radius")
#     plt.ylabel("Ripley's L")
#     plt.title(name)
#     plt.legend()
    
#     plt.tight_layout()
#     plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/interactions/{core}_{name}_curve.png", bbox_inches='tight')
#     plt.close()

# # integrals
# signed = adata_s.obsm["ripley_signed"]
# absolute = adata_s.obsm["ripley_abs"]

# signed_mean = np.nanmean(signed, axis=0)
# absolute_mean = np.nanmean(absolute, axis=0)

# signed_mat = interaction_to_matrix(signed_mean, interactions, celltypes)
# abs_mat    = interaction_to_matrix(absolute_mean, interactions, celltypes)

# plt.figure(figsize=(8,6))

# sns.heatmap(
#     signed_mat,
#     xticklabels=celltypes,
#     yticklabels=celltypes,
#     cmap="RdBu_r",
#     vmax = 40000,
#     center=0,
#     square=True,
#     cbar_kws={"label": "Signed integral"}
# )

# plt.title(f"{core} - Signed interaction")
# plt.xticks(rotation=45, ha="right")
# plt.yticks(rotation=0)
# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/{core}_signed_heatmap.png", bbox_inches='tight')
# plt.close()

# plt.figure(figsize=(8,6))

# sns.heatmap(
#     abs_mat,
#     xticklabels=celltypes,
#     yticklabels=celltypes,
#     cmap="Reds",
#     square=True,
#     cbar_kws={"label": "Absolute integral"}
# )

# plt.title(f"{core} - Absolute interaction")
# plt.xticks(rotation=45, ha="right")
# plt.yticks(rotation=0)
# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/{core}_absolute_heatmap.png", bbox_inches='tight')
# plt.close()

# # z-score
# z_mean_inter = np.nanmean(z_mean, axis=1) #!!!!!!! cannot be mean because of irregular steps in radii, max radius is 250
# z_mat = interaction_to_matrix(z_mean_inter, interactions, celltypes)


# plt.figure(figsize=(10,8))
# sns.heatmap(
#     z_mat,
#     xticklabels=celltypes,
#     yticklabels=celltypes,
#     cmap="RdBu_r",
#     center=0,
#     vmin=-200, vmax=200,   # IMPORTANT: consistent scaling
#     square=True,
#     cbar_kws={"label": "Z-score"}
# )

# plt.title(f"{core} - Z-score heatmap")
# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/{core}_z_heatmap.png", bbox_inches='tight')
# plt.close()

# n_r = z_mean.shape[1]

# fig, axes = plt.subplots(1, n_r, figsize=(6*n_r, 6))

# if n_r == 1:
#     axes = [axes]

# for k in range(n_r):
#     z_k = z_mean[:, k]  # z-scores for this radius

#     # convert to matrix (same function you already use)
#     z_mat = interaction_to_matrix(z_k, interactions, celltypes)

#     sns.heatmap(
#         z_mat,
#         xticklabels=celltypes,
#         yticklabels=celltypes,
#         cmap="RdBu_r",
#         center=0,
#         vmin=-200, vmax=200,  # keep consistent scale
#         square=True,
#         cbar=(k == n_r - 1),  # only show one colorbar
#         ax=axes[k]
#     )

#     axes[k].set_title(f"r = {radii[k]}")

# plt.suptitle(f"{core} - Z-score per radius", y=1.02)
# plt.tight_layout()

# plt.savefig(
#     f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/{core}_z_heatmap_per_radius.png",
#     bbox_inches='tight'
# )
# plt.close()

# # z-score curves
# for i, name in enumerate(interactions):
#     plt.figure(figsize=(5,4))
    
#     plt.plot(radii, z_mean[i])
#     plt.axhline(0, linestyle="--")
    
#     plt.xlabel("Radius")
#     plt.ylabel("Z-score")
#     plt.title(name)
    
#     plt.tight_layout()
#     plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/interactions/zcurve_{core}_{name}.png", bbox_inches='tight')
#     plt.close()

# # multiple interactions in one plot
# plt.figure(figsize=(6,5))

# for i in range(len(interactions)):
#     plt.plot(radii, z_mean[i], alpha=0.3)

# plt.axhline(0, linestyle="--")
# plt.xlabel("Radius")
# plt.ylabel("Z-score")
# plt.title(f"{core} - all interactions")

# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/interactions/{core}_z_all.png", bbox_inches='tight')
# plt.close()




# # print("Radii:", radii)
# # print("Interactions:", interactions)
# # print(adata_s.uns["ripley_sim_mean"])
# # print(adata_s.uns["ripley_sim_std"])

# n_r = len(adata_s.uns["ripley_params"]["radii"])
# radii = adata_s.uns["ripley_params"]["radii"]

# fig, axes = plt.subplots(nrows=1, ncols=n_r, figsize=(15, 10))
# axes = axes.flatten()

# for k in range(n_r):
#     std_k = []
    
#     for ss in adata_s.uns["ripley_sim_std"].values():
#         std_k.append(ss[:, k])   # all interactions for radius k
    
#     std_k = np.concatenate(std_k)
#     std_k = std_k[np.isfinite(std_k)]
    
#     axes[k].hist(np.log10(std_k), bins=20)
#     axes[k].set_title(f"r = {radii[k]}")
#     axes[k].set_xlabel("log10(std)")
#     axes[k].set_ylabel("count")

# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/cross_ripley_sim_log_std_per_radius.png", bbox_inches='tight')
# plt.close()


# fig, axes = plt.subplots(nrows=1, ncols=n_r, figsize=(15, 10))
# axes = axes.flatten()

# for k in range(n_r):
#     std_k = []
    
#     for ss in adata_s.uns["ripley_sim_std"].values():
#         std_k.append(ss[:, k])   # all interactions for radius k
    
#     std_k = np.concatenate(std_k)
#     std_k = std_k[np.isfinite(std_k)]
    
#     axes[k].hist(std_k, bins=20)
#     axes[k].set_title(f"r = {radii[k]}")
#     axes[k].set_xlabel("std")
#     axes[k].set_ylabel("count")

# plt.tight_layout()
# plt.savefig(f"/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}_old/spatial/patching/{patch_size}um_{overlap}um_gpt/cross_ripley_sim_std_per_radius.png",bbox_inches='tight')
# plt.close()


# # adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/adata_per_sample/Neutro_Epi_extImm_pooled_A_EM_N/T23_004535_110005_1.h5ad')

# # adata.obs['celltypenew'] = adata.obs[phenotyping_level].copy()

# # adata.obs['celltypenew'] = adata.obs['celltypenew'].cat.add_categories(['non_malignant'])

# # adata.obs['celltypenew'] = adata.obs['celltypenew'].where(
# #     adata.obs['celltypenew'].isin(['Tumor_cells', 'Stromal']),
# #     'non_malignant'
# # )

# # sq.pl.spatial_scatter(adata, color="celltypenew", size=1,shape=None,figsize=(10,10))
# # plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/Neutro_Epi_extImm_pooled_A_EM_N/spatial/per_sample/T23_004535_110005_1/Stromal_Tumor_cells_spatial_scatter.png', dpi=300, bbox_inches='tight')



