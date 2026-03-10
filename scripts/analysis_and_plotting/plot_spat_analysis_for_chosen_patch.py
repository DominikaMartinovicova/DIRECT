import matplotlib.pyplot as plt
import pickle
import os
import seaborn as sns
import pandas as pd
import squidpy as sq
import numpy as np

patch = 'T23_004535_110005_1_window_0'
overlap=50
patch_size=5000
phenotyping_level = 'Neutro_Epi_extImm_pooled_A_EM_N'
input_dir = f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/analysis/{phenotyping_level}/spatial/per_patch/{patch_size}um_{overlap}um/{patch}/'

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Prepare data and run analysis
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
centrality_scores=pd.read_csv(os.path.join(input_dir, 'centrality_scores.csv'), index_col=0)
print(centrality_scores)
# Plot centrality scores
y='closeness_centrality'
plt.figure(figsize=(8,6))
sns.barplot(data=centrality_scores, x=centrality_scores.index, y=y)
plt.title(f'Centrality Scores for Patch {patch}')
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_{y}.png')
plt.close()

nhood_enrichment_path = os.path.join(input_dir, 'neighbors_enrichment.pkl')
with open(nhood_enrichment_path, 'rb') as f:
    nhood_enrichment = pickle.load(f)

plt.figure(figsize=(8,8))
sns.heatmap(nhood_enrichment['zscore'], cmap='vlag', center=0, xticklabels=centrality_scores.index, yticklabels=centrality_scores.index)
plt.title(f'Neighborhood enrichment for Patch {patch}')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_nhood_enrichment.png')
plt.close()

interaction_matrix_path = os.path.join(input_dir, 'interaction_matrix.pkl')
with open(interaction_matrix_path, 'rb') as f:
    interaction_matrices = pickle.load(f)
plt.figure(figsize=(8,8))
sns.heatmap(interaction_matrices, cmap='vlag', center=0, xticklabels=centrality_scores.index, yticklabels=centrality_scores.index)
plt.title(f'Interaction matrix for Patch {patch}')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_interaction_matrix.png')
plt.close()


coocurrence_path = os.path.join(input_dir, 'co_occurrence_probabilities.pkl')
with open(coocurrence_path, 'rb') as f:
    cooccurrence_matrices = pickle.load(f)

print(cooccurrence_matrices)
occ = cooccurrence_matrices["occ"]                 # (cell, cell, distance)
distances = cooccurrence_matrices["interval"]
cell_types = cooccurrence_matrices["cell_type_key"]

b_idx = cell_types.index("B cell")

dist_edges = cooccurrence_matrices["interval"]              # length 50
dist_centers = 0.5 * (dist_edges[:-1] + dist_edges[1:])  # length 49

# choose a colormap (tab20 has 20 distinct colors)
cmap = plt.get_cmap("tab20")
colors = [cmap(i / len(cell_types)) for i in range(len(cell_types))]

fig, ax = plt.subplots(figsize=(8, 6))

for j, ct in enumerate(cell_types):
    ax.plot(
        dist_centers,
        occ[b_idx, j, :],
        label=ct,
        linewidth=1,
        color=colors[j]
    )

# reference line (random expectation)
ax.axhline(1, linestyle="--", linewidth=1)

ax.set_xlabel("Distance (µm)")
ax.set_ylabel("Co-occurrence score")
ax.set_title("B cell co-occurrence vs distance (all cell types)")

ax.legend(
    bbox_to_anchor=(1.05, 1),
    loc="upper left",
    fontsize=8,
    frameon=False
)

plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_cooccurrence_mean.png')
plt.close()



# Plot ripley's L statistics
ripleys_L_path = os.path.join(input_dir, 'dict_ripleys_L.pkl')
with open(ripleys_L_path, 'rb') as f:
    dict_ripley = pickle.load(f)

df = dict_ripley['L_stat']

plt.figure(figsize=(12, 7))
sns.lineplot(data=df,x='bins',y='stats',hue='Neutro_Epi_extImm_pooled_A_EM_N',estimator=None, palette="tab20")  # ensures raw lines are drawn

plt.legend(title='Cell type', bbox_to_anchor=(1.02, 1), loc='upper left')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/{patch}_ripleys_L.png', bbox_inches='tight')


# Plot cross Ripley's statistics for cells of interest
type_i = 'B_cell'
type_j = 'T_cell_regulatory'
ripleys_L_path = os.path.join(input_dir, 'dict_ripleys_L.pkl')
with open(os.path.join(ripleys_L_path, f"cross_ripley_{type_i}_vs_{type_j}.pkl"), "rb") as f:
    res = pickle.load(f)


r = res["r"]
observed = res["L_observed"]
sims = res["L_simulations"]
csr = res["csr_expectation"]

# simulation envelope (5–95%)
lower = np.percentile(sims, 5, axis=0)
upper = np.percentile(sims, 95, axis=0)

plt.figure(figsize=(6,5))

# simulation envelope
plt.fill_between(r, lower, upper, alpha=0.3, label="Permutation envelope (5–95%)")

# observed L
plt.plot(r, observed, linewidth=2, label="Observed L(r)")

# CSR expectation
plt.plot(r, csr, linestyle="--", label="CSR expectation")

plt.xlabel("Radius (r)")
plt.ylabel("L(r)")
plt.title(f"Cross-type Ripley's L: {res['type_i']} vs {res['type_j']}")
plt.legend()
plt.tight_layout()

plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/analysis/{phenotyping_level}/spatial/patching/{patch_size}um_{overlap}um/cross_ripley_{type_i}_vs_{type_j}.png', bbox_inches='tight')




