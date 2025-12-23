#import pandas as pd
from matplotlib.gridspec import GridSpec
import scanpy as sc
import spapros as sp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle
import holoviews as hv

pheno_level = "pooled_A_EM"

with open(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/evaluator_results_{pheno_level}.pkl', 'rb') as f:
    results = pickle.load(f)

print("Spapros Results:")
print(results['results'])


# Plot
summary = results['summary']
print("Summary:")
print(summary)
default_cmaps = {
    "cluster_similarity nmi_5_20": "Greens",
    "cluster_similarity nmi_21_60": "Greens",
    "knn_overlap mean_overlap_AUC": "Oranges",
    "forest_clfs accuracy": "Reds",
    "forest_clfs perct acc > 0.8": "Reds",
    "marker_corr": "Purples",
    "gene_corr 1 - mean": "Blues",
    "gene_corr perct max < 0.8": "Blues",
    "other": "Greys"}

gs = GridSpec(1, len(summary.columns))
plt.figure(figsize=(1 * len(summary.columns), 1.5))
for i, col in enumerate(summary.columns):
    yticklabels = bool(i == 0)
    ax = plt.subplot(gs[i])
    sns.heatmap(summary[[col]], annot=True, fmt='.2f', cmap=default_cmaps[col], ax=ax, cbar=False, vmin=0, vmax=1, yticklabels=yticklabels)
    plt.tick_params(axis="x", which="major", labelbottom=True, bottom=True, top=False, labeltop=False)
    plt.tick_params(axis="y", which="major")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    ax.set_yticklabels(ax.get_yticklabels(), rotation=90)
plt.suptitle(f'Spapros Evaluator Summary Metrics {pheno_level}')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/plot_summary_{pheno_level}.svg', format='svg', bbox_inches='tight')
plt.close()

class_matrix = results['results']['forest_clfs']['xenium_io']
class_matrix.to_csv(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/confusion_matrix_{pheno_level}.csv')
plt.figure(figsize=(18,18))
sns.heatmap(class_matrix, annot=True, fmt='.2f')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/conf_matrix_{pheno_level}.svg', format='svg', bbox_inches='tight')
plt.title('Confusion Matrix Classifiation Accuracy')
plt.tight_layout()
plt.close()

# Subset Tcells
tcell_class_matrix = class_matrix.loc[class_matrix.index.str.contains('T cell'), class_matrix.columns.str.contains('T cell')]
tcell_class_matrix.to_csv(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/tcell_confusion_matrix_{pheno_level}.csv')
plt.figure(figsize=(8,8))
sns.heatmap(tcell_class_matrix, annot=True, fmt='.2f')
plt.title(f'Confusion Matrix Classification Accuracy - T Cells {pheno_level}')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/tcell_conf_matrix_{pheno_level}.svg', format='svg',bbox_inches='tight')
plt.close()

genecorr_matrix = results['results']['gene_corr']['xenium_io']
plt.figure(figsize=(12,10))
sns.heatmap(genecorr_matrix, cmap='bwr')
plt.title(f'Gene Correlation Matrix {pheno_level}')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/genecorr_matrix_{pheno_level}.svg', format='svg', bbox_inches='tight')
plt.close()

clust_sim = results['results']['cluster_similarity']['xenium_io']
plt.figure()
sns.lineplot(data=clust_sim)
plt.title(f'Cluster Similarity {pheno_level}')
plt.xlabel('Number of clusters')
plt.ylabel('NMI')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/cluster_similarity_{pheno_level}.svg', format='svg', bbox_inches='tight')
plt.close()

knn_overlap = results['results']['knn_overlap']['xenium_io']
plt.figure()
sns.lineplot(data=knn_overlap)
plt.title(f'KNN Overlap {pheno_level}')
plt.xlabel('Number of clusters')
plt.ylabel('Mean kNN overlap')
plt.tight_layout()
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/knn_overlap_{pheno_level}.svg', format='svg', bbox_inches='tight')
plt.close()

# Plot a Sankey flow diagram for cell type mapping
# print("Plotting Sankey diagrams...")
# class_matrix_melted = class_matrix.reset_index().melt(id_vars='index')
# class_matrix_melted.columns = ['True Cell Type', 'Predicted Cell Type', 'Accuracy']
# print(class_matrix_melted.head())
# sankey = hv.Sankey(class_matrix_melted,label='Celltype Mapping')
# plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/sankey_diagram.svg', format='svg', bbox_inches='tight')
# plt.close()

# tcell_class_matrix_melted = tcell_class_matrix.reset_index().melt(id_vars='index')
# tcell_class_matrix_melted.columns = ['True Cell Type', 'Predicted Cell Type', 'Accuracy']
# tcell_sankey = hv.Sankey(tcell_class_matrix_melted, label='T Celltype Mapping')
# plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/tcell_sankey_diagram.svg', format='svg', bbox_inches='tight')
# plt.close()
