#import pandas as pd
from matplotlib.gridspec import GridSpec
import scanpy as sc
import spapros as sp
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import pickle

with open('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/evaluator_results.pkl', 'rb') as f:
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
plt.suptitle('Spapros Evaluator Summary Metrics')
plt.savefig(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/plot_summary.svg', format='svg', bbox_inches='tight')
plt.close()

class_matrix = results['results']['forest_clfs']['xenium_io']
plt.figure(figsize=(18,18))
sns.heatmap(class_matrix, annot=True, fmt='.2f')
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/conf_matrix.svg', format='svg', bbox_inches='tight')
plt.title('Confusion Matrix Classifiation Accuracy')
plt.tight_layout()
plt.close()

# Subset Tcells
tcell_class_matrix = class_matrix.loc[class_matrix.index.str.contains('T cell'), class_matrix.columns.str.contains('T cell')]
plt.figure(figsize=(8,8))
sns.heatmap(tcell_class_matrix, annot=True, fmt='.2f')
plt.title('Confusion Matrix Classification Accuracy - T Cells')
plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/tcell_conf_matrix.svg', format='svg',bbox_inches='tight')
plt.close()

genecorr_matrix = results['results']['gene_corr']['xenium_io']
plt.figure(figsize=(12,10))
sns.heatmap(genecorr_matrix, cmap='bwr')
plt.title('Gene Correlation Matrix')
plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/genecorr_matrix.svg', format='svg', bbox_inches='tight')
plt.close()

clust_sim = results['results']['cluster_similarity']['xenium_io']
plt.figure()
sns.lineplot(data=clust_sim)
plt.title('Cluster Similarity')
plt.xlabel('Number of clusters')
plt.ylabel('NMI')
plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/cluster_similarity.svg', format='svg', bbox_inches='tight')
plt.close()

knn_overlap = results['results']['knn_overlap']['xenium_io']
plt.figure()
sns.lineplot(data=knn_overlap)
plt.title('KNN Overlap')
plt.xlabel('Number of clusters')
plt.ylabel('Mean kNN overlap')
plt.tight_layout()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/plots/spapros/knn_overlap.svg', format='svg', bbox_inches='tight')
plt.close()