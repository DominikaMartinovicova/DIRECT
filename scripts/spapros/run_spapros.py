#
# Evaluate Xenium IO probeset panel on NSCLC scRNA-seq dataset
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        python run_spapros.py

print('Running run_spapros.py...')

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import Libraries
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import pandas as pd
import scanpy as sc
import spapros as sp
import matplotlib.pyplot as plt
import numpy as np
import pickle

sc.settings.verbosity = 1
sc.logging.print_header()
print(f"spapros=={sp.__version__}")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Load scRNAseq dataset and probeset
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Load data
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/NSCLC/scRNAseq/data/final_scRNAseq_atlas_Salcher.h5ad')
print(adata)
gene_list = pd.read_csv("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/raw/slide_3/cell_feature_matrix/features.tsv.gz", 
                 sep="\t", 
                 header=None, 
                 compression="gzip")

# Filter out control probes
controls = ['Codeword', 'Probe']
controls = "|".join(controls)
gene_list = gene_list[~gene_list[1].str.contains(controls)][1].to_list()
print(gene_list)

# Check which genes are present in the dataset
valid_genes = [g for g in gene_list if g in adata.var_names]
missing = set(gene_list) - set(valid_genes)
print("Missing genes:", missing)
print("Using", len(valid_genes), "valid genes")

# Detect genes with zero variance
stds = np.std(adata.X.toarray(), axis=0)
zero_var_genes = np.where(stds == 0)[0]
zero_var_gene_names = adata.var_names[zero_var_genes]
print("Zero-variance genes:", zero_var_gene_names)
print("Count:", len(zero_var_genes))

zero_var_in_list = [g for g in zero_var_gene_names if g in gene_list]
print("Zero-variance genes in the gene list:", zero_var_in_list)

# Remove zero-variance genes from valid_genes
valid_genes = [g for g in valid_genes if g not in zero_var_gene_names]
print("Using", len(valid_genes), "valid genes after removing zero-variance genes")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Run spapros evaluation
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set up an Evaluator
print('Initiating evaluator...')
evaluator = sp.ev.ProbesetEvaluator(adata, celltype_key='Neutro_Epi_extImm', scheme="full", verbosity=2, results_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros')

print('Evaluating...')
evaluator.evaluate_probeset(valid_genes, set_id='xenium_io')

print(evaluator.summary_results)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Save results
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot
plt.figure()
evaluator.plot_summary()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/plot_summary.png', dpi=300, bbox_inches='tight')
plt.close()

plt.figure()
evaluator.plot_confusion_matrix()
plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/plot_conf_matrix.png', dpi=300, bbox_inches='tight')
plt.close()

# plt.figure()
# evaluator.plot_marker_corr()
# plt.savefig('/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros_5k/plot_marker_corr_5k.png', dpi=300, bbox_inches='tight')
# plt.close()

pickleable_data = {
    "results": evaluator.results,
    "summary": evaluator.summary_results,
}

with open("/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros/evaluator_results.pkl", "wb") as f:
    pickle.dump(pickleable_data, f)


