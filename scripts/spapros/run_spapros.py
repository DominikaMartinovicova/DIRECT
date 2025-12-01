print('Running run_spapros.py...')

import pandas as pd
import scanpy as sc
import spapros as sp
import matplotlib.pyplot as plt

sc.settings.verbosity = 1
sc.logging.print_header()
print(f"spapros=={sp.__version__}")


# Load dataset
adata = sc.read_h5ad('/net/beegfs/groups/tgac/dmartinovicova_new/NSCLC/scRNAseq/data/final_scRNAseq_atlas_Salcher_5k.h5ad')
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

# Set up an Evaluator
print('Initiating evaluator...')
evaluator = sp.ev.ProbesetEvaluator(adata, celltype_key='Neutro_Epi_extImm', scheme="full", verbosity=2, results_dir='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros_5k')

print('Evaluating...')
evaluator.evaluate_probeset(valid_genes, set_id='xenium_io')

print(evaluator.summary_results)

# Plot
evaluator.plot_summary(save='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros_5k/plot_summary_5k.png', dpi=300)
evaluator.plot_confusion_matrix(save='/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/results/spapros_5k/plot_conf_matrix_5k.png', dpi=300)