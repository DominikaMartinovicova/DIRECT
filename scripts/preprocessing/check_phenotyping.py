#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# check_phenotyping.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#
#
#
# Adapted by: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/python/check_phenotyping.py \
        -i {input.combined_adatas} \
        -o {output.checked_adatas} \
        --output_plot {params.out_plot_dir}
"""


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
import anndata as ad
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import csr_matrix
import argparse
import os
import leidenalg

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 check_phenotyping.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Create celltype specific signature matrices  ')
    parser.add_argument('-i', help='path to combined Xenium dirs metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('-threads', help='n threads to use',
                        dest='threads',
                        type=int)
    parser.add_argument('--phen_level', help='phenotyping level',
                        dest='phen_level',
                        type=str)
    parser.add_argument('-o', help='path to output checked xenium dirs metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_plot', help='path to plot of UMAP',
                        dest='output_plot',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
os.makedirs(args.output_plot, exist_ok=True)
os.makedirs(args.output_plot + '/violinplots/', exist_ok=True)
os.makedirs(args.output_plot + '/dotplots/', exist_ok=True)
os.makedirs(args.output_plot + '/umaps/', exist_ok=True)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
print('Reading adata...')
adata_combined= sc.read_h5ad(args.input)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Plotting
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Plot UMAP with celltypes and clinical labels
#-------------------------------------------------------------------------------
print('Plotting...')
colors = [args.phen_level, 'slide','sample_type', 'pt_id']

for color in colors:
    sc.pl.umap(adata_combined, color=color, show = False, size=5, layer='raw_counts')
    plt.savefig(args.output_plot + '/UMAP_' + color + '.png', bbox_inches='tight')
    plt.close()

# Plot UMAP with gene expression
#-------------------------------------------------------------------------------
markers = [
    # B cell
    'ADAM28', 'BANK1', 'CD19', 'CD27', 'CD37', 'CD79A', 'CD79B', 'CD83', 'FCMR', 'IGHM', 'MS4A1', 'SELL', 'TNFRSF13C',
    # cDC1 (includes cDC2, DC mature)
    'CD1A', 'CD1C', 'CD40', 'CD80', 'CHIT1', 'CLEC10A', 'CXCL14', 'FABP3', 'FCER1A', 'FLT3LG', 'GPR34', 'IL23A', 'ITGAE', 'MMP12', 'PLA2G7', 
    'PLD4', 'S100B', 'SLC1A3', 'TREM2',
    # Endothelial cell
    'PLVAP',
    # Stromal
    'ACTA2', 'DCN', 'FN1', 'LUM', 'NOTCH3', 'PDGFRA', 'RGS5',
    # Macrophage (includes Macrophage alveolar)
    'AIF1','CD44', 'CD68', 'CD74', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D', 'CXCL1', 'CXCL5', 'FCGR1A', 'FLT1', 'GLIPR2', 'IL4', 'MARCO', 'VCAN',
    # Mast cell
    'AREG', 'CPA3', 'HPGDS', 'MS4A2', 'P2RX1', 'PTGS1',
    # Monocyte classical
    'CD14', 'CSF1', 'CSF1R', 'FCGR3A', 'ITGAM', 'LILRA5', 'LILRB2', 'C1QB', 'CLEC12A',

    # Monocyte non-classical
    # (no unique markers in the given gene list, left blank)

    # Neutrophils (includes NAN and TAN pooled)
    'S100A9',
    # NK cell (includes T cell NK-like)
    'CST7', 'CD247', 'IQGAP2', 'KLRD1', 'KLRK1', 'NKG7', 'RUNX3', 'SAMD3', 'SMAD3', 'UBE2C',
    # pDC
    'IRF8', 'LILRA4', 'LILRB4', 'MPEG1', 'SPIB', 'SYK',
    # Plasma cell
    'CD38', 'DERL3', 'FKBP11', 'IGHG1', 'IGHG2', 'IGHG3', 'IGHGP', 'JCHAIN', 'PRDX4', 'SDC1', 'SSR4', 'STAT3', 'TNFRSF17',
    # T cell (general)
    'BATF3', 'CD2', 'CD28', 'CD3E', 'CORO1A', 'FAS', 'FASLG', 'GPR171', 'IL2', 'IL4R', 'JAK1', 'STAT4', 'CX3CR1',
    # T cell CD4
    'CD4', 'TCF7',
    # T cell CD8 activated
    'CD3D', 'CD8A', 'EOMES', 'GNLY', 'GZMB', 'GZMK', 'KLRB1', 'CXCR6',
    # T cell CD8 terminally exhausted
    'DGKA', 'HAVCR2', 'LAG3', 'PDCD1', 'TOX', 'TIGIT',
    # T cell regulatory
    'CTLA4', 'FOXP3', 'IL2RA', 'IRF4', 'CXCL13',
    # Tumor cells
    'AKT1', 'ANPEP', 'CFC1', 'CEACAM8', 'CD47', 'CD274', 'EGFR', 'KRAS', 'RNF43', 'SOX9', 'CXCL2'
]


sc.pl.dotplot(adata_combined, markers, groupby = args.phen_level)
plt.savefig(args.output_plot + '/dotplots/dotplot_celltype.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

#for marker in markers:
#    if marker in adata_combined.var_names:
#        sc.pl.umap(adata_combined, color=marker, show = False, size=5, layer='raw_counts')
#        plt.savefig(args.output_plot + '/umaps/rUMAP_' + marker + '.png', bbox_inches='tight')
#        plt.close()
#        sc.pl.umap(adata_combined, color=marker, show = False, size=5)
#        plt.savefig(args.output_plot + '/umaps/nUMAP_' + marker + '.png', bbox_inches='tight')
#        plt.close()
#    else:
#        print(f'{marker} is NOT in adata')

# # Plot selected markers as violin plots
# #-------------------------------------------------------------------------------
# markers = ['CD19','CD79A','MS4A1', 'BANK1', # B cells
#            'PLVAP','RGS5', # Endothelial
#            'AIF1','MARCO', 'CD68','CD163', # Macrophage
#            'CPA3','HPGDS','MS4A2', # Mast cell
#            'VCAN','FCGR3A','ITGAM','CSF1R','LILRB2', # Monocyte
#            'GNLY','KLRD1','KLRK1','NKG7','CD247', # NK
#            'S100A9', # Neutrophil
#            'IGHG1','IGHG2','IGHG3','IGHGP','JCHAIN', # Plasma cell
#            'LUM','DCN','FN1','ACTA2','PDGFRA', # Stromal
#            'CD3E','CD3D','CD2','CD28','STAT4', # T cell
#            'CD4','TCF7', # T cell CD4
#            'CD8A','CD8B','CXCR6','GZMK', # T cell CD8
#            'CTLA4','FOXP3','IL2RA','CXCL13', # T reg
#            'CEACAM8','AKT1','EGFR','KRAS', # Tumor cell
#            'CLEC10A','CD1A','CD1C','CCL19', # DC
#            'LILRA4','IRF8','MPEG1','SPIB', # pDC
#            ]

# print('Plotting dotplot...')
# sc.pl.dotplot(adata_combined, markers, groupby = 'celltype')
# plt.savefig(args.output_plot + '/dotplots/selected_dotplot_celltype.png', dpi=300, bbox_inches='tight')
# plt.close()

# print('Plotting raw violinplot...')
# sc.pl.stacked_violin(adata_combined, markers, groupby='celltype', show = False, layer='raw_counts', swap_axes=True)
# plt.savefig(args.output_plot + '/violinplots/raw_selected_violin_celltype.png', bbox_inches='tight')
# plt.close()

# print('Plotting nlog1p violinplot...')
# sc.pl.stacked_violin(adata_combined, markers, groupby='celltype', show = False, swap_axes=True)
# plt.savefig(args.output_plot + '/violinplots/nlog1p_selected_violin_celltype.png', bbox_inches='tight')
# plt.close()


# Plot subselected markers
#-------------------------------------------------------------------------------
markers = [
    'MS4A1','TNFRSF13C', # B cell
    'PLVAP', 'FLT1','RGS5', # Endothelial cell
    'CD74', 'CD163', 'CD68', 'MARCO',  # Macrophage (includes Macrophage alveolar)
    'CPA3','MS4A2', # Mast cell
    'CSF1R', 'LILRB2', 'CLEC12A', 'VCAN', # Monocyte classical
    'S100A9','SELL', # Neutrophils (includes NAN and TAN pooled)
    'GNLY', 'NKG7','KLRD1',  # NK cell (includes T cell NK-like)
    'IGHGP', 'JCHAIN','IGHG2', # Plasma cell
    'DCN', 'LUM', 'FN1', # Stromal
    'CD2', 'CD3E',  # T cell (general)
    'CD4', 'TCF7', # T cell CD4
    'CD8A', 'GZMB', 'GZMK', 'CXCR6', # T cell CD8 activated
    'DGKA', 'HAVCR2', 'LAG3', 'TOX', 'TIGIT', # T cell CD8 terminally exhausted
    'CTLA4', 'FOXP3', 'IL2RA', 'CXCL13', # T cell regulatory
    'EGFR', 'KRAS','ANPEP',     # Tumor cells
    'ITGAE','S100B','CLEC10A','SLC1A3', # cDC1 (includes cDC2, DC mature)
    'IRF8', 'MPEG1', 'LILRA4', 'SPIB','PLD4'  # pDC
]

print('Plotting dotplot...')
sc.pl.dotplot(adata_combined, markers, groupby = args.phen_level)
plt.savefig(args.output_plot + '/dotplots/subselected_dotplot_celltype.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting raw violinplot...')
sc.pl.stacked_violin(adata_combined, markers, groupby=args.phen_level, show = False, layer='raw_counts', swap_axes=True)
plt.savefig(args.output_plot + '/violinplots/raw_subselected_violin_celltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting nlog1p violinplot...')
sc.pl.stacked_violin(adata_combined, markers, groupby=args.phen_level, show = False, swap_axes=True)
plt.savefig(args.output_plot + '/violinplots/nlog1p_subselected_violin_celltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()


# List of non immune cell types you want to exclude
non_immune_cells = [
    "Alveolar cell type 1",
    "Alveolar cell type 2",
    "Ciliated",
    "Club",
    "Endothelial cell",
    "Epithelial cell",
    "Stromal",
    "transitional club/AT2",
    "Tumor cells",
    "other"
]

# Plot subselected markers
#-------------------------------------------------------------------------------
markers = [
    'MS4A1','TNFRSF13C', # B cell
    'CD74', 'CD163', 'CD68', 'MARCO',  # Macrophage (includes Macrophage alveolar)
    'CPA3','MS4A2', # Mast cell
    'CSF1R', 'LILRB2', 'CLEC12A', 'VCAN', # Monocyte classical
    'S100A9','SELL', # Neutrophils (includes NAN and TAN pooled)
    'GNLY', 'NKG7','KLRD1',  # NK cell (includes T cell NK-like)
    'IGHGP', 'JCHAIN','IGHG2', # Plasma cell
    'CD2', 'CD3E',  # T cell (general)
    'CD4', 'TCF7', # T cell CD4
    'CD8A', 'GZMB', 'GZMK', 'CXCR6', # T cell CD8 activated
    'DGKA', 'HAVCR2', 'LAG3', 'TOX', 'TIGIT', # T cell CD8 terminally exhausted
    'CTLA4', 'FOXP3', 'IL2RA', 'CXCL13', # T cell regulatory
    'ITGAE','S100B','CLEC10A','SLC1A3', # cDC1 (includes cDC2, DC mature)
    'IRF8', 'MPEG1', 'LILRA4', 'SPIB','PLD4'  # pDC
]


# Subset the AnnData object to only those cell types
adata_sub = adata_combined[~adata_combined.obs[args.phen_level].isin(non_immune_cells)].copy()
print('Plotting dotplot...')
sc.pl.dotplot(adata_sub, markers, groupby = args.phen_level)
plt.savefig(args.output_plot + '/dotplots/subselected_dotplot_subcelltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting raw violinplot...')
sc.pl.stacked_violin(adata_sub, markers, groupby=args.phen_level, show = False, layer='raw_counts')
plt.savefig(args.output_plot + '/violinplots/raw_subselected_violin_subcelltype_notswaped.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting nlog1p violinplot...')
sc.pl.stacked_violin(adata_sub, markers, groupby=args.phen_level, show = False)
plt.savefig(args.output_plot + '/violinplots/nlog1p_subselected_violin_subcelltype_notswaped.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()


print('Done!')