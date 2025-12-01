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
colors = ['celltype', 'slide','sample_type', 'pt_id']

for color in colors:
    sc.pl.umap(adata_combined, color=color, show = False, size=5, layer='raw_counts')
    plt.savefig(args.output_plot + '/UMAP_' + color + '.png', bbox_inches='tight')
    plt.close()

# Plot UMAP with gene expression
#-------------------------------------------------------------------------------
markers = ['CD19','CD79A','MS4A1','ADAM28', 'BANK1','FCMR','SELL','CD79B','IGHM', # B cells
           'CHIT1','TREM2','FCER1A','FABP3','ITGAE','SLC1A3','PLD4','CLEC10A','CD1A','CD1C', # DC
           'PLVAP', # Endothelial
           'LUM','DCN','FN1','RGS5','NOTCH3','ACTA2','PDGFRA', # Stromal
           'AIF1','CXCL1','FCGR1A','MARCO','VCAN','EGF', 'CDKN1B','CD74','CD68','CD163','IL4', # Macrophage
           'AREG','CPA3','HPGDS','MS4A2','P2RX1','PTGS1', # Mast cell
           'FCGR3A','ITGAM','CSF1','CSF1R','CD14','LILRA5','LILRB2','C1QB','CLEC12A', # Monocyte
           'S100A9', # Neutrophil
           'IQGAP2','KLRD1','KLRK1','NKG7','RUNX3','SAMD3','SMAD3','UBE2C','CD247','CST7', # NK
           'LILRA4','LILRB4','IRF8','MPEG1','SPIB','SYK','IL23A',# pDC
           'DERL3','FKBP11','IGHG1','IGHG2','IGHG3','IGHGP','PRDX4','SDC1','SSR4','CD38','STAT3','JCHAIN','PRDX4','TNFRSF17', # Plasma cell
           'CD3E','BATF3','FAS','FASLG','GPR171','JAK1','IL4R','IL2','CX3CR1','CD2','CD28','CD8B','STAT4', # T cell
           'CD4','TCF7', # T cell CD4
           'CD8A','GNLY','KLRB1','CD3D','CXCR6','EOMES','GZMB','GZMK', # T cell CD8
           'LAG3','TOX','PDCD1','TIGIT','HAVCR2','DGKA',# T cell exhausted
           'FOXP3','IL2RA','CTLA4','IRF4','CXCL13', # T reg
           'ANPEP', 'SOX9','RNF43','CEACAM8','AKT1','CFC1','CD47','CD274','EGFR','KRAS','CXCL2' # Tumor cell
           ]

sc.pl.dotplot(adata_combined, markers, groupby = 'celltype')
plt.savefig(args.output_plot + '/dotplots/dotplot_celltype.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

for marker in markers:
    if marker in adata_combined.var_names:
        sc.pl.umap(adata_combined, color=marker, show = False, size=5, layer='raw_counts')
        plt.savefig(args.output_plot + '/umaps/rUMAP_' + marker + '.png', bbox_inches='tight')
        plt.close()
        sc.pl.umap(adata_combined, color=marker, show = False, size=5)
        plt.savefig(args.output_plot + '/umaps/nUMAP_' + marker + '.png', bbox_inches='tight')
        plt.close()
    else:
        print(f'{marker} is NOT in adata')

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
markers = ['MS4A1', 'BANK1', # B cells
           'CD3E','CD2', # T cell
           'CD4','TCF7', # T cell CD4
           'CD8A', # T cell CD8
           'CLEC10A','CD1C', # DC
           'LILRA4','IRF8','MPEG1', # pDC
           'PLVAP','RGS5', # Endothelial
           'DCN','PDGFRA', # Stromal
           'CD163','CD68', # Macrophage
           'CPA3','MS4A2', # Mast cell
           'VCAN', # Monocyte
           'GNLY','NKG7', # NK
           'S100A9', # Neutrophil
           'IGHG1','IGHGP','JCHAIN', # Plasma cell
           'CTLA4','FOXP3', # T reg
           'EGFR','KRAS' # Tumor cell
           ]

print('Plotting dotplot...')
sc.pl.dotplot(adata_combined, markers, groupby = 'celltype')
plt.savefig(args.output_plot + '/dotplots/subselected_dotplot_celltype.svg',format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting raw violinplot...')
sc.pl.stacked_violin(adata_combined, markers, groupby='celltype', show = False, layer='raw_counts', swap_axes=True)
plt.savefig(args.output_plot + '/violinplots/raw_subselected_violin_celltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting nlog1p violinplot...')
sc.pl.stacked_violin(adata_combined, markers, groupby='celltype', show = False, swap_axes=True)
plt.savefig(args.output_plot + '/violinplots/nlog1p_subselected_violin_celltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()


# List of cell types you want to include
selected_celltypes = ['B_cell', 'Dendritic_cell', 'Endothelial_cell', 'Macrophage', 'Mast_cell', 'Monocyte','NK_cell', 
                      'Neutrophil', 'Plasma_cell', 'Fibroblast', 'CD4_T_cell', 'CD8_T_cell', 'Regulatory_T_cell']  # change to your desired cell types

# Subset the AnnData object to only those cell types
adata_sub = adata_combined[adata_combined.obs['celltype'].isin(selected_celltypes)].copy()
print('Plotting dotplot...')
sc.pl.dotplot(adata_sub, markers, groupby = 'celltype')
plt.savefig(args.output_plot + '/dotplots/subselected_dotplot_subcelltype.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting raw violinplot...')
sc.pl.stacked_violin(adata_sub, markers, groupby='celltype', show = False, layer='raw_counts')
plt.savefig(args.output_plot + '/violinplots/raw_subselected_violin_subcelltype_notswaped.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()

print('Plotting nlog1p violinplot...')
sc.pl.stacked_violin(adata_sub, markers, groupby='celltype', show = False)
plt.savefig(args.output_plot + '/violinplots/nlog1p_subselected_violin_subcelltype_notswaped.svg', format='svg', dpi=300, bbox_inches='tight')
plt.close()


print('Done!')