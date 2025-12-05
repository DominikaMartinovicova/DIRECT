#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  preprocess_Xenium.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Filter cells and genes based on QC metrics.
#
#   0 Import libraries and parse arguments
#   1 Assign cells to their corresponding samples based on coordinates
#   2 Add metadata about samples to adata.obs
#   3 Calculate QC metrics
#   4 Filter genes and cells
#   5 Check QC after filtering
#   6 Save the filtered data
#
# Author: Dominika Martinovicova (d.martiovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/preprocessing/preprocess_Xenium.py \
        -i {input.compressed_Xenium} \
        --input_coords {input.coordinates} \
        --input_meta {input.metadata} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        -o {output.preprocessed_Xenium} \
        --output_plot_QA {output.quality_plots} \
        --output_plot_spatial {output.spatial_plots}
"""
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import spatialdata as sd
from spatialdata_io import xenium
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import squidpy as sq
import argparse
import os
import numpy as np
import pandas as pd
from shapely.geometry import Point, Polygon

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 preprocess_Xenium.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Filter cells and genes based on QC metrics.  ')
    parser.add_argument('-i', help='path to compressed Xenium metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('--input_coords', help='path to coordinates of samples',
                        dest='input_coordinates',
                        type=str)
    parser.add_argument('--input_meta', help='path to metadata of samples',
                        dest='input_metadata',
                        type=str)
    parser.add_argument('--input_dir', help='path to compressed Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('--output_dir', help='path to preprocessed Xenium output dir',
                        dest='output_dir',
                        type=str)
    parser.add_argument('-o', help='path to output preprocessed Xenium metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_plot_QA', help='path to  quality plots',
                        dest='output_plot_QA',
                        type=str)
    parser.add_argument('--output_plot_spatial', help='path to spatial plots',
                        dest='output_plot_spatial',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

# Read data and create directories
#-------------------------------------------------------------------------------
os.makedirs(args.output_plot_QA, exist_ok=True) # Make directory for QC plots
os.makedirs(args.output_plot_spatial, exist_ok=True) # Make directory for spatial plots
sdata = sd.read_zarr(args.input_dir) # Read spatial data
adata = sdata.tables["table"] # Subset sdata on adata
print(adata)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Assign cells to their corresponding samples based on coordinates
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read in coordinates of each sample
df_coordinates = pd.read_csv(args.input_coordinates, comment='#')

# Create a dictionary with polygon coordinates for each sample
polygons_dict = {}
for sample_name, group in df_coordinates.groupby('Selection'):
    coordinates = list(zip(group['X'], group['Y']))     # Create list of (x, y) tuples
    polygons_dict[sample_name] = Polygon(coordinates)   # Save coordinates to the corresponding sample name

# Assign cells to samples based on polygon coordinates
sample_assignment=[]
for point_coords in adata.obsm['spatial']:
    point = Point(point_coords)                          # Save x, y coordinates into point variable
    assigned_sample = None 
    for sample_name, polygon in polygons_dict.items():
        if polygon.contains(point):                      # Check if point is within polygon
            assigned_sample = sample_name
            break                                        # Exit loop once the sample is found
    sample_assignment.append(assigned_sample)

# Add column to adata with specified sample assignment
adata.obs['sample'] = sample_assignment


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Add metadata about samples to adata.obs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
metadata = pd.read_csv(args.input_metadata, dtype={'pt_id':str, 'disease_stage':str})  # Patient ids are numbers, read as strings

obs = adata.obs.copy()
obs_merged = obs.merge(metadata, how='left', on = 'sample', sort=False)  # merge metadata to adata.obs based on sample names
obs_merged.index = obs.index  # keep original index
adata.obs=obs_merged

# Visualize if data looks correct
#-------------------------------------------------------------------------------
# List features you want to display
info = ['sample', 'T_number', 'sample_type', 'pt_id', 'treatment_scheme', 'regression']

# Loop over each feature and create spatial plot
for element in info:
    sq.pl.spatial_scatter(adata,library_id="spatial",shape=None,color=element,size = 1)
    output_file = args.output_plot_spatial + f'/spatial_{element}.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

# Remove cells that are NaN in sample assignment (not in any polygon)
#-------------------------------------------------------------------------------
print(f"Adata shape before removing cells not in any sample in {args.input_dir}: " + str(adata.shape))
adata = adata[~adata.obs['sample'].isna(), :]
print(f"Adata shape after removing cells not in any sample in {args.input_dir}: " + str(adata.shape))

sc.pl.highest_expr_genes(adata, n_top=20)
plt.savefig(args.output_plot_QA + "/before_filtering_highest_expr_genes.svg", format='svg')

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 3 Quality control
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sc.pp.calculate_qc_metrics(adata, percent_top=(10, 20, 50, 150), inplace=True)

cprobes = (adata.obs["control_probe_counts"].sum() / adata.obs["total_counts"].sum() * 100)
cwords = (adata.obs["control_codeword_counts"].sum() / adata.obs["total_counts"].sum() * 100)
print(f"Negative DNA probe count % in {args.input_dir}: {cprobes}")
print(f"Negative decoding count % in {args.input_dir}: {cwords}")

# Detection statistics
#-------------------------------------------------------------------------------
# Histograms
fig, axs = plt.subplots(1, 2, figsize=(10, 4))

axs[0].set_title("Total transcripts counts per cell")
axs[0].set_xlim(0, 2000)
sns.histplot(adata.obs["total_counts"], kde=True,bins=50, ax=axs[0])

axs[1].set_title("Unique transcripts per cell")
axs[1].set_xlim(0, 250)
sns.histplot(adata.obs["n_genes_by_counts"],kde=True,bins=50,ax=axs[1])

plt.subplots_adjust(wspace=0.4)  # Increase horizontal space between plots
plt.savefig(args.output_plot_QA + "/before_filtering_histQuality_plots_detection.svg", format='svg', bbox_inches='tight')

# Violinplot
sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts'], jitter=0.4, multi_panel=True, show = False)
plt.savefig(args.output_plot_QA + "/before_filtering_violinQuality_plots_detection.png")

# Scatterplot
print("Plotting and saving n_genes_by_counts vs total_counts")
sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts', show = False)
plt.savefig(args.output_plot_QA + "/before_filtering_total_counts_vs_ngenes.png")

# Spatial distribution of these metrics
fig, axes = plt.subplots(1, 2, figsize=(16, 5))
sc.pl.spatial(adata, color='total_counts', spot_size=10, title="Spatial Distribution of Total Transcript Counts", cmap='viridis_r', show=False, ax=axes[0])
sc.pl.spatial(adata, color='n_genes_by_counts', spot_size=10, title="Spatial Distribution of Number of Genes Detected per Cell", cmap='viridis_r', show=False, ax=axes[1])
plt.savefig(args.output_plot_spatial + "/before_filtering_spatialQuality_plots_detection.png")


# Segmentation statistics
#-------------------------------------------------------------------------------
# Histograms
fig, axs = plt.subplots(1, 3, figsize=(16, 5))

axs[0].set_title("Area of segmented cells")
axs[0].set_xlim(0, 500)
sns.histplot(adata.obs["cell_area"],kde=True,bins=50,ax=axs[0])

axs[1].set_title('Distribution of Nucleus Areas')
sns.histplot(adata.obs['nucleus_area'], bins=50, ax=axs[1], kde=True)

axs[2].set_title("Nucleus ratio")
plt.xlim(0,1)
sns.histplot(adata.obs["nucleus_area"] / adata.obs["cell_area"],bins=50, kde=True,ax=axs[2])

plt.savefig(args.output_plot_QA + "/before_filtering_histQuality_plots_segmentation.svg", format='svg')

# Scatterplot
fig, axs = plt.subplots(1, 2, figsize=(16, 5))

sns.scatterplot(x=adata.obs['cell_area'], y=adata.obs['nucleus_area'], alpha=0.5, ax=axs[0])
axs[0].set_title("Cell Area vs Nucleus Area")

sns.scatterplot(x=adata.obs['cell_area'], y=adata.obs['total_counts'], alpha=0.5, ax=axs[1])
axs[1].set_title("Cell Area vs Total Transcript Counts")

plt.savefig(args.output_plot_QA + "/before_filtering_scatterQuality_plots_segmentation.png")

# Spatial distribution of these metrics
fig, axes = plt.subplots(1, 2, figsize=(16, 5))
sc.pl.spatial(adata, color='cell_area', spot_size=10, title="Spatial Distribution of Cell Area", cmap='viridis_r', show=False, ax=axes[0])
sc.pl.spatial(adata, color='nucleus_area', spot_size=10, title="Spatial Distribution of Nucleus Area", cmap='viridis_r', show=False, ax=axes[1])
plt.savefig(args.output_plot_spatial + "/before_filtering_spatialQuality_plots_segmentation.png")


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 4 Filter cells and genes (no gene filtering due to limited gene panel)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Define cut offs
#--------------------------------------------------------------------------------
min_transcripts_per_cell = 20       # at least 20 transcripts per cell
#max_transcripts_per_cell = np.quantile(adata.obs['total_counts'], 0.98)       # at most 98 percentile
#min_cells_per_gene = 100            # gene expressed in at least 100 cells

# Check how many cells and genes would be filtered out
adata.obs['filtered'] = (adata.obs['transcript_counts'] < min_transcripts_per_cell) 

# Print the number of cells that would be filtered out
print(f"Number of cells to be filtered out from slide {args.input_dir}: {adata.obs['filtered'].sum()} out of {adata.n_obs}")

# Plot where the filtered cells appear on the slide
sq.pl.spatial_scatter(adata,library_id="spatial",shape=None,color='filtered',size = 0.5)
plt.savefig(args.output_plot_spatial + f'/spatial_filtered.png', dpi=300, bbox_inches='tight')
plt.close()

# # Calculate the number of cells each gene is expressed in
# gene_cell_counts = np.array((adata.X > 0).sum(axis=0)).flatten()
# genes_to_keep = gene_cell_counts >= min_cells_per_gene
# filtered_genes = adata.var_names[~genes_to_keep]
# print(f"Number of genes to be filtered out from {args.input_dir}: {len(filtered_genes)} out of {adata.n_vars}")
# print(f"Genes to be filtered out from {args.input_dir}:")
# print(filtered_genes)

# Apply filtering
sc.pp.filter_cells(adata, min_counts=min_transcripts_per_cell)
#sc.pp.filter_cells(adata, max_counts=max_transcripts_per_cell)
#sc.pp.filter_genes(adata, min_cells=min_cells_per_gene) no filtering of genes as quality seems good of all genes

# Remove the 'filtered' column as it's no longer needed
# adata.obs.drop(columns=['filtered'], inplace=True)


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 6 Save sdata with new preprocessed adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
sdata.tables['table'] = adata   # Change adata in sdata for preprocessed version
os.rmdir(args.output_dir) # Remove empty .zarr dir snakemake created, else it won't (over)write it
sdata.write(args.output_dir) # Save Xenium data in zarr obj for quick reading later















