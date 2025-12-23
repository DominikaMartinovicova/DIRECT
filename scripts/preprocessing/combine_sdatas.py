
#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# combine_sdatas.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine raw counts sdata from all the slides to one dataset and preprocess together
# to analyze spatial patterns.
#
#   0 Import libraries and parse arguments
#   1 Read data and check if raw layers contain raw counts (integers)
#   2 Preprocessing
#       a. Normalization
#       b. PCA
#       c. kNN
#       d. PAGA
#       e. UMAP
#       f. Leiden clustering
#   3 Save
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/preprocessing/combine_sdatas.py \
	    --input_dir {params.in_dir} \
        --phen_level {params.phen_level} \
        -o {output.combined_sdatas} \
        --output_plot {output.output_plots}
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
from spatialdata.transformations import Translation
import spatialdata_plot
from shapely.affinity import translate

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 combine_adatas.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  # Combine raw counts adata from all the slides to one sdata dataset and preprocess together to analyze spatial patterns.  ')
    parser.add_argument('-i', help='path to phenotyped Xenium dirs metadata file',
                        dest='input',
                        type=str)
    parser.add_argument('--input_dir', help='path to phenotyped Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('--phen_level', help='phenotyping level',
                        dest='phen_level')
    parser.add_argument('-o', help='path to output combined xenium dirs metadata file',
                        dest='output',
                        type=str)
    parser.add_argument('--output_plot', help='path to output combined xenium plots dir',
                        dest='output_plot',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
#os.makedirs(args.output_plot, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# From input directory read all zarr files and combine them into one adata
#--------------------------------------------------------------------------------
#input_dir = args.input_dir
phen_level = 'Neutro_Epi_extImm'
input_dir = f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/phenotyped/{phen_level}/'

sdatas_dict={}
for folder in os.listdir(input_dir):    # Loop over all folders and files in the directory
    if folder.endswith(".zarr"):    # Only select .zarr folders 
        slide = folder.replace('.zarr','')
        print("Reading " + folder)
        sdata_tmp = sd.read_zarr(input_dir + folder)    # Read the spatial data
        sdata_tmp.tables['table'].obs_names = sdata_tmp.tables['table'].obs_names + '_' + slide  # Make obs names unique by adding slide name
        sdata_tmp.tables["table"].obs['slide'] = slide  # Add slide name to obs
        sdatas_dict[slide] = sdata_tmp     # Create a dictionary with all sdatas
print(sdatas_dict)
sdata_combined = sd.concatenate(sdatas_dict)  # Combine all sdatas into one      
print("Combined sdata:")
print(sdata_combined)











#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# X Save data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#os.rmdir(args.output_dir) # Remove empty .zarr dir snakemake created, else it won't (over)write it
#sdata_combined.write(args.output)
sdata_combined.write(f'/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/data/combined/{phen_level}_combined_sdatas.zarr')

