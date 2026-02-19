#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# patching.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Save each core or tissue piece as a separate adata object. 
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Create adata objects for each core or tissue piece
#       *Save cores individually as adata 
#
#
# Author: Dominika Martinovicova (d.martinovicova@amsterdamumc.nl)
#
# Usage:
#        """
#        """

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0 Import libraries and parse arguments
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
import scanpy as sc
import numpy as np
import pandas as pd
import os
import argparse
import warnings
warnings.filterwarnings("ignore")

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 core_level_adata.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Save each core or tissue piece as a separate adata object.  ')
    parser.add_argument('-i', help='path to adata sample subset',
                        dest='input',
                        type=str)
    parser.add_argument('-o', '--output_adata', help='path to output dir with adata per core or tissue piece',
                        dest='output_dir_adata',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()
os.makedirs(args.output_dir_adata, exist_ok=True)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
adata = sc.read_h5ad(args.input)
print(adata)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Create adata objects for each core or tissue piece
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cores = adata.obs['sample'].unique().tolist()

for core in cores:
    adata_core = adata[adata.obs['sample'] == core].copy()

    # Save adata object for each core or tissue piece
    output_path = os.path.join(args.output_dir_adata, f'{core}.h5ad')
    adata_core.write_h5ad(output_path)



