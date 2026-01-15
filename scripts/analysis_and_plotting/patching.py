#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# patching.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
#   Divide each piece of tissue into patches. 
#
#   0 Import libraries and parse arguments
#   1 Read data
#   2 Divide each piece of tissue into patches
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
import squidpy as sq
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os
import pickle
import argparse

# Parse arguments from commandline
#--------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 patching.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Divide each piece of tissue into patches.  ')
    parser.add_argument('-i', help='path to adata sample subset',
                        dest='input',
                        type=str)
    parser.add_argument('-o_patches', help='path to output dir with patches per sample',
                        dest='output_dir_patches',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1 Read  data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Read adata
print('Reading data...')
adata = sc.read_h5ad(args.input)

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 2 Divide each piece of tissue into patches
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++