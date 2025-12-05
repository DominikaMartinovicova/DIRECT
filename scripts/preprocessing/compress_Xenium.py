#!/usr/bin/python3
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# compress_Xenium.py
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Convert raw Xenium data to compressed .zarr
#
# Author: Mischa Steketee (m.f.b.steketee@amsterdamumc.nl)
#
# Usage:
"""
        python3 scripts/preprocessing/compress_Xenium.py \
        -i {input.raw_Xenium} \
        --input_dir {params.in_dir} \
        --output_dir {params.out_dir} \
        --threads {threads} \
        -o {output.compressed_Xenium} 
"""
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Import Libraries
#-------------------------------------------------------------------------------
import spatialdata as sd
import spatialdata_io
from spatialdata_io import xenium
import argparse
import os

#-------------------------------------------------------------------------------
# 1.1 Parse command line arguments
#-------------------------------------------------------------------------------
def parse_args():
    "Parse inputs from commandline and returns them as a Namespace object."
    parser = argparse.ArgumentParser(prog = 'python3 compress_Xenium.py',
        formatter_class = argparse.RawTextHelpFormatter, description =
        '  Convert raw Xenium data to compressed .zarr  ')
    parser.add_argument('-i', help='path to cell feature',
                        dest='input',
                        type=str)
    parser.add_argument('--input_dir', help='path to raw Xenium dir',
                        dest='input_dir',
                        type=str)
    parser.add_argument('--output_dir', help='path to compressed output dir',
                        dest='output_dir',
                        type=str)
    parser.add_argument('--threads', help='path to config file',
                        dest='threads',
                        type=int)
    parser.add_argument('-o', help='path to output metadata file',
                        dest='output',
                        type=str)
    args = parser.parse_args()
    return args

args = parse_args()

print(args)
#-------------------------------------------------------------------------------
# 2.1 Read data
#-------------------------------------------------------------------------------
## Read spatial data from Xenium experiment
sdata = xenium(args.input_dir, cells_as_circles = False, n_jobs = args.threads)

#-------------------------------------------------------------------------------
# 3.1 Write to compressed dir
#-------------------------------------------------------------------------------
## Remove empty .zarr dir snakemake created, else it won't (over)write it
os.rmdir(args.output_dir)

## Save Xenium data in zarr obj for quick reading later
sdata.write(args.output_dir)
