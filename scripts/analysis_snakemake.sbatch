#!/bin/bash
#SBATCH --job-name=Xenium_spatial_analysis
#SBATCH --partition=defq                 # adjust for your cluster 
#SBATCH --cpus-per-task=20              # total CPU cores
#SBATCH --mem=150G  
#SBATCH --time=72:00:00
#SBATCH --output=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/analysis/%x_%j.out
#SBATCH --error=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/analysis/%x_%j.err
#SBATCH --mail-user=d.martinovicova@amsterdamumc.nl
#SBATCH --mail-type=END,FAIL

# Activate environment 
source /home/P083831/miniconda3/etc/profile.d/conda.sh 

# activate your env 
conda activate /appdata/users/P083831/conda_envs/snakemake

# Go to project root
cd /net/beegfs/groups/tgac/dmartinovicova_new/DIRECT

# Run (CPU)
snakemake -s Snakefile_statial -c 20 --use-conda