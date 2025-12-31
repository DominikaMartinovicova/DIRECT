#!/bin/bash
#SBATCH --job-name=Xenium_spatial_analysis
#SBATCH --partition=defq                 # adjust for your cluster 
#SBATCH --cpus-per-task=32              # total CPU cores
#SBATCH --mem=80G  
#SBATCH --time=168:00:00
#SBATCH --output=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/analysis/%x_%j.out
#SBATCH --error=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/analysis/%x_%j.err
#SBATCH --mail-user=d.martinovicova@amsterdamumc.nl
#SBATCH --mail-type=END,FAIL

# Activate environment 
source /home/P083831/miniconda3/etc/profile.d/conda.sh 

# activate your env 
conda activate /appdata/users/P083831/conda_envs/xenium_py310

# Run (CPU)
srun python /net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/scripts/analysis_and_plotting/spatial_analysis.py