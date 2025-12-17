#!/bin/bash
#SBATCH --job-name=spapros
#SBATCH --partition=defq                 # adjust for your cluster 
#SBATCH --cpus-per-task=32              # total CPU cores
#SBATCH --mem=400G  
#SBATCH --time=168:00:00
#SBATCH --output=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/spapros/%x_%j.out
#SBATCH --error=/net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/logs/spapros/%x_%j.err
#SBATCH --mail-user=d.martinovicova@amsterdamumc.nl
#SBATCH --mail-type=END,FAIL

# Activate environment 
source /home/P083831/miniconda3/etc/profile.d/conda.sh 

# activate your env 
conda activate /appdata/users/P083831/conda_envs/spapros

# Run (CPU)
srun python /net/beegfs/groups/tgac/dmartinovicova_new/DIRECT/scripts/spapros/run_spapros_2.py