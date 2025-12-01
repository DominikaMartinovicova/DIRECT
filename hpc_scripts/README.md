# HPC Scripts

This directory contains backup copies of scripts from the HPC cluster for the DIRECT project analysis.

## Purpose

- **Backup**: Prevent accidental removal of important scripts
- **Version Control**: Track changes to scripts over time
- **Sync**: Keep local and HPC cluster scripts synchronized

## Usage

### Backing up scripts from HPC

To back up scripts from the HPC cluster to this repository:

```bash
# Using rsync (recommended - handles subdirectories and special characters)
rsync -avz user@hpc:/path/to/scripts/ ./hpc_scripts/

# Or using scp with recursive flag for directories
scp -r user@hpc:/path/to/scripts/* ./hpc_scripts/
```

### Syncing scripts to HPC

To sync scripts from this repository back to the HPC cluster:

```bash
# Using rsync (recommended)
# Add --delete flag to remove files on destination that don't exist in source
rsync -avz --delete ./hpc_scripts/ user@hpc:/path/to/scripts/

# Or using scp with recursive flag
scp -r ./hpc_scripts/* user@hpc:/path/to/scripts/
```

**Note**: Be cautious with the `--delete` flag as it will remove files on the destination that don't exist in the source.

## Directory Structure

Add your scripts here organized by purpose, for example:
- `preprocessing/` - Data preprocessing scripts
- `analysis/` - Analysis scripts
- `jobs/` - SLURM/PBS job submission scripts
- `utils/` - Utility scripts

## Notes

- Always commit and push changes after adding new scripts
- Use meaningful commit messages to track what changed
- Consider adding script documentation in comments
