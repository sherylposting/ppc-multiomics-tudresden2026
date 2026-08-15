#!/bin/bash
#SBATCH --job-name=v1.3homer_dmrs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/v1.3homer_%j.out
#SBATCH --error=logs/v1.3homer_%j.err

source /software/rapids/r24.04/Miniconda3/23.10.0-1/etc/profile.d/conda.sh
conda activate "$HOME/conda-envs/homer"

bash "$HOME/p_dna15/DMRichR/HOMER.sh" \
    mm39 \
    "$SLURM_CPUS_PER_TASK" \
    "/home/shli842i/p_dna15/DMRichR/v1.3"
