#!/bin/bash
#SBATCH --job-name=v2.1homer_dmrs
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=logs/v2.1homer_%j.out
#SBATCH --error=logs/v2.1homer_%j.err

module purge
module load release/23.04
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

export R_LIBS_USER=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR2

# global variables - check these
export VERSION=v2.1

export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR2
export WORKDIR=/home/shli842i/p_dna15/DMRichR/$VERSION
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

echo "Starting run $VERSION at $(date)"

Rscript --vanilla /home/shli842i/p_dna15/DMRichR/code/HOMER_v1.R

echo "Finished at $(date)"

source /software/rapids/r24.04/Miniconda3/23.10.0-1/etc/profile.d/conda.sh
conda activate "$HOME/conda-envs/homer"

bash "$HOME/p_dna15/DMRichR/HOMER.sh" \
    mm39 \
    "$SLURM_CPUS_PER_TASK" \
    "$WORKDIR"
