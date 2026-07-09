#!/bin/bash
#SBATCH --job-name=v1.1_DMRichR
#SBATCH --output=logs/v1.1_DMRichR_%j.out
#SBATCH --error=logs/v1.1_DMRichR_%j.err
#SBATCH --time=72:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

module purge
module load release/23.04
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

export R_LIBS_USER=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1

# global variables - check these
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1
export WORKDIR=/home/shli842i/p_dna15/DMRichR
mkdir -p logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

# global variables - edit per run
export VERSION=v1.1

echo "Starting run $VERSION at $(date)"

Rscript /home/shli842i/p_dna15/DMRichR/code/DMRpreliminary_v1.R

echo "Finished at $(date)"
