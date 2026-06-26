#!/bin/bash
#SBATCH --job-name=granges-annot
#SBATCH --output=logs/granges-annot_%j.out
#SBATCH --error=logs/granges-annot_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=1

module purge
module load release/25.06 GCC/13.3.0 R

cd /home/shli842i/p_dna15

mkdir -p scripts/logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

echo "TMPDIR=$TMPDIR"
echo "Starting at $(date)"

# check these
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks
export WORKDIR=/home/shli842i/p_dna15

Rscript code/granges-annot_v1.0.R

echo "Finished at $(date)"
