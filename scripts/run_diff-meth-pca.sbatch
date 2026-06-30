#!/bin/bash
#SBATCH --job-name=diff-meth-pca_v2.4
#SBATCH --output=logs/diff-meth-pca_v2.4_%j.out
#SBATCH --error=logs/diff-meth-pca_v2.4_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=64G
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

Rscript code/diff-meth-pca_v2.0.R

echo "Finished at $(date)"
