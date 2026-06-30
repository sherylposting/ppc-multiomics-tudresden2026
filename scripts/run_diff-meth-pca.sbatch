#!/bin/bash
#SBATCH --job-name=diff-meth-pca_v3
#SBATCH --output=logs/diff-meth-pca_v3_%j.out
#SBATCH --error=logs/diff-meth-pca_v3_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

cd /home/shli842i/p_dna15

# global variables - check these
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks
export WORKDIR=/home/shli842i/p_dna15
mkdir -p scripts/logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

# global variables - edit per run
export VERSION=v3.1

echo "Starting run $VERSION at $(date)"

Rscript code/diff-meth-pca_v3.R

echo "Finished at $(date)"
