#!/bin/bash
#SBATCH --job-name=diffmeth-methraw
#SBATCH --output=logs/diffmeth-methraw_%j.out
#SBATCH --error=logs/diffmeth-methraw_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=2

module purge
module load release/25.06 GCC/13.3.0 R

cd /home/shli842i/p_dna15

mkdir -p scripts/logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

echo "TMPDIR=$TMPDIR"
echo "Starting at $(date)"

Rscript code/diff-meth-volcanoplot_v1.0.R

echo "Finished at $(date)"
