#!/bin/bash
#SBATCH --job-name=diff-meth_v2.0
#SBATCH --output=logs/diff-meth_v2.0_%j.out
#SBATCH --error=logs/diff-meth_v2.0_%j.err
#SBATCH --time=6:00:00
#SBATCH --mem=400G
#SBATCH --cpus-per-task=8

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

Rscript code/diff-meth_v2.0.R

echo "Finished at $(date)"
