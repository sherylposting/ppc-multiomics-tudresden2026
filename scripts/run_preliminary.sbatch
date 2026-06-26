#!/bin/bash
#SBATCH --job-name=preliminary2.1-withnorm-plotboth
#SBATCH --output=logs/preliminary2.1-withnorm-plotboth_%j.out
#SBATCH --error=logs/preliminary2.1-withnorm-plotboth_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=1

module purge
module load release/25.06 GCC/13.3.0 R

cd /home/shli842i/p_dna15

mkdir -p scripts/logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

echo "TMPDIR=$TMPDIR"
echo "Starting at $(date)"

Rscript code/preliminary_v2.1.R

echo "Finished at $(date)"
