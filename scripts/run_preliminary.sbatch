#!/bin/bash
#SBATCH --job-name=preliminary_v3
#SBATCH --output=logs/preliminary_v3_%j.out
#SBATCH --error=logs/preliminary_v3_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=256G
#SBATCH --cpus-per-task=1

module purge
module load release/25.06 GCC/13.3.0 R

cd /home/shli842i/p_dna15

# global variables - check these
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks
export WORKDIR=/home/shli842i/p_dna15
mkdir -p scripts/logs /data/horse/ws/shli842i-p_dna15_1/tmp
export TMPDIR=/data/horse/ws/shli842i-p_dna15_1/tmp

# global variables - edit per run
export VERSION=v3.2
export FILT_LO_COUNT=3
export OMITTED=1

echo "Starting run $VERSION at $(date) with $FILT_LO_COUNT x cutoff and sample $OMITTED omitted"

Rscript code/preliminary_v3.R

echo "Finished at $(date)"
