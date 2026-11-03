#!/bin/bash
#SBATCH --job-name=v2.2_DMRichR
#SBATCH --output=logs/v2.2_DMRichR_%j.out
#SBATCH --error=logs/v2.2_DMRichR_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1

# loading modules to run R
# note - you may need to check module avail R/4.2.1 to ensure that these prerequisite modules are still correct
module purge
module load release/23.04
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load R/4.2.1

# global variables - check these
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR2
export VERSION=v2.2 # arbitrary, for your own record-keeping
export WORKDIR=/home/shli842i/p_dna15/DMRichR/$VERSION
export SAMPLEINFO_PATH=/home/shli842i/p_dna15/DMRichR/sample_info.xlsx
export DATADIR=/projects/p_dna15/data/EM_seq_files/CpG_reports
export RESPATH=/projects/p_dna15/data

echo "Starting run $VERSION at $(date)"

export R_LIBS_USER=$LIBPATH # not sure if necessary, maybe needed by R

Rscript /home/shli842i/p_dna15/DMRichR/code/run_DMR.R

echo "Finished at $(date)"
