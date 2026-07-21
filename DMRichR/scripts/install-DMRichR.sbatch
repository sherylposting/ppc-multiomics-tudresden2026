#!/bin/bash
#SBATCH --job-name=install_DMRichR
#SBATCH --output=logs/install_DMRichR_%j.out
#SBATCH --error=logs/install_DMRichR_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

module purge
module load release/23.04
module load GCC/11.3.0
module load OpenMPI/4.1.4
module load HDF5/1.13.2
module load R/4.2.1
module list

export CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2023-04-25 # approx the last time DMRichR was updated
export LIBPATH=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1_DMRichR

# set path to the custom DMRichR package. important!
export DMRICHR_PACKAGE_PATH=/home/shli842i/DMRichR

# important to load it --vanilla! something about the Rprofile maybe? causes the incorrect system library to be loaded
Rscript --vanilla install-DMRichR.R