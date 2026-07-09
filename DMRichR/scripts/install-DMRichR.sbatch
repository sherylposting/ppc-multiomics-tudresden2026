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
module load HDF5/1.14.3
module load R/4.2.1
module list

export R_LIBS_USER=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1

R