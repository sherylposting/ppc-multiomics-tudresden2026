#!/bin/bash
#SBATCH --job-name=install_methylkit
#SBATCH --output=install_methylkit_%j.out
#SBATCH --error=install_methylkit_%j.err
#SBATCH --time=02:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module purge
module load release/23.04  GCC/11.3.0  OpenMPI/4.1.4 R/4.2.1

export R_LIBS_USER=/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1
mkdir -p "$R_LIBS_USER"

# Install BiocManager if needed
Rscript -e '
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager", repos="https://cloud.r-project.org")
'

# Install methylKit and dependencies
R -e '
BiocManager::install(
    "methylKit",
    ask = FALSE,
    lib = "/data/horse/ws/shli842i-p_dna15_1/rpacks_4.2.1",
    update = TRUE
)
'
