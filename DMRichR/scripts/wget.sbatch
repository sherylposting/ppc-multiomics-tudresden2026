#!/bin/bash
#SBATCH --job-name=wget
#SBATCH --output=logs/wget_%j.out
#SBATCH --error=logs/wget_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

cd /home/shli842i/p_dna15

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz
gunzip mm39.fa.gz

module load release/24.04  GCC/13.2.0  OpenMPI/4.1.6 SAMtools/1.19.2

samtools faidx mm39.fa
