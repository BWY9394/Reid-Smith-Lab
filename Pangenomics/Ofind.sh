#!/bin/bash
#SBATCH --job-name=orthofinder_51
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH -t 7-00:00 # time (D-HH:MM)
#SBATCH --partition=week


# Load required modules
module load OrthoFinder/2.5.4-foss-2020b

# Input directory containing protein FASTA files (51 genomes)
INPUT_DIR="panfaba/Longest_AA/"

# Run OrthoFinder with DIAMOND (much faster for large datasets)
orthofinder -f ${INPUT_DIR}


