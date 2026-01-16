#!/bin/bash
#SBATCH --job-name=orthofinder_51
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH -t 5-00:00 # time (D-HH:MM)
#SBATCH --partition=week


# Load required modules
module purge
module load OrthoFinder/2.5.4-foss-2020b

# Set number of threads
THREADS=32

# Input directory containing protein FASTA files (51 genomes)
INPUT_DIR="./panfaba/Longest_AA"

# Output directory (OrthoFinder will create this)
OUTPUT_DIR="./panfaba/Longest_AA/orthofinder_results" #having a folder here already created with the folder name will throw an error

# Run OrthoFinder with DIAMOND (much faster for large datasets)
orthofinder -f ${INPUT_DIR} \
            -t ${THREADS} \
            -a ${THREADS} \
            -S diamond \
            -o ${OUTPUT_DIR}

echo "OrthoFinder analysis complete!"
echo "Results in: ${OUTPUT_DIR}"
