#!/bin/bash
#SBATCH --time=1-3:00:00
#SBATCH --nodes=1
#SBATCH --job-name=squant26
#SBATCH --partition=week
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH --output=STAR_aln26.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.weeyang@latrobe.edu.au

module load STAR/2.7.10b-GCC-11.3.0


# -- Configure these variables --------------------
SAMPLES_DIR="Raw/060826_cowpea_xx_3_5mM_N"
STAR_INDEX="STAR_index"
OUTPUT_DIR="STAR_xx_2026"
THREADS=12
# -------------------------------------------------
mkdir -p "$OUTPUT_DIR"

echo "Only process samples that actually exist"

shopt -s nullglob
files=( "$SAMPLES_DIR"/*.fastq.gz )

echo "Found ${#files[@]} FASTQ files"

for R1 in "${files[@]}"; do

    filename=$(basename "$R1")
    samp="${filename%.fastq.gz}"

    echo "Aligning ${samp}"
    echo "Input: ${R1}"

    STAR \
        --runThreadN "$THREADS" \
        --genomeDir "$STAR_INDEX" \
        --readFilesIn "$R1"  \
        --readFilesCommand zcat \
        --twopassMode Basic \
        --outFileNamePrefix "${OUTPUT_DIR}/${samp}_" \
        --outSAMtype BAM SortedByCoordinate

done

echo "Finished STAR alignment for each file in main directory."
