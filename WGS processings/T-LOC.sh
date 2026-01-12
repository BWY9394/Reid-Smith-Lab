#!/bin/bash
#SBATCH -t 5-00:00 # time (D-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --job-name="TLC17PY3"
#SBATCH --partition=week

# Load required modules
#python2
#module load BWA/0.7.17-GCC-10.2.0
#module load SAMtools/1.11-GCC-10.2.0
#module load R/4.2.1-foss-2022a
#module load Python/2.7.18-GCCcore-10.2.0

#python3
module load BWA/0.7.17-GCC-10.2.0
module load SAMtools/1.11-GCC-10.2.0
module load R/4.2.1-foss-2022a
module load Python/3.8.6-GCCcore-10.2.0
module load Pysam/0.17.0-GCC-11.2.0



# Define paths to your files
FASTQ1="/data/group/medaglab/project/BWee/T-LOC/Samples/S17_DKDN250032138-1A_2373VCLT4_L3_1.fq"
FASTQ2="/data/group/medaglab/project/BWee/T-LOC/Samples/S17_DKDN250032138-1A_2373VCLT4_L3_2.fq"
GENOME="/data/group/medaglab/project/BWee/T-LOC/Samples/VuGenome.fasta"
TDNA="/data/group/medaglab/project/BWee/T-LOC/Samples/L2_FUN_exon_guides_vector_with_RUBY_wrapped.fasta"

# Extract sample ID (e.g., "S1") from FASTQ1 filename
SAMPLE=$(basename "$FASTQ1" | cut -d'_' -f1)

# Define output folder using the sample ID
OUTPUT="/data/group/medaglab/project/BWee/T-LOC/Samples/${SAMPLE}_output_morecores_rawreadscontrol"

# Run T-LOC
python /data/group/medaglab/project/BWee/T-LOC/T-LOC_py3.py \
    --fastq ${FASTQ1},${FASTQ2} \
    --genome $GENOME \
    --TDNA $TDNA \
    --output $OUTPUT
