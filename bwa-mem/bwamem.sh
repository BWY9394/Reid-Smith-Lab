#!/bin/bash
#SBATCH --job-name=bwaS15 #edit job name to whatever works for you. 6-10 char max recommended
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --partition=day

# =====================
# CREATE OUTPUT FOLDER
# =====================
mkdir -p results #edit this to whatever makes sense for you.

# =====================
# INPUTS
# =====================
REF=WGS/VuGenome.fasta
R1=WGS/S15_DKDN250032136-1A_2373VCLT4_L3_1_trimmed.fq
R2=WGS/S15_DKDN250032136-1A_2373VCLT4_L3_2_trimmed.fq

OUT=results/sample15 #So your samples will called sample15.sorted.bam or sample15.idx.stats It's the prefix to your output files OK????
THREADS=8

# =====================
# LOAD MODULES
# =====================
module load BWA/0.7.17-GCCcore-11.2.0
module load SAMtools/1.16.1-GCC-11.3.0

# =====================
# INDEX REFERENCE
# =====================
bwa index $REF
samtools faidx $REF

# =====================
# MAP TO REFERENCE
# =====================
bwa mem -t $THREADS -M $REF $R1 $R2 \
| samtools view -b -q 10 \
| samtools sort -@ $THREADS -o ${OUT}.sorted.bam

# =====================
# INDEX BAM
# =====================
samtools index ${OUT}.sorted.bam

# =====================
# QC / MAPPING REPORT
# =====================
samtools flagstat ${OUT}.sorted.bam > ${OUT}.flagstat.txt
samtools idxstats ${OUT}.sorted.bam > ${OUT}.idxstats.txt

echo "Illumina mapping complete"
