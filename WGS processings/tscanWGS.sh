#!/bin/bash
#SBATCH -t 1-00:00 # time (D-HH:MM)
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --job-name="TSNS04"
#SBATCH --partition=week

module load Python/3.9.6-GCCcore-11.2.0
module load BWA/0.7.17-GCCcore-11.2.0
module load SAMtools/1.15.1-GCC-11.2.0

python3 tdnascan2.py -1 WGS/S4_DKDN250032125-1A_2373VCLT4_1_trimmed.fq -2 WGS/S4_DKDN250032125-1A_2373VCLT4_2_trimmed.fq -t WGS/L2_FUN_exon_guides_vector_with_RUBY_wrapped.fasta -g WGS/VuGenome.fasta -p WGS_S04dummy
