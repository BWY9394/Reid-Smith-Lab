#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name="BBdukS16"
#SBATCH --partition=day

module load BBMap/38.90-GCC-10.2.0   # Load BBTools (includes bbduk.sh)
# Input files
RAW_R1="Raw/S16/S16_DKDN250032137-1A_2373VCLT4_L3_1.fq.gz"   # forward (R1)
RAW_R2="Raw/S16/S16_DKDN250032137-1A_2373VCLT4_L3_2.fq.gz"   # reverse (R2)

# Output files
OUT_DIR="Trimmed"
BASE="S16_DKDN250032137-1A_2373VCLT4_L3"
mkdir -p "${OUT_DIR}"

bbduk.sh \
  in1="${RAW_R1}" in2="${RAW_R2}" \
  out1="${OUT_DIR}/${BASE}_1_trimmed.fq" \
  out2="${OUT_DIR}/${BASE}_2_trimmed.fq" \
  ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 \
  qtrim=rl trimq=10 tpe=t \
  minlen=50 \

