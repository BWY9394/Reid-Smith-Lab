#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=8192
#SBATCH --job-name="BBduk_all"
#SBATCH --partition=day

module load BBMap/38.90-GCC-10.2.0

mkdir=TrimmedTestSepFold

RAW_ROOT="RawTestSepFold/"
OUT_ROOT="TrimmedTestSepFold/"

for SAMPLE_DIR in ${RAW_ROOT}/*/; do
    SAMPLE=$(basename "${SAMPLE_DIR}")
    OUT_DIR="${OUT_ROOT}/${SAMPLE}"
    mkdir -p "${OUT_DIR}"

    for R1 in ${SAMPLE_DIR}/*_1.fq.gz; do
        R2="${R1/_1.fq.gz/_2.fq.gz}"

        if [[ ! -f "${R2}" ]]; then
            echo "WARNING: Missing R2 for ${R1}, skipping"
            continue
        fi

        BASE=$(basename "${R1}" _1.fq.gz)

        echo "Processing ${SAMPLE} / ${BASE}"

        bbduk.sh \
          in1="${R1}" in2="${R2}" \
          out1="${OUT_DIR}/${BASE}_1_trimmed.fq.gz" \
          out2="${OUT_DIR}/${BASE}_2_trimmed.fq.gz" \
          ref=adapters.fa ktrim=r k=23 mink=11 hdist=1 \
          qtrim=rl trimq=10 \
          minlen=50 \
          threads=${SLURM_CPUS_PER_TASK}
    done
done
