#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name="C1_24NCBI"
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au # send-to address

#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b
#/data/group/medaglab/project/BWee/Salmon/Raw/FUN_TC_Vu/

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Quantifying reads via Salmon
# -- Configure these variables only ------------------------------------------
SAMPLES_DIR="Raw/FUN_TC_Vu" #Won't work if you don't have the right directory
SAMPLE_PATTERN="C{1..24}"  #Won't work if your sample pattern is wrong
SALMON_INDEX="salmon_indexVuNCBI" #Won't work if you don't call the same/exisiting indexed transcriptome+genome created before/relevant for your study
OUTPUT_DIR="quantsC1_C24_NCBI"  #Whatever you like
SUMMARY_FILE="salmon_read_summary_C1_C24_NCBI.tsv" #Whatever you like
THREADS=8 #scale up if needed, e.g. to 12, but also scale up ram.
# ----------------------------------------------------------------------------

module load Salmon/1.4.0-GCC-11.2.0

for fn in "${SAMPLES_DIR}"/${SAMPLE_PATTERN}; do
  samp=$(basename "${fn}")
  echo "Processing sample ${samp}"
  salmon quant -i "${SALMON_INDEX}" -l A \
    -1 "${fn}/${samp}_1.fq.gz" \
    -2 "${fn}/${samp}_2.fq.gz" \
    -p "${THREADS}" --validateMappings -o "${OUTPUT_DIR}/${samp}_quant"
done

echo -e "Sample\tNumReads\tNumMapped\tPercentMapped" > "${OUTPUT_DIR}/${SUMMARY_FILE}"

for d in "${OUTPUT_DIR}"/*_quant; do
  samp=$(basename "$d" _quant)
  numReads=$(jq  '.num_processed'   "$d/aux_info/meta_info.json")
  numMapped=$(jq '.num_mapped'      "$d/aux_info/meta_info.json")
  pctMapped=$(jq '.percent_mapped'  "$d/aux_info/meta_info.json")
  echo -e "${samp}\t${numReads}\t${numMapped}\t${pctMapped}" >> "${OUTPUT_DIR}/${SUMMARY_FILE}"
done