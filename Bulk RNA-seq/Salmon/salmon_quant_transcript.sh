#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64
#SBATCH --job-name="salC1|24"
#SBATCH --partition=day


#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Run salmon quant

for fn in Raw/C{1..24}; do
  samp=$(basename "${fn}")
  echo "Processing sample ${samp}"

  # -l A: salmon auto-detects library type (stranded/unstranded)
  # -1: forward reads
  # -2: reverse reads (ignore if single-end)
  # -p 8: use 8 CPU threads
  # --validateMappings: selective alignment mode
  # -o: output directory for quant results

  salmon quant -i salmon_indexVu -l A \
    -1 "${fn}/${samp}_1.fq.gz" \
    -2 "${fn}/${samp}_2.fq.gz" \
    -p 32 --validateMappings -o "quantsC1_C24/${samp}_quant"
done
