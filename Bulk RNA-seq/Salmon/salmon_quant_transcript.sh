#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name="salZBF1to28"
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au # send-to address

#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b
#/data/group/medaglab/project/BWee/Salmon/Raw/FUN_TC_Vu/

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Run salmon quant

for fn in Raw/ZBF/ZBF{1..28}; do
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
    -p 8 --validateMappings -o "quantsZBF1_ZBF28/${samp}_quant"
done


# Now make an output file
echo -e "Sample\tNumReads\tNumMapped\tPercentMapped" > quantsZBF1_ZBF28/salmon_read_summary_ZBF1_ZBF28.tsv



#Extract mapping stats out into tab-separated table 
for d in quantsZBF1_ZBF28/*_quant; do
   samp=$(basename "$d" _quant)
   numReads=$(jq '.num_processed' "$d/aux_info/meta_info.json")
   numMapped=$(jq '.num_mapped' "$d/aux_info/meta_info.json")
   pctMapped=$(jq '.percent_mapped' "$d/aux_info/meta_info.json")

   echo -e "${samp}\t${numReads}\t${numMapped}\t${pctMapped}" \
     >> quantsZBF1_ZBF28/salmon_read_summary_ZBF1_ZBF28.tsv
done


# Optional: email summary to yourself
mail -s "Salmon Quant Finished" -a quantsZBF1_ZBF2/salmon_read_summary_ZBF1_ZBF28.tsv BWeeYang@ltu.edu.au <<< "All samples finished. See attached summary."

