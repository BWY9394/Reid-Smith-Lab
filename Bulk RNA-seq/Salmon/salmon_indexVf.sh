#!/bin/bash
#SBATCH --time=1-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --job-name="indAmbD"
#SBATCH --partition=day
#send-to address
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au 


#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b

#Reference for pre-processings steps
#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

#pre-processing steps that you can just run in HPC terminal. Just contain here for easy reference
#raw files
#1. Genome:
#Vfaba_824_v1.0.softmasked.fa.gz
#250504_PBA_Amberley_pseudomolecules_v1_chr1_parts+unanchored_contigs.fasta.gz
#2  Transcripts (All splice variants),UTRs and exons:
#Vfaba_824_v1.1.transcript.fa.gz
#VFABA.AMBERLEY.pgsb.r1.Oct2025.transcripts.fa

#Salmon indexing requires the names of the genome targets, so that with:
grep "^>" <(gunzip -c ./Genomics/250504_PBA_Amberley_pseudomolecules_v1_chr1_parts+unanchored_contigs.fasta.gz) | cut -d " " -f 1 > decoysVf_Amberly_HCplusLC.txt
sed -i.bak -e 's/>//g' decoysVf_Amberly_HCplusLC.txt

#Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
cat \
  ./Genomics/VFABA.AMBERLEY.pgsb.r1.Oct2025.salmon.transcriptome.fa \
  <(gzip -dc "./Genomics/250504_PBA_Amberley_pseudomolecules_v1_chr1_parts+unanchored_contigs.fasta.gz") \
  | gzip > PBA_Amberley_v1_HC_LC.gentrome.fa.gz


#module load
module load Salmon/1.4.0-GCC-11.2.0

#Index
#-t = transcripts (your cDNA file) #-i = index directory   --gencode
#you don't need the gencode:  zcat ./Genomics/Vfaba_824_v1.1.transcript.fa.gz | head -n 5  will display entry of your gene. If there is no "|" pipe, you're good e.g. >ENST00000335137.4|ENSG000001234| means need gencode
#output
salmon index -t PBA_Amberley_v1_HC_LC.gentrome.fa.gz -d decoysVf_Amberly_HCplusLC.txt -p 12 -i salmon_indexVf_Amberly_HCplusLC -k 31
