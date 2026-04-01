#!/bin/bash
#SBATCH --time=5-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=256G
#SBATCH --job-name="salidxVf"
#SBATCH --partition=week
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au # send-to address


#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b

#Reference for pre-processings steps
#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

#pre-processing steps that you can just run in HPC terminal. Just contain here for easy reference
#raw files
#1. Genome:
#Vfaba_824_v1.0.softmasked.fa.gz
#2  Transcripts (All splice variants),UTRs and exons:
#Vfaba_824_v1.1.transcript.fa.gz

#Salmon indexing requires the names of the genome targets, so that with:
#grep "^>" <(gunzip -c ./Genomics/Vfaba_824_v1.0.softmasked.fa.gz) | cut -d " " -f 1 > decoysVf.txt
#sed -i.bak -e 's/>//g' decoysVf.txt

#Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
#cat ./Genomics/Vfaba_824_v1.1.transcript.fa.gz \    ./Genomics/Vfaba_824_v1.0.softmasked.fa.gz \    > Vfaba_824_v1.1.gentrome.fa.gz
#zcat Vfaba_824_v1.1.gentrome.fa.gz | head -n 5

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Index
#-t = transcripts (your cDNA file) #-i = index directory   --gencode
#you don't need the gencode:  zcat ./Genomics/Vfaba_824_v1.1.transcript.fa.gz | head -n 5  will display entry of your gene. If there is no "|" pipe, you're good e.g. >ENST00000335137.4|ENSG000001234| means need gencode
#output
salmon index -t Vfaba_824_v1.1.gentrome.fa.gz -d decoysVf.txt -p 12 -i salmon_indexVf -k 31