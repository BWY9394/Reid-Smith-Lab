#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --job-name="salidx"
#SBATCH --partition=day

#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b

#Reference for pre-processings steps
#https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/

#pre-processing steps that you can just run in HPC terminal. Just contain here for easy reference
#raw files
#1. Genome:
#Vunguiculata_540_v1.0.softmasked.fa.gz
#2  Transcripts (All splice variants),UTRs and exons:
#Vunguiculata_540_v1.2.transcript.fa.gz

#Salmon indexing requires the names of the genome targets, so that with:
#grep "^>" <(gunzip -c ./Genomics/Vunguiculata_540_v1.0.softmasked.fa.gz) | cut -d " " -f 1 > decoysVu.txt
#sed -i.bak -e 's/>//g' decoysVu.txt

#Along with the list of decoys salmon also needs the concatenated transcriptome and genome reference file for index. NOTE: the genome targets (decoys) should come after the transcriptome targets in the reference
#cat ./Genomics/Vunguiculata_540_v1.2.transcript.fa.gz ./Genomics/Vunguiculata_540_v1.0.softmasked.fa.gz > Vunguiculata_540_v1.0.gentrome.fa.gz

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Index
#-t = transcripts (your cDNA file) #-i = index directory   --gencode
#y ou don't need the gencode:  zcat ./Genomics/Vunguiculata_540_v1.2.transcript.fa.gz | head -n 5  will display entry of your gene. If there is no "|" pipe, you're good e.g. >ENST00000335137.4|ENSG000001234| means need gencode
#output
salmon index -t Vunguiculata_540_v1.0.gentrome.fa.gz -d decoysVu.txt -p 12 -i salmon_index
