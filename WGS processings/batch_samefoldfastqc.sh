#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=4096
#SBATCH --job-name="fastqc"
#SBATCH --partition=day


module load FastQC/0.11.9-Java-11
module load Python/3.10.4-GCCcore-11.3.0  
module load MultiQC/1.12-foss-2021b  

mkdir FastQC_out/ #ensure you have a directory as precedes


for sample in *.fastq.gz 
 do  
   echo $sample  
   describer=$(echo ${sample} | sed 's/.fastq.gz//')  
   echo $describer  
   
   
   fastqc -t 24 $describer.fastq.gz --outdir=./FastQC_out/

 done 
 
cd FastQC_out
multiqc . #that's really all there is to it.
 
