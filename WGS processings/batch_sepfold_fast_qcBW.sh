#!/bin/bash
#SBATCH --time=23:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=12G
#SBATCH --job-name="fastqc"
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au # send-to address

module load FastQC/0.11.9-Java-11
module load Python/3.10.4-GCCcore-11.3.0  
module load MultiQC/1.12-foss-2021b  

mkdir -p FastQC_out/

for sample in */*.fq.gz
do
  echo "$sample"
  
  describer=$(basename "$sample" .fq.gz)
  echo "$describer"
  
  fastqc -t 24 "$sample" --outdir=./FastQC_out/
done

cd FastQC_out
multiqc .
