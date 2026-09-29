#!/bin/bash
#SBATCH --job-name=Ofind6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 5-00:00:00 # time (D-HH:MM:SS)
#SBATCH --partition=week

#check that you are in the directory just before your primary_transcripts folder where your longest primary protein transcripts

# Load required modules
module load OrthoFinder/2.5.4-foss-2020b

# Input directory containing protein FASTA files (51 genomes)
INPUT_DIR="./8splegumes/primary_transcripts/" 

# Run OrthoFinder with DIAMOND (much faster for large datasets)
orthofinder -f ${INPUT_DIR} -t 32 -a 4

#-t : dictates no. total for DIAMOND/search-related work
#-a : number of parallel threads used for OrthoFinder’s internal analysis tasks, which are often RAM-intensive.

