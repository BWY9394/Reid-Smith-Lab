# Quick walkthrough for bwa-mem files

#### _Disclaimer_
Currently, have only written script for single file bwa-mem output. Will tackle batch-processing at later point in time i.e indexing just once so that when running a bunch of alignments to reference in parallel for different samples, no chance of errors being thrown.

## Necessary files
Your raw fastq files:

Single:
> yourfiles_L1.fq.gz
>

Paired:
> yourfiles_L1_1.fq.gz
> 
> yourfiles_L1_2.fq.gz

Script:
> bwamem.sh

## To run

First change directories in HPC to wherever your raw files are. In my case my raw files are in a subfolder that I will call directly later
>cd  /data/group/medaglab/project/BWee/TDNAscan/TDNAscan

Open up bwamem.sh inside WinSCP's native text editor. Under inputs, path and file name directory to as you have set up in your own HPC folder.

Example:
> REF=WGS/VuGenome.fasta
> 
> R1=WGS/S15_DKDN250032136-1A_2373VCLT4_L3_1_trimmed.fq
> 
> R2=WGS/S15_DKDN250032136-1A_2373VCLT4_L3_2_trimmed.fq
> 
> #The WGS/ suffix indicates that my samples are in a subfolder inside the current working directory.

Now run the script.

>sbatch  /data/group/medaglab/project/BWee/scripts/bwamem.sh

Notice that my bash script is in a different directory- if you have decided to put the bash script in the same working directory you have set initially, feel free, just omit the file directories preceding bwamem.sh

The script currently saves the output in a new folder. This new folder is by default set/named to results, but you can edit it under

> mkdir -p results15 #edit this to whatever makes sense for you.

## Run-time
For a genome of 519 Mb, mapping of paired end (PE) WGS reads approximating 20x coverage,  150bp read length, 6 GB of raw data per sample, 78 million reads per sample (~36 per "end" of PE), took approximately 4 hours with 8 cores/threads, 32GB ram assigned on HPC. Doing this in Geneious would take 3-4 days per sample conmputer with 16GB ram and 6 cores (12 threads).

## Viewing your data

The output files are:
> sample1.sorted.bam.bai
>
> sample1.sorted.bam
>
> sample1.flagstat.txt
>
> sample1.idxstats.txt
>

Download these to a folder on your PC with Genenious/IGV.

Transfer these to a folder of your choosing in Geneious. It will ask for a reference, which just means where is your reference genome in Geneious, so just select where this is in Geneious.




