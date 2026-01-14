# WGS pipelines for T-LOC and TDNA-scan

## First set-up HPC access
Don't bother trying to do this on your local especially if we you have multiple samples.
Let your PI know that you are requesting access (i.e. AskICT on latrobe intranet)
Also set-up
> 1. PuTTy (HPC terminal access)
> 2. WinSCP (File transfer and browsing with UX interactivity to HPC and any other remote access location)


## Raw files
Ok, so what do we actually get from Novogene?
1. Raw files (zipped folder)
> Example
> 1. S1_DKDN250032122-1A_2373VCLT4_L1_1.fq
> 2. S1_DKDN250032122-1A_2373VCLT4_L1_2.fq

FASTQC doesnt care it's zipped, so just port it into HPC folder to begin with

2. QC report
> Example
> X201SC25083173-Z01-F001_Report.html in a folder like "02.Report_X201SC25083173-Z01-F001"

The QC report is useful, but we would also like to generate some other QC metrics.

The QC report also has adapter information, which will be important for trimming later:
> Sequences of adapter
> P5 adapter：
> P5→P7’(5’→3’)
> AATGATACGGCGACCACCGAGATCTACAC[i5*]ACACTCTTTCCCTACACGACGCTCTTCCGATCT

> P7 adapter：
> P5→P7’(5’→3’)
> GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7*]ATCTCGTATGCCGTCTTCTGCTTG

> *i5/i7 represent index sequences. They are provided here to show the complete adapter structure. Requirements for including or excluding indices in adapter trimming software may vary. Please refer to the specific documentation of the relevant software.



## QC report
Technically Novogene has already provided this, but may want to also look at other quality metrics i.e. adapter content
Run FastQC in the folder where your files are
Script required:
> batch_fastq.sh

The command(s) would be:
> cd /data/group/medaglab/project/BWee/hemp_RNAseq/01-fastq
>
> sbatch /data/group/medaglab/project/BWee/scripts/batch_fastqc.sh

Notice that my bash script is in a different directory- if you have decided to put the bash script in the same working directory you have set initially, feel free, just omit the file directories preceding batch_fastqc.sh


You will get:
> 1. fastqc PER sample PER paired-end file.
> 2. multiqc html file that combines all the fastqc reports into just one html file

Open up your multiqc files (multiqc_report.html) and inspect generated output.

## Trimming
You can do this on Geneious, but (MUCH)faster to just run it on HPC. We're talking couple of minutes vs. an hour on your regular University-issue machine.
Will be using BBDukTrim, but feel free to use whatever pleases you e.g. Trimmomatic etc. 
A *bit* of debate as to whether is necessary, but meh, seems to be gold standard so will do so until have someone says otherwise

Script required:
> BBduk_single.sh
>
> BBduk_single.sh for single samples at one time
>
> BBduk_batch.sh for batch processing if all in one folder
>
> BBduk_sep_fold.sh if in subfolders

Next, navigate to where you have stored the folder containg yourr raw files. In my case its in a subfolder inside Trimming. Also create a folder called "Trimmed".

The command(s) would be :
> cd /data/group/medaglab/project/BWee/Trimming
> 
> mkdir -p Trimmed #makes a folder called Trimmed

Now, if open up BBduk-single.sh and edit in your file names for raw reads to be trimmed
> RAW_R1="Raw/S16/S16_DKDN250032137-1A_2373VCLT4_L3_1.fq.gz"   # forward (R1)
>
> RAW_R2="Raw/S16/S16_DKDN250032137-1A_2373VCLT4_L3_2.fq.gz"   # reverse (R2)
>
> BASE="S16_DKDN250032137-1A_2373VCLT4_L3" #just the file name, no _1.fq.gz or _2.fq.gz

Then run the script:
> sbatch /data/group/medaglab/project/BWee/scripts/BBduk_single.sh  

You will get:
Trimmed reads in your designated output folder under "Trimmed"

Once familiar, switch over to either the BBduk_batch.sh or BBduk_sep_fold.sh scripts for bulk sample processing, just do:
> sbatch /data/group/medaglab/project/BWee/scripts/BBduk_batch.sh #Or BBduk_batch_sep_fold.sh if your samples are in subfolders in the Trimming subfolder

No need to edit files names as it just scans the folder containing your raw reads.

## Running T-LOC
What is T-LOC

## Running TDNA-scan
What is TDNA-scan

