# WGS pipelines for T-LOC and TDNA-scan

## First set-up HPC access
Don't bother trying to do this on your local especially if we you have multiple samples.
Let your PI know that you are requesting access (i.e. AskICT on latrobe intranet)
Also set-up
> 1. PuTty (HPC terminal access)
> 2. WinSCP (File transfer and browsing with UX interactivity to HPC and any other remote access location


## Raw files
Ok, so what do we actually get from Novogene?
1. Raw files (zipped folder)
> Example
> 1. S1_DKDN250032122-1A_2373VCLT4_L1_1.fq
> 2. S1_DKDN250032122-1A_2373VCLT4_L1_2.fq
FASTQC doesnt care it's zipped, so just port it into HPC folder to begin with
2. QC report
> Example
> X201SC25083173-Z01-F001_Report.html in a folder like e.g. 02.Report_X201SC25083173-Z01-F001 (Folder)

The QC report is useful, but we would also like to generate some other QC metrics.
The QC report also has adapter information:
> Sequences of adapter
> P5 adapter：
> P5→P7’(5’→3’)
> AATGATACGGCGACCACCGAGATCTACAC[i5*]ACACTCTTTCCCTACACGACGCTCTTCCGATCT

> P7 adapter：
> P5→P7’(5’→3’)
> GATCGGAAGAGCACACGTCTGAACTCCAGTCAC[i7*]ATCTCGTATGCCGTCTTCTGCTTG

> *i5/i7 represent index sequences. They are provided here to show the complete adapter structure. Requirements for including or excluding indices in adapter trimming software may vary. Please refer to the specific documentation of the relevant software.



## QC report
Technically Novogene has already provided this, but want to also look at other quality metrics i.e. adapter content
Run FastQC
> batch_fastq.sh

You will get:
> 1. fastqc PER sample PER paired-end file.
> 2. multiqc html file that combines all the fastqc reports into just one html file

Open up your multiqc files (multiqc_report.html) and inspect generated output.

## Trimming
You can do this on Geneious, but (MUCH)faster to just run it on HPC. We're talking couple of minutes vs. an hour on your regular University-issue machine.
Will be using BBDukTrim, but feel free to use whatever pleases you e.g. Trimmomatic etc. 
