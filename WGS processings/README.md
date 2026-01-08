# WGS pipelines for T-LOC and TDNA-scan

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

## Trimming
You can do this on Geneious, but (MUCH)faster to just run it on HPC. We're talking couple of minutes vs. an hour on your regular University-issue machine.
