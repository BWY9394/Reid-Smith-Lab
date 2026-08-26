# OrthoFinder walkthrough

Based off https://davidemms.github.io/menu/tutorials.html

Adapted for LTU systems and Reid lab infrastructure

### Running OrthoFinder

2 options: 
1) Run on HPC **
2) Run on Reid lab HPC
3) Run on your own computer *

\* For Linux: follow the instructions below.
\* For Macs: it’s easiest to use Bioconda to install and then follow the instructions below.
\* For Windows it is best to install the Windows Subsystem for Linux and then continue as for Linux below.
\*\* Would recommend this option if you are running a huge number of species and/or the species have very large genomes e.g. Faba which has a genome of 13Gb.

You will also need python.

### 1. Setting it up(Skip to Step 2 if picking option 1 or 2) ###


### 2. Running an analysis

OrthoFinder uses the high-performance protein sequence aligner DIAMOND by default to perform all-versus-all sequence similarity searches, providing a much faster alternative to BLAST for large-scale phylogenomic analyses. These sequence similarity results are used to infer orthogroups, which are groups of genes descended from a single ancestral gene in the last common ancestor of the species being compared. OrthoFinder subsequently uses phylogenetic analysis to infer orthologous and paralogous relationships within these groups.


What it requires:
1. Proteomes formatted as either .fa”, “.faa”, “.fasta”, “.fas” or “.pep”., that only contain the primary transcript.

It is often the case that you will need to generate a fasta file containing just the primary transcript of each protein. In this case(and most), we define this as the longest representative transcript.

To do this we run:
```bash
for f in *fa ; do python ~/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done
```
This command processes every .fa file in the current directory and runs OrthoFinder's primary_transcript.py script on each one.
Syntax:
for f in *fa — finds files in the current directory whose names end in fa and assigns each filename to the variable f, one at a time.

The output is saved in a new primary_transcripts/ folder inside the same directory as each input FASTA file. The OrthoFinder script explicitly creates that folder automatically.
```
/home/user/proteomes/
├── species1.fa
├── species2.fa
├── species3.fa
└── primary_transcripts/
    ├── species1.fa
    ├── species2.fa
    └── species3.fa
```
