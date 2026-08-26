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

### 1. (Skip to Step 2 if picking option 1 or 2)


### 2. Running an analysis

OrthoFinder uses the high-performance protein sequence aligner DIAMOND by default to perform all-versus-all sequence similarity searches, providing a much faster alternative to BLAST for large-scale phylogenomic analyses. These sequence similarity results are used to infer orthogroups, which are groups of genes descended from a single ancestral gene in the last common ancestor of the species being compared. OrthoFinder subsequently uses phylogenetic analysis to infer orthologous and paralogous relationships within these groups.
