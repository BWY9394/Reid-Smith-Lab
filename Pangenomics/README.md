# Pangenomics 101

What is pangenomics? Essentially, think of it as a means of representing all of the genomic variation found in a group of organisms, in contrast to a traditional reference that is (usually) assembled from a single individual.

## How do we display pangenomics infomation?

A great starting link: https://doi.org/10.1093/bib/bbae588

## First work out what kind of pangenome..
...that you are working on/want to generate.

1. Presence-abscence variation (PAV) pangenome: composed of a ‘core genome’ and an ‘accessory genome’. The core genome is the full set of genes that are present in every member of the population, while the accessory genome is composed of genes present in a subset of the population.
- focuses on gene presence and absence, does not account for gene location, allelic diversity, or intergenic sequence
2. Representative sequence pangenome: composed of carefully selected genomic sequences so that as much genomic variation from the population as possible is represented using as little sequence as possible.
  - same structure as traditional ref genome, but with aditional contigs containing supplementary genomic sequences.
3. Pangenome Graph/Graphical Pangenome: May be sequence or gene-oriented.
- Sequence-oriented pangenome graph models genomic sequence variation as well as its location relative to other genomic sequences from the population
- Gene-oriented pangenome graph models the genes found within the population and their order relative to other members of the population


## Constructing a pangenome

Two methods:
1. Homolog-based strategy: Widely used for bacteria pangenomes
2. Map-to-pan strategy

_Will finish later_

## How to analyze abscence/presence with OrthoFinder
https://www.cd-genomics.com/resource-pan-genome-analysis-bioinformatics.html#guide4
Useful resource

# Using OrthoFinder for Presence/Absence Analysis
What files to give OrthoFinder

For each genome, use the protein FASTA file containing the longest protein per gene:

> *.longest.aa.fa

That means you will have 51 FASTA files, each with one protein per gene. Put all of them together in a single directory, like this:

> orthofinder_input/
> 
> ├── genome01.longest.aa.fa 
>
> ├── genome02.longest.aa.fa
> 
> ├── ...
>
> ├── genome51.longest.aa.fa


Make sure nothing else is in that folder.

---

## Running OrthoFinder

To start OrthoFinder, first create a main folder to contain input and output files:
> mkdir -p /data/group/medaglab/project/BWee/orthofinder_tutorial
> cd /data/group/medaglab/project/BWee/orthofinder_tutorial

If you are running from local, first download via:
>  wget https://github.com/davidemms/OrthoFinder/releases/latest/download/OrthoFinder.tar.gz

If running it from HPC, check where it's already been installed
> module spider Ortho

Results:

    Description:
      OrthoFinder is a fast, accurate and comprehensive platform for comparative genomics


    This module can be loaded directly: module load OrthoFinder/2.5.4-foss-2020b

    Help:
      Description
      ===========
      OrthoFinder is a fast, accurate and comprehensive platform for comparative genomics


      More information
      ================
       - Homepage: https://github.com/davidemms/OrthoFinder


If you do life is easy and just run the OFind.sh script

> sbatch /data/group/medaglab/project/BWee/scripts/Ofind.sh

However, before running, please change your directories in the script as follows.

    # Input directory containing protein FASTA files (51 genomes)
      INPUT_DIR="./panfaba/Longest_AA"
      orthofinder -f orthofinder_input -t 24

    # Output directory (OrthoFinder will create this)
      OUTPUT_DIR="./panfaba/Longest_AA/orthofinder_results" #having a folder here already created with the folder name will throw an error

No extra configuration is needed.

---

## What OrthoFinder does for you

1. Performs all-vs-all protein similarity search (using DIAMOND)  
2. Clusters proteins into orthogroups (gene families)  
3. Differentiates paralogs from orthologs  
4. Generates a presence/absence matrix  

You don’t have to manually pre-cluster, rename genes, or filter paralogs yourself.

---

## Finding the presence/absence matrix

When OrthoFinder finishes, look in this folder inside the output directory:

> Results_/
> 
>  └── Orthogroups/
> 
>    ├── Orthogroups.tsv
> 
>    ├── Orthogroups.GeneCount.tsv
> 
>    ├── Orthogroups_SingleCopyOrthologues.tsv


The file `Orthogroups.GeneCount.tsv` is the one you want.  
- Rows = gene families  
- Columns = genomes  
- Values = number of genes from each family in each genome  

To get presence/absence data:  
- `>0` → present  
- `0` → absent  

---

## Two important notes

### 1. File naming

OrthoFinder uses the FASTA filenames as genome IDs, so keep them clean and simple:

> Genome01.fa 
> Genome02.fa

You *can* use longer filenames like `VFABA.29H.pgsb.r1.Oct2025.longest.aa.fa`, but the output tables can become harder to read.

### 2. Paralogs are expected

Some gene families may have multiple copies in a genome. For example:

GenomeA = 2 GenomeB = 1 GenomeC = 0


For presence/absence analysis, consider any count ≥1 as present, and 0 as absent.

---

## When OrthoFinder might not be the best choice

| Situation                     | Recommended tool        |
|-------------------------------|-----------------------|
| Read-based presence/absence   | Coverage mapping      |
| Bacterial pangenomes          | Panaroo or Roary      |
| Large datasets (500+ genomes) | MMseqs2               |


---



