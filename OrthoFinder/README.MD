# OrthoFinder walkthrough

Based off https://davidemms.github.io/menu/tutorials.html

Adapted for LTU systems and Reid lab infrastructure

### Running OrthoFinder

2 options: 
1) Run on HPC **
2) Run on Reid lab HPC *
3) Run on your own computer *

\* For Linux: follow the instructions below.  
\* For Macs: it’s easiest to use Bioconda to install and then follow the instructions below.  
\* For Windows it is best to install the Windows Subsystem for Linux and then continue as for Linux below.  
\*\* Would recommend this option if you are running a huge number of species and/or the species have very large genomes e.g. Faba which has a genome of 13Gb. Also if you don't want to muck around too much with with installing python + cloning the OrthoFinder git repo 

### 1. Setting it up (Skip to Step 2 if picking option 1, Step 2-alternate if running from Reid lab HPC) ###

TBC ( or just follow here https://davidemms.github.io/orthofinder_tutorials/downloading-and-running-orthofinder.html)

### 2. Running an analysis

OrthoFinder uses the high-performance protein sequence aligner DIAMOND by default to perform all-versus-all sequence similarity searches, providing a much faster alternative to BLAST for large-scale phylogenomic analyses. These sequence similarity results are used to infer orthogroups, which are groups of genes descended from a single ancestral gene in the last common ancestor of the species being compared. OrthoFinder subsequently uses phylogenetic analysis to infer orthologous and paralogous relationships within these groups.

The assumed folder structure:

```
user/Ortho/
├── 8splegumes #whatever you have name your project
│   ├── species1_fullproteome.fa
│   ├── species2_fullproteome.fa
│   ├── species3_fullproteome.fa
│   ├── ...
└── tools
    └── primary_transcripts.py

```

We assume you are in <user/Ortho/>

What OrthoFinder requires:
1. Proteomes formatted as either .fa”, “.faa”, “.fasta”, “.fas” or “.pep”., that only contain the primary transcript.

It is often the case that you will need to generate a fasta file containing just the primary transcript of each protein. In this case (and most), we define this as the longest representative transcript.

To check multiple fastas at once, run relative to where you have placed your proteomes:
```bash
# cd to whereever your proteomes are
for f in ./8splegumes/*.fa; do
    echo "===== $f ====="
    head -n 10 "$f"
    echo
done
```


Here is what an example subset of what each proteome fasta will contain:
```
===== ./8splegumes/Araport11_pep_20220914_representative_gene_model.fa =====
>AT1G06190.5 | Symbols: RHON1 | Rho-N 1 | chr1:1892467-1893681 REVERSE LENGTH=404
MFYLSICLLYVVYLLGNTGYTLSDSRCFFNSAVSRRTLAILPCSSCLDHKNGRLKSVPNRSSFVCRASSGGYRRNPDFSR
LNKHGYRGNNRQSGGREDFDIENSDMLSSRNGPLFNLSSSPKFQATSSPGPREKEIVELFRKVQAQLRARAAAKKEEKKI
EEASKGQGKESETVDSLLKLLRKHSGEQSKRQVSKFSSQGEVQGDTVDKQDRTGNLVTSGNKDNNASSFTRPTSSFRRKS
PVPRSQSPPAYSSEATFDQSSSYSVTWTQKKDTVELHDEPEHEPAYEHEHEPENESEPGPVTTMLEPDSELKPESSSFYQ
EEEDDDVTFDVLSQDDGILDVLSDDDESLDDADEDSDEAEEEAVKDLSELKLVELRGIAKSRGLKGLSKMKKAELVELLG
SDSS*
>AT1G06620.1 | Symbols: no symbol available | no full name available | chr1:2025618-2027094 FORWARD LENGTH=365
MESSLPQVAALDRSTLLKAFDETKTGVKGLIDAGITEIPSIFRAPPATLTSPKPPSSSDFSIPTIDLKGGGTDSITRRSL
VEKIGDAAEKWGFFQVINHGIPMDVLEKMIDGIREFHEQDTEVKKGFYSRDPASKMVYSSNFDLFSSPAANWRDTLGCYT

===== ./8splegumes/Carietinum_492_v1.0.protein_primaryTranscriptOnly.fa =====
>Ca_25697 pacid=40134984 transcript=Ca_25697 locus=Ca_25697 ID=Ca_25697.v1.0.492 annot-version=v1.0
MQLKRLHGKWIRPQAAVEPAQPQPQPQPQPVTAPATPDFVGKIMDFIEAQMAHNAKVLESQSRLEASLRTLQESLNSVVK
DVDAIQAHLVIKLSFEDIANIGR*
>Ca_25703 pacid=40134985 transcript=Ca_25703 locus=Ca_25703 ID=Ca_25703.v1.0.492 annot-version=v1.0
MERKRTQTGSKFVRQWYEIGTLTVDYRLLHYFLCYILVPKAGNYSQITKFELQIMYFMFEGCPLNWAYFVRSVMFASKDA
VGLPYAKEVTRILDHFKLDFSSEVMYTQTKDNLLDLSVLLSMKIIWNEQQHAYVHKGHMGQSSTTAQAADEVEDTASDDD
HVSNCMVMDKLDSLQTFITT*
>Ca_25699 pacid=40134986 transcript=Ca_25699 locus=Ca_25699 ID=Ca_25699.v1.0.492 annot-version=v1.0
MPRHLISDAHEWINEIPTVLVYYPAKPQPRERAWRNQWEKKTLIHQVLDCSPTNRERELGLDRRETGKARGSLELNCVSV
TVHDLDLFSVGLGHAHDTRTYDTRLRWKHSNNPRTLKKLVVGPWVGSIGPPLVCTGWLVPSAGDALLALIGRSDMKFLNY

===== ./8splegumes/CBA_Captain.aa.fa =====
>Ca_CA_31830.1
MWNLLSWVMEVVALMTIALANDSGRPPDWQDFIDIIALLVINSTIRDDEEVSELGPIVEKPRNKKFELRESVGNQESSWF
KC*
>Ca_CA_30584.1
MKTFPTCGDISNVDGFNVPMSLTPQGGSGDCKTSSCPANINDVCPSELQVVGSDGSVIACKNACLALNQDQYCCRGDYNT
EEKCPPTNYSIIFENQSPDAYSYAYDDKSSTTFTCFTNPSYAITFCPAT*
>Ca_CA_30710.1
MEVVLPIEVEIPSLKVLMEAKLEESEWIQTPYEKKIRPREFQEGDLVLKKILPIQKDYLI*
>Ca_CA_31903.1
MSKFDIDPGSISNRLRWLANGPSLHVFSYTGYVINGYTFYTKEQDDQTTMQNSGVTLVAEAMHVSTNGPSLHVFSYTGYV

...

```
To get the primary transcript only for the proteome of each species, we run the following command in the same folder where you have placed your proteomes:

```bash
 for f in ./8splegumes/*.fa; do
  python3 ./tools/primary_transcripts.py  $f ; 
  done
```

This command processes every .fa file in the current directory and runs OrthoFinder's primary_transcript.py script on each one.  

Syntax:  
- for f in *fa — finds files in the current directory whose names end in fa and assigns each filename to the variable f, one at a time.

The output is saved in a new primary_transcripts/ folder inside the same directory as each input FASTA file. The OrthoFinder script explicitly creates that folder automatically.
```
/user/Ortho/project/
├── 8splegumes
│   ├── species1_fullproteome.fa
│   ├── species2_fullproteome.fa
│   ├── species3_fullproteome.fa
│   ├── ...
└── tools
    └── primary_transcripts.py
├── primary_transcripts
│   ├── species1_primaryproteome.fa
│   ├── species2_primaryproteome.fa
│   ├── species3_primaryproteome.fa
|   ├── ..


```
Fantastic. Now we are ready to run the OrthoFinder bash script with:

```bash
sbatch <filepath of where you have placed your script Ofind script>/Ofind.sh
#in my case
# sbatch /data/group/reidlab/project/Ben/scripts/OFind.sh
```
But first check that your input directory is correct:  

```bash
#!/bin/bash
#SBATCH --job-name=Ofind6
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH -t 1-00:00:00 # time (D-HH:MM:SS)
#SBATCH --partition=day

#check that you are in the directory just before your primary_transcripts folder where your longest primary protein transcripts

# Load required modules
module load OrthoFinder/2.5.4-foss-2020b

# Input directory containing protein FASTA files (51 genomes)
INPUT_DIR="./8splegumes/primary_transcripts/" 
#assume you are in the directory one level up from ./8splegumes, so in Ortho for me

# Run OrthoFinder with DIAMOND (much faster for large datasets)
orthofinder -f ${INPUT_DIR} -t 32 -a 4

#-t : dictates no. total for DIAMOND/search-related work
#-a : number of parallel threads used for OrthoFinder’s internal analysis tasks, which are often RAM-intensive.

```

And that's it!

The results will look like:
```
/user/Ortho/
├── 8splegumes/
│   ├── species1_fullproteome.fa
│   ├── species2_fullproteome.fa
│   ├── species3_fullproteome.fa
│   └── ...
│
├── primary_transcripts/
│   ├── species1_primaryproteome.fa
│   ├── species2_primaryproteome.fa
│   ├── species3_primaryproteome.fa
│   └── ...
│
├── tools/
│   └── primary_transcripts.py
│
└── OrthoFinder/
    └── Results_Aug26/
        ├── Comparative_Genomics_Statistics/
        │   ├── Orthogroups_SpeciesOverlaps.tsv
        │   ├── Statistics_Overall.tsv
        │   └── Statistics_PerSpecies.tsv
        │
        ├── Gene_Duplication_Events/
        │   └── Duplications.tsv
        │
        ├── Orthogroup_Sequences/
        │   ├── OG0000000.fa
        │   ├── OG0000001.fa
        │   └── ...
        │
        ├── Orthogroups/
        │   ├── Orthogroups.tsv
        │   ├── Orthogroups.GeneCount.tsv
        │   ├── Orthogroups_SingleCopyOrthologues.txt
        │   └── Orthogroups_UnassignedGenes.tsv
        │
        ├── Orthologues/
        │   ├── species1/
        │   │   ├── species1__v__species2.tsv
        │   │   └── species1__v__species3.tsv
        │   └── ...
        │
        ├── Phylogenetic_Hierarchical_Orthogroups/
        │   ├── N0.tsv
        │   ├── N1.tsv
        │   └── ...
        │
        ├── Resolved_Gene_Trees/
        │   └── Resolved_Gene_Trees.txt
        │
        ├── Single_Copy_Orthologue_Sequences/
        │   ├── OG0000000.fa
        │   ├── OG0000001.fa
        │   └── ...
        │
        ├── Species_Tree/
        │   └── SpeciesTree_rooted.txt
        │
        ├── WorkingDirectory/
        │   └── ...
        │
        ├── Log.txt
        └── Citation.txt
```

What's we need from here:
- The orthogroup file we use for downstream applications to other analysis is in <Orthogroups/Orthogroups.tsv>
- The unassigned genes (single member orthogroups) in Orthogroups/<Orthogroups_UnassignedGenes.tsv>
- The fasta files that are generated for each orthogroup in <Orthogroup_Sequences/>
-  Species Tree generated to be helpful as well for inferring duplication/speciation event in <Species_Tree/>, you can open this in Dendroscope
- Lastly, for some stats, <Comparative_Genomics_Statistics/>

For an explanation of this and other outputs, refer to:(https://davidemms.github.io/orthofinder_tutorials/exploring-orthofinders-results.html)



### 2-alternate. Running an analysis on Reid-Lab HPC/Your personal computer.. 

TBC ( or just follow here https://davidemms.github.io/orthofinder_tutorials/downloading-and-running-orthofinder.html)

