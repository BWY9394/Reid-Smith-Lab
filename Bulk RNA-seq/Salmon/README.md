# RNA-seq pipeline using Salmon for mapping, tximport for export of counts & abundance, DESeq2 for DEG analysis

## References
https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/ #For indexing

https://combine-lab.github.io/salmon/getting_started/#indexing-txome #For running Salmon after indexing

https://bioconductor.statistik.tu-dortmund.de/packages/3.18/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#why-un-normalized-counts #DESeq2 tute

https://cloud.wikis.utexas.edu/wiki/spaces/bioiteam/pages/47731451/Mapping+tutorial #General information bank for bioinfomatics

## First set-up HPC access

You could run this on your local machine, the reason not to is because you have space issues (like me), and you also have access to the HPC (every student can request access with PI approval).
Make sure your PI has granted you access (AskICT on La Trobe intranet).

Also set up:

> 1. **PuTTy** — terminal access to HPC
> 2. **WinSCP** — file transfer and browsing with GUI

---

## Raw files for indexing

We need two main reference files for Salmon indexing:

1. **Genome (soft-masked)**

> Example: `Vunguiculata_540_v1.0.softmasked.fa.gz`

2. **Transcriptome (all splice variants, UTRs and exons)**

> Example: `Vunguiculata_540_v1.2.transcript.fa.gz`

> Note: Transcriptome and genome versions mismatch (v1.2 vs v1.0) is fine, all comes from v1.2 folder hosted on Phytozome.

---

## Pre-processing for indexing 

### 1. Generate decoy list

Salmon requires a list of genome sequences to act as **decoys**.

Run the below code (adjust genome/transcriptome names to your use-case), note that it's 2 lines of separate code:
```bash
grep "^>" <(gunzip -c ./Genomics/Vunguiculata_540_v1.0.softmasked.fa.gz) | cut -d " " -f 1 > decoysVu.txt
sed -i.bak -e 's/>//g' decoysVu.txt
```

This produces `decoysVu.txt` containing all chromosome/scaffold names. Change Vu to whatever works for you, for me I'm starting with cowpea. So. Yea. More for housekeeping so you know what is what


### 2. Concatenate transcripts + genome

Salmon’s **decoy-aware mode** requires the genome sequences to come **after** transcripts in a single FASTA (`gentrome`).

Run the below code (adjust genome/transcriptome names to your use-case):

```bash
cat ./Genomics/Vunguiculata_540_v1.2.transcript.fa.gz \
    ./Genomics/Vunguiculata_540_v1.0.softmasked.fa.gz \
    > Vunguiculata_540_v1.0.gentrome.fa.gz
```

---

## Salmon Indexing

### Module load

```bash
module load Salmon/1.4.0-GCC-11.2.0
```

### Index command

```bash
salmon index \
  -t Vunguiculata_540_v1.0.gentrome.fa.gz \
  -d decoysVu.txt \
  -p 12 \
  -i salmon_index
```

**Explanation of flags**:

* `-t`: transcript + genome FASTA (`gentrome`)
* `-d`: decoy list (genome sequences)
* `-p`: number of CPU threads
* `-i`: output index directory

> **Important:** `--gencode` **not used** — headers do **not** contain GENCODE-style pipe-separated IDs (e.g., `>ENST00000335137.4|ENSG000001234`), so Salmon will correctly use the first whitespace-delimited token as transcript ID.

You can check with:

```bash
zcat ./Genomics/Vunguiculata_540_v1.2.transcript.fa.gz | head -n 5
```

---

## Notes

* Full transcripts with UTRs improve isoform quantification and mapping accuracy.
* Decoy-aware indexing reduces false-positive mapping to transcripts.
* Indexing can take a few minutes to hours depending on genome size and HPC resources.
* Why do we create an index?
  * Creating an index for a computer database (which is basically what any reference genome/transcriptome is), allows for quick access to any "record" (gene+gene metadata), given a short "key" (usually ID of splice variant/isoform for a specific gene locus).
  * In other words, creating an index for a reference sequence allows it to more rapidly place a read on that sequence at a location where it knows at least a piece of the read matches perfectly or with only a few mismatches.
  * By jumping right to these spots in the genome, rather than trying to fully align the read to every place in the genome, it saves a ton of time.
  * Also, once you index once, you shouldn't need to index it again for other mapping runs for the same species with the same pipeline.

---

## Salmon Quantification

Once the Salmon index is ready, the next step is to quantify transcript expression for all samples. Here, we assume **paired-end FASTQ files** are stored in individual folders named `C1` → `C24` inside the `Raw` directory.

### SLURM Batch Script

Run salmon_quant_transcript.sh script on HPC:

```bash
sbatch /data/group/medaglab/project/BWee/scripts/salmon_quant_transcripts.sh
```
Which looks like below:
```bash
#!/bin/bash
#SBATCH --time=20:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --job-name="salZBF1to28"
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@ltu.edu.au # send-to address

#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b
#/data/group/medaglab/project/BWee/Salmon/Raw/FUN_TC_Vu/

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Run salmon quant

for fn in Raw/ZBF/ZBF{1..28}; do
  samp=$(basename "${fn}")
  echo "Processing sample ${samp}"

  # -l A: salmon auto-detects library type (stranded/unstranded)
  # -1: forward reads
  # -2: reverse reads (ignore if single-end)
  # -p 8: use 8 CPU threads
  # --validateMappings: selective alignment mode
  # -o: output directory for quant results

  salmon quant -i salmon_indexVu -l A \
    -1 "${fn}/${samp}_1.fq.gz" \
    -2 "${fn}/${samp}_2.fq.gz" \
    -p 8 --validateMappings -o "quantsZBF1_ZBF28/${samp}_quant"
done


# Now make an output file
echo -e "Sample\tNumReads\tNumMapped\tPercentMapped" > quantsZBF1_ZBF28/salmon_read_summary_ZBF1_ZBF28.tsv


#Extract mapping stats out into tab-separated table 
for d in quantsZBF1_ZBF28/*_quant; do
   samp=$(basename "$d" _quant)
   numReads=$(jq '.num_processed' "$d/aux_info/meta_info.json")
   numMapped=$(jq '.num_mapped' "$d/aux_info/meta_info.json")
   pctMapped=$(jq '.percent_mapped' "$d/aux_info/meta_info.json")

   echo -e "${samp}\t${numReads}\t${numMapped}\t${pctMapped}" \
     >> quantsZBF1_ZBF28/salmon_read_summary_ZBF1_ZBF28.tsv
done


# Optional: email summary to yourself
mail -s "Salmon Quant Finished" -a quantsZBF1_ZBF2/salmon_read_summary_ZBF1_ZBF28.tsv BWeeYang@ltu.edu.au <<< "All samples finished. See attached summary."

```

 **Explanation of flags**

* `-l A` — automatically detect library type (stranded/unstranded)
* `-1` — forward reads
* `-2` — reverse reads (omit if single-end)
* `-p 32` — number of CPU threads
* `--validateMappings` — selective alignment mode for improved accuracy
* `-i` — input Salmon index
* `-o` — output directory for per-sample quantification results

---

### Notes

* Each sample’s output includes **transcript-level quantification** (TPM, counts).
* The script processes samples sequentially in a **single SLURM job**. For very large datasets, consider submitting each sample as a separate job to parallelize.
* Make sure your folder structure matches the script and FASTQ files are named correctly (`C1_1.fq.gz`, `C1_2.fq.gz`, … `C24_1.fq.gz`, `C24_2.fq.gz`).

---


## DESeq2 (DEG analysis) and tximport (extracting out your TPMs for gene abundance comparison, raw counts for DEG)

First, we run tximport to export files out. This will be done in R locally, so please move your RNAseq quantification .sf containing output folders out to a local location on your PC.

## References
https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html#salmon
 
## Raw files for quantification
* List of isoform-level reference of transcript ID to gene ID. It is important that the transcript ID matches your reference transcriptome used for the initial mapping steps in salmon
** In this example it's Vunguiculata_540_v1.2.transcript.fa from Phytozome, download it from https://phytozome-next.jgi.doe.gov/info/Vunguiculata_v1_2
* .sf files per sample

```plaintext
#Example of file structure
Root ([P:\LAB - OPV - Dugald_Reid\Whereeveryourfilesare])/
├── quantsZBF1_ZBF28/
│ ├── ZBF1_quant
│   ├── aux_info
│   ├── libParams
│   ├── logs
│   ├── cm_info.JSON
│   └── lib_format_counts.JSON
│   └── quant.sf (this is where your tpm, counts, read length per sample sits in)
│ ├── ZBF2_quant
...
```

## Script needed
* DESeq2_salmon-ZBF1_28.R

## Packages needed in R
```r
library(tidyverse)
library(tximport)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
```
## Set working directory (i.e. where your raw files actually are.)
```r
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/")
getwd() #check working directory set correctly
```

### 1. Set up files names and folder paths
```r
# Vector of folder names containing Salmon quantification outputs
folders <-  data.frame(folder_names=list.files("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/"))
write.csv(folders,file="filenames.csv",row.names=FALSE)
#already generated the folder path so let's use that
folders = read.csv("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/filenames.csv", header = TRUE)
# Converted explicitly to character to avoid factor-related issues
folders = as.character(folders$folder_names)
# Sort folders in natural (human-readable) order
# Ensures ZBF2 comes before ZBF10, etc.
folders <- gtools::mixedsort(folders) # Inspect sorted folder names
folders
# Generate full file paths to each Salmon quant.sf file
# Assumes standard Salmon output structure: <sample>_quant/quant.sf
files = file.path( folders, "quant.sf")
# Assign sample names to each file
# These names must match rownames in the sample metadata (sampleTable)
names(files) <- paste0("ZBF", 1:28)
```
This will give you the below object structure, linking the quant.sf file location with a specific sample ID that makes sense for that sample.+
```r
files
> head(files)
                 ZBF1                  ZBF2                  ZBF3                  ZBF4                  ZBF5                  ZBF6 
"ZBF1_quant/quant.sf" "ZBF2_quant/quant.sf" "ZBF3_quant/quant.sf" "ZBF4_quant/quant.sf" "ZBF5_quant/quant.sf" "ZBF6_quant/quant.sf" 
```
You then run the below:
```r
# Returns TRUE if every expected file is present. So each rowname will be what is named to your salmon output abundnace/count/readlength per sample, as well as where that is stored, so should not have ANY  mismatches. 
all(file.exists(files)) #If FALSE, do not proceed until the file -> filename link is fixed and returns TRUE
all(file.exists(files))
```
This checks if your file paths to sample ID for naming samples is correct i.e. For ZBF1 at the supposed path of ZBF1_quant/quant.sf, does it exist.


### 2. Set up tx2gene annotation file 
This is necessary to map transcript ID used in salmon to gene-level IDs that tximport will summarize transcript isoforms to.

In this case, it is the transcript.fa files, which can be downloaded from:

https://phytozome-next.jgi.doe.gov/info/Vunguiculata_v1_2

Click downloaad and look for Vunguiculata_540_v1.2.transcript.fa

The annotation file will also be present in their database, will be a txt file with annotation_info as a suffix. This file is provided as an example in this repo (ExampleData).

We'll use the annotation file to create the tx2gene file, bnot the .fa file, but that's because I already know the 2 are (and shoudld by default) be synced. We'll explore how to check that later.
```r
annotations <-  read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) #54484 annotated transcripts
data.table::uniqueN(annotations$locusName) #Number of unique transcripts, discards repeated rows per transcript #31948
tx2gene <- annotations[,c(3,2)] #54484
head(tx2gene)
```
This will give you 2 columns, first being the transcript name, the second being gene IDs
```r
> head(tx2gene)
    transcriptName      locusName
1 Vigun01g000100.1 Vigun01g000100
2 Vigun01g000200.1 Vigun01g000200
3 Vigun01g000200.2 Vigun01g000200
4 Vigun01g000250.1 Vigun01g000250
5 Vigun01g000300.1 Vigun01g000300
6 Vigun01g000300.2 Vigun01g000300
```

You can also check that there are no duplicates in the transcriptID with:
```r
any(duplicated(tx2gene[,1]))
```
which will give
```r
> any(duplicated(tx2gene[,1]))
[1] FALSE
```

Great! Technically that's all we need to then start tximport. However, at this stage it's good to do some sanity checks as below. For instance, first we want to check that what we used for the annotation file to actual .fa file link.
```r
###Optional
#now, if you're actually not sure whether your annotations files is set up correctly, run the below
transcript_fasta_data <- read.fasta(file = "P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.transcript.fa.gz", as.string = TRUE) #54484 annotated transcripts

transcript_fasta_data_df <- do.call(rbind, lapply(transcript_fasta_data, function(seq_obj) {
  data.frame(
    transcriptName = attr(seq_obj, "name"),
    Vuannotation = attr(seq_obj, "Annot"),
    ft_ni_sequence = as.character(seq_obj),
    stringsAsFactors = FALSE
  ) 
})) #54484 annotated transcripts
```
What we really want there is the number of rows, which will correspond with the number of unique transcript IDs provided to salmon for mapping of reads to the reference transcriptome. We can see that is matches (54484 unqiue transcript IDs) with the the tx2gene (54484 unqiue transcript IDs) dataframe created from the annotations file.

I also provide below some additional code for counting unique isoforms per gene and saving it into a table for export.
```r
#you may also be interested in:
#For accounting, first let's determine how many isoforms per transcript
duplicate_counts_df_counts <- as.data.frame(table(annotations$locusName)) |> #table() counts how many times each unique locusName appears in your annotations data frame.
  subset(Freq > 1) |>       # keep only loci with >1 transcript
  setNames(c("locusName", "Transcript_Isoforms")) #gives no. isoforms per transcript if there are >1
#Ok, so let's summarize to gene-level, make sure it matches above number of 31948
annotations_GL <- annotations %>%
  dplyr::distinct(locusName, .keep_all = TRUE) #31948 more accurate because some genetranscripts might not have .1 read. Some start with .2

#Could be useful to know also how many isoforms per gene
annotations_GL <- left_join(annotations_GL,duplicate_counts_df_counts, by="locusName", relationship="one-to-one")
write.csv(annotations_GL,
          file="P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info_GL.csv",
          row.names = FALSE)
          
```

### 3. Generating TPMs, raw counts and other summarizations with tximport and tidyverse
Great! Now we are ready to create our tximport-mediated object for some downstream analysis.
```r
txi.salmon.tsv = tximport(files,
                          type = "salmon",
                          tx2gene = tx2gene,
                          dropInfReps = TRUE,
                          countsFromAbundance = "no") # don't worry about this too much, if you are using DESeq2, it will do its own thing with count normalization. And you will use TPMs for gene abundance anyway.

#save it into a dataframe
final <- data.frame(txi.salmon.tsv)
head(final)
```

And that gives you a dataframe to extract TPMs and the like from!
```r
> head(final)
               abundance.ZBF1 abundance.ZBF2 abundance.ZBF3 abundance.ZBF4 abundance.ZBF5 abundance.ZBF6 abundance.ZBF7 abundance.ZBF8 abundance.ZBF9 abundance.ZBF10
Vigun01g000100      24.582169      23.796251      27.654636      26.291820      24.318085      29.925334      23.540728      21.786672      24.270277       24.711348
Vigun01g000200     652.507434     691.774213     593.485490     596.899794     583.840262     458.590237     616.001666     637.918191     542.073237      475.285512
```
You'll also observe that the samples are correctly labelled with what you provided with the "files" argument in tximport(files,...). So everything works

Next, you will also want to set up some metadata for 1. renaming your files from 


```
### 3. Now set up filepaths to tell tximport where to look in to get your .sf files
```r
folders <-  data.frame(folder_names=list.files("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/"))
getwd()
write.csv(folders,file="filenames.csv",row.names=FALSE)
folders = read.csv("filenames.csv", header = TRUE)
folders = as.character(folders$folder_names)
folders <- gtools::mixedsort(folders)
folders
#generate a data frame with all filesnames and paths
files = file.path( folders, "quant.sf")
files
```

Done correctly, it'll look like:
```r
> files
 [1] "ZBF1_quant/quant.sf"  "ZBF2_quant/quant.sf"  "ZBF3_quant/quant.sf"  "ZBF4_quant/quant.sf"  "ZBF5_quant/quant.sf"  "ZBF6_quant/quant.sf"  "ZBF7_quant/quant.sf" 
 [8] "ZBF8_quant/quant.sf"  "ZBF9_quant/quant.sf"  "ZBF10_quant/quant.sf" "ZBF11_quant/quant.sf" "ZBF12_quant/quant.sf" "ZBF13_quant/quant.sf" "ZBF14_quant/quant.sf"
[15] "ZBF15_quant/quant.sf" "ZBF16_quant/quant.sf" "ZBF17_quant/quant.sf" "ZBF18_quant/quant.sf" "ZBF19_quant/quant.sf" "ZBF20_quant/quant.sf" "ZBF21_quant/quant.sf"
[22] "ZBF22_quant/quant.sf" "ZBF23_quant/quant.sf" "ZBF24_quant/quant.sf" "ZBF25_quant/quant.sf" "ZBF26_quant/quant.sf" "ZBF27_quant/quant.sf" "ZBF28_quant/quant.sf"
```

### 4. Now you're ready to extract out the necessary info from your sf files
```r
#now do tximport
txi.salmon.tsv = tximport(files,
                          type = "salmon",
                          tx2gene = tx2gene,
                          dropInfReps = TRUE,
                          countsFromAbundance = "no")


#save into dataframe
final <- data.frame(txi.salmon.tsv)
```
An issue is that your sample names may not be informative and correct for downstream processing (i.e. they are named ZBF1, ZBF2 etc. which does not result in correct groupings later)

### 5. Adding some metadata information into saved headers for sanity
```r
# create headers
header = c(paste("TPM.", folders,sep=""), paste("counts.", folders,sep=""), paste("length.", folders,sep=""), "countsFromAbundance")
print(header)
datacheckheader <- data.frame(final=colnames(final),header=header)

#save now so that you can add metadata in later easily and align with current sample names
write.csv( datacheckheader, file= "ZBF1_28_metadata.csv", row.names=F, quote=F)
#just check that things line up

#now change header names
colnames(final) = header
colnames(final)

#add in the locus
final$locusName <- rownames(final)
#View(final)
final <- final[,c(86,1:85)]

#load metadata
metadata <- readxl::read_xlsx("Z:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/01.RawData/Copy of Sample Register.xlsx")
colnames(metadata)

ZBFmetadata <- (dplyr::filter(metadata,str_detect(`Exp Sample ID`, "ZBF")))

#extract TPM
header=(c("locusName", c(paste("TPM",ZBFmetadata$Description, sep=".")))
)
header <- str_remove(header, pattern= "VuFUN ")
header
TPM = data.frame(final[,grep('TPM.', colnames(final))])
TPM <- cbind(locusName = rownames(TPM), TPM)
colnames(TPM)=header
View(TPM)
write.csv( TPM, file= "ZBFTPM.csv", row.names=F, quote=F)


#extract only counts and save
header=(c("locusName", c(paste("counts",ZBFmetadata$Description, sep=".")))
)
header <- str_remove(header, pattern= "VuFUN ")
header
counts = data.frame(final[,grep('counts.', colnames(final))])
counts <- cbind(locusName = rownames(counts), counts)
counts
colnames(counts)=header
View(counts)
write.csv( counts, file= "ZBFcounts.csv", row.names=F, quote=F)
```



