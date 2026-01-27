# RNA-seq pipelines for Vigna unguiculata

## References
https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/ #For indexing
https://combine-lab.github.io/salmon/getting_started/#indexing-txome #For running Salmon after indexing

## First set-up HPC access

You could run this on your local machine, the reason not to is because you have space issues (like me), and you also have access to the HPC (every student can request access with PI approval).
Make sure your PI has granted you access (AskICT on La Trobe intranet).

Also set up:

> 1. **PuTTy** — terminal access to HPC
> 2. **WinSCP** — file transfer and browsing with GUI

---

## Raw files

We need two main reference files for Salmon indexing:

1. **Genome (soft-masked)**

> Example: `Vunguiculata_540_v1.0.softmasked.fa.gz`

2. **Transcriptome (all splice variants, UTRs and exons)**

> Example: `Vunguiculata_540_v1.2.transcript.fa.gz`

> Note: Transcriptome and genome versions mismatch (v1.2 vs v1.0) — .

---

## Pre-processing

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

---

## Notes

* Full transcripts with UTRs improve isoform quantification and mapping accuracy.
* Decoy-aware indexing reduces false-positive mapping to transcripts.
* Indexing can take a few minutes to hours depending on genome size and HPC resources.

---

Absolutely! Here’s your **Salmon quantification step** rewritten in the **same GitHub Markdown tutorial style** as your previous sections:

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
#SBATCH --job-name="salC1|24"
#SBATCH --partition=day
#SBATCH --mail-type=ALL
#SBATCH --mail-user=BWeeYang@latrobe.edu.au # send-to address

#Available salmon versions
#Salmon/1.4.0-GCC-11.2.0
#Salmon/1.4.0-gompi-2020b

#module load
module load Salmon/1.4.0-GCC-11.2.0

#Run salmon quant

for fn in Raw/C{1..24}; do
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
    -p 8 --validateMappings -o "quantsC1_C24/${samp}_quant"
done

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




