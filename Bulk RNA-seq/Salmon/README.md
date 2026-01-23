# RNA-seq pipelines for Vigna unguiculata

## References
https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
https://combine-lab.github.io/salmon/getting_started/#indexing-txome

## First set-up HPC access

Don’t bother trying this on your local machine, especially if you have multiple samples.
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

> Note: Transcriptome and genome versions mismatch (v1.2 vs v1.0) — acceptable for this build.

---

## Pre-processing

### 1. Generate decoy list

Salmon requires a list of genome sequences to act as **decoys**.

```bash
grep "^>" <(gunzip -c ./Genomics/Vunguiculata_540_v1.0.softmasked.fa.gz) | cut -d " " -f 1 > decoysVu.txt
sed -i.bak -e 's/>//g' decoysVu.txt
```

This produces `decoysVu.txt` containing all chromosome/scaffold names.

---

### 2. Concatenate transcripts + genome

Salmon’s **decoy-aware mode** requires the genome sequences to come **after** transcripts in a single FASTA (`gentrome`):

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



Do you want me to add that too?
