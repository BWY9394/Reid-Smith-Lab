### Walkthrough for STAR to QC RNA-seq reads 
#### Disclaimer: I wrote this in 20 minutes so still needs some cleaning up. Also specific to LTU HPC and WSL2.


First we index our genome with
1. Genome
2. gtf/gff3 file

and we then also align in the same script with the produced

3. STAR_INDEX

and your raw reads  

4. Paired raw reads


```bash
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --job-name=STAR_index
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=STRindAlign.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.weeyang@latrobe.edu.au #change to your email.

module load STAR/2.7.10b-GCC-11.3.0

echo "start indexing"

STAR \
  --runMode genomeGenerate \
  --runThreadN 12 \
  --genomeDir STAR_index \
  --genomeFastaFiles Genomics/Vunguiculata_540_v1.0.softmasked.fa \
  --sjdbGTFfile Genomics/Vunguiculata_540_v1.2.gene_exons.gtf \
  --sjdbOverhang 149

echo "finish indexing"

echo "start alignment"

STAR \
  --runThreadN 12 \
  --genomeDir ./STAR_index \
  --readFilesIn \
    ./Raw/C5/C5_1.fq.gz \
    ./Raw/C5/C5_2.fq.gz \
  --readFilesCommand zcat \
  --twopassMode Basic \
  --outFileNamePrefix STAR_test_ \
  --outSAMtype BAM SortedByCoordinate

echo "finished elignment"
```



This will give a bunch of stuff, but most importantly will produce your bam file e.g. STAR_test_Aligned.sortedByCoord.out.bam

Index it with (run directly in HPC terminal line by line)
```bash
module load SAMtools/1.16.1-GCC-11.3.0
samtools index /where/your/bam/file/is/STAR_test_Aligned.sortedByCoord.out.bam
module purge
```

Next, we convert the GFF3 annotation to a BED12 transcript annotation, which is the format required by the RSeQC tools we're using.

This step is best done on a linux machine, the HPC doesn't have the necssary tools (that I know off) to do this. This facilitated (for me) via [WSL2](https://learn.microsoft.com/en-us/windows/wsl/install).

We will require gffread.

```bash
#in WSL2
#install gffread
sudo apt install gffread
gffread     Genomics/Vunguiculata_540_v1.2.gene_exons.gff3     --bed     -o Genomics/Vunguiculata_540_v1.2.RSeQC.bed
head Genomics/Vunguiculata_540_v1.2.RSeQC.bed | column -t

# the above produced 13 columns for me as I have some thing like
# 13 contig_106      11433   11965   VigunL059001.1.v1.2     100     +       11548   11965   0,0,0   1       532,   0,       CDS=11548:11965;CDSphase=0;geneID=VigunL059001.v1.2
# the 13th column is from phytozome gff3, not necessary. 
# so we remove with
awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}'     Genomics/Vunguiculata_540_v1.2.RSeQC.bed     > Genomics/Vunguiculata_540_v1.2.RSeQC.BED12.bed
 
#check with the following command, should be nothing that pops up. if you do find Ben or chat
awk 'NF != 12 {print "BAD:", NF, $0; exit}'    Genomics/Vunguiculata_540_v1.2.RSeQC.BED12.bed

#you can also check your exon start stop per genes make sense.
head -20000 Genomics/Vunguiculata_540_v1.2.RSeQC.BED12.bed | column -t
```

Next we run RSeQC using:
1. The aligned BAM file
2. The corresponding BAM index (.bai)
3. The BED12 reference annotation

And aggregate results with MULTQC.

```bash
#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --job-name=SRRSeQC
#SBATCH --partition=day
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=SRRSeQC.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=b.weeyang@latrobe.edu.au

# If you did not index your bam file, uncomment this step
#module load SAMtools/1.16.1-GCC-11.3.0

#echo "Start Indexing your aligned reads against the genome.bam file"

#samtools index /data/group/reidlab/project/Ben/STAR/AlignedBam/STAR_test_Aligned.sortedByCoord.out.bam

#echo "Finished Indexing"


# Clear environment and load the exact cluster module
module purge
module load RSeQC/4.0.0-foss-2020a-Python-3.8.2


BAM_DIR="AlignedBam"
BED_REF="Genomics/Vunguiculata_540_v1.2.RSeQC.BED12.bed"
OUT_DIR="RSeQC/multqc_results/"


echo "Starting job on $(hostname)"

for bam in "$BAM_DIR"/*.bam; do
    sample=$(basename "$bam"  .bam)
    echo "Processing QC metrics for: $sample"
    
    # Run Read Distribution
    read_distribution.py -i "$bam" -r "$BED_REF" > "$OUT_DIR/${sample}.read_distribution.txt"
    
    # Run Infer Experiment
    infer_experiment.py -i "$bam" -r "$BED_REF" > "$OUT_DIR/${sample}.infer_experiment.txt"
    
    # Run BAM Statistics
    bam_stat.py -i "$bam" > "$OUT_DIR/${sample}.bam_stat.txt"
    
    # Run Inner Distance (For Paired-End data)
    inner_distance.py -i "$bam" -r "$BED_REF" -o "$OUT_DIR/${sample}"
    
    # Run TIN (Note: Can be time consuming on large BAMs, comment out if running too slow)
    tin.py -i "$bam" -r "$BED_REF" > "$OUT_DIR/${sample}.tin.txt"
done

echo "RSeQC job finished."


module purge
module load MultiQC/1.12-foss-2021b


echo "Aggregate with MultQC"

multiqc "$OUT_DIR" -o "$OUT_DIR" --title "RSeQC_QC_Report"

echo "RSeqQC and MultQC pipeline complete"







```

Inspect your multqc report to see (most importantly for us), the read_distribution.txt output
