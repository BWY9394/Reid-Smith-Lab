### Walkthrough for STAR to QC RNA-seq reads 
#### Disclaimer: I wrote this in 20 minutes so still needs some cleaning up. Also specific to HPC and WSL2.


First we index our genome with
1. Genome
2. gtf/gff3 file

```bash


```



This will give a bunch of stuff, but most importantly the folder called something like STAR_INDEX


Then we align with STAR
1. STAR_INDEX
2. Raw Reads

This will produce your bam file e.g. STAR_test_Aligned.sortedByCoord.out.bam

Index it with (run directly in HPC terminal)
```bash
module load SAMtools/1.16.1-GCC-11.3.0
samtools index /data/group/reidlab/project/Ben/STAR/AlignedBam/STAR_test_Aligned.sortedByCoord.out.bam


Next, we obtain our bed file.

This step is best done on a linux machine, the HPC doesn't have the necssary tools (that I know off) to do this.

We will require gffread.

```bash
#in WSL2

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

Next we do RSeQC on:
1. The 
