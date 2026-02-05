BiocManager::install("apeglm")
library(tidyverse)
library(tximport)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(seqinr)
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/")
getwd()

#generate filenames/folderpaths
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
files
# Sanity check: confirm that all quant.sf files exist
# Returns TRUE if every expected file is present. So each rowname will be what is named to your salmon output abundnace/count/readlength per sample, as well as where that is stored, so should not have ANY  mismatches. 
all(file.exists(files)) #If FALSE, do not proceed until the file -> filename link is fixed and returns TRUE
all(file.exists(files))


#annotations
annotations <-  read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) #54484 annotated transcripts
data.table::uniqueN(annotations$locusName) #Number of unique transcripts, discards repeated rows per transcript #31948
tx2gene <- annotations[,c(3,2)]
str(tx2gene)
any(duplicated(tx2gene[,1]))

###Optional
#now, if you're actually not sure whether your annotations files is set up correctly, run the below
transcript_fasta_data <- read.fasta(file = "P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.transcript.fa.gz", as.string = TRUE) #54484 annotated transcripts
#
transcript_fasta_data_df <- do.call(rbind, lapply(transcript_fasta_data, function(seq_obj) {
  data.frame(
    transcriptName = attr(seq_obj, "name"),
    Vuannotation = attr(seq_obj, "Annot"),
    ft_ni_sequence = as.character(seq_obj),
    stringsAsFactors = FALSE
  ) 
})) #54484 annotated transcripts
str(transcript_fasta_data_df)


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


#now do tximport
txi.salmon.tsv = tximport(files,
                          type = "salmon",
                          tx2gene = tx2gene,
                          dropInfReps = TRUE,
                          countsFromAbundance = "no")


#save into dataframe
final <- data.frame(txi.salmon.tsv)
str(final)

#some metadata matching
filenames <-  file.path( folders, "quant.sf")
filenames
sampleNames <- sub("_quant/quant.sf", "", filenames)
sampleNames #

metadata <- readxl::read_xlsx("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/01.RawData/Copy of Sample Register.xlsx") #load meta data in
ZBFmetadata <- metadata %>%
  filter(Project == "Zn Binding and Filamentation") #filter out for just ZBF experiment
head(ZBFmetadata)



## ---------------------------
## Save Output
## ---------------------------
#1. Look at saving raw counts + TPM.

#1.1 Extract only counts and save
header=(c("locusName", c(paste("counts",ZBFmetadata$Description, sep=".")))) #create headers based on metadata, adding counts/TPM to start of it, with a "." to separate.
header <- str_remove(header, pattern= "VuFUN ") #also remove VuFUN, unnecessary 
header
counts = data.frame(final[,grep('counts.', colnames(final))])
colnames(counts) #you may have an extra column called countsfromAbundance to remove
counts <- cbind(locusName = rownames(counts), counts[,-29]) #generate the same locusName column that from rownames, matching the saved locusName column in header, also remove last column as pointless
colnames(counts)
colnames(counts)=header 
View(counts)
write.csv( counts, file= "ZBFcounts.csv", row.names=F, quote=F)

#1.2 Extract only TPM and save
header=(c("locusName", c(paste("TPM",ZBFmetadata$Description, sep=".")))
)
header <- str_remove(header, pattern= "VuFUN ")
header
colnames(final)
TPM <-  data.frame(final[,grep('abundance', colnames(final))])
colnames(TPM)
TPM <- cbind(locusName = rownames(TPM), TPM)
colnames(TPM)
colnames(TPM)=header
View(TPM)
write.csv( TPM, file= "ZBFTPM.csv", row.names=F, quote=F)


#1.3 Perhaps you would like mean +SE out of TPM?

TPM_long = pivot_longer(TPM,
                        cols = -locusName,
                        names_to = 'Sample',
                        values_to = 'TPM')

head(TPM_long$Sample)

#1.4 Add back in metadata
colnames(TPM)
TPMheaders <- colnames(TPM[,-1])
ZBFmetadata_tpm <- data.frame(Sample=TPMheaders,ZBFmetadata[,c(4,5)])
ZBFmetadata_tpm <- left_join(TPM_long,ZBFmetadata_tpm,by="Sample")
View(ZBFmetadata_tpm)

colnames(ZBFmetadata_tpm)

#1.5 Do stats for TPM file.
TPM_mean_comb <- plyr::ddply(
  ZBFmetadata_tpm,
  c("locusName", "Identifying.Allele"),
  summarise,
  mean = mean(TPM, na.rm = TRUE),
  sd   = sd(TPM, na.rm = TRUE),
  n    = sum(!is.na(TPM)),
  sem  = sd(TPM, na.rm = TRUE) / sqrt(sum(!is.na(TPM)))
)

writexl:: write_xlsx(TPM_mean_comb,
                      "ZBF1_28_TPMstats2.xlsx")




#2. Construct DESeqDataset object from previously generated txi object
# sampleTable contains sample-level metadata (conditions, genotypes, etc.)
colnames(ZBFmetadata)
sampleTable <- ZBFmetadata[,c(5)]
sampleTable
sampleTable <- data.frame(sampleTable)
rownames(sampleTable) <- sampleNames
sampleNames
sampleTable$Identifying.Allele <- as.factor(sampleTable$Identifying.Allele)
levels(sampleTable$Identifying.Allele)
sampleTable #so your rownames will be sampleID that was assigned to the sample name when building your txi salmon object, i.e sampleNames

# Design formula specifies Identifying.Allele as the main factor to model
# (i.e., differences between constructs / lines / replicates)
dds <- DESeqDataSetFromTximport(txi.salmon.tsv, sampleTable, ~ Identifying.Allele) #~ this is the factor variable, in this case it's the construct representing each line/rep
#refactor for comparisons
dds$Identifying.Allele <- relevel(dds$Identifying.Allele, ref = "WT13")


# Run the DESeq2 pipeline:
# - estimate size factors
# - estimate dispersion
# - fit negative binomial models
dds <- DESeq(dds)

# View available result contrasts based on the design formula
resultsNames(dds)
#
res <- results(dds)


vsd <- vst(dds, blind=TRUE)


plotPCA(vsd, intgroup=c("Identifying.Allele"),returnData = FALSE)

#build a better PCA
pca_data <- plotPCA(
  vsd,
  intgroup = "Identifying.Allele",
  returnData = TRUE
)
pca_data <- ren
colnames(ZBFmetadata)
pca_data <- dplyr::rename(pca_data, `RNA Sample Code`=name)
pca_data <- left_join(pca_data,ZBFmetadata[c(3,4)],by="RNA Sample Code")
pca_data

percent_var <- round(100 * attr(pca_data, "percentVar"), 1)

levels(pca_data$Identifying.Allele)

#subset(pca_data, !Identifying.Allele %in% c("WT13", "WT15"))
ggplot(pca_data, aes(PC1, PC2, color = Identifying.Allele)) +
  geom_point(
    size = 2.8,
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = Description),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey70"
  )+ 
  theme_bw(base_size = 12)+
  labs(
    title = "Principal Component Analysis of RNA-seq Samples",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Allele"
  ) +
  scale_color_brewer(palette = "Set2") +
  #theme_classic(base_size = 12) 
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )

#Get normalized raw counts for each gene per sample
resnormcounts <- counts(dds, normalized=T)

write.csv(resnormcounts, file = "norm_counts_results.csv")


#WT13 vs DBL21
res <- results(dds, contrast=c("Identifying.Allele","WT13","DBL21"),alpha=0.1)
res
?results
resultsNames(dds)
res_tableOE <- lfcShrink(dds, coef="Identifying.Allele_DBL21_vs_WT13", type="apeglm")
res_tableOE <- lfcShrink(
  dds,
  contrast = c("Identifying.Allele", "WT13", "DBL21"),
  type = "apeglm"
)
res_tableOE
summary(res)
res
sum(res$padj < 0.05, na.rm=TRUE)
resOrdered <- res[order(res$padj),]
resOrderedDF <- as.data.frame(resOrdered)
resOrderedDF
write.csv(resOrderedDF, file = "results_pvals_WT13vsDBL21.csv")

shrink_plot <- plotMA(res_tableOE,main="LFC Shrinkage")  # if using lfcShrink (ashr or apeglm)
noshrink_plot <- plotMA(res)

dev.off()
pdf(file="test.pdf",height=6)
par(mfrow=c(1,2))
plotMA(res,main="Identifying.Allele WT13 vs DBL21:\nNo Shrinkage")
plotMA(res_tableOE,main="Identifying.Allele WT13 vs DBL21:\nLFC Shrinkage")
par(mfrow=c(1,1))  # reset
dev.off()


# pivot_wider(names_from = sample, values_from = SE)
               
  #pivot_wider(names_from = Identifying.Allele, values_from = TPM) #increases number of columns, decreasingrows. specify column to get name of output column from, as well as column to get cell values from.



any(duplicated(counts[,1])) #sho

ZBFmetadata <- data.frame(counthead=coounthead,ZBFmetadata[,c(4,5)])
str(counts)
str(ZBFmetadata)




## ---------------------------
## PCA of RNA-seq counts
## ---------------------------

library(ggplot2)

# 1. Prepare count matrix
count_mat <- counts[, sapply(counts, is.numeric)]
rownames(count_mat) <- counts$locusName

# sanity check
stopifnot(ncol(count_mat) == nrow(ZBFmetadata))

# 2. Log-transform
log_counts <- log2(count_mat + 1)

# 3. PCA (samples as rows)
pca <- prcomp(t(log_counts), scale. = TRUE)
keep <- rowSums(count_mat > 0) >= 5
count_mat <- count_mat[keep, ]

log_counts <- log2(count_mat + 1)
pca <- prcomp(t(log_counts), scale. = TRUE)
# variance explained
percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

# 4. PCA dataframe with metadata
pca_df <- data.frame(
  Sample = colnames(count_mat),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Genotype = ZBFmetadata$`Identifying Allele`)

pca_df$Sample

# 5. Plot PCA
ggplot(pca_df, aes(PC1, PC2,group=Genotype, color=Genotype)) +
  geom_text(aes(label=Sample, color=Genotype),check_overlap = FALSE, size = 3, vjust = -0.5)+
  geom_point(size = 3) +
  #stat_ellipse(level = 0.95, linewidth = 1) +
  labs(
    title = "PCA of RNA-seq samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_minimal()

ggplot(pca_df, aes(PC1, PC2, shape = Genotype)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95, linewidth = 1) +
  labs(
    title = "PCA of RNA-seq samples",
    x = paste0("PC1 (", percent_var[1], "%)"),
    y = paste0("PC2 (", percent_var[2], "%)")
  ) +
  theme_minimal()
