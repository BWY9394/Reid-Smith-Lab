BiocManager::install("apeglm")
BiocManager::install("ashr")
install.packages("styler")
# Author: Benjamin Wee
# Purpose: tximport and DeSeq2 on vufun mutants, along with heatmap and venns, from Salmon workflow.
# Prelim data comparisons: For each timepoint, per genotype, what is the difference with +N
# Github: https://github.com/BWY9394/Reid-Smith-Lab/tree/main/Bulk%20RNA-seq/Salmon

# Step 1: Load packages and set working directory
library(tidyverse)
library(tximport)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(seqinr)
library(ashr)
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/C1_C24_VuFUN_timecourse_2025/quantsC1_C24/")
getwd()



# Step 2: generate filenames/folderpaths
folders <-  data.frame(folder_names=list.files("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/C1_C24_VuFUN_timecourse_2025/quantsC1_C24"))
folders
getwd()
write.csv(folders,file="filenames.csv",row.names=FALSE)
#Step 2: Start here if rerunnign analayis
folders = read.csv("filenames.csv", header = TRUE)
folders = as.character(folders$folder_names)
folders <- gtools::mixedsort(folders)
folders
#generate a data frame with all filesnames and paths
files = file.path(folders, "quant.sf")
files
names(files) <- paste0("C", 1:24) # Assign sample names to each file. These names MUST match rownames in the sample metadata (sampleTable)
files
getwd()

# Step 4: Sanity check to confirm that all quant.sf files exist
# Returns TRUE if every expected file is present. So each rowname will be what is named to your salmon output abundnace/count/readlength per sample, as well as where that is stored, so should not have ANY  mismatches.
# If FALSE, do not proceed until the file -> filename link is fixed and returns TRUE
all(file.exists(files))
head(files, n = 24)

# Step 5: Load and clean annotations for tx2gene dataframe that tximport requires
annotations <- read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) # 54484 annotated transcripts
data.table::uniqueN(annotations$locusName) # Number of unique transcripts, discards repeated rows per transcript #31948
tx2gene <- annotations[, c(3, 2)] # save annotations into df for tximport to summarize gene-level from transcript isoforms
str(tx2gene)
any(duplicated(tx2gene[, 1])) # if not false, means you have MULTIPLE column 1 rownames, which means raw annot is NOT clean

# Step 5 Optional: Getting additional metrics out such as no. transcript isoforms per gene
# now, if you're actually not sure whether your annotations files is set up correctly, run the below
transcript_fasta_data <- read.fasta(file = "P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.transcript.fa.gz", as.string = TRUE) # 54484 annotated transcripts
#
transcript_fasta_data_df <- do.call(rbind, lapply(transcript_fasta_data, function(seq_obj) {
  data.frame(
    transcriptName = attr(seq_obj, "name"),
    Vuannotation = attr(seq_obj, "Annot"),
    ft_ni_sequence = as.character(seq_obj),
    stringsAsFactors = FALSE
  )
})) # 54484 annotated transcripts
str(transcript_fasta_data_df)

duplicate_counts_df_counts <- as.data.frame(table(annotations$locusName)) |> # table() counts how many times each unique locusName appears in your annotations data frame.
  subset(Freq > 1) |> # keep only loci with >1 transcript
  setNames(c("locusName", "Transcript_Isoforms")) # gives no. isoforms per transcript if there are >1
# Ok, so let's summarize to gene-level, make sure it matches above number of 31948
annotations_GL <- annotations %>%
  dplyr::distinct(locusName, .keep_all = TRUE) # 31948 more accurate because some genetranscripts might not have .1 read. Some start with .2
annotations_GL <- left_join(annotations_GL, duplicate_counts_df_counts, by = "locusName", relationship = "one-to-one") # Could be useful to know also how many isoforms per gene
write.csv(annotations_GL,
          file = "P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info_GL.csv",
          row.names = FALSE
)


# Step 6: Tximport to summarize counts,TPM, avgereadlngth on a gene-level basis
#navigate to folders path
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/C1_C24_VuFUN_timecourse_2025/quantsC1_C24/")
txi.salmon.tsv = tximport(files,
                            type = "salmon",
                            tx2gene = tx2gene,
                            dropInfReps = TRUE,
                            countsFromAbundance = "no") # don't worry about this too much, if you are using DESeq2, it will do its own thing with count normalization. And you will use TPMs for gene abundance anyway.

final <- data.frame(txi.salmon.tsv) # save into dataframe. Has TPM, counts, avge read length all in one df. Will extract into single dataframes
str(final)
colnames(final)
str(files)

# Step 7: Load in relevant metadata for this experiment per sample, so as to match with colnames of saved tximport object
metadata <- readxl::read_xlsx("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/01.RawData/Copy of Sample Register.xlsx") # load meta data in
str(metadata)
metadata$Project <- factor(metadata$Project)
levels(metadata$Project)
vufunmetadata <- metadata %>%
  filter(Project == "vufun vs WT +/- N timecourse" ) # filter out for just vufun N tc exp experiment
colnames(vufunmetadata)


# REMINDER. At this step make sure that your metadata is 100% correct or you will be running incorrect DEGs and outputting incorrect data

# Step 8: Save Counts
header <- c("locusName", c(paste("counts", vufunmetadata$CombID,vufunmetadata$Rep, sep = "."))) # Create column headers by combining "counts" prefix with sample descriptions from metadata. Output will be "counts.Description1", "counts.Description2", etc.
header
counts <- data.frame(final[, grep("counts.", colnames(final))]) # Extract only the count columns from the final dataframe
colnames(counts) # Check current column names (may include an unwanted "countsFromAbundance" column)
counts <- cbind(locusName = rownames(counts), counts[, -25]) # counts[,-25] removes the 25th column which is likely the unwanted countsFromAbundance column. cbind() combines the new locusName column with the existing data
colnames(counts) # Verify the column names after adding locusName and removing the extra column
colnames(counts) <- header # Rename all columns using the cleaned header vector created earlier. This replaces the default/messy column names with the standardized format. And we know
head(counts)
str(counts)

write.csv(counts, file = "vufun_TC_Ntrt_rawcounts.csv", row.names = FALSE) # Save output

# Step 9: Save TPMs
header <- c("locusName", c(paste("TPM", vufunmetadata$CombID,vufunmetadata$Rep, sep = ".")))
header
colnames(final)
TPM <- data.frame(final[, grep("abundance", colnames(final))])
colnames(TPM)
TPM <- cbind(locusName = rownames(TPM), TPM)
colnames(TPM)
colnames(TPM) <- header
View(TPM)
write.csv(TPM, file = "vufun_TC_Ntrt_TPM.csv", row.names = F, quote = F)
TPM <- read.csv("vufun_TC_Ntrt_TPM.csv")
# Step 10: Generate some stats for TPM
TPM_long <- pivot_longer(TPM,
                         cols = -locusName,
                         names_to = "Sample",
                         values_to = "TPM"
)
head(TPM_long$Sample)
colnames(TPM)
TPMheaders <- colnames(TPM[, -1])

vufunmetadata_tpm <- data.frame(Sample = TPMheaders, vufunmetadata[, c("Treatment","Day","Rep","Identifying Allele","Description","CombID")])
vufunmetadata_tpm <- left_join(TPM_long, vufunmetadata_tpm, by = "Sample")
View(vufunmetadata_tpm)
colnames(vufunmetadata_tpm)
write.csv(TPM, file = "vufunTPM_longformat_metadata.csv", row.names = F, quote = F)

TPM_mean_comb <- plyr::ddply(
  vufunmetadata_tpm,
  c("locusName", "Identifying.Allele","Day","Treatment","CombID"),
  summarise,
  mean = mean(TPM, na.rm = TRUE),
  sd   = sd(TPM, na.rm = TRUE),
  n    = sum(!is.na(TPM)),
  sem  = sd(TPM, na.rm = TRUE) / sqrt(sum(!is.na(TPM)))
) # Generates mean, sd, sem, n for each gene per sample.

TPM_mean_comb <- readxl::read_xlsx("vufun_TPMstats.xlsx")

writexl::write_xlsx(
  TPM_mean_comb,
  "vufun_TPMstats.xlsx"
)


TPM_mean_comb <- readxl::read_xlsx("vufun_TPMstats.xlsx")


# Step 11: Last component needed for DESeq analysis, which is to provide experimental setup details in relation to experiment treatment factors
filenames <- file.path(folders, "quant.sf") # we created folders earlier, object containing folder names which hold .sf files
filenames # so now we have filepaths that are correct for each sample
sampleNames <- sub("_quant/quant.sf", "", filenames)
sampleNames # just one of many ways to get the sampleNames(literal) sample names)
sampleNames

colnames(vufunmetadata_tpm)
sampleTable <- vufunmetadata[, c("CombID")]
sampleTable
sampleTable <- data.frame(sampleTable)
sampleTable$`RNA Sample Code` <- colnames(txi.salmon.tsv$counts)
sampleTable_full <- sampleTable %>%
  left_join(vufunmetadata[,c("Treatment","Day","Identifying Allele","RNA Sample Code")],by="RNA Sample Code") %>%
  column_to_rownames("RNA Sample Code")
all(rownames(sampleTable_full) == colnames(txi.salmon.tsv$counts))
sampleTable_full # #  sampleTable, which contains a data frame of sample ID assigned per sample when building txi.salmon.tsv object in the rownames, while the first (and only) columns consist of the experiment treatment factor
sampleTable_full[] <- lapply(sampleTable_full, as.factor) #make everything a factor
str(sampleTable_full)
sampleTable_full <- dplyr::rename(sampleTable_full, `Genotype` = "Identifying Allele")



# Step 12.1: Carry out DEG analysis
# Use CombID as the single factor (it already encodes all condition info)
colnames(sampleTable_full)
ddsfullexpfactors <- DESeqDataSetFromTximport(txi.salmon.tsv, 
                                sampleTable_full, 
                                ~ Treatment + Day + Genotype + Genotype:Treatment+Day:Treatment)

ddscompl <- DESeq(ddsfullexpfactors)
resultsNames(ddscompl)

rescomp <- results(ddscompl, name =c("Genotype_WT..IT86D.1010._vs_vufun2","Genotype_WT..IT86D.1010._vs_vufun2"),alpha=0.05)
)

summary(rescomp)



#Step 12 workaround: Carry out DEG analysis or pre-grouped samples.

dds <- DESeqDataSetFromTximport(txi.salmon.tsv, 
                                              sampleTable_full, 
                                              ~ CombID)

dds <- DESeq(dds)
resultsNames(dds)
levels(dds$CombID)






# Step 12.1 test 1. Can I do results extraction despite not specifying ref level when doing DE analysis?
resD1 <- results(dds, contrast=list(c("Treatment_No.Nitrate_vs_Nitrate..","WT_D3_Nplus"),alpha=0.05))
summary(resD1)


# Get all levels
levels_CombID <- levels(sampleTable_full$CombID)

# Generate all pairwise combinations
combos <- combn(levels_CombID, 2, simplify = FALSE)
combos

?results

# Run all contrasts and store results
results_list <- lapply(combos, function(pair) {
  res <- results(dds, contrast = c("CombID", pair[1], pair[2]))
  res <- lfcShrink(dds, contrast = c("CombID", pair[1], pair[2]), res = res, type = "ashr")
  return(res)
})

results_list <- lapply(combos, function(pair) {
  res <- results(dds, contrast = c("CombID", pair[1], pair[2]),alpha = 0.05)
  return(res)
})

summary(results_list[[28]])

















# Step 12: Carry out DEG analysis, with WT13 AND WT15 as a reference.
# Use CombID as the single factor (it already encodes all condition info)
dds <- DESeqDataSetFromTximport(txi.salmon.tsv, 
                                sampleTable_full, 
                                ~CombID)

dds <- DESeq(dds, contrasts=c)


resD1 <- results(dds, contrast=c("CombID","WT_D3_Nplus","vufun2_D3_Nminus"))


summary(resD1)

levels(sampleTable_full$CombID)

# Step 13: Global Transcriptomic response
vsd <- vst(dds, blind = TRUE) # apply variance stabilizing transformation (VSD) to reduce dependence of variance on mean experession (more homoscedastic). Corrects for natural high variance at high experssion levels, low variance at low experience levels, making data visualization (PCA/Heatmaps) clearer.
plotPCA(vsd, intgroup = c("CombID"), returnData = FALSE, pcsToUse = c(1,2)) # raw PCA. Good, but want more control over certain variables
pca_data <- plotPCA(
  vsd,
  intgroup = "CombID",
  pcsToUse = c(1,2),
  returnData = TRUE
) # build a better PCA, so first save pca data into dataframe for PC1 nad PC2

pca_data


pca_data <- dplyr::rename(pca_data, `RNA Sample Code` = name) # need to to do this to merge some metadta backinto this PCA summary of PC1 and PC2
pca_data
colnames(pca_data)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
levels(pca_data$Identifying.Allele)
# subset(pca_data, !Identifying.Allele %in% c("WT13", "WT15"))

PCA3 <- ggplot(pca_data, aes(PC1, PC3, color = Day ,shape = Treatment)) + # now create custom PCA
  geom_point(
    size = 2.8,
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = `Identifying.Allele`), # assign labels for each dotpoint
    size = 4, # text size
    max.overlaps = 20,
    box.padding = 0.3, # toggle higher to make line visible
    point.padding = 0.3, # line to dot space ratio. Higher = more distance from line to dot
    segment.color = "grey70" # line to text color
  ) +
  theme_bw() +
  labs(
    title = "PCA of vufun timecourse experiment: PC1 vs PC3",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Allele"
  ) +
  scale_color_brewer(palette = "Set2") +
  # theme_classic(base_size = 12)
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, vjust = 0.5, face = "bold"),
    strip.text.x = element_text(size = 8, face = "bold", margin = margin(0.1, 0, 0.1, 0, "cm")),
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

PCA1 <- ggplot(pca_data, aes(PC1, PC2, color = Day ,shape = Treatment)) + # now create custom PCA
  geom_point(
    size = 2.8,
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = `Identifying.Allele`), # assign labels for each dotpoint
    size = 4, # text size
    max.overlaps = 20,
    box.padding = 0.3, # toggle higher to make line visible
    point.padding = 0.3, # line to dot space ratio. Higher = more distance from line to dot
    segment.color = "grey70" # line to text color
  ) +
  theme_bw() +
  labs(
    title = "PCA of vufun timecourse experiment: PC1 vs PC2",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC3 (", percent_var[2], "% variance)"),
    color = "Allele"
  ) +
  scale_color_brewer(palette = "Set2") +
  # theme_classic(base_size = 12)
  theme(
    axis.title.x = element_text(color = "black", size = 12, face = "bold"),
    axis.title.y = element_text(color = "black", size = 12, face = "bold"),
    axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10, vjust = 0.5, face = "bold"),
    strip.text.x = element_text(size = 8, face = "bold", margin = margin(0.1, 0, 0.1, 0, "cm")),
    # panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank()
  )

library(gridExtra)
combinedPCA <- grid.arrange(PCA3,PCA1,nrow=1)
ggsave(combinedPCA, file="ZBF_pca_comps.pdf",height=25,width=45, units="cm")

# Step 14: Get normalized counts for each gene per sample. Can use for gene abundance comparisons
resnormcounts <- counts(dds, normalized = T)
resnormcounts <- as.data.frame(resnormcounts)
resnormcounts
head(vufunmetadata$CombID, n = 24)
colnames(resnormcounts) <- vufunmetadata$CombID
colnames(resnormcounts)
header <-  c(c(paste(vufunmetadata$CombID,vufunmetadata$Rep, sep = ".R"))) # Create column headers by combining "counts" prefix with sample descriptions from metadata. Output will be "counts.Description1", "counts.Description2", etc.
colnames(resnormcounts) <- header
resnormcounts <- resnormcounts %>%
  rownames_to_column("locusName")
resnormcounts
write.csv(resnormcounts, file = "vufun_norm_counts_results.csv", row.names = FALSE)

#Ok now what if we want to plot some shizzzz
RRgenes <- readxl::read_xlsx("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab/Panfaba work/RobRoy_Vuadded.xlsx")
colnames(RRgenes)
RRgenes <- dplyr::rename(RRgenes, locusName = Vunguiculata)

vuresnormcountslong <- pivot_longer(resnormcounts,
                                  cols = -locusName,
                                  names_to = "Sample",
                                  values_to = "normcounts"
)

vunormcountheaders <- colnames(resnormcounts[,-1])
vunormcountheaders
vumetadata_counts <- data.frame(Sample = vunormcountheaders, vufunmetadata[, c("RNA Sample Code","Identifying Allele","Treatment","Day","CombID","Rep")])

vumetadata_counts <- left_join(vuresnormcountslong, vumetadata_counts, by = "Sample")
vumetadata_counts$CombID <- factor(vumetadata_counts$CombID,levels=c("vufun2_D1_Nminus","vufun2_D1_Nplus","WT_D1_Nminus","WT_D1_Nplus","vufun2_D3_Nminus","vufun2_D3_Nplus","WT_D3_Nminus","WT_D3_Nplus" ) )

vutoplot <- vumetadata_counts %>%
  dplyr::filter(locusName %in% c("Vigun04g109900", "Vigun01g085000", "Vigun11g064200", "Vigun05g197300", "Vigun02g057800", "Vigun03g206900"))
vutoplot <- left_join(vutoplot, RRgenes[, c("Gene Symbol", "locusName")], by = "locusName")
vutoplot$Rep <- factor(vutoplot$Rep)
ggplot(vutoplot, aes(x = CombID, y = normcounts)) +
  geom_boxplot(fill = NA, color = "gray90", alpha = 0.05, outlier.shape = NA) +
  geom_point(aes(color = Rep)) +
  # geom_bar(stat="")+
  #  geom_text(aes(label=Rep))+
  facet_wrap(~`locusName`, nrow = 2, scales = "free") +
  labs(
    title = "Normalized counts of select genes\nfor ZBF lines", y = "normalized counts", x = "ZBF line"
  ) +
  theme_bw() +
  scale_color_manual(values = c("1" = "#4477AA", "2" = "#EE6677", "3" = "#228833", "4" = "#CCBB44")) +
  theme(
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = -40, hjust = 0, face = "bold"),
    axis.text.y = element_text(size = 15, vjust = 0.5, face = "bold"),
    strip.text.x = element_text(size = 15, face = "bold", margin = margin(0.1, 0, 0.1, 0, "cm")),
    strip.background = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    theme(
      legend.key.size = unit(15, "cm"), # change legend key size
      legend.key.height = unit(15, "cm"), # change legend key height
      legend.key.width = unit(15, "cm"), # change legend key width
      legend.title = element_text(size = 14), # change legend title font size
      legend.text = element_text(size = 10)
    )
  )


#Heatmaps to plot for lfc


# Step 15: Get DEG results (fyi log2foldchange > 0.0585 = fold change of 1.5)
# extra annotations table for nodulation genes
RRgenes <- readxl::read_xlsx("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab/Panfaba work/RobRoy_Vuadded.xlsx")
colnames(RRgenes)
RRgenes <- dplyr::rename(RRgenes, locusName = Vunguiculata)


# D1, vufun +N
resultsNames(dds) # these are all the generated comparisons.

resD1_vufun <- results(dds, contrast=c("CombID","vufun2_D1_Nplus","vufun2_D1_Nminus"),alpha=0.05)
summary(resD1_vufun)
OEresD1_vufun <- lfcShrink(dds, contrast = c("CombID", "vufun2_D1_Nplus", "vufun2_D1_Nminus"), type = "ashr", res = resD1_vufun)
resD1_vufundf <- as.data.frame(resD1_vufun)
resD1_vufundf <- resD1_vufundf %>%
  rownames_to_column("locusName")
OEresD1_vufundf <- as.data.frame(OEresD1_vufun)
OEresD1_vufundf <- OEresD1_vufundf %>%
  rownames_to_column("locusName")
OEresD1_vufundf <- OEresD1_vufundf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD1_vufun <- left_join(resD1_vufundf, OEresD1_vufundf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD1_vufun <- left_join(finresD1_vufun, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD1_vufun <- left_join(finresD1_vufun, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD1_vufun <- finresD1_vufun[order(finresD1_vufun$padj), ]
finresD1_vufun$DEG <- ifelse(finresD1_vufun$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD1_vufun, path = "results_pvals_D1_vufun_plusN_vs_minusN.xlsx")
getwd()

# D1, WT +N
resultsNames(dds) # these are all the generated comparisons.

resD1_WT <- results(dds, contrast=c("CombID","WT_D1_Nplus","WT_D1_Nminus"),alpha=0.05)
summary(resD1_WT)
OEresD1_WT <- lfcShrink(dds, contrast = c("CombID", "WT_D1_Nplus", "WT_D1_Nminus"), type = "ashr", res = resD1_WT)
summary(OEresD1_WT)
resD1_WTdf <- as.data.frame(resD1_WT)
resD1_WTdf <- resD1_WTdf %>%
  rownames_to_column("locusName")
OEresD1_WTdf <- as.data.frame(OEresD1_WT)
OEresD1_WTdf <- OEresD1_WTdf %>%
  rownames_to_column("locusName")
OEresD1_WTdf <- OEresD1_WTdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD1_WT <- left_join(resD1_WTdf, OEresD1_WTdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD1_WT <- left_join(finresD1_WT, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD1_WT <- left_join(finresD1_WT, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD1_WT <- finresD1_WT[order(finresD1_WT$padj), ]
finresD1_WT$DEG <- ifelse(finresD1_WT$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD1_WT, path = "results_pvals_D1_WT_plusN_vs_minusN.xlsx")



# D3, vufun +N
resultsNames(dds) # these are all the generated comparisons.

resD3_vufun <- results(dds, contrast=c("CombID","vufun2_D3_Nplus","vufun2_D3_Nminus"),alpha=0.05)
summary(resD3_vufun)
OEresD3_vufun <- lfcShrink(dds, contrast = c("CombID", "vufun2_D3_Nplus", "vufun2_D3_Nminus"), type = "ashr", res = resD3_vufun)
summary(OEresD3_vufun)
resD3_vufundf <- as.data.frame(resD3_vufun)
resD3_vufundf <- resD3_vufundf %>%
  rownames_to_column("locusName")
OEresD3_vufundf <- as.data.frame(OEresD3_vufun)
OEresD3_vufundf <- OEresD3_vufundf %>%
  rownames_to_column("locusName")
OEresD3_vufundf <- OEresD3_vufundf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD3_vufun <- left_join(resD3_vufundf, OEresD3_vufundf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD3_vufun <- left_join(finresD3_vufun, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD3_vufun <- left_join(finresD3_vufun, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD3_vufun <- finresD3_vufun[order(finresD3_vufun$padj), ]
finresD3_vufun$DEG <- ifelse(finresD3_vufun$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD3_vufun, path = "results_pvals_D3_vufun_plusN_vs_minusN.xlsx")
getwd()




# D3, WT +N
resultsNames(dds) # these are all the generated comparisons.

resD3_WT <- results(dds, contrast=c("CombID","WT_D3_Nplus","WT_D3_Nminus"),alpha=0.05)
summary(resD3_WT)
OEresD3_WT <- lfcShrink(dds, contrast = c("CombID", "WT_D3_Nplus", "WT_D3_Nminus"), type = "ashr", res = resD3_WT)
summary(OEresD3_WT)
resD3_WTdf <- as.data.frame(resD3_WT)
resD3_WTdf <- resD3_WTdf %>%
  rownames_to_column("locusName")
OEresD3_WTdf <- as.data.frame(OEresD3_WT)
OEresD3_WTdf <- OEresD3_WTdf %>%
  rownames_to_column("locusName")
OEresD3_WTdf <- OEresD3_WTdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD3_WT <- left_join(resD3_WTdf, OEresD3_WTdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD3_WT <- left_join(finresD3_WT, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD3_WT <- left_join(finresD3_WT, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD3_WT <- finresD3_WT[order(finresD3_WT$padj), ]
finresD3_WT$DEG <- ifelse(finresD3_WT$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD3_WT, path = "results_pvals_D3_WT_plusN_vs_minusN.xlsx")


# MA plots
dev.off()
pdf(file = "vufun_singleGeno_singleTP_plusN.pdf", height = 8, width = 12)
par(mfrow = c(3, 4))
plotMA(resD1_vufun, main = "D1,vufun +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_vufun, main = "D1,vufun +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD1_WT, main = "D1,WT +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_WT, main = "D1,WT +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_vufun, main = "D3,vufun +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_vufun, main = "D3,vufun +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_WT, main = "D3,WT +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_WT, main = "D3,WT +N:\nWith Shrinkage", ylim = c(-4, 4))

dev.off()


#Step 16 Collect all dataframes for single genotype, per day, Treatment effects
# Ok let's start pulling shit in
SingleGenoSingleDay_plusN <- list(
  finresD1_vufun %>%
    select(locusName,
           lfc_D1_vufun_plusN = log2FoldChange,
           lfcs_D1_vufun_plusN = l2FCshrunk,
           padj_D1_vufun_plusN = padj
    ),
  finresD1_WT %>%
    select(locusName,
           lfc_D1_WT_plusN = log2FoldChange,
           lfcs_D1_WT_plusN = l2FCshrunk,
           padj_D1_WT = padj
    ),
  finresD3_vufun %>%
    select(locusName,
           lfc_D3_vufun_plusN = log2FoldChange,
           lfcs_D3_vufun_plusN = l2FCshrunk,
           padj_D3_vufun_plusN = padj
    ),
  finresD3_WT %>%
    select(locusName,
           lfc_D3_WT_plusN = log2FoldChange,
           lfcs_D3_WT_plusN = l2FCshrunk,
           padj_D3_WT_plusN = padj
    )
)


library(purrr)

SingleGenoSingleDay_plusN_merged <- purrr::reduce(SingleGenoSingleDay_plusN, inner_join, by = "locusName")
SingleGenoSingleDay_plusN_merged <- SingleGenoSingleDay_plusN_merged %>%
  dplyr::distinct(locusName, .keep_all = TRUE) # 31337

writexl::write_xlsx(SingleGenoSingleDay_plusN_merged, path = "lfc_pval_vufun_singleGeno_singleTP_plusN_combined.xlsx")

#for padj<0.05
colnames(SingleGenoSingleDay_plusN_merged)
SingleGenoSingleDay_plusN_merged_DEG <- SingleGenoSingleDay_plusN_merged %>%
  filter(padj_D1_vufun_plusN < 0.05 | padj_D1_WT < 0.05 | padj_D3_vufun_plusN < 0.05 | padj_D3_WT_plusN < 0.05)

writexl::write_xlsx(SingleGenoSingleDay_plusN_merged_DEG, path = "lfc_pval_vufun_singleGeno_singleTP_plusN_combined_DEG.xlsx")

#Step16 HM
# Create base dataframes to work with
SG_SD_plusN_rronly <- left_join(RRgenes[, c("Gene Symbol", "Function", "locusName")],
                                SingleGenoSingleDay_plusN_merged,
                         by = "locusName"
) # All Rob Roy genes #343

SG_SD_plusN_noNA <- SG_SD_plusN_rronly %>%
  drop_na(starts_with("lfcs_")) # All DEG Rob Roy genes #321



# Make Groupings
row_orderRRall <- data.frame(
  locusName = SG_SD_plusN_noNA$locusName,
  Grouping = SG_SD_plusN_noNA$Function,
  geneSym = SG_SD_plusN_noNA$`Gene Symbol`
)

row_orderRRall$Grouping <- factor(row_orderRRall$Grouping,
                                  levels = c(
                                    "Early Signalling",
                                    "Host Range Restriction",
                                    "Rhizobial Infection",
                                    "Nodule Organogenesis",
                                    "Autoregulation of Nodule Number",
                                    "Bacterial Maturation",
                                    "Symbiosome Maturation",
                                    "Nodulation Metabolism and Transport",
                                    "Senescence",
                                    "Defense",
                                    "DAP-seq target"
                                  )
)

row_orderRRall <- row_orderRRall %>%
  arrange(Grouping)

SG_SD_plusN_noNA <- left_join(row_orderRRall[, c("locusName", "Grouping")], SG_SD_plusN_noNA, by = "locusName")
colnames(SG_SD_plusN_noNA)
SG_SD_plusN_noNA <- SG_SD_plusN_noNA[, -c(2, 3)]

# Make matrix
SG_SD_plusN_noNA_rr_lfc_mat <- SG_SD_plusN_noNA %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

SingleGenoSingleDay_plusN_merged_DEG_mat <- SingleGenoSingleDay_plusN_merged_DEG %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

padj_mat <- SG_SD_plusN_noNA %>%
  select(locusName, starts_with("padj_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

colnames(padj_mat) <- sub("padj_", "", colnames(padj_mat))

# Create empty character matrix
sig_mat <- matrix("", nrow = nrow(padj_mat), ncol = ncol(padj_mat))
rownames(sig_mat) <- rownames(padj_mat)
colnames(sig_mat) <- colnames(padj_mat)

# Fill with significance levels
sig_mat[padj_mat < 0.05] <- "*"
sig_mat[padj_mat < 0.01] <- "**"
sig_mat[padj_mat < 0.001] <- "***"


# Create column order
colnames(SG_SD_plusN_noNA_rr_lfc_mat)
col_order <- c("lfcs_DBL21", "lfcs_DBL19", "lfcs_HA11", "lfcs_HA3", "lfcs_ADMD17", "lfcs_WT13")

library(circlize)

# set up colors
colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#2166AC", "white", "#B2182B")
)

colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#3B4CC0", "#F7F7F7", "#B40426")
) # best for white backgrounds

colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#0072B2", "#F0E442", "#D55E00")
)
colfun <- colorRamp2(
  c(-4, -1, 0, 1, 4),
  c("#313695", "#74ADD1", "#F7F7F7", "#FDAE61", "#A50026")
)


colnames(SG_SD_plusN_noNA_rr_lfc_mat)
col_order <- c("lfcs_D1_vufun_plusN", "lfcs_D1_WT_plusN",    "lfcs_D3_vufun_plusN", "lfcs_D3_WT_plusN" )
column_split_factor <- factor(col_order, levels = col_order)
column_split_factor
# Heatmap

dev.off()
pdf("SG_SD_plusN_RR.pdf", height = 12)
Heatmap(SG_SD_plusN_noNA_rr_lfc_mat,
              col = colfun,
              use_raster = TRUE,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              row_labels = row_orderRRall$geneSym,
              # Show row names (from the matrix) on the right
              show_row_names = TRUE,
              row_names_side = "right",
              # put them on the right
              cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
              # cluster_rows=TRUE,
              # clustering_method_columns="ward.D",
              # clustering_method_rows = "ward.D",
              row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
              row_split = factor(row_orderRRall$Grouping),
              row_gap = unit(2, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (sig_mat[i, j] != "") {
                  grid.text(
                    sig_mat[i, j],
                    x, y,
                    gp = gpar(fontsize = 2, col = "red")
                  )
                }
              },
              # row_km = 12,
              # column_km=9,
              # column_km_repeats=100,
              # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
              row_title_gp = gpar(font = 2, fontsize = 4),
              row_title_rot = 0,
              row_names_gp = gpar(fontsize = 3),
              cluster_column_slices = FALSE,
              column_split = column_split_factor, # <- visually separates minusP vs plusP
              column_gap = unit(2, "mm"),
              column_title = "vufun TC +N exp: Within genotype differences per timepoint to +N",
              # column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
              heatmap_legend_param = list(
                title = "shrunklog2FC",
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 8),
                at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
                labels = seq(-4, 4, by = 2) # what to print next to these ticke
              )
)
dev.off()

Heatmap(SingleGenoSingleDay_plusN_merged_DEG_mat,
        col = colfun,
        use_raster = TRUE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        #row_labels = row_orderRRall$geneSym,
        # Show row names (from the matrix) on the right
        show_row_names = FALSE,
        #row_names_side = "right",
        # put them on the right
        #cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
        # cluster_rows=TRUE,
        # clustering_method_columns="ward.D",
        # clustering_method_rows = "ward.D",
       # row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
      #  row_split = factor(row_orderRRall$Grouping),
      #  row_gap = unit(2, "mm"),
        row_km = 22,
        # column_km=9,
        # column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        row_title_gp = gpar(font = 2, fontsize = 6),
        row_title_rot = 0,
        row_names_gp = gpar(fontsize = 6),
        cluster_column_slices = FALSE,
        column_split = column_split_factor, # <- visually separates minusP vs plusP
        column_gap = unit(2, "mm"),
       # column_title = "vufun TC +N exp: Within genotype differences per timepoint to +N:While Transcriptome",
        # column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(
          title = "shrunklog2FC",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
          labels = seq(-4, 4, by = 2) # what to print next to these ticke
        )
)


# What if we want DAP-seq list?
DAPseqVu <- read.csv("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab/CowpeaDAPseq.csv", header = TRUE)
SingleGenoSingleDay_plusN_merged_DEG_annot <- SingleGenoSingleDay_plusN_merged_DEG %>%
  left_join(annotations_GL[, c("locusName", "Best.hit.arabi.name", "Best.hit.arabi.defline", "Best.hit.rice.name", "Best.hit.rice.defline", "Transcript_Isoforms")], by = "locusName") %>%
  left_join(RRgenes, by = "locusName") %>%
  left_join(DAPseqVu, by = "locusName") %>%
  left_join(resnormcounts, by = "locusName")

writexl::write_xlsx(SingleGenoSingleDay_plusN_merged_DEG_annot, path = "vufun_SG_SD_plusN_intersects_fullannot.xlsx")
SingleGenoSingleDay_plusN_merged_DEG


#Going back to Step 15 for Between Genotype differences
# D1, vufun vs WT, -N
resultsNames(dds) # these are all the generated comparisons.

resD1_minusN <- results(dds, contrast=c("CombID","vufun2_D1_Nminus","WT_D1_Nminus"),alpha=0.05)
summary(resD1_minusN)
OEresD1_minusN <- lfcShrink(dds, contrast = c("CombID", "vufun2_D1_Nminus", "WT_D1_Nminus"), type = "ashr", res = resD1_minusN)
resD1_minusNdf <- as.data.frame(resD1_minusN)
resD1_minusNdf <- resD1_minusNdf %>%
  rownames_to_column("locusName")
OEresD1_minusNdf <- as.data.frame(OEresD1_minusN)
OEresD1_minusNdf <- OEresD1_minusNdf %>%
  rownames_to_column("locusName")
OEresD1_minusNdf <- OEresD1_minusNdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD1_minusN <- left_join(resD1_minusNdf, OEresD1_minusNdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD1_minusN <- left_join(finresD1_minusN, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD1_minusN <- left_join(finresD1_minusN, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD1_minusN <- finresD1_minusN[order(finresD1_minusN$padj), ]
finresD1_minusN$DEG <- ifelse(finresD1_minusN$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD1_minusN, path = "results_pvals_D1_minusN_vufun_vs_WT.xlsx")
getwd()

# D1, vufun vs WT, +N
resultsNames(dds) # these are all the generated comparisons.

resD1_plusN <- results(dds, contrast=c("CombID","vufun2_D1_Nplus","WT_D1_Nplus"),alpha=0.05)
summary(resD1_plusN)
OEresD1_plusN <- lfcShrink(dds, contrast = c("CombID", "vufun2_D1_Nplus", "WT_D1_Nplus"), type = "ashr", res = resD1_plusN)
resD1_plusNdf <- as.data.frame(resD1_plusN)
resD1_plusNdf <- resD1_plusNdf %>%
  rownames_to_column("locusName")
OEresD1_plusNdf <- as.data.frame(OEresD1_plusN)
OEresD1_plusNdf <- OEresD1_plusNdf %>%
  rownames_to_column("locusName")
OEresD1_plusNdf <- OEresD1_plusNdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD1_plusN <- left_join(resD1_plusNdf, OEresD1_plusNdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD1_plusN <- left_join(finresD1_plusN, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD1_plusN <- left_join(finresD1_plusN, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD1_plusN <- finresD1_plusN[order(finresD1_plusN$padj), ]
finresD1_plusN$DEG <- ifelse(finresD1_plusN$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD1_plusN, path = "results_pvals_D1_plussN_vufun_vs_WT.xlsx")
getwd()



# D3, vufun vs WT, -N
resultsNames(dds) # these are all the generated comparisons.

resD3_minusN <- results(dds, contrast=c("CombID","vufun2_D3_Nminus","WT_D3_Nminus"),alpha=0.05)
summary(resD3_minusN)
OEresD3_minusN <- lfcShrink(dds, contrast = c("CombID", "vufun2_D3_Nminus", "WT_D3_Nminus"), type = "ashr", res = resD3_minusN)
resD3_minusNdf <- as.data.frame(resD3_minusN)
resD3_minusNdf <- resD3_minusNdf %>%
  rownames_to_column("locusName")
OEresD3_minusNdf <- as.data.frame(OEresD3_minusN)
OEresD3_minusNdf <- OEresD3_minusNdf %>%
  rownames_to_column("locusName")
OEresD3_minusNdf <- OEresD3_minusNdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD3_minusN <- left_join(resD3_minusNdf, OEresD3_minusNdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD3_minusN <- left_join(finresD3_minusN, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD3_minusN <- left_join(finresD3_minusN, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD3_minusN <- finresD3_minusN[order(finresD3_minusN$padj), ]
finresD3_minusN$DEG <- ifelse(finresD3_minusN$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD3_minusN, path = "results_pvals_D3_minussN_vufun_vs_WT.xlsx")
getwd()




# D3, vufun vs WT, +N
resultsNames(dds) # these are all the generated comparisons.

resD3_plusN <- results(dds, contrast=c("CombID","vufun2_D3_Nplus","WT_D3_Nplus"),alpha=0.05)
summary(resD3_plusN)
OEresD3_plusN <- lfcShrink(dds, contrast = c("CombID", "vufun2_D3_Nplus", "WT_D3_Nplus"), type = "ashr", res = resD3_plusN)
resD3_plusNdf <- as.data.frame(resD3_plusN)
resD3_plusNdf <- resD3_plusNdf %>%
  rownames_to_column("locusName")
OEresD3_plusNdf <- as.data.frame(OEresD3_plusN)
OEresD3_plusNdf <- OEresD3_plusNdf %>%
  rownames_to_column("locusName")
OEresD3_plusNdf <- OEresD3_plusNdf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD3_plusN <- left_join(resD3_plusNdf, OEresD3_plusNdf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD3_plusN <- left_join(finresD3_plusN, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD3_plusN <- left_join(finresD3_plusN, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD3_plusN <- finresD3_plusN[order(finresD3_plusN$padj), ]
finresD3_plusN$DEG <- ifelse(finresD3_plusN$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD3_plusN, path = "results_pvals_D3_plussN_vufun_vs_WT.xlsx")
getwd()


# MA plots
dev.off()
pdf(file = "vufun_betweengenotypes_singleTP_singleTrt.pdf", height = 8, width = 12)
par(mfrow = c(3, 4))
plotMA(resD1_minusN, main = "D1,WT vs vufun2 -N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_minusN, main = "D1,WT vs vufun2 -N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD1_plusN, main = "D1,WT vs vufun2 +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_plusN, main = "D1,WT vs vufun2 +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_minusN, main = "D3,WT vs vufun2 -N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_minusN, main = "D3,WT vs vufun2 -N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_plusN, main = "D3,WT vs vufun2 +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_plusN, main = "D3,WT vs vufun2 +N:\nWith Shrinkage", ylim = c(-4, 4))

dev.off()

#Step 16:revising to consolidate 
SingleNSingleDay_vufunvsWT <- list(
  finresD1_minusN %>%
    select(locusName,
           lfc_finresD1_minusN = log2FoldChange,
           lfcs_finresD1_minusN = l2FCshrunk,
           padj_finresD1_minusN = padj
    ),
  finresD1_plusN %>%
    select(locusName,
           lfc_finresD1_plusN = log2FoldChange,
           lfcs_finresD1_plusN = l2FCshrunk,
           padj_finresD1_plusN = padj
    ),
  finresD3_minusN %>%
    select(locusName,
           lfc_finresD3_minusN = log2FoldChange,
           lfcs_finresD3_minusN = l2FCshrunk,
           padj_finresD3_minusN = padj
    ),
  finresD3_plusN %>%
    select(locusName,
           lfc_finresD3_plusN = log2FoldChange,
           lfcs_finresD3_plusN = l2FCshrunk,
           padj_finresD3_plusN = padj
    )
)
library(purrr)

SingleNSingleDay_vufunvsWT_merged <- purrr::reduce(SingleNSingleDay_vufunvsWT, inner_join, by = "locusName")
SingleNSingleDay_vufunvsWT_merged <- SingleNSingleDay_vufunvsWT_merged %>%
  dplyr::distinct(locusName, .keep_all = TRUE) # 31337

writexl::write_xlsx(SingleNSingleDay_vufunvsWT_merged, path = "lfc_pval_vufun_VS_WT_combined.xlsx")

#for padj<0.05
colnames(SingleNSingleDay_vufunvsWT_merged)
SingleNSingleDay_vufunvsWT_merged_DEG <- SingleNSingleDay_vufunvsWT_merged %>%
  filter(padj_finresD1_minusN < 0.05 | padj_finresD1_plusN < 0.05 | padj_finresD3_minusN < 0.05 | padj_finresD3_plusN < 0.05)

writexl::write_xlsx(SingleNSingleDay_vufunvsWT_merged_DEG, path = "lfc_pval_vufun_VS_WT_combined_DEG.xlsx")


#Step16 HM
# Create base dataframes to work with
ST_SD_vufun2vsWT_rronly <- left_join(RRgenes[, c("Gene Symbol", "Function", "locusName")],
                                     SingleNSingleDay_vufunvsWT_merged_DEG,
                                by = "locusName"
) # All Rob Roy genes #343
colnames(ST_SD_vufun2vsWT_rronly)

ST_SD_vufun2vsWT_rronly_no_NA <- ST_SD_vufun2vsWT_rronly %>%
  drop_na(starts_with("lfcs_")) # All DEG Rob Roy genes #321

colnames(ST_SD_vufun2vsWT_rronly_no_NA)

# Make Groupings
row_orderRRall <- data.frame(
  locusName = ST_SD_vufun2vsWT_rronly_no_NA$locusName,
  Grouping = ST_SD_vufun2vsWT_rronly_no_NA$Function,
  geneSym = ST_SD_vufun2vsWT_rronly_no_NA$`Gene Symbol`
)

row_orderRRall$Grouping <- factor(row_orderRRall$Grouping,
                                  levels = c(
                                    "Early Signalling",
                                    "Host Range Restriction",
                                    "Rhizobial Infection",
                                    "Nodule Organogenesis",
                                    "Autoregulation of Nodule Number",
                                    "Bacterial Maturation",
                                    "Symbiosome Maturation",
                                    "Nodulation Metabolism and Transport",
                                    "Senescence",
                                    "Defense",
                                    "DAP-seq target"
                                  )
)

row_orderRRall <- row_orderRRall %>%
  arrange(Grouping)

ST_SD_vufun2vsWT_rronly_no_NA <- left_join(row_orderRRall[, c("locusName", "Grouping")], ST_SD_vufun2vsWT_rronly_no_NA, by = "locusName")
View(ST_SD_vufun2vsWT_rronly_no_NA)
SG_SD_plusN_noNA <- SG_SD_plusN_noNA[, -c(2, 3)]

# Make matrix
ST_SD_vufun2vsWT_rronly_no_NA_rr_lfc_mat <- ST_SD_vufun2vsWT_rronly_no_NA %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

ST_SD_vufun2vsWT_merged_DEG_mat <- SingleNSingleDay_vufunvsWT_merged_DEG %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

vufunvsWT_padj_mat <- ST_SD_vufun2vsWT_rronly_no_NA %>%
  select(locusName, starts_with("padj_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

colnames(vufunvsWT_padj_mat) <- sub("padj_", "", colnames(vufunvsWT_padj_mat))

# Create empty character matrix
vufunvsWT_sig_mat <- matrix("", nrow = nrow(vufunvsWT_padj_mat), ncol = ncol(vufunvsWT_padj_mat))
rownames(vufunvsWT_sig_mat) <- rownames(vufunvsWT_padj_mat)
colnames(vufunvsWT_sig_mat) <- colnames(vufunvsWT_padj_mat)

# Fill with significance levels
vufunvsWT_sig_mat[vufunvsWT_padj_mat < 0.05] <- "*"
vufunvsWT_sig_mat[vufunvsWT_padj_mat < 0.01] <- "**"
vufunvsWT_sig_mat[vufunvsWT_padj_mat < 0.001] <- "***"


# Create column order
colnames(ST_SD_vufun2vsWT_rronly_no_NA_rr_lfc_mat)
col_order <- c("lfcs_DBL21", "lfcs_DBL19", "lfcs_HA11", "lfcs_HA3", "lfcs_ADMD17", "lfcs_WT13")

library(circlize)

# set up colors
colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#2166AC", "white", "#B2182B")
)

colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#3B4CC0", "#F7F7F7", "#B40426")
) # best for white backgrounds

colfun <- colorRamp2(
  c(-4, 0, 4),
  c("#0072B2", "#F0E442", "#D55E00")
)
colfun <- colorRamp2(
  c(-4, -1, 0, 1, 4),
  c("#313695", "#74ADD1", "#F7F7F7", "#FDAE61", "#A50026")
)


colnames(ST_SD_vufun2vsWT_rronly_no_NA_rr_lfc_mat)
col_order <- c("lfcs_finresD1_minusN", "lfcs_finresD3_minusN" ,   "lfcs_finresD1_plusN", "lfcs_finresD3_plusN"  )
column_split_factor <- factor(col_order, levels = col_order)
column_split_factor
# Heatmap

dev.off()
pdf("ST_SD_vufun_RR_locusName.pdf", height = 12)
Heatmap(ST_SD_vufun2vsWT_rronly_no_NA_rr_lfc_mat,
        col = colfun,
        use_raster = TRUE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        row_labels = row_orderRRall$locusName,
        # Show row names (from the matrix) on the right
        show_row_names = TRUE,
        row_names_side = "right",
        # put them on the right
        cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
        # cluster_rows=TRUE,
        # clustering_method_columns="ward.D",
        # clustering_method_rows = "ward.D",
        row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
        row_split = factor(row_orderRRall$Grouping),
        row_gap = unit(2, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          if (vufunvsWT_sig_mat[i, j] != "") {
            grid.text(
              vufunvsWT_sig_mat[i, j],
              x, y,
              gp = gpar(fontsize = 8, col = "red")
            )
          }
        },
        # row_km = 12,
        # column_km=9,
        # column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        row_title_gp = gpar(font = 2, fontsize = 8),
        row_title_rot = 0,
        row_names_gp = gpar(fontsize = 8),
        cluster_column_slices = FALSE,
        column_order = order(col_order),
        column_split = col_order, 
        column_gap = unit(2, "mm"),
        column_title = "vufun TC +N exp: Within genotype differences per timepoint to +N",
        # column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(
          title = "shrunklog2FC",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
          labels = seq(-4, 4, by = 2) # what to print next to these ticke
        )
)
dev.off()

Heatmap(SingleGenoSingleDay_plusN_merged_DEG_mat,
        col = colfun,
        use_raster = TRUE,
        cluster_columns = FALSE,
        cluster_rows = TRUE,
        #row_labels = row_orderRRall$geneSym,
        # Show row names (from the matrix) on the right
        show_row_names = FALSE,
        #row_names_side = "right",
        # put them on the right
        #cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
        # cluster_rows=TRUE,
        # clustering_method_columns="ward.D",
        # clustering_method_rows = "ward.D",
        #row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
        #row_split = factor(row_orderRRall$Grouping),
        #row_gap = unit(2, "mm"),
        #row_km = 22,
        # column_km=9,
        # column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        row_title_gp = gpar(font = 2, fontsize = 6),
        row_title_rot = 0,
        row_names_gp = gpar(fontsize = 6),
        cluster_column_slices = FALSE,
        #column_split = column_split_factor, # <- visually separates minusP vs plusP
        column_gap = unit(2, "mm"),
       # column_title = "vufun TC +N exp: Within genotype differences per timepoint to +N:While Transcriptome",
        # column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(
          title = "shrunklog2FC",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
          labels = seq(-4, 4, by = 2) # what to print next to these ticke
        )
)


dev.off()















# What if we want DAP-seq list?
DAPseqVu <- read.csv("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab/CowpeaDAPseq.csv", header = TRUE)
SingleNSingleDay_vufunvsWT_merged_DEG_annot <- SingleNSingleDay_vufunvsWT_merged_DEG %>%
  left_join(annotations_GL[, c("locusName", "Best.hit.arabi.name", "Best.hit.arabi.defline", "Best.hit.rice.name", "Best.hit.rice.defline", "Transcript_Isoforms")], by = "locusName") %>%
  left_join(RRgenes, by = "locusName") %>%
  left_join(DAPseqVu, by = "locusName") %>%
  left_join(resnormcounts, by = "locusName")

writexl::write_xlsx(SingleNSingleDay_vufunvsWT_merged_DEG_annot, path = "vufunvsWT_ST_SD_fullannot.xlsx")
SingleGenoSingleDay_plusN_merged_DEG


#here's another thing. map back to ZBF mutants!
#Dataset with all ZBF mutants:
getwd()
ZBFbothrronly_DEG <- readxl::read_xlsx("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/ZBF1315_onlyDEG.xlsx")
#Dataset with all vufun TC (within genotype diffs)
SingleGenoSingleDay_plusN_merged <- readxl::read_xlsx("lfc_pval_vufun_singleGeno_singleTP_plusN_combined.xlsx")
#Dataset with  all vufun TC (between genotype diffs)
SingleNSingleDay_vufunvsWT_merged <- readxl::read_xlsx("lfc_pval_vufun_VS_WT_combined.xlsx")

ZBF_vufunSG_vufunBG <- ZBFbothrronly_DEG %>%
  left_join(SingleGenoSingleDay_plusN_merged, by = "locusName") %>%
  left_join(SingleNSingleDay_vufunvsWT_merged, by = "locusName") 

colnames(ZBF_vufunSG_vufunBG)

ZBF_vufunSG_vufunBG_DEG <- ZBF_vufunSG_vufunBG[,-c(1,2)]  %>%
  filter(padj_ADMD17_ref15 < 0.05 | padj_DBL19_ref15 < 0.05 | padj_DBL21_ref15 < 0.05 | padj_HA11_ref15 < 0.05 | padj_HA3_ref15 < 0.05 | padj_WT13_ref15 < 0.05 |
           padj_ADMD17_ref13 < 0.05 | padj_DBL19_ref13 < 0.05 | padj_DBL21_ref13 < 0.05 | padj_HA11_ref13 < 0.05 | padj_HA3_ref13 < 0.05 | padj_WT15_ref13 < 0.05|
           padj_D1_vufun_plusN < 0.05 | padj_D1_WT < 0.05 | padj_D3_vufun_plusN < 0.05 | padj_D3_WT_plusN < 0.05|
           padj_finresD1_minusN < 0.05 | padj_finresD1_plusN < 0.05 | padj_finresD3_minusN < 0.05 | padj_finresD3_plusN < 0.05
           ) 


ZBF_vufunSG_vufunBG_DEG <- ZBF_vufunSG_vufunBG_DEG %>% 
  filter(!if_all(starts_with("lfcs_"), is.na))

ZBF_vufunSG_vufunBG_DEGrronly <- left_join(RRgenes[, c("Gene Symbol", "Function", "locusName")],
                                           ZBF_vufunSG_vufunBG_DEG,
                         by = "locusName"
) # All Rob Roy genes

ZBF_vufunSG_vufunBG_DEGrronly_noNA <- ZBF_vufunSG_vufunBG_DEGrronly %>% 
  filter(!if_all(starts_with("lfcs_"), is.na)) #remove rows if NA in ALL LFCs-containg columns

# Make Groupings
row_orderRRall <- data.frame(
  locusName = ZBF_vufunSG_vufunBG_DEGrronly_noNA$locusName,
  Grouping = ZBF_vufunSG_vufunBG_DEGrronly_noNA$Function,
  geneSym = ZBF_vufunSG_vufunBG_DEGrronly_noNA$`Gene Symbol`
)

row_orderRRall$Grouping <- factor(row_orderRRall$Grouping,
                                    levels = c(
                                      "Early Signalling",
                                      "Host Range Restriction",
                                      "Rhizobial Infection",
                                      "Nodule Organogenesis",
                                      "Autoregulation of Nodule Number",
                                      "Bacterial Maturation",
                                      "Symbiosome Maturation",
                                      "Nodulation Metabolism and Transport",
                                      "Senescence",
                                      "Defense",
                                      "DAP-seq target"
                                    )
)

ALLrr_lfc_mat <- ZBF_vufunSG_vufunBG_DEGrronly_noNA %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

all_padj_mat <- ZBF_vufunSG_vufunBG_DEGrronly_noNA %>%
  select(locusName, starts_with("padj_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()

# colnames(padj_mat) <- sub("padj_", "", colnames(padj_mat))

# Create empty character matrix
all_sig_mat <- matrix("", nrow = nrow(all_padj_mat), ncol = ncol(all_padj_mat))
rownames(all_sig_mat) <- rownames(all_padj_mat)
colnames(all_sig_mat) <- colnames(all_padj_mat)

# Fill with significance levels
all_sig_mat[all_padj_mat < 0.05] <- "*"
all_sig_mat[all_padj_mat < 0.01] <- "**"
all_sig_mat[all_padj_mat < 0.001] <- "***"

colnames(ALLrr_lfc_mat)
col_order <- colnames(ALLrr_lfc_mat)
column_split_factor <- factor(col_order, levels = col_order)
dev.off()
pdf("ZBF_vufunall_genesym.pdf", height = 24,width=12)
Heatmap(ALLrr_lfc_mat,
              col = colfun,
              use_raster = TRUE,
              cluster_columns = FALSE,
              cluster_rows = TRUE,
              row_labels = row_orderRRall$geneSym,
              # Show row names (from the matrix) on the right
              show_row_names = TRUE,
              row_names_side = "right",
              # put them on the right
              cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
              # cluster_rows=TRUE,
              # clustering_method_columns="ward.D",
              # clustering_method_rows = "ward.D",
              row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
              row_split = factor(row_orderRRall$Grouping),
              row_gap = unit(2, "mm"),
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (all_sig_mat[i, j] != "") {
                  grid.text(
                    all_sig_mat[i, j],
                    x, y,
                    gp = gpar(fontsize = 6, col = "red")
                  )
                }
              },
              # row_km = 12,
              # column_km=9,
              # column_km_repeats=100,
              # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
              row_title_gp = gpar(font = 2, fontsize = 6),
              row_title_rot = 0,
              row_names_gp = gpar(fontsize = 6),
              cluster_column_slices = FALSE,
              column_order = column_split_factorall,
              #column_split = column_split_factorall, # <- visually separates minusP vs plusP
              column_gap = unit(2, "mm"),
              column_title = "ZBF +vuFUN TC",
              # column_title_side = "bottom",
              column_title_gp = gpar(fontsize = 15, fontface = "bold"),
              heatmap_legend_param = list(
                title = "shrunklog2FC",
                title_gp = gpar(fontsize = 10, fontface = "bold"),
                labels_gp = gpar(fontsize = 8),
                at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
                labels = seq(-4, 4, by = 2) # what to print next to these ticke
              )
)
dev.off()


#More comparisons

SingleNSingleDay_vufunvsWT_merged_DEG <- readxl::read_xlsx(path = "lfc_pval_vufun_VS_WT_combined_DEG.xlsx")
colnames(SingleNSingleDay_vufunvsWT_merged_DEG)

#Leave out if not useful
SingleNSingleDay_vufunvsWT_merged_DEG <- SingleNSingleDay_vufunvsWT_merged_DEG %>%
  left_join(SingleGenoSingleDay_plusN_merged,by="locusName")
vufun2_vs_WT_plusN_D3_DEG <- SingleNSingleDay_vufunvsWT_merged_DEG %>%
  filter(padj_finresD3_plusN < 0.05)

super_annot <- left_join(annotations_GL,RRgenes,by="locusName")
vufun2_vs_WT_plusN_D3_DEG_annot <- left_join(vufun2_vs_WT_plusN_D3_DEG,super_annot[,c("locusName","Function","Gene Name")],by="locusName")
colnames(vufun2_vs_WT_plusN_D3_DEG)
ST_SD_vufun2vsWT_merged_DEG_mat <- vufun2_vs_WT_plusN_D3_DEG[,-c(14:19)] %>%
  select(locusName, starts_with("lfcs_")) %>%
  column_to_rownames("locusName") %>%
  as.matrix()
colnames(ST_SD_vufun2vsWT_merged_DEG_mat)
col_order <- colnames(ST_SD_vufun2vsWT_merged_DEG_mat)
column_split_factor <- factor(col_order, levels = col_order)

# Fantastic, now we do some heatmaps
dev.off()
pdf("vufun2vsWTD3_plusN_locusName.pdf", height = 12)
set.seed(123)
ht=Heatmap(ST_SD_vufun2vsWT_merged_DEG_mat,
        col = colfun,
        use_raster = TRUE,
        cluster_columns = FALSE,
        #row_labels = row_orderRRall$geneSym,
        # Show row names (from the matrix) on the right
        show_row_names = FALSE,
        #row_names_side = "right",
        # put them on the right
        #cluster_row_slices = FALSE, # IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
        cluster_rows=TRUE,
        cluster_row_slices = TRUE,
        clustering_method_columns="ward.D",
        clustering_method_rows = "ward.D",
        #row_dend_reorder = FALSE, # turned off because cant see difference due to magnitude of genes (5000+)
        #row_split = factor(row_orderRRall$Grouping),
        row_gap = unit(2, "mm"),
        row_km = 26,
        # column_km=9,
        # column_km_repeats=100,
         row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        row_title_gp = gpar(font = 2, fontsize = 6),
        row_title_rot = 0,
        row_names_gp = gpar(fontsize = 6),
        cluster_column_slices = FALSE,
        column_split = column_split_factor, # <- visually separates minusP vs plusP
        column_gap = unit(2, "mm"),
        #column_names_rot=90,
        # column_title = "vufun TC +N exp: Within genotype differences per timepoint to +N:While Transcriptome",
        # column_title_side = "bottom",
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        heatmap_legend_param = list(
          title = "shrunklog2FC",
          title_gp = gpar(fontsize = 10, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          at = seq(-4, 4, by = 2), # where to pick the ticks at, by denotes max
          labels = seq(-4, 4, by = 2) # what to print next to these ticke
        )
)
HM=draw(ht)
dev.off()
#extract clusters
library(magrittr)
r.dend <- row_dend(HM)  #If needed, extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x)) #check/confirm size gene clusters
# Extract clusters as a data.frame
clusterlist = row_order(HM)

clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(locusName = rownames(ST_SD_vufun2vsWT_merged_DEG_mat[rcl.list[[i]],]), #rownames(tf.log)[clusterlist[[i]]], #if cluster somehow only has one gene. https://www.biostars.org/p/465304/
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)


clu_dfannot <- clu_df %>%
  left_join(super_annot,by= "locusName") %>%
  left_join(DAPseqVu,by="locusName")

View(clu_df)
writexl::write_xlsx(
  clu_dfannot,
  "extractedclusters_DEG_vufunvsWT_plusN_D3.xlsx"
)

#
# Great now, lets do some general stats
str(SingleGenoSingleDay_plusN_merged)
str(SingleNSingleDay_vufunvsWT_merged)

# Define your lines
lines <- c("D1_vufun_plusN", "D1_WT_plusN", "D3_vufun_plusN","D3_WT_plusN")
lines13 <- c("DBL21", "DBL19", "HA11", "HA3", "ADMD17", "WT15")
# Loop through each line and calculate up/down DEGs
lines <- c("D1_vufun_plusN", "D1_WT_plusN", "D3_vufun_plusN","D3_WT_plusN")

deg_summary <- lapply(lines, function(line) {
  
  lfc_col <- paste0("lfc_", line)
  padj_col <- paste0("padj_", line)
  
  df <- data.frame(
    lfc = SingleGenoSingleDay_plusN_merged[[lfc_col]],
    padj = SingleGenoSingleDay_plusN_merged[[padj_col]]
  ) %>%
    filter(!is.na(padj), padj < 0.05)
  
  data.frame(
    line = line,
    up = sum(df$lfc > 0),
    down = sum(df$lfc < 0),
    total_DEG = nrow(df)
  )
}) %>% bind_rows()

deg_summary

df_long <- deg_summary %>%
  pivot_longer(cols = c(up, down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Direction = ifelse(Direction == "up", "Up", "Down"),
    Treatment = line
  )

a=ggplot(df_long, aes(x = Treatment, y = ifelse(Direction == "Up", Count, -Count), fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4.5
  ) +
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) +
  #scale_x_discrete(limits = c("DBL21", "DBL19", "HA11","HA3", "ADMD17" , "WT13")) + # custom order
  labs(
    title = "Upregulated and Downregulated DEGs for within genotype differences to +N",
    x = "ZBF Line",
    y = "Gene Count",
    fill = "Regulation"
  ) +
  ylim(-7000, 7000) +
  theme_bw() +
  theme(
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 15, vjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )

a
+#df_longbetween

linesbetween <- c(
  "finresD1_minusN",
  "finresD1_plusN",
  "finresD3_minusN",
  "finresD3_plusN"
)

deg_summarybetween <- lapply(linesbetween, function(line) {
  
  lfc_col  <- paste0("lfc_", line)
  padj_col <- paste0("padj_", line)
  
  df <- data.frame(
    lfc  = SingleNSingleDay_vufunvsWT_merged[[lfc_col]],
    padj = SingleNSingleDay_vufunvsWT_merged[[padj_col]]
  ) %>%
    filter(!is.na(padj), padj < 0.05)
  
  data.frame(
    contrast = line,
    up = sum(df$lfc > 0),
    down = sum(df$lfc < 0),
    total_DEG = nrow(df)
  )
  
}) %>% bind_rows()

deg_summarybetween



df_longbetween <- deg_summarybetween %>%
  pivot_longer(cols = c(up, down), names_to = "Direction", values_to = "Count") %>%
  mutate(
    Direction = ifelse(Direction == "up", "Up", "Down"),
    Treatment = contrast
  )

df_longbetween

b=ggplot(df_longbetween, aes(x = Treatment, y = ifelse(Direction == "Up", Count, -Count), fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4.5
  ) +
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) +
  #scale_x_discrete(limits = c("DBL21", "DBL19", "HA11","HA3", "ADMD17" , "WT13")) + # custom order
  labs(
    title = "Upregulated and Downregulated DEGs for genotypic differences per Treatment and Day",
    x = "ZBF Line",
    y = "Gene Count",
    fill = "Regulation"
  ) +
  ylim(-1000, 1000) +
  theme_bw() +
  theme(
    axis.title.x = element_text(color = "black", size = 15, face = "bold"),
    axis.title.y = element_text(color = "black", size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 15, vjust = 0.5, face = "bold"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.key.size = unit(1, "cm"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 10)
  )
b

library(gridExtra)
combined <- grid.arrange(a,b,ncol=1)
ggsave(combined, file = "vufun_DEGbarplotsboth.pdf", height = 25, width = 25, units = "cm")
