library


BiocManager::install("apeglm")
BiocManager::install("ashr")
install.packages("styler")
# Author: Benjamin Wee
# Purpose: tximport and DeSeq2 on ljfun mutants, along with heatmap and venns, from Salmon workflow.
# Prelim data comparisons: For each timepoint, per genotype, what is the difference with +N
# Github: https://github.com/BWY9394/Reid-Smith-Lab/tree/main/Bulk%20RNA-seq/Salmon

##### Step 1: Load packages and set working directory #####
library(tidyverse)
library(tximport)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)
library(ggrepel)
library(seqinr)
library(ashr)
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/L1_L24_LjFUN_timecourse2025/quantsL1_L24/")
getwd()
save.image(file=".RData")
load(file="P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/L1_L24_LjFUN_timecourse2025/.RData")
setwd("Z:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/L1_L24_LjFUN_timecourse2025/quantsL1_L24")


#Step 2 generate filenames/folderpaths ####
folders <-  data.frame(folder_names=list.files("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/L1_L24_LjFUN_timecourse2025/quantsL1_L24"))
folders
getwd()
write.csv(folders,file="filenames.csv",row.names=FALSE)
#Step 2: Start here if rerunnign analayis
folders = read.csv("filenames.csv", header = TRUE)
folders = as.character(folders$folder_names)
folders <- gtools::mixedsort(folders)
folders
#generate a data frame with all filesnames and paths
files = file.path(head(folders,-1), "quant.sf")
files
names(files) <- paste0("L", 1:24) # Assign sample names to each file. These names MUST match rownames in the sample metadata (sampleTable)
files
getwd()


# Step 3: Sanity check to confirm that all quant.sf files exist  #####
# Returns TRUE if every expected file is present. So each rowname will be what is named to your salmon output abundnace/count/readlength per sample, as well as where that is stored, so should not have ANY  mismatches.
# If FALSE, do not proceed until the file -> filename link is fixed and returns TRUE
all(file.exists(files))
head(files, n = 24)

# Step 4: Load and clean annotations for tx2gene dataframe that tximport requires  #####
#annotations <- read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) # 54484 annotated transcripts

annotations <- read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Lotus/Phytozome/Gifu v1.3/download.20260424.121032/Phytozome/PhytozomeV14/LjaponicusGifu/v1.3/annotation/LjaponicusGifu_1046_v1.3.P14.annotation_info.txt.gz", header = TRUE, sep = "\t", stringsAsFactors = FALSE) # 54484 annotated transcripts

View(annotations)

transcript_fasta_data <-  read.fasta(file = "P:/LAB - OPV - Dugald_Reid/Genomics/Lotus/transcriptsgifu.fa", as.string = TRUE) # 54484 transcripts that were used for mapping

View(transcript_fasta_data_df)
#
transcript_fasta_data_df <- do.call(rbind, lapply(transcript_fasta_data, function(seq_obj) {
  data.frame(
    transcriptName = attr(seq_obj, "name"),
    Vuannotation = attr(seq_obj, "Annot"),
    ft_ni_sequence = as.character(seq_obj),
    stringsAsFactors = FALSE
  )
}))
transcript_fasta_data_df$locusName <- transcript_fasta_data_df$transcriptName
test_transcript_fasta_data_df_removed <-   transcript_fasta_data_df %>%
  mutate(
    locusName = str_remove(locusName, "_.*$")
  ) %>%
  mutate(locusName = str_remove(locusName, "\\.\\d+$")
  ) 
View(annotations2)
colnames(annotations)
annotations <- rename(annotations, "Protein-Accession" = "transcriptName")
colnames(annotations)
annotations2 <- left_join(annotations,test_transcript_fasta_data_df_removed[,c(1,4)],by="transcriptName") #but then lose a lot, better if we do from the official annotations file.




colnames(test_transcript_fasta_data_df_removed)
data.table::uniqueN(test_transcript_fasta_data_df_removed$locusName) # Number of unique transcripts, discards repeated rows per transcript #30586 #extracted from transcriptome used for mapping
data.table::uniqueN(annotations2$locusName) # Number of unique transcripts, discards repeated rows per transcript #29419 #extracted from Dugald's list of annotations. Need to check discrepency.
tx2geneLj <- tibble::as_tibble(test_transcript_fasta_data_df_removed[,c(1,4),])
head(tx2geneLj)# save annotations into df for tximport to summarize gene-level from transcript isoforms


any(duplicated(tx2geneLj[, 1])) # if not false, means you have MULTIPLE column 1 rownames, which means you have multiple rows with the same transcript isoform! You can have multipled column 2 matches, (locus of the transcript, but not thetranscripts themselves!)

## Step 6: Tximport to summarize counts,TPM, avgereadlngth on a gene-level basis ######
# navigate to folders path
txi.salmon.tsv <- tximport(files,
                           type = "salmon",
                           tx2gene = tx2geneLj,
                           dropInfReps = TRUE,
                           countsFromAbundance = "no"
) # don't worry about this too much, if you are using DESeq2, it will do its own thing with count normalization. And you will use TPMs for gene abundance anyway.

final <- data.frame(txi.salmon.tsv) # save into dataframe. Has TPM, counts, avge read length all in one df. Will extract into single dataframes
str(final)
colnames(final)
str(files)


# Step 7: Load in relevant metadata for this experiment per sample, so as to match with colnames of saved tximport object #####
metadata <- readxl::read_xlsx("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/01.RawData/Copy of Sample Register.xlsx") # load meta data in
str(metadata)
metadata$Project <- factor(metadata$Project)
levels(metadata$Project)
ljfunmetadata <- metadata %>%
  filter(Project == "ljfun vs WT +/- N timecourse") # filter out for just ljfun N tc exp experiment
ljfunmetadata$CombID <- factor(ljfunmetadata$CombID)
#levels(ljfunmetadata$CombID) if you make th
ljfunmetadata$CountsID <- paste("counts",ljfunmetadata$`RNA Sample Code`, sep="." )
ljfunmetadata$TPMID <- paste("abundance",ljfunmetadata$`RNA Sample Code`, sep="." )

# Step 8: Save Counts #####
#header <- c("locusName", c(paste("counts", ljfunmetadata$CombID, ljfunmetadata$Rep, sep = "."))) # Create column headers by combining "counts" prefix with sample descriptions from metadata. Output will be "counts.Description1", "counts.Description2", etc.
header <- setNames(ljfunmetadata$CombID, ljfunmetadata$CountsID)
header
counts <- data.frame(final[, grep("counts.", colnames(final))]) # Extract only the count columns from the final dataframe
colnames(counts) # Check current column names (may include an unwanted "countsFromAbundance" column)
counts <- cbind(locusName = rownames(counts), counts[, -25]) # counts[,-25] removes the 25th column which is likely the unwanted countsFromAbundance column. cbind() combines the new locusName column with the existing data
colnames(counts) # Verify the column names after adding locusName and removing the extra column
#colnames(counts) <- header # Rename all columns using the cleaned header vector created earlier. This replaces the default/messy column names with the standardized format. And we know
#
#str(counts)
#Apply the transformation
colnames(counts)[-1] <- header[colnames(counts)[-1]] #won't work if the cols in header are factors.

head(counts)
View(counts)
write.csv(counts, file = "ljfun_TC_Nitrate_Treatments_rawcounts.csv", row.names = FALSE) # Save output


# Step 9: Save TPM #####
#header <- c("locusName", c(paste("TPM", ljfunmetadata$CombID, ljfunmetadata$Rep, sep = "."))) # Create column headers by combining "counts" prefix with sample descriptions from metadata. Output will be "counts.Description1", "counts.Description2", etc.
header <- setNames(ljfunmetadata$CombID, ljfunmetadata$TPMID)
header
TPM <- data.frame(final[, grep("abundance", colnames(final))]) # Extract only the count columns from the final dataframe
TPM <- cbind(locusName = rownames(TPM),TPM)
head(TPM) # Check current column names (may include an unwanted "countsFromAbundance" column)
head(TPM)
str(TPM)
#Apply the transformation
colnames(TPM)[-1] <- header[colnames(TPM)[-1]]
head(TPM)
write.csv(TPM, file = "ljfun_TC_Nitrate_Treatments_TPM.csv", row.names = FALSE) # Save output
View(TPM)


# Step 10: Generate some stats for TPM #####
TPM_long <- pivot_longer(TPM,
                         cols = -locusName,
                         names_to = "Sample",
                         values_to = "TPM"
)
head(TPM_long$Sample)
colnames(TPM)
TPMheaders <- colnames(TPM[, -1])

ljfunmetadata_tpm <- data.frame(Sample = TPMheaders, ljfunmetadata[, c("Treatment", "Day", "Rep", "Identifying Allele", "Description", "CombID")])
ljfunmetadata_tpm <- left_join(TPM_long, ljfunmetadata_tpm, by = "Sample")
View(ljfunmetadata_tpm)
colnames(ljfunmetadata_tpm)
write.csv(ljfunmetadata_tpm, file = "ljfunTPM_longformat_metadata.csv", row.names = F, quote = F)

TPM_mean_comb <- plyr::ddply(
  ljfunmetadata_tpm,
  c("locusName", "Identifying.Allele", "Day", "Treatment", "CombID"),
  summarise,
  mean = mean(TPM, na.rm = TRUE),
  sd   = sd(TPM, na.rm = TRUE),
  n    = sum(!is.na(TPM)),
  sem  = sd(TPM, na.rm = TRUE) / sqrt(sum(!is.na(TPM)))
) # Generates mean, sd, sem, n for each gene per sample.


writexl::write_xlsx(
  TPM_mean_comb,
  "ljfun_TPMstats.xlsx"
)

TPM_mean_comb <- readxl::read_xlsx("ljfun_TPMstats.xlsx")

##### Step 11: Last component needed for DESeq analysis, which is to provide experimental setup details in relation to experiment treatment factors #####
filenames <- file.path(folders, "quant.sf") # we created folders earlier, object containing folder names which hold .sf files
filenames # so now we have filepaths that are correct for each sample
sampleNames <- sub("_quant/quant.sf", "", filenames)
sampleNames # just one of many ways to get the sampleNames(literal) sample names)
sampleNames

colnames(ljfunmetadata_tpm)
sampleTable <- ljfunmetadata[, c("CombID")]
sampleTable
sampleTable <- data.frame(sampleTable)
sampleTable$`RNA Sample Code` <- colnames(txi.salmon.tsv$counts)
sampleTable_full <- sampleTable %>%
  left_join(ljfunmetadata[, c("Treatment", "Day", "Identifying Allele", "RNA Sample Code")], by = "RNA Sample Code") %>%
  column_to_rownames("RNA Sample Code")
all(rownames(sampleTable_full) == colnames(txi.salmon.tsv$counts))
sampleTable_full # #  sampleTable, which contains a data frame of sample ID assigned per sample when building txi.salmon.tsv object in the rownames, while the first (and only) columns consist of the experiment treatment factor
sampleTable_full[] <- lapply(sampleTable_full, as.factor) # make everything a factor
str(sampleTable_full)
sampleTable_full <- dplyr::rename(sampleTable_full, `Genotype` = "Identifying Allele")


# Step 12 workaround: Carry out DEG analysis or pre-grouped samples.

dds <- DESeqDataSetFromTximport(
  txi.salmon.tsv,
  sampleTable_full,
  ~CombID
)

dds <- DESeq(dds)
resultsNames(dds)
levels(dds$CombID)

 
# Step 13: Global Transcriptomic response PCAs from normalized counts ####
vsd <- vst(dds, blind = TRUE) # apply variance stabilizing transformation (VSD) to reduce dependence of variance on mean experession (more homoscedastic). Corrects for natural high variance at high experssion levels, low variance at low experience levels, making data visualization (PCA/Heatmaps) clearer.
plotPCA(vsd, intgroup = c("CombID"), returnData = FALSE, pcsToUse = c(1, 2)) # raw PCA. Good, but want more control over certain variables
pca_data <- plotPCA(
  vsd,
  intgroup = "CombID",
  pcsToUse = c(1, 2),
  returnData = TRUE
) # build a better PCA, so first save pca data into dataframe for PC1 nad PC2

pca_data


pca_data <- dplyr::rename(pca_data, `RNA Sample Code` = name) # need to to do this to merge some metadta backinto this PCA summary of PC1 and PC2
pca_data
colnames(pca_data)
percent_var <- round(100 * attr(pca_data, "percentVar"), 1)
levels(pca_data$Genotype)
# subset(pca_data, !Identifying.Allele %in% c("WT13", "WT15"))

PCA1 <- ggplot(pca_data, aes(PC1, PC2, color = Day, shape = Treatment)) + # now create custom PCA
  geom_point(
    size = 2.8,
    alpha = 0.85
  ) +
  geom_text_repel(
    aes(label = `Genotype`), # assign labels for each dotpoint
    size = 4, # text size
    max.overlaps = 20,
    box.padding = 0.3, # toggle higher to make line visible
    point.padding = 0.3, # line to dot space ratio. Higher = more distance from line to dot
    segment.color = "grey70" # line to text color
  ) +
 # xlim(-75,75)+
  theme_bw() +
  labs(
    title = "PCA of ljfun timecourse experiment: PC1 vs PC2",
    x = paste0("PC1 (", percent_var[1], "% variance)"),
    y = paste0("PC2 (", percent_var[2], "% variance)"),
    color = "Timepoint"
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
PCA1
ggsave(PCA1, file = "Global_PCA_lj.pdf", height = 25, width = 22.25, units = "cm")
ggsave(PCA1, file = "Global_PCA_lj.jpg", height = 25, width = 22.25, units = "cm")

##### Step 14: Get normalized counts for each gene per sample. Can use for gene abundance comparisons ####
resnormcounts <- counts(dds, normalized = T)
resnormcounts <- as.data.frame(resnormcounts)
resnormcounts
head(ljfunmetadata$CombID, n = 24)
colnames(resnormcounts) <- ljfunmetadata$CombID
colnames(resnormcounts)
header <- c(c(paste(ljfunmetadata$CombID, ljfunmetadata$Rep, sep = ".R"))) # Create column headers by combining "counts" prefix with sample descriptions from metadata. Output will be "counts.Description1", "counts.Description2", etc.
colnames(resnormcounts) <- header
resnormcounts <- resnormcounts %>%
  rownames_to_column("locusName")
resnormcounts
write.csv(resnormcounts, file = "ljfun_norm_counts_results.csv", row.names = FALSE)

ljresnormcountslong <- pivot_longer(resnormcounts,
                                    cols = -locusName,
                                    names_to = "Sample",
                                    values_to = "normcounts"
)
ljnormcountheaders <- colnames(resnormcounts[, -1])
ljfunmetadata
ljmetadata_counts <- data.frame(Sample = ljnormcountheaders, ljfunmetadata[, c("RNA Sample Code", "Identifying Allele", "Treatment", "Day", "CombID", "Rep")])
ljmetadata_counts
ljmetadata_counts <- left_join(ljresnormcountslong, ljmetadata_counts, by = "Sample")
ljmetadata_counts$CombID <- factor(ljmetadata_counts$CombID)
levels(ljmetadata_counts$CombID)
ljmetadata_counts$Rep <- factor(ljmetadata_counts$Rep)
ljmetadata_counts$CombID <- factor(ljmetadata_counts$CombID, levels = c("ljfun3_D1_Nminus", "ljfun3_D1_Nplus" , "ljfun3_D3_Nminus", "ljfun3_D3_Nplus" , "WT_D1_Nminus"  ,   "WT_D1_Nplus"  ,    "WT_D3_Nminus"  ,   "WT_D3_Nplus"))
ljtoplot <- ljmetadata_counts %>%
  dplyr::filter(locusName %in% c("LotjaGi2g1v0279100", "LotjaGi2g1v0259200", "LotjaGi1g1v0199700","LotjaGi3g1v0487600"))
View(ljmetadata_counts)
head(ljmetadata_counts,n=20)

#Plots
ggplot(ljtoplot, aes(x = CombID, y = normcounts)) +
  geom_boxplot(fill = NA, color = "gray90", alpha = 0.05, outlier.shape = NA) +
  geom_point(aes(color = Rep)) +
  # geom_bar(stat="")+
  #  geom_text(aes(label=Rep))+
  facet_wrap(~`locusName`, nrow = 1, scales = "free") +
  labs(
    title = "Normalized counts of select genes\nfor ljFun lines", y = "normalized counts", x = "ZBF line"
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

# Step 15: Get DEG results for all comparisons (between genotypes and within genotypes). log-fold changes and adj p-vals ####



# Step 15.1 Within Genotype Response to +N ####
resultsNames(dds)
print(ljfunmetadata$CombID)


##### Step 15.1, within genotype response to +N #####
# D1, ljfun +N
resultsNames(dds) # these are all the generated comparisons.

resD1_ljfun <- results(dds, contrast = c("CombID", "ljfun3_D1_Nplus", "ljfun3_D1_Nminus"), alpha = 0.05) # the last field is what you are contrasting to i.e. fold change compared to N-minus conditions at D1 for ljfun3.
summary(resD1_ljfun)
OEresD1_ljfun <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D1_Nplus", "ljfun3_D1_Nminus"), type = "ashr", res = resD1_ljfun)
resD1_ljfundf <- as.data.frame(resD1_ljfun)
resD1_ljfundf <- resD1_ljfundf %>%
  rownames_to_column("locusName")
OEresD1_ljfundf <- as.data.frame(OEresD1_ljfun)
OEresD1_ljfundf <- OEresD1_ljfundf %>%
  rownames_to_column("locusName")
OEresD1_ljfundf <- OEresD1_ljfundf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD1_ljfun <- left_join(resD1_ljfundf, OEresD1_ljfundf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
#finresD1_ljfun <- left_join(finresD1_ljfun, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
#finresD1_ljfun <- left_join(finresD1_ljfun, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD1_ljfun <- finresD1_ljfun[order(finresD1_ljfun$padj), ]
finresD1_ljfun$DEG <- ifelse(finresD1_ljfun$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD1_ljfun, path = "results_pvals_D1_ljfun_plusN_vs_minusN.xlsx")
view(finresD1_ljfun)
getwd()

# D1, WT +N
resultsNames(dds) # these are all the generated comparisons.

resD1_WT <- results(dds, contrast = c("CombID", "WT_D1_Nplus", "WT_D1_Nminus"), alpha = 0.05)
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


# D3, ljfun +N
resultsNames(dds) # these are all the generated comparisons.

resD3_ljfun <- results(dds, contrast = c("CombID", "ljfun3_D3_Nplus", "ljfun3_D3_Nminus"), alpha = 0.05)
summary(resD3_ljfun)
OEresD3_ljfun <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D3_Nplus", "ljfun3_D3_Nminus"), type = "ashr", res = resD3_ljfun)
summary(OEresD3_ljfun)
resD3_ljfundf <- as.data.frame(resD3_ljfun)
resD3_ljfundf <- resD3_ljfundf %>%
  rownames_to_column("locusName")
OEresD3_ljfundf <- as.data.frame(OEresD3_ljfun)
OEresD3_ljfundf <- OEresD3_ljfundf %>%
  rownames_to_column("locusName")
OEresD3_ljfundf <- OEresD3_ljfundf %>%
  dplyr::rename(l2FCshrunk = log2FoldChange, lfcshrunkSE = lfcSE)
finresD3_ljfun <- left_join(resD3_ljfundf, OEresD3_ljfundf[, c("locusName", "l2FCshrunk", "lfcshrunkSE")], by = "locusName")
finresD3_ljfun <- left_join(finresD3_ljfun, annotations_GL[, c("locusName", "Best.hit.arabi.defline")], by = "locusName")
finresD3_ljfun <- left_join(finresD3_ljfun, RRgenes[, c("locusName", "Gene Symbol", "Orthogroup", "Function", "Protein class/Molecular function")], by = "locusName")
finresD3_ljfun <- finresD3_ljfun[order(finresD3_ljfun$padj), ]
finresD3_ljfun$DEG <- ifelse(finresD3_ljfun$padj < 0.05, "Yes", "No")
writexl::write_xlsx(finresD3_ljfun, path = "results_pvals_D3_ljfun_plusN_vs_minusN.xlsx")
getwd()


# D3, WT +N
resultsNames(dds) # these are all the generated comparisons.

resD3_WT <- results(dds, contrast = c("CombID", "WT_D3_Nplus", "WT_D3_Nminus"), alpha = 0.05)
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
pdf(file = "ljfun_singleGeno_singleTP_plusN.pdf", height = 8, width = 12)
par(mfrow = c(3, 4))
plotMA(resD1_ljfun, main = "D1,ljfun +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_ljfun, main = "D1,ljfun +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD1_WT, main = "D1,WT +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_WT, main = "D1,WT +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_ljfun, main = "D3,ljfun +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_ljfun, main = "D3,ljfun +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_WT, main = "D3,WT +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_WT, main = "D3,WT +N:\nWith Shrinkage", ylim = c(-4, 4))

dev.off()


##### Step 15.2, between genotype responses, ljfun vs. WT ####
# D1, ljfun vs WT, -N
resultsNames(dds) # these are all the generated comparisons.

resD1_minusN <- results(dds, contrast = c("CombID", "ljfun3_D1_Nminus", "WT_D1_Nminus"), alpha = 0.05)
summary(resD1_minusN)
OEresD1_minusN <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D1_Nminus", "WT_D1_Nminus"), type = "ashr", res = resD1_minusN)
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
writexl::write_xlsx(finresD1_minusN, path = "results_pvals_D1_minusN_ljfun_vs_WT.xlsx")
getwd()

# D1, ljfun vs WT, +N
resultsNames(dds) # these are all the generated comparisons.

resD1_plusN <- results(dds, contrast = c("CombID", "ljfun3_D1_Nplus", "WT_D1_Nplus"), alpha = 0.05)
summary(resD1_plusN)
OEresD1_plusN <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D1_Nplus", "WT_D1_Nplus"), type = "ashr", res = resD1_plusN)
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
writexl::write_xlsx(finresD1_plusN, path = "results_pvals_D1_plussN_ljfun_vs_WT.xlsx")
getwd()


# D3, ljfun vs WT, -N
resultsNames(dds) # these are all the generated comparisons.

resD3_minusN <- results(dds, contrast = c("CombID", "ljfun3_D3_Nminus", "WT_D3_Nminus"), alpha = 0.05)
summary(resD3_minusN)
OEresD3_minusN <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D3_Nminus", "WT_D3_Nminus"), type = "ashr", res = resD3_minusN)
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
writexl::write_xlsx(finresD3_minusN, path = "results_pvals_D3_minussN_ljfun_vs_WT.xlsx")
getwd()


# D3, ljfun vs WT, +N
resultsNames(dds) # these are all the generated comparisons.

resD3_plusN <- results(dds, contrast = c("CombID", "ljfun3_D3_Nplus", "WT_D3_Nplus"), alpha = 0.05)
summary(resD3_plusN)
OEresD3_plusN <- lfcShrink(dds, contrast = c("CombID", "ljfun3_D3_Nplus", "WT_D3_Nplus"), type = "ashr", res = resD3_plusN)
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
writexl::write_xlsx(finresD3_plusN, path = "results_pvals_D3_plussN_ljfun_vs_WT.xlsx")
getwd()


# MA plots
dev.off()
pdf(file = "ljfun_betweengenotypes_singleTP_singleTrt.pdf", height = 8, width = 12)
par(mfrow = c(3, 4))
plotMA(resD1_minusN, main = "D1,WT vs ljfun3 -N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_minusN, main = "D1,WT vs ljfun3 -N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD1_plusN, main = "D1,WT vs ljfun3 +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD1_plusN, main = "D1,WT vs ljfun3 +N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_minusN, main = "D3,WT vs ljfun3 -N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_minusN, main = "D3,WT vs ljfun3 -N:\nWith Shrinkage", ylim = c(-4, 4))
plotMA(resD3_plusN, main = "D3,WT vs ljfun3 +N:\nNo Shrinkage", ylim = c(-4, 4))
plotMA(OEresD3_plusN, main = "D3,WT vs ljfun3 +N:\nWith Shrinkage", ylim = c(-4, 4))

dev.off()


##### Step 16 Collect all dataframes  #####
##### Step 16.1 Collect for within Genotype differences (for single genotype, per day, Treatment effects) #####
# Ok let's start pulling shit in
SingleGenoSingleDay_plusN <- list(
  finresD1_ljfun %>%
    select(locusName,
           lfc_D1_ljfun_plusN = log2FoldChange,
           lfcs_D1_ljfun_plusN = l2FCshrunk,
           padj_D1_ljfun_plusN = padj
    ),
  finresD1_WT %>%
    select(locusName,
           lfc_D1_WT_plusN = log2FoldChange,
           lfcs_D1_WT_plusN = l2FCshrunk,
           padj_D1_WT_plusN = padj
    ),
  finresD3_ljfun %>%
    select(locusName,
           lfc_D3_ljfun_plusN = log2FoldChange,
           lfcs_D3_ljfun_plusN = l2FCshrunk,
           padj_D3_ljfun_plusN = padj
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

writexl::write_xlsx(SingleGenoSingleDay_plusN_merged, path = "lfc_pval_ljfun_singleGeno_singleTP_plusN_combined.xlsx")
SingleGenoSingleDay_plusN_merged <- readxl::read_xlsx(path = "lfc_pval_ljfun_singleGeno_singleTP_plusN_combined.xlsx")

# for padj<0.05
colnames(SingleGenoSingleDay_plusN_merged)
SingleGenoSingleDay_plusN_merged_DEG <- SingleGenoSingleDay_plusN_merged %>%
  filter(padj_D1_ljfun_plusN < 0.05 | padj_D1_WT_plusN < 0.05 | padj_D3_ljfun_plusN < 0.05 | padj_D3_WT_plusN < 0.05)

writexl::write_xlsx(SingleGenoSingleDay_plusN_merged_DEG, path = "lfc_pval_vufun_singleGeno_singleTP_plusN_combined_DEG.xlsx")
SingleGenoSingleDay_plusN_merged_DEG <- readxl::read_xlsx(path = "lfc_pval_vufun_singleGeno_singleTP_plusN_combined_DEG.xlsx")


##### Step 16.2 Collect for between Genotype differences (between genotypes, what is the difference per Day and Treatment condition) ####
SingleNSingleDay_ljfunvsWT <- list(
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

SingleNSingleDay_ljfunvsWT_merged <- purrr::reduce(SingleNSingleDay_ljfunvsWT, inner_join, by = "locusName")
SingleNSingleDay_ljfunvsWT_merged <- SingleNSingleDay_ljfunvsWT_merged %>%
  dplyr::distinct(locusName, .keep_all = TRUE) # 31337

writexl::write_xlsx(SingleNSingleDay_ljfunvsWT_merged, path = "lfc_pval_ljfun_VS_WT_combined.xlsx")
SingleNSingleDay_ljfunvsWT_merged <- readxl::read_xlsx(path = "lfc_pval_ljfun_VS_WT_combined.xlsx")
# for padj<0.05
colnames(SingleNSingleDay_ljfunvsWT_merged)
SingleNSingleDay_ljfunvsWT_merged_DEG <- SingleNSingleDay_ljfunvsWT_merged %>%
  filter(padj_finresD1_minusN < 0.05 | padj_finresD1_plusN < 0.05 | padj_finresD3_minusN < 0.05 | padj_finresD3_plusN < 0.05)

writexl::write_xlsx(SingleNSingleDay_ljfunvsWT_merged_DEG, path = "lfc_pval_ljfun_VS_WT_combined_DEG.xlsx")
SingleNSingleDay_ljfunvsWT_merged_DEG <- readxl::read_xlsx(path = "lfc_pval_ljfun_VS_WT_combined_DEG.xlsx")

# What if we want DAP-seq list?
DAPseqVu <- read.csv("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab/CowpeaDAPseq.csv", header = TRUE)
SingleNSingleDay_ljfunvsWT_merged_DEG_annot <- SingleNSingleDay_ljfunvsWT_merged_DEG %>%
  left_join(annotations_GL[, c("locusName", "Best.hit.arabi.name", "Best.hit.arabi.defline", "Best.hit.rice.name", "Best.hit.rice.defline", "Transcript_Isoforms")], by = "locusName") %>%
  left_join(RRgenes, by = "locusName") %>%
  left_join(DAPseqVu, by = "locusName") %>%
  left_join(resnormcounts, by = "locusName")

writexl::write_xlsx(SingleNSingleDay_ljfunvsWT_merged_DEG_annot, path = "ljfunvsWT_ST_SD_fullannot.xlsx")
SingleGenoSingleDay_plusN_merged_DEG


##### This bit is for DEG barplots #####
# Great now, lets do some general stats
View(SingleGenoSingleDay_plusN_merged)
str(SingleGenoSingleDay_plusN_merged)
str(SingleNSingleDay_ljfunvsWT_merged)
SingleGenoSingleDay_plusN_merged <- readxl::read_xlsx(path = "lfc_pval_ljfun_singleGeno_singleTP_plusN_combined.xlsx") # all genes
# Define your lines
lines <- c("D1_ljfun_plusN", "D1_WT_plusN", "D3_ljfun_plusN", "D3_WT_plusN")
# Loop through each line and calculate up/down DEGs
lines <- c("D1_ljfun_plusN", "D1_WT_plusN", "D3_ljfun_plusN", "D3_WT_plusN")
colnames(SingleGenoSingleDay_plusN_merged)
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

a <- ggplot(df_long, aes(x = Treatment, y = ifelse(Direction == "Up", Count, -Count), fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4.5
  ) +
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) +
  # scale_x_discrete(limits = c("DBL21", "DBL19", "HA11","HA3", "ADMD17" , "WT13")) + # custom order
  labs(
    title = "Upregulated and Downregulated DEGs for within genotype differences to +N",
    x = "Genotypes (ljfun, WT)",
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
+# df_longbetween
  
  SingleNSingleDay_ljfunvsWT_merged <- readxl::read_xlsx(path = "lfc_pval_ljfun_VS_WT_combined.xlsx") # all genes
colnames(SingleNSingleDay_ljfunvsWT_merged)
linesbetween <- c(
  "finresD1_minusN",
  "finresD1_plusN",
  "finresD3_minusN",
  "finresD3_plusN"
)

deg_summarybetween <- lapply(linesbetween, function(line) {
  lfc_col <- paste0("lfc_", line)
  padj_col <- paste0("padj_", line)
  
  df <- data.frame(
    lfc  = SingleNSingleDay_ljfunvsWT_merged[[lfc_col]],
    padj = SingleNSingleDay_ljfunvsWT_merged[[padj_col]]
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

b <- ggplot(df_longbetween, aes(x = Treatment, y = ifelse(Direction == "Up", Count, -Count), fill = Direction)) +
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5),
            fontface = "bold", size = 4.5
  ) +
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) +
  # scale_x_discrete(limits = c("DBL21", "DBL19", "HA11","HA3", "ADMD17" , "WT13")) + # custom order
  labs(
    title = "Upregulated and Downregulated DEGs for genotypic differences per Treatment and Day",
    x = "Genotypes (ljfun, WT)",
    y = "Gene Count",
    fill = "Regulation"
  ) +
  ylim(-1100, 1000) +
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
combined <- grid.arrange(a, b, ncol = 1)
ggsave(combined, file = "ljfun_DEGbarplotsboth.pdf", height = 25, width = 25, units = "cm")



#Step X: Venns #####
library(dplyr)
library(VennDiagram)
#install.packages("ggvenn")
library(ggvenn)
SingleNSingleDay_ljfunvsWT_merged
SingleNSingleDay_vufunvsWT_merged <- readxl::read_xlsx(path = "P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/C1_C24_VuFUN_timecourse_2025/quantsC1_C24/lfc_pval_vufun_VS_WT_combined.xlsx") # all genes



SingleNSingleDay_ljfunvsWT_merged_DEG_D1PlusN <- SingleNSingleDay_ljfunvsWT_merged %>%
  filter(padj_finresD1_plusN  < 0.05 ) %>% 
  pull(locusName) 
SingleNSingleDay_ljfunvsWT_merged_DEG_D3PlusN <- SingleNSingleDay_ljfunvsWT_merged %>%
  filter(padj_finresD3_plusN  < 0.05 ) %>%
  pull(locusName)

SingleNSingleDay_vufunvsWT_merged_DEG_D1PlusN <- SingleNSingleDay_vufunvsWT_merged %>%
  filter(padj_finresD1_plusN  < 0.05 ) %>% 
  pull(locusName) 
SingleNSingleDay_vufunvsWT_merged_DEG_D3PlusN <- SingleNSingleDay_vufunvsWT_merged %>%
  filter(padj_finresD3_plusN  < 0.05 ) %>%
  pull(locusName)




ljfun <- ggvenn(
  list(
    "D1"   = SingleNSingleDay_ljfunvsWT_merged_DEG_D1PlusN,
    "D3" = SingleNSingleDay_ljfunvsWT_merged_DEG_D3PlusN
  ),
  fill_color = c("#1B9E77", "#D95F02"),  # pastel complexvenn-like palette
  fill_alpha = 0.35,                     # soft transparency
  stroke_size = 0.4,                     # thin clean outlines
  stroke_color = "grey30",               # darker, crisp edge
  set_name_size = 5,                     # larger set labels
  #set_name_fontface = "bold",
  text_size = 5.2,                       # slightly larger counts
  text_color = "grey15",
  show_percentage = FALSE,
) +
  labs(title="ljfun-3 vs WT DEG overlap between D1 and D3")+
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ))
ljfun
vufun <- ggvenn(
  list(
    "D1"   = SingleNSingleDay_vufunvsWT_merged_DEG_D1PlusN,
    "D3" = SingleNSingleDay_vufunvsWT_merged_DEG_D3PlusN
  ),
  fill_color = c("#1B9E77", "#D95F02"),  # pastel complexvenn-like palette
  fill_alpha = 0.35,                     # soft transparency
  stroke_size = 0.4,                     # thin clean outlines
  stroke_color = "grey30",               # darker, crisp edge
  set_name_size = 5,                     # larger set labels
  #set_name_fontface = "bold",
  text_size = 5.2,                       # slightly larger counts
  text_color = "grey15",
  show_percentage = FALSE,
) +
  labs(title="vufun-2 vs WT DEG overlap between D1 and D3")+
  theme(
    plot.title = element_text(
      size = 14,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 8)
    ))


vufun

View(SingleNSingleDay_ljfunvsWT_merged_DEG_D3PlusN)


CombinedVenns <- grid.arrange(ljfun,vufun,nrow=2,ncol=1)
ggsave(CombinedVenns, file="CombinedVenns_ljfun3_vufun2.jpg",units="cm",height=15,width=15)
