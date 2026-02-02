#setwd
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/")
getwd()
#load libraries
library(tximport)
library(tidyverse)
library("seqinr")




#just want to check what we actually provided to salmon for mapping
transcript_fasta_data <- read.fasta(file = "P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.transcript.fa.gz", as.string = TRUE) #54484 annotated transcripts

transcript_fasta_data_df <- do.call(rbind, lapply(transcript_fasta_data, function(seq_obj) {
  data.frame(
    transcriptName = attr(seq_obj, "name"),
    Vuannotation = attr(seq_obj, "Annot"),
    ft_ni_sequence = as.character(seq_obj),
    stringsAsFactors = FALSE
  ) 
})) #54484 annotated transcripts

#let's get some annotation files
annotations <-  read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE) #54484 annotated transcripts

#so everything checks out in terms of numbers

#Now let's up the tximport rawfile
View(annotations)

#For accounting, first let's determine how many isoforms per transcript
duplicate_counts_df_counts <- as.data.frame(table(annotations$locusName)) |> #table() counts how many times each unique locusName appears in your annotations data frame.
  subset(Freq > 1) |>       # keep only loci with >1 transcript
  setNames(c("locusName", "Transcript_Isoforms")) #gives no. isoforms per transcript if there are >1

#Ok, so then let's count how many actual transcripts we have if summarizing to gene-level
data.table::uniqueN(annotations$locusName) #Number of unique transcripts, discards repeated rows per transcript #31948

#Ok, so let's summarize to gene-level, make sure it matches above number of 31948
annotations_GL <- annotations %>%
  dplyr::distinct(locusName, .keep_all = TRUE) #31948 more accurate because some genetranscripts might not have .1 read. Some start with .2

#Could be useful to know also how many isoforms per gene
annotations_GL <- left_join(annotations_GL,duplicate_counts_df_counts, by="locusName", relationship="one-to-one")
write.csv(annotations_GL,file="P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info_GL.csv",row.names = FALSE)

#Ok, let's generate the tx2gene annotation file
colnames(annotations)
tx2gene <- annotations[,c(3,2)]
head(tx2gene)
#Does any transcript ID map to more than one gene/locus?
any(duplicated(tx2gene[,1])) #should be false otherwise ur code is doodoo
locusName = sort(unique(tx2gene$locusName)) 
#read a samples file with all folder names in one column
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


#navigate to folders path
setwd("P:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/03.Mapped/ZBF1_ZBF28_ZBF_2025/quantsZBF1_ZBF28/")

#now do tximport
txi.salmon.tsv = tximport(files,
                          type = "salmon",
                          tx2gene = tx2gene,
                          dropInfReps = TRUE,
                          countsFromAbundance = "no")


#save into dataframe
final <- data.frame(txi.salmon.tsv)

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


metadata <- readxl::read_xlsx("Z:/LAB - OPV - Dugald_Reid/Genomics/RNAseq datasets/220126 cowpea lotus faba Novogene/X201SC25100870-Z01-F001/01.RawData/Copy of Sample Register.xlsx")
colnames(metadata)

ZBFmetadata <- (dplyr::filter(metadata,str_detect(`Exp Sample ID`, "ZBF")))
colnames(TPM)
print(ZBFmetadata$``)
print(ZBFmetadata$Description)

header=(c("locusName", c(paste("TPM",ZBFmetadata$`Identifying Allele`, sep=".")))
)

write.csv( TPM, file= "ZBFTPM.csv", row.names=F, quote=F)
getwd()
print(header)
colnames(TPM)=header


TPM_long = pivot_longer(TPM,
                        cols = -locusName,
                        names_to = 'sample',
                        values_to = 'TPM')

TPM_long$sample = str_remove(TPM_long$sample, pattern = ".rep1|.rep2|.rep3")

TPM_mean = TPM_long %>% 
  group_by(LOC, sample) %>% 
  summarise(TPM = mean(TPM)) %>% 
  pivot_wider(names_from = sample, values_from = TPM) #increases number of columns, decreasingrows. specify column to get name of output column from, as well as column to get cell values from.

write.csv(TPM_mean, 'TPM_mean.csv',row.names=F)

TPM_SE = TPM_long %>% 
  group_by(LOC, sample) %>% 
  summarise(SE = sd(TPM) / sqrt(n())) %>% 
  pivot_wider(names_from = sample, values_from = SE)

write.csv(TPM_SE, 'TPM_SE.csv',row.names=F)

