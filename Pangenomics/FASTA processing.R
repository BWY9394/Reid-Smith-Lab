# Script to just load and check fasta sequences

setwd("P:/LAB - OPV - Dugald_Reid/Genomics/panfaba/annotations")
getwd()


install.packages("seqinr")
library('Biostrings')
library("seqinr")
library(dplyr)
prot_fasta_data <-readAAStringSet("HEDIN_TMP5.check.renamed.phase.complete_proteins.fa") # Use readAAStringSet for protein sequences,  readDNAStringSet for DNA
cds_fasta_data <- readDNAStringSet("HEDIN_TMP5.check.renamed.phase.complete_CDS.fa")


#seqinr approach
species3_fasta <- read.fasta(file = "./Hedin2/VFABA.HEDIN2.pgsb.r3.Oct2025.longest.aa.fa", as.string = TRUE)

species7_fasta <- read.fasta(file = "./NPZ3/VFABA.NPZ3.pgsb.r1.Oct2025.longest.aa.fa", as.string = TRUE)


species3_AA_df <- data.frame(
  id  = names(species3_fasta),
  seq = sapply(species3_fasta, paste0, collapse = ""),
  stringsAsFactors = FALSE
)

species7_AA_df <- data.frame(
  id  = names(species7_fasta),
  seq = sapply(species7_fasta, paste0, collapse = ""),
  stringsAsFactors = FALSE
)


#Use below if need to filter for longest transcript per gene
cds_longest <- cds_df %>% #now create helper columns gene_id and seq_len
mutate(
  gene_id = sub("\\.[0-9]+$", "", id), #create gene_id column and then replace whatever comes after \\. with a space ("")
  seq_len = nchar(seq) #counts the number of characters in nchar(seq)
) %>%
  group_by(gene_id) %>% #now group by gene_id
  slice_max(seq_len, n = 1, with_ties = FALSE) %>%
  #slice_max selec row swith smallest or largest values of a variable. n= no. rows to keep per group. if 2 transcripts have same length, keep only one.
  ungroup() %>% #removes grouping metadata
  select(id, seq) %>% #drops helper columns created.
  mutate(seq = toupper(seq)) #converts seq column entries to caps
#conver back to fasta

str(cds_longest)

#seqs_split <- lapply(species7_AA_df$seq, function(x) strsplit(x, "")[[1]]) # now defunct

#splits seq column entries into separated character values for wrapping functionality when saving back into fasta #substitute the df for whatever applies.
seqs_split <- setNames(
  lapply(species7_AA_df$seq, function(x) strsplit(x, "")[[1]]),
  species7_AA_df$id
) 


write.fasta(
  sequences = seqs_split,
  names     = names(seqs_split),
  nbchar    = 70,
  file.out  = "youredittedfastafile.fa"
)



getwd()
#now defunct
#write.fasta(
 # sequences = seqs_split,
  #names     = species7_AA_df$id,
  #nbchar    = 70,
  #file.out  = "HEDIN_TMP5.check.renamed.phase.longest_CDS.fa",
#)

