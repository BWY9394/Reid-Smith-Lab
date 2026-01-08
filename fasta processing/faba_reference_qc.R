setwd("Z:/LAB - OPV - Dugald_Reid/Genomics/faba")
getwd() #confirm is same as the above. Change the directory address to whereever your raw file locations are i.e. the fasta files

install.packages("seqinr") #same thing to install dplyr if n
library("seqinr")
library(dplyr)

#Alternative step if using Biostrings, but honestly just use seqinr if possible
#prot_fasta_data <-readAAStringSet("HEDIN_TMP5.check.renamed.phase.complete_proteins.fa") # Use readAAStringSet for protein sequences,  readDNAStringSet for DNA
#cds_fasta_data <- readDNAStringSet("HEDIN_TMP5.check.renamed.phase.complete_CDS.fa")


#seqinr approach
##protein loaded in as well just out of interest 
prot_fasta_data <- read.fasta(file = "HEDIN_TMP5.check.renamed.phase.complete_proteins.fa", as.string = TRUE)
##CDS is the real meat here
###load CDS file
cds_fasta_data <- read.fasta(file = "HEDIN_TMP5.check.renamed.phase.complete_CDS.fa", as.string = TRUE)
###Convert to data frame
cds_df <- data.frame(
  id  = names(cds_fasta_data),
  seq = sapply(cds_fasta_data, paste0, collapse = ""),
  stringsAsFactors = FALSE
)
###always good to check if your fasta file has repeats i.e. is fucked from the start
duplicated_ids <- cds_df$id[data.table::duplicated(cds_df$id)]

###save if you want, otherwise skip step
write.csv(cds_df,file="HEDIN_TMP5_check_renamed_phase_complete_CDS.csv",row.names=FALSE)

#now ready to filter
##Keep only longest CDS per gene
cds_longest <- cds_df %>%
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

##convert back to fasta
seqs_split <- lapply(cds_longest$seq, function(x) strsplit(x, "")[[1]]) #splits seq column entries into separated character values for wrapping functionality when saving back into fasta

##save into fasta
write.fasta(
  sequences = seqs_split,
  names     = cds_longest$id,
  nbchar    = 70,
  file.out  = "HEDIN_TMP5.check.renamed.phase.longest_CDS.fa",
  )
##sanity check
nrow(cds_df)
length(read.fasta("HEDIN_TMP5.check.renamed.phase.longest_CDS.fa"))

