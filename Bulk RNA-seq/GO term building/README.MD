# Quick walkthrough for building GO term database

#### _Disclaimer_
This was written in 10 minutes there will be errors until I have time to proofread/test further.

## Necessary files
Your annotation file from same genome release you used for mapping

Annotation file:
- Vunguiculata_540_v1.2.P14.annotation_info

Packages:
- Building the database
```
 library(dplyr)
 library(tidyr)
 library(GO.db)
```

- Running the GO term analysis:
```
 library(clusterProfiler)
 library(GOSemSim)
```

## To run
```r
#Load in annotation fille
annotations <- read.delim("P:/LAB - OPV - Dugald_Reid/Genomics/Reference Genomes/Vunguiculata/v1.2/annotation/Vunguiculata_540_v1.2.P14.annotation_info.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#start cleaning to prepare dataset for building the GO term database for cowpea
df_filtered2 <- read.csv("Vunguiculata_540_v1.2.P14.annotation_info.csv")
gene_info <- df_filtered2[,c(2,10)]
colnames(gene_info)
gene_info <- dplyr::rename(gene_info,GID = geneID) #"toreplace"= "originalcolumn"
gene_info %>%
  filter(is.na(GID)) %>%
  distinct(GID)
gene_info$CHROMOSOME <- sub("Vigun([0-9]+)g.*", "\\1", gene_info$GID)
fsym <- gene_info[,c(1,3)]

#create GO table
safe_ontology <- function(goid) {
  term <- GOTERM[[goid]]
  if (is.null(term)) return(NA_character_)
  Ontology(term)
}



go_table <- df_filtered2 %>%
  dplyr::select(GID = geneID, GO) %>%
  dplyr::filter(GO != "") %>%
  tidyr::separate_rows(GO, sep = "[ ]+") %>%
    dplyr::mutate(
    EVIDENCE = "IEA",
    ontology = vapply(GO, safe_ontology, FUN.VALUE = character(1))
  )


go_table %>%
  filter(is.na(ontology)) %>%
  distinct(GO)

go_table <- go_table %>%
  filter(!GO %in% c("GO:0016573","GO:0016571","GO:0034968",
                       "GO:0016575","GO:0018142","GO:0035065","GO:0016570"))

colnames(go_table)

#table created, now make GO term db
library(AnnotationForge)

dir.create("output4")

file.remove("output3/org.Vigna.unguiculata.eg.sqlite")
unlink("output3/org.Vunguiculata.eg.sqlite", force = TRUE)
file.access("output3", 2)  # returns 0 if writable
colnames(go_table)
go_table <- dplyr::rename(go_table,EVIDENCE = evidence) #"toreplace"= "originalcolumn"
go_table2 <- go_table[,c("GID","GO","EVIDENCE")] #requirements for GO db

makeOrgPackage(
  gene_info = fsym,
  go        = go_table2,
  version   = "0.1",
  maintainer= "B Wee Yang <b.weeyang@latrobe.edu.au>",  # <- must be First Last <email>
  author    = "B Wee Yang",                              # <- no comma
  outputDir = "output4",                                 # <-output folder it'll be saved to
  tax_id    = "3917",                                    #taxonomy ID fo cowpea
  genus     = "Vigna",
  species   = "unguiculata",
  goTable   = "go",
  packageName = "org.Vunguiculata.eg.db"
)
  
#Now install the packages.
install.packages("output4/org.Vunguiculata.eg.db", repos = NULL, type = "source")
#install.packages(file.path(getwd(), "output4/org.Vigna.unguiculata.eg.db_0.1.tar.gz"),
#                 repos = NULL, type = "source")
#install.packages("org.Vunguiculata.eg.db") 
library(org.Vunguiculata.eg.db) # this needs to then be moved over to the same folder your R packages are saved in. 


# Quick enrichment test

# example gene list - replace with your interesting genes. Don't do this manually like below, just load in a list form your DEGs for example.
geneList <- c("Vigun01g000100", "Vigun02g012345", "Vigun03g000800")

gsea_res <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Vunguiculata.eg.db,
  keyType      = "ENTREZID",
  ont          = "BP",   # or "MF", "CC", "ALL"
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.05,
  verbose      = FALSE
)

#various visualizations
dotplot(gsea_res)
ridgeplot(gsea_res)
gseaplot2(gsea_res, geneSetID = 1)

gsea_simplified <- simplify(
  gsea_res,
  cutoff      = 0.7,
  by          = "p.adjust",
  select_fun  = min,
  measure     = "Wang"
)

dotplot(gsea_simplified)
ridgeplot(gsea_simplified)
gseaplot2(gsea_simplified, geneSetID = 1)










#Random shit
str(go_table)

go_table <- df_filtered2 %>%
  dplyr::select(gene_id = geneID, GO) %>%
  dplyr::filter(GO != "") %>%
  tidyr::separate_rows(GO, sep = "[ ]+") %>%
  dplyr::rename(go_id = GO) %>%
  dplyr::mutate(
    evidence = "IEA",
    ontology = vapply(go_id, function(x) {
      out <- Ontology(GOTERM[[x]])
      if (is.null(out)) NA_character_ else out
    }, FUN.VALUE = character(1))
  )

  
  str(gene_info)

gene2go <- gene_info %>%
  select(geneID, GO) %>%
  epl64;;filter(GO != "") %>%                 # remove empty GO entries
  separate_rows(GO, sep = "\\s+") %>%  # split multiple GO terms into separate rows
  distinct()                           # remove duplicates just in case

library(clusterProfiler)

enrich_res <- enricher(
  gene = your_gene_vector,   # your list of interesting genes
  TERM2GENE = gene2go
)
```