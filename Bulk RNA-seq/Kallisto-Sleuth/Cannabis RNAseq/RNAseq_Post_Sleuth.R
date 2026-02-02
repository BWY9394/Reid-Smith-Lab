BiocManager::install("clusterProfiler")



#Required Packages
library("tidyverse")
library("ggfortify")
library("ggrepel")
library("gridExtra")
library("dplyr")
library("reshape2")

library("ggforce")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("clusterProfiler")
library(GOSemSim)
library(enrichplot)
library(viridis)



#Working Directory
setwd("C:/Users/BWeeY/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
setwd("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
setwd("C:/Users/OBerkowtiz/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/")
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/RNAseq/HP0001/Kallisto Output HP01/")

#Computed TPMs from Kalisto workflow and gene-level annotation file
tpmmeans =read.csv("TPM_mean.csv",header=T, stringsAsFactors = T)
tpmraw=read.csv("TPM.csv",header=T,stringsAsFactors = T)
tpmse=read.csv("TPM_SE.csv",header=T, stringsAsFactors = T)
tpmselong=pivot_longer(tpmse,cols=2:21,names_to="ID",values_to="tpmSE")
countsraw=read.csv("counts.csv",header=T, stringsAsFactors = T)
annot=read.csv("annot_genelevel.csv",header=T)


#Making a base reference table for TPMs coming out of Kallisto workflow
TPMmeans_reftable <- left_join(tpmmeans,annot,by="LOC")
write.csv(TPMmeans_reftable,file="TPMmeans_reftable.csv",row.names=FALSE) #note to self that I rearrange it manually in excel +renamed the anotation column

TPMmeans_reftable <- read.csv("TPMmeans_reftable.csv")

feature2add <- read.csv("HP01_feature_table.csv")
colnames(feature2add)
feature2add=feature2add[,c(8,6)]
colnames(feature2add)
feature2add <- feature2add %>% distinct(LOC,.keep_all = TRUE)
addedfeature=inner_join(TPMmeans_reftable,feature2add,by="LOC")

write.csv(addedfeature,file="addedfeature.CSV",row.names=FALSE) #so this is your tpms with no filters.
addedfeature <- read.csv("addedfeature.CSV")


annot_with_chromosome <-left_join(annot,feature2add,by="LOC")  #Useful to have the gene level annot also have chromosome number on it
write.csv(annot_with_chromosome,file="annot_genelevel_with_chromosome.csv",row.names=FALSE) # THIS HAS gene level annotations with chromosome attached. useful for making ref tables.

#Applying the TPM>10 in at least one tissue filter to DEG genes. Criteria: FDR<0.05, log2FC>1, <-1, per TISSUE TYPE across treatments
DEGtable <- read.csv("filt_ALL.csv",stringsAsFactors = TRUE,header=TRUE) #11k DEG (counts) 
DEGtable_mod <- read.csv("filt_ALL_no_0P.csv",stringsAsFactors = TRUE,header=TRUE) #5.3k DEG (counts)

mod=DEGtable[,c(1,2,3)]
mod_without0POL=DEGtable_mod[,c(1,2,3)]

DEG_TPM_all <- left_join(mod,y=TPMmeans_reftable[,c(1,4:23)],by="LOC")
write.csv(DEG_TPM_all,"DEG_TPM_all.csv",row.names=FALSE)

DEG_TPM_all_no0POL <- left_join(mod_without0POL,y=TPMmeans_reftable[,c(1,4:23)],by="LOC")
write.csv(DEG_TPM_all_no0POL,"DEG_TPM_all_no0POL.CSV",row.names=FALSE)


DEGtable <-  read.csv("DEG_TPM_all.csv",stringsAsFactors = TRUE,header=TRUE) #11224 DEG
DEGtable_mod <- read.csv("DEG_TPM_all_no0POL.CSV",stringsAsFactors = TRUE,header=TRUE) #5341 DEG

colnames(DEGtable)

DEGtable$LRtotal=rowSums(DEGtable[,4:7])
DEGtable$OLtotal=rowSums(DEGtable[,8:11])
DEGtable$MLtotal=rowSums(DEGtable[,12:15])
DEGtable$YStotal=rowSums(DEGtable[,16:19])
DEGtable$TFtotal=rowSums(DEGtable[,20:23])
colnames(DEGtable)


colnames(DEGtable_mod)
DEGtable_mod$LRtotal=rowSums(DEGtable_mod[,4:7])
DEGtable_mod$OLtotal=rowSums(DEGtable_mod[,8:11])
DEGtable_mod$MLtotal=rowSums(DEGtable_mod[,12:15])
DEGtable_mod$YStotal=rowSums(DEGtable_mod[,16:19])
DEGtable_mod$TFtotal=rowSums(DEGtable_mod[,20:23])
colnames(DEGtable_mod)


colorder <- c("LOC","Feature","Annotation","LR_0","LR_0.25","LR_1","LR_2","OL_0","OL_0.25","OL_1","OL_2",
              "ML_0","ML_0.25","ML_1","ML_2","YS_0","YS_0.25","YS_1","YS_2","TF_0","TF_0.25","TF_1","TF_2")


tester2 <- DEGtable[rowSums(DEGtable[, c(24:28)] > 10) > 0, ] #10327, matches output below
tester2 = tester2%>% dplyr::rename_with(~colorder,.cols=c(1:23))
colnames(tester2)
write.csv(tester2[,1:23],file="DEG_TPM_all_atleast10TPM.CSV",row.names=FALSE)

tester <- DEGtable_mod[rowSums(DEGtable_mod[, c(24:28)] > 10) > 0, ] #4871 genes, matches output below
tester = tester%>% dplyr::rename_with(~colorder,.cols=c(1:23))
write.csv(tester[,1:23],file="DEG_TPM_all_no_0P_atleast10TPM.CSV",row.names=FALSE)


#this bit is to make the list for the write-up.
tester3 <- left_join(tester,DEGtable, by="LOC",relationship = "one-to-one")
tester3 <- tester3[,-c(3:29)]
tester <- read.csv("DEG_TPM_all_atleast10TPM.CSV",stringsAsFactors = TRUE)

ref_pi<- read.csv("phosphate_genes.csv")
ref_pi <-dplyr:: rename(ref_pi,"LOC"="gene")
ref_pi$list <-"phosphate"
ref_n <- read.csv("nitrogen_genes.csv")
ref_n <- dplyr::rename(ref_n,"LOC"="gene")
ref_n$list <-"nitrate"
ref_s <- read.csv("sulfur_genes.csv")
ref_s <- dplyr::rename(ref_s,"LOC"="gene")
ref_s$list <-"sulfate"
ref_fe <- read.csv("iron_genes.csv")
ref_fe <- dplyr::rename(ref_fe,"LOC"="gene")
ref_fe$list <-"iron"
ref_cann <- read.csv("secondarymetabolite_genes.csv")
ref_cann <- dplyr::rename(ref_cann,"LOC"="gene")
ref_cann$list <-"2ndMet"
ref_all <- full_join(ref_pi,ref_n) %>% full_join(ref_s) %>% full_join(ref_fe) %>% full_join(ref_cann)

tester4 <- left_join(tester3,ref_all,by="LOC",relationship="one-to-many")
write.csv(tester4,file="DEG_endpipe.csv",row.names = FALSE)


#sanity check for the above codes by going from raw TPM files first before merging with DEG list
DEGTPMmeans=read.csv("TPM_mean.csv")

DEGTPMmeans$LRtotal=rowSums(DEGTPMmeans[,2:5])
DEGTPMmeans$MLtotal=rowSums(DEGTPMmeans[,6:9])
DEGTPMmeans$OLtotal=rowSums(DEGTPMmeans[,10:13])
DEGTPMmeans$TFtotal=rowSums(DEGTPMmeans[,14:17])
DEGTPMmeans$YStotal=rowSums(DEGTPMmeans[,18:21])
colnames(DEGTPMmeans)

DEGTPMmeansfilt <- DEGTPMmeans[,c(1,22,23,24,25,26)] #29283 genes
DEGTPMmeansfilt <- DEGTPMmeansfilt[rowSums(DEGTPMmeansfilt[, -1] > 10) > 0, ] #18,665 genes NOT DEG yet but filtered to keep genes only >10 in at least 1 organ type

DEGTPMfiltered <- DEGTPMmeansfilt[,1]
DEGTPMfiltered <- as.data.frame(DEGTPMfiltered)
DEGTPMfiltered <- dplyr::rename(DEGTPMfiltered,"LOC"="DEGTPMfiltered") #18,665 genes NOT DEG yet.
write.csv(DEGTPMfiltered,file="Countsmorethan10TPM.csv",row.names=FALSE) #This file is NOT DEG.  Filtered for expression TPM>10 across all samples.

testingfilter2 <- inner_join(DEG_TPM_all,DEGTPMfiltered, by="LOC", relationship = "one-to-one") #10327 genes.Needs to match the above pipe

testingfilter <- inner_join(DEG_TPM_all_no0POL,DEGTPMfiltered, by="LOC",relationship = "one-to-one") #4871 genes.Needs to match the above pipe



#Reference List from above workflow. Criteria: FDR<0.05, log2FC>1, <-1, per TISSUE TYPE across treatments. TPM>10 per TISSUE TYPE across treatments as well.
DEG <- read.csv("DEG_TPM_all_no_0P_atleast10TPM.CSV")
DEG2 <- read.csv("DEG_TPM_all_atleast10TPM.csv") 

#modified gene lists to include respective list identifier for merging with DEG2 or "DEG_TPM_all_atleast10TPM.csv"
# you will merge this output of ref_all with an appropriate file which contains all the lfc values, e.g. filt_all
# so this will net you a list that has a) genes that are DEG and meet TPM filter criteria, while also containing LFC values, as well as gene list specific sorting capabilities.
ref_pi<- read.csv("phosphate_genes.csv")
ref_pi$gene.list <- "phosphate"
ref_pi <-dplyr:: rename(ref_pi,"LOC"="gene")
ref_n <- read.csv("nitrogen_genes.csv")
ref_n <- dplyr::rename(ref_n,"LOC"="gene")
ref_n$gene.list <- "nitrate"
ref_s <- read.csv("sulfur_genes.csv")
ref_s <- dplyr::rename(ref_s,"LOC"="gene")
ref_s$gene.list <- "sulfate"
ref_fe <- read.csv("iron_genes.csv")
ref_fe <- dplyr::rename(ref_fe,"LOC"="gene")
ref_fe$gene.list <- "fe"
ref_cann <- read.csv("secondarymetabolite_genes.csv")
ref_cann <- dplyr::rename(ref_cann,"LOC"="gene")
ref_cann$gene.list <- "cann"
ref_all <- full_join(ref_pi,ref_n) %>% full_join(ref_s) %>% full_join(ref_fe) %>% full_join(ref_cann)



#Gene lists
ref_pi<- read.csv("phosphate_genes.csv")
ref_pi <-dplyr:: rename(ref_pi,"LOC"="gene")
ref_n <- read.csv("nitrogen_genes.csv")
ref_n <- dplyr::rename(ref_n,"LOC"="gene")
ref_s <- read.csv("sulfur_genes.csv")
ref_s <- dplyr::rename(ref_s,"LOC"="gene")
ref_fe <- read.csv("iron_genes.csv")
ref_fe <- dplyr::rename(ref_fe,"LOC"="gene")
ref_cann <- read.csv("secondarymetabolite_genes.csv")
ref_cann <- dplyr::rename(ref_cann,"LOC"="gene")
ref_all <- full_join(ref_pi,ref_n) %>% full_join(ref_s) %>% full_join(ref_fe) %>% full_join(ref_cann)

#Making a gene lists with transformed TPM values.These will go into the master summary xls file. "HP0001_RNAseq_summary_tables 
ref_pi_TPMsattached <- inner_join(ref_pi,TPMmeans_reftable[,c(1,4:23)],by="LOC") #so you start of with a list of genes that are not DE.
view(ref_pi_TPMsattached)
DEGs_pi <- left_join(ref_pi,DEG,by="LOC") # now you have another specific list where genes are DE. NAs across rows represent genes whereby it was not DE.
ref_pi_TPMsattached$DE <-DEGs_pi$name  #Using the name column from that DE list to identify which ones are DE. same order, so can just paste.
view(ref_pi_TPMsattached)
ref_pi_TPMsattached[, c(6:25)] <- apply(ref_pi_TPMsattached[, c(6:25)], 2, function(x) log2(x + 1))
write.csv(ref_pi_TPMsattached,file = "ref_pi_TPMsattached.csv",row.names=FALSE)

ref_n_TPMsattached <- inner_join(ref_n,TPMmeans_reftable[,c(1,4:23)],by="LOC")
view(ref_n_TPMsattached)
DEGs_n <- left_join(ref_n,DEG,by="LOC")
ref_n_TPMsattached$DE <-DEGs_n$name
ref_n_TPMsattached[, c(5:24)] <- apply(ref_n_TPMsattached[, c(5:24)], 2, function(x) log2(x + 1))
write.csv(ref_n_TPMsattached,file = "ref_n_TPMsattached.csv",row.names=FALSE)

ref_s_TPMsattached <- inner_join(ref_s,TPMmeans_reftable[,c(1,4:23)],by="LOC")
view(ref_s_TPMsattached)
DEGs_s <- left_join(ref_s,DEG,by="LOC")
ref_s_TPMsattached$DE <-DEGs_s$name
ref_s_TPMsattached[, c(5:24)] <- apply(ref_s_TPMsattached[, c(5:24)], 2, function(x) log2(x + 1))
write.csv(ref_s_TPMsattached,file = "ref_s_TPMsattached.csv",row.names=FALSE)

ref_fe_TPMsattached <- inner_join(ref_fe,TPMmeans_reftable[,c(1,4:23)],by="LOC")
view(ref_fe_TPMsattached)
DEGs_fe <- left_join(ref_fe,DEG,by="LOC")
ref_fe_TPMsattached$DE <-DEGs_fe$name
ref_fe_TPMsattached[, c(5:24)] <- apply(ref_fe_TPMsattached[, c(5:24)], 2, function(x) log2(x + 1))
write.csv(ref_fe_TPMsattached,file = "ref_fe_TPMsattached.csv",row.names=FALSE)

ref_cann_TPMsattached <- inner_join(ref_cann,TPMmeans_reftable[,c(1,4:23)],by="LOC")
view(ref_cann_TPMsattached)
DEGs_cann <- left_join(ref_cann,DEG,by="LOC")
ref_cann_TPMsattached$DE <-DEGs_cann$name
ref_cann_TPMsattached[, c(6:25)] <- apply(ref_cann_TPMsattached[, c(6:25)], 2, function(x) log2(x + 1))
write.csv(ref_cann_TPMsattached,file = "ref_cann_TPMsattached.csv",row.names=FALSE)

#Now to structure the DEG tables with TPM cutoffs. Criteria: FDR<0.05, log2FC>1, <-1, per TISSUE TYPE across treatments. TPM>10 per TISSUE TYPE across treatments as well. Also want to add in reference genes.
DEG <- read.csv("DEG_TPM_all_no_0P_atleast10TPM.CSV")
DEG2 <- read.csv("DEG_TPM_all_atleast10TPM.csv") #why is there a tpm value for genes that are specifically expressed in this list

colnames(DEG)
colnames(DEG2)

TPMmeans_reftable <- read.csv("TPMmeans_reftable.csv")
#TPMmeans_reftable[, 4:ncol(TPMmeans_reftable)] <- apply(TPMmeans_reftable[, 4:ncol(TPMmeans_reftable)], 2, function(x) log2(x + 1))


row_to_append <- TPMmeans_reftable[c(20394,17078,18772), ]   # Select the row from df1
colnamestorename <- colnames(DEG)

colnames(row_to_append) <- colnamestorename # Get the column names of df1
colnames(row_to_append)

df1 <- rbind(DEG, row_to_append) #paste in the reference genes


colnames(DEG)
df1[, 4:ncol(df1)] <- apply(df1[, 4:ncol(df1)], 2, function(x) log2(x + 1))  #log transformations
write.csv(df1, file="data_heatmap_base_TPM10.csv",row.names=FALSE)


df2 <- df1
df2 <- cbind(df2[,c(1:3)],t(apply(df2[,-c(1:3)], 1, function(x) (x - mean(x)) / sd(x)))) #zscoring per gene across all samples
write.csv(df2, file="data_heatmap_zscored.csv_TPM10.csv",row.names=FALSE)

#now also want to to generate a list of the genes that are in OL0P only. should have 5457 genes only expressed in OL0P.Let's ID them.
df3 <- full_join(DEG2,DEG[,1:2],by="LOC")
colnames(df3)

df3_with_na <- subset(df3, is.na(Feature.y))
df3_with_na <- df3_with_na[,-24]
colnames(df3_with_na)
df3_with_na[, 4:ncol(df3_with_na)] <- apply(df3_with_na[, 4:ncol(df3_with_na)], 2, function(x) log2(x + 1))
rownames(df3_with_na) <- NULL

df3_with_na_zscore <- df3_with_na
df3_with_na_zscore <- cbind(df3_with_na_zscore[,c(1:3)],t(apply(df3_with_na_zscore[,-c(1:3)], 1, function(x) (x - mean(x)) / sd(x))))


#PCA Analysis
library("gridExtra")
library("PCAtools")
library("scales")
library("pcaExplorer")
#library("sleuth")
library("plyr")


##sampleinfo file for grouping 
sampleinfo2=read.csv("sampleinfo2.csv",header=T,stringsAsFactors = T) #for tpm_mean.csv
sampleinfo2$Replicate <- as.factor(sampleinfo2$Replicate)
sampleinfo2$Treatment <- as.factor(sampleinfo2$Treatment)
sampleinfo2$Organ <- as.factor(sampleinfo2$Organ)
sampleinfo=read.csv("sampleinfo.csv",header=T,stringsAsFactors = T) #for tpm.csv
sampleinfo$Replicate <- as.factor(sampleinfo$Replicate)
sampleinfo$Treatment <- as.factor(sampleinfo$Treatment)
sampleinfo$Organ <- as.factor(sampleinfo$Organ)
sampleinfo3=read.csv("sampleinfo3.csv",header=T, stringsAsFactors = T) #tpm_mean.csv but reordered for corrected treatment levels.
#PCA analysis
#remove non-tpm values
tpmraw_pca <- tpmraw[,-1]

#normalisation
logtpm <- log2(tpmraw_pca+1)
all(colnames(logtpm) == rownames(sampleinfo$ID))#see if order matches that of sample info

#subsetting into differnt organ types
LR_logtpm <- logtpm[,c(1:12)]
ML_logtpm <- logtpm[,c(13:24)]
YS_logtpm <- logtpm[,c(25:36)]
TF_logtpm <- logtpm[,c(37:48)]
YS_logtpm <- logtpm[,c(49:60)]


#Calculating PCs for all and per organ type
pcDat <- prcomp(t(logtpm))
LRpcDat <- prcomp(t(LR_logtpm))
MLpcDat <- prcomp(t(ML_logtpm))
YSpcDat <- prcomp(t(YS_logtpm))
TFpcDat <- prcomp(t(TF_logtpm))
YSpcDat <- prcomp(t(YS_logtpm)) 

#trying

#calculate % total variance as explained by each principal component
percentVar <- pcDat$sdev^2 / sum(pcDat$sdev^2 )*10

#plot pca based on all data 
a=autoplot(pcDat,
         data=sampleinfo2,
         shape="Treatment", size=1.5,
         color="Organ",
         #label=TRUE, label.label="ID", label.repel=T, #label.size= 1.5, label.repel=T,
         frame=TRUE, 
         #frame.type="norm"
)+ #ggforce::geom_mark_ellipse(aes(fill = Organ,color=Organ))+
  ggtitle("PCA_all")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )
a
ggsave(a,file="pca_ALL.pdf",width=28, height =15, units="cm")


#plot out PCA select data frame based on target organ.
autoplot(MLpcDat)
OL <- autoplot(autoplot(OLpcDat,
         data = subset(sampleinfo2,Organ == "OL"), #change to fit specific tissue.
         colour="Treatment", 
         shape="Replicate", size=2.5,#if looking at organ specific level then change to replicate
         #label = TRUE, label.label = "ID",label.size=4,label.repel = T,
         #frame = TRUE, frame.type = "norm", remove hashtag when figure how to make frame with with just 3 datapoints
         ))+ ggforce::geom_mark_ellipse(aes(fill = Treatment,
                                          color = Treatment))+
        ggtitle("OL")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )


LR <- autoplot(autoplot(LRpcDat,
                        data = subset(sampleinfo2,Organ == "LR"), #change to fit specific tissue.
                        colour="Treatment", 
                        shape="Replicate", size=2.5,#if looking at organ specific level then change to replicate
                        #label = TRUE, label.label = "ID",label.size=4,label.repel = T,
                        #frame = TRUE, frame.type = "norm", remove hashtag when figure how to make frame with with just 3 datapoints
))+ ggforce::geom_mark_ellipse(aes(fill = Treatment,color = Treatment))+
  ggtitle("LR")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )

LR
YS <- autoplot(autoplot(YSpcDat,
                        data = subset(sampleinfo2,Organ == "YS"), #change to fit specific tissue.
                        colour="Treatment", 
                        shape="Replicate", size=2.5,#if looking at organ specific level then change to replicate
                        #label = TRUE, label.label = "ID",label.size=4,label.repel = T,
                        #frame = TRUE, frame.type = "norm", remove hashtag when figure how to make frame with with just 3 datapoints
)+ ggforce::geom_mark_ellipse(aes(fill = Treatment,
                                  color = Treatment)))+ggtitle("YS")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )
TF <- autoplot(autoplot(TFpcDat,
                        data = subset(sampleinfo2,Organ == "TF"), #change to fit specific tissue.
                        colour="Treatment", 
                        shape="Replicate", size=2.5,#if looking at organ specific level then change to replicate
                        #label = TRUE, label.label = "ID",label.size=4,label.repel = T,
                        #frame = TRUE, frame.type = "norm", remove hashtag when figure how to make frame with with just 3 datapoints
))+ ggforce::geom_mark_ellipse(aes(fill = Treatment,
                                  color = Treatment))+ggtitle("TF")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )
ML <- autoplot(autoplot(MLpcDat,
                        data = subset(sampleinfo2,Organ == "ML"), #change to fit specific tissue.
                        colour="Treatment", 
                        shape="Replicate", size=2.5,#if looking at organ specific level then change to replicate
                        #label = TRUE, label.label = "ID",label.size=4,label.repel = T,
                        #frame = TRUE, frame.type = "norm", #remove hashtag when figure how to make frame with with just 3 datapoints
))+ ggforce::geom_mark_ellipse(aes(fill = Treatment,
                                  color = Treatment))+ggtitle("ML")+theme_bw()+
  theme(axis.title.x = element_text(color="black", size=10, face="bold"),
        axis.title.y = element_text(color="black", size=10, face="bold"),
        axis.text.x = element_text(size = 8, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 8, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 8, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        #panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x=element_blank()
  )
ML


d <- grid.arrange(LR,OL,ML,YS,TF,ncol=3,nrow=2)

ggsave(d,file="pca_perorgan2.pdf",height =25, width=30, units="cm")



#Heatmaps####
#START HERE FOR HEATMAP PLOTS IF USING THE TPM10 FILTER.
dataHM_base <- read.csv("data_heatmap_base_TPM10.csv") 
dataHM_ZS <- read.csv("data_heatmap_zscored.csv_TPM10.csv")
dataHM_TAU <- read.csv("tau.csv")
ref_pi<- read.csv("phosphate_genes.csv")
ref_pi <-dplyr:: rename(ref_pi,"LOC"="gene") #remove LOC115712693 
View(ref_pi)
ref_pi <- ref_pi[-79,]
ref_n <- read.csv("nitrogen_genes.csv")
ref_n <- dplyr::rename(ref_n,"LOC"="gene") #alreadyy added
ref_s <- read.csv("sulfur_genes.csv")
ref_s <- dplyr::rename(ref_s,"LOC"="gene")
ref_fe <- read.csv("iron_genes.csv")
ref_fe <- dplyr::rename(ref_fe,"LOC"="gene")
ref_cann <- read.csv("secondarymetabolite_genes.csv")
ref_cann <- dplyr::rename(ref_cann,"LOC"="gene")
ref_all <- read.csv("ref_all.csv")
write.csv(ref_all,"ref_all.csv",row.names=FALSE)

dataHM_ZS_pi <- left_join(ref_pi,dataHM_ZS,by="LOC")
dataHM_ZS_s <-  left_join(ref_s,dataHM_ZS,by="LOC")
dataHM_ZS_n <- left_join(ref_n,dataHM_ZS,by="LOC") 
dataHM_ZS_fe <- left_join(ref_fe,dataHM_ZS,by="LOC")
dataHM_ZS_cann <- left_join(ref_cann,dataHM_ZS,by="LOC")

dataHM_pi <- left_join(ref_pi,dataHM_base,by="LOC")
dataHM_n <- left_join(ref_n,dataHM_base,by="LOC")
dataHM_s <- left_join(ref_s,dataHM_base,by="LOC")
dataHM_fe <- left_join(ref_fe,dataHM_base,by="LOC")
dataHM_cann <- left_join(ref_cann,dataHM_base,by="LOC")

dataHMtau <- left_join(ref_pi,dataHM_TAU,by="LOC")
dataHMtau2 <- left_join(ref_n,dataHM_TAU,by="LOC")


dataHM_ZS_feat <-  full_join(dataHM_ZS_pi, dataHM_ZS_s)%>%
  full_join(dataHM_ZS_n)

#remove NAs
dataHM_ZS_feat <- dataHM_ZS_feat %>% filter_at(vars(8:ncol(.)), any_vars(!is.na(.))) 
dataHM_ZS_pi <- dataHM_ZS_pi%>%  filter_at(vars(8:ncol(.)), any_vars(!is.na(.)))
dataHM_pi<- dataHM_pi%>%  filter_at(vars(8:ncol(.)), any_vars(!is.na(.)))
dataHM_ZS_fe <- dataHM_ZS_fe%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_fe<- dataHM_fe%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_ZS_n <- dataHM_ZS_n%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_n<- dataHM_n%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_ZS_s <- dataHM_ZS_s%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_s <- dataHM_s%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_ZS_cann <- dataHM_ZS_cann%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
dataHM_cann <- dataHM_cann%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))

dataHMtau <- dataHMtau%>%  filter_at(vars(8:ncol(.)), any_vars(!is.na(.)))
dataHMtau2 <- dataHMtau2%>%  filter_at(vars(7:ncol(.)), any_vars(!is.na(.)))
#Ok, now is where you are selecting what to plot. So choose your biological question
#setting row grouping order
row_order=data.frame(ID=dataHM_pi$LOC,
                     Grouping=dataHM_pi$gene.family) #pick what you are plotting

row_order2=data.frame(ID=dataHM_n$LOC,
                     Grouping=dataHM_n$gene.family) #pick what you are plotting


print(row_order$Grouping)
#$Grouping=factor(row_order$Grouping)
#row_order$ID=as.factor(row_order$ID)                    

row_order$Grouping=factor(row_order$Grouping,levels=c("PHT1","PHT3","PHT4","PHT5/VPT","PHF","VPE","PHO1","PHO2","NLA","ANR","WRKY","MYB","PHR1","SPX","PP-InsP","MGD","DGD","SQD","NMT/PEAMT","GDPD","PLD","PECP","RNS","PAP","TAP46","RNF4")) #P

row_order2$Grouping=factor(row_order2$Grouping,levels=c("AMT","NPF","NRT3","AAT","NIA","NIR","GLN","GDH","CEP","CIPK","ANR1","NIGT","NLP","Ref")) #N


row_order$Grouping=factor(row_order$Grouping,levels=c("SULTR","APS","SAT","SDI","LSU1","SQD2","Ref")) #S


row_order$Grouping=factor(row_order$Grouping,levels=c("PRI1","IDF1","bHLH","FRO2","NAS","YSL","FRD","FPN","TCR","VTL","AHA","PDR9","IRT1","Ref")) #Fe


row_order$Grouping=factor(row_order$Grouping,levels=c("CBD","THC","IPP","Ref")) #Cann

levels(row_order$Grouping)

#ok time to make col_order.
colnames(df_pi)
col_order=data.frame(ID=c("LR:P0","LR:P0.25","LR:P1","LR:P2",
                          "OLP:0","OLP:0.25","OLP:1","OLP:2",
                          "MLP:0","MLP:0.25","MLP:1","MLP:2",
                          "YSP:0","YSP:0.25","YSP:1","YSP:2",
                          "TFP:0","TFP:0.25","TFP:1","TFP:2"),
                     Grouping=c(rep("LR",4),rep("OL",4),rep("ML",4),rep("YS",4),rep("TF",4))
)
col_order$Grouping=factor(col_order$Grouping,levels=c("LR","OL","ML","YS","TF"))
levels(col_order$Grouping)

col_ordertau=data.frame(Grouping=c("tau0P","tau025P","tau1P","tau2P"))
col_ordertau$Grouping=factor(col_ordertau$Grouping,levels=c("tau0P","tau025P","tau1P","tau2P"))


#melt your gene symbol into row numbers as row names.

justOL0P <- df3_with_na[,c(1,2,3,8)]

rownames(justOL0P) <- NULL

senescinggenesTPM <- df3_with_na %>% 
  tibble:: column_to_rownames(var = "LOC")

senescinggenesjustOL0P <- justOL0P %>% 
  tibble:: column_to_rownames(var = "LOC")



senescinggenesZscore <- df3_with_na_zscore %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_base <- dataHM_base %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_ZS <- dataHM_ZS %>% 
  tibble:: column_to_rownames(var = "LOC")

dataHM_ZS_pi <- dataHM_ZS_pi %>% 
  tibble:: column_to_rownames(var = "LOC")

dataHM_pi <- dataHM_pi %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_ZS_n <- dataHM_ZS_n %>% 
  tibble:: column_to_rownames(var = "LOC")

dataHM_n <- dataHM_n %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_ZS_s <- dataHM_ZS_s %>% 
  tibble:: column_to_rownames(var = "LOC")
dataHM_s <- dataHM_s %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_ZS_fe <- dataHM_ZS_fe %>% 
  tibble:: column_to_rownames(var = "LOC")
dataHM_fe <- dataHM_fe %>% 
  tibble:: column_to_rownames(var = "LOC")



dataHM_ZS_cann <- dataHM_ZS_cann %>% 
  tibble:: column_to_rownames(var = "LOC")
dataHM_cann <- dataHM_cann %>% 
  tibble:: column_to_rownames(var = "LOC")


dataHM_ZS_feat <- dataHM_ZS_feat %>% 
  tibble:: column_to_rownames(var = "LOC")

dataHMtau <- dataHMtau %>% 
  tibble:: column_to_rownames(var = "LOC")

dataHMtau2 <- dataHMtau2 %>% 
  tibble:: column_to_rownames(var = "LOC")

#Selecting columns, changing to data matrix as required by package

df_pi <- dataHM_ZS_pi[,-c(1:6)]
df_pi <- data.matrix(df_pi)

df_pi_log2tpm <- dataHM_pi[,-c(1:6)]
df_pi_log2tpm <- data.matrix(df_pi_log2tpm)



df_fe <- dataHM_ZS_fe[,-c(1:5)]
df_fe <- data.matrix(df_fe)
df_fe_log2tpm <- dataHM_fe[,-c(1:5)]
df_fe_log2tpm <- data.matrix(df_fe_log2tpm)


df_n <- dataHM_ZS_n[,-c(1:5)]
df_n <- data.matrix(df_n)
df_n_log2tpm <- dataHM_n[,-c(1:5)]
df_n_log2tpm <- data.matrix(df_n_log2tpm)



df_s <- dataHM_ZS_s[,-c(1:5)]
df_s <- data.matrix(df_s)
df_s_log2tpm <- dataHM_s[,-c(1:5)]
df_s_log2tpm <- data.matrix(df_s_log2tpm)

df_cann <- dataHM_ZS_cann[,-c(1:6)]
df_cann <- data.matrix(df_cann)
df_cann_log2tpm <- dataHM_cann[,-c(1:6)]
df_cann_log2tpm <- data.matrix(df_cann_log2tpm)

df_base<- dataHM_base[,-c(1,2)]
df_base <- data.matrix(df_base)

df_ZS<- dataHM_ZS[,-c(1,2)]
df_ZS <- data.matrix(df_ZS)

head(senescinggenesTPM)


df_senescing<- senescinggenesZscore[,-c(1,2)]
df_senescing <- data.matrix(df_senescing)

view(df_senescing_log2tpm)
df_senescing_log2tpm<- senescinggenesTPM[,-c(1,2)]
df_senescing_log2tpm <- data.matrix(df_senescing_log2tpm)


df_senescingjustOL0P_log2tpm<- df_senescing_log2tpm[,5]
df_senescingjustOL0P_log2tpm <- data.matrix(df_senescingjustOL0P_log2tpm,header=TRUE)

df_tau_pi <- dataHMtau[,-c(1:6)]
df_tau_pi <- data.matrix(df_tau_pi)

df_tau_nitrate <- dataHMtau2[,-c(1:5)]
df_tau_nitrate <- data.matrix(df_tau_nitrate)

#
#setting colors
colfun = circlize::colorRamp2(c(0,2.5,5,7.5,10,12.5,15), #gives range of the scale
                              c("#4b3991","#408ea4","#fcffad","#fba453", "#c5242a","#8a0033","#1f000c")) #for log2tpm+1

colfun2 = circlize::colorRamp2(c(-2,0,2), #gives range of the scale
                               c("#408ea4","lightyellow1","#c5242a")) #for zscoring

colfun3 = circlize::colorRamp2(c(0,0.5,1), #gives range of the scale
                               c("#408ea4","lightyellow1","#c5242a")) #for zscoring


colfun4 = circlize::colorRamp2(c(0,2.5,5,7.5,10,12.5,15), #gives range of the scale
                              c("#250f51","#4b3991","#408ea4","#fcffad","#fba453", "#c5242a","#8a0033")) #for log2tpm+1


pal <- colorRampPalette(c("skyblue", "yellow", "orange", "darkred"))(4)
colors = structure(pal, names = c("0.00", "0.20","0.85", "1.00"))
colors = structure(pal, names = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

colors<-colorRamp2(c(0,0.24,0.25,0.83,0.84,0.85), c("darkblue","skyblue","white","#c5242a","#fcffad", "darkred"))
colors<-colorRamp2(c(0,0.84,0.85), c("white","#8a0033", "black"))
pdf("Pi_HM_v4_alt.pdf",height=14)
pdf("N_HM_v3_alt.pdf",height=8)

levels(factor)
logtpmHM=Heatmap(df_pi_log2tpm,
                 col =colfun4,
                 #use_raster=FALSE,
                 cluster_columns=FALSE,
                 cluster_rows=TRUE,
                 show_row_names = TRUE,
                 cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT ON IF YOU WANT TO REODER PER YOUR LIST.
                 # cluster_rows=TRUE,
                 #clustering_method_columns="ward.D",
                 #clustering_method_rows = "ward.D",
                 row_dend_reorder = TRUE, #turned off because cant see difference due to magnitude of genes (5000+)
                 # row_km = 12,
                 #column_km=9,
                 #column_km_repeats=100,
                 # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus 
                 row_gap = unit(2, "mm"),
                 row_split=factor(row_order$Grouping), #row_split=paste(row_order$Grouping)
                 row_title_gp = gpar(font = 2,fontsize=6),
                 row_names_gp = gpar(fontsize = 6),
                 column_split = factor(col_order$Grouping),
                 cluster_column_slices = FALSE,
                 column_gap=unit(2,"mm"),
                 column_title="N Regulatory genes", 
                 #column_title_side = "bottom",
                 column_title_gp=gpar(fontsize=15,fontface="bold")
)

draw(logtpmHM)

dev.off()

zscoreHM=Heatmap(df_pi,
                 col =colfun2,
                 #use_raster=FALSE,
                 cluster_columns=FALSE,
                 #cluster_rows=TRUE,
                 show_row_names = TRUE,
                 cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
                 cluster_rows=TRUE,
                 #clustering_method_columns="ward.D",
                 #clustering_method_rows = "ward.D",
                 row_dend_reorder = TRUE,
                 # row_km = 12,
                 #column_km=9,
                 #column_km_repeats=100,
                 # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
                 row_gap = unit(2, "mm"),
                 row_split=factor(row_order$Grouping), #row_split=paste(row_order$Grouping)
                 row_title_gp = gpar(font = 2,fontsize=6),
                 row_names_gp = gpar(fontsize = 6),
                 column_split = factor(col_order$Grouping),
                 cluster_column_slices = FALSE,
                 column_gap=unit(2,"mm"),
                 column_title="Pi Regulatory genes", 
                 #column_title_side = "bottom",
                 column_title_gp=gpar(fontsize=15,fontface="bold")
)



draw(zscoreHM)
#this bit is for tau scoring
pdf("tau_scoring.pdf",height=14)
tauscorepi=Heatmap(df_tau_pi,
                 col =colors,
                 #use_raster=FALSE,
                 cluster_columns=FALSE,
                 #cluster_rows=TRUE,
                 show_row_names = TRUE,
                 cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
                 cluster_rows=TRUE,
                 #clustering_method_columns="ward.D",
                 #clustering_method_rows = "ward.D",
                 row_dend_reorder = TRUE,
                 # row_km = 12,
                 #column_km=9,
                 #column_km_repeats=100,
                 # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
                 row_gap = unit(2, "mm"),
                 row_split=factor(row_order$Grouping), #row_split=paste(row_order$Grouping)
                 row_title_gp = gpar(font = 2,fontsize=6),
                 row_names_gp = gpar(fontsize = 6),
                 column_split = factor(col_ordertau$Grouping),
                 cluster_column_slices = FALSE,
                 column_gap=unit(1,"mm"),
                 column_title="Phosphate Gene Specificity", 
                 #column_title_side = "bottom",
                 column_title_gp=gpar(fontsize=15,fontface="bold"),
                 width=unit(2,"cm")
)
draw(tauscorepi)
logtpmHM+tauscorepi

dev.off()
tauscorenitrate=Heatmap(df_tau_nitrate,
                   col =colors,
                   #use_raster=FALSE,
                   cluster_columns=FALSE,
                   #cluster_rows=TRUE,
                   show_row_names = TRUE,
                   cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
                   cluster_rows=TRUE,
                   #clustering_method_columns="ward.D",
                   #clustering_method_rows = "ward.D",
                   row_dend_reorder = TRUE,
                   # row_km = 12,
                   #column_km=9,
                   #column_km_repeats=100,
                   # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
                   row_gap = unit(2, "mm"),
                   row_split=factor(row_order2$Grouping), #row_split=paste(row_order$Grouping)
                   row_title_gp = gpar(font = 2,fontsize=6),
                   row_names_gp = gpar(fontsize = 6),
                   column_split = factor(col_ordertau$Grouping),
                   cluster_column_slices = FALSE,
                   column_gap=unit(1,"mm"),
                   column_title="Nitrate Gene Specificity", 
                   #column_title_side = "bottom",
                   column_title_gp=gpar(fontsize=15,fontface="bold"),
                   width=unit(2,"cm")
)

logtpmHM+tauscorenitrate

draw(tauscorenitrate)
dev.off()


#heatmaps for ALL DEGs.
pdf("AllDEG_kmeans_HM_forjim.pdf",height=7)

hfkmeans=Heatmap(df_base,
                 col =colfun,
                 use_raster=FALSE,
                 cluster_columns=FALSE,
                 show_row_names = FALSE,
                 cluster_rows=TRUE,
                 #clustering_method_columns="ward.D",
                 clustering_method_rows = "ward.D",
                 row_dend_reorder = TRUE, #turned off because cant see difference due to magnitude of genes (5000+)
                 row_km = 8,
                 #column_km=9,
                 #column_km_repeats=100,
                 row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
                 row_gap = unit(2, "mm"),
                 row_names_gp = gpar(fontsize = 6),
                 column_split = factor(col_order$Grouping),
                 #cluster_column_slices = FALSE
                 column_gap=unit(2,"mm"),
                 column_title="DEGS -OL0P kmeans",
                 column_title_gp=gpar(fontsize=15,fontface="bold")
)

draw(hfkmeans)

hf=Heatmap(df_ZS,
           col =colfun2,
           use_raster=FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           cluster_rows=TRUE,
           #clustering_method_columns="ward.D",
           clustering_method_rows = "ward.D",
           row_dend_reorder = TRUE, #turned off because cant see difference due to magnitude of genes (5000+)
           row_km = 8,
           #column_km=9,
           #column_km_repeats=100,
           #row_km_repeats = 1000,#divide into clusters by kmeans, repeat 100 use consensus
           row_gap = unit(2, "mm"),
           row_names_gp = gpar(fontsize = 6),
           column_split = factor(col_order$Grouping),
           #cluster_column_slices = FALSE
           column_gap=unit(2,"mm"),
           column_title="DEGS -OL0P kmeans",
           column_title_gp=gpar(fontsize=15,fontface="bold")
)


library("magick")
draw(hf)
dev.off()

###DO NOT DELETE ANYTHING AFTER THIS. USED FOR GO TERM ANALYSIS.###
pdf("AllDEG_kmeans_HM_test3.pdf",height=7)
set.seed(124)
class(df_ZS)
hf=Heatmap(df_ZS,
           col =colfun2,
           use_raster=TRUE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           cluster_rows=TRUE,
           cluster_row_slices = TRUE,
           #clustering_method_columns="ward.D",
           clustering_method_rows = "ward.D",
           show_row_dend = FALSE,
           row_dend_reorder = TRUE,
           row_km = 9,
           row_km_repeats = 1000,
           #column_km=9,
           #column_km_repeats=100,
           row_gap = unit(2, "mm"),
           row_names_gp = gpar(fontsize = 6),
           column_split = factor(col_order$Grouping),
           #cluster_column_slices = FALSE
           column_gap=unit(2,"mm"),
           column_title="DEGS -OL0P kmeans",
           column_title_gp=gpar(fontsize=15,fontface="bold")
)
HM=draw(hf)
#testing to see if I get the same clusters from a tpm heatmap
pdf("AllDEG_kmeans_HM_forjim2.pdf",height=7)
set.seed(124)
hb=Heatmap(df_base,
           col =colfun,
           use_raster=FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           cluster_rows=TRUE,
           cluster_row_slices = TRUE,
           #clustering_method_columns="ward.D",
           clustering_method_rows = "ward.D",
           row_dend_reorder = TRUE,
           row_km = 9,
           row_km_repeats = 1000,
           #column_km=9,
           #column_km_repeats=100,
           row_gap = unit(2, "mm"),
           row_names_gp = gpar(fontsize = 6),
           column_split = factor(col_order$Grouping),
           #cluster_column_slices = FALSE
           column_gap=unit(2,"mm"),
           column_title="DEGS -OL0P kmeans",
           column_title_gp=gpar(fontsize=15,fontface="bold")
)
HM2=draw(hb)

dev.off()




r.dend <- row_dend(HM)  #If needed, extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x)) #check/confirm size gene clusters

#testing for the tpm heatmaps
r.dend2 <- row_dend(HM2)  #If needed, extract row dendrogram
rcl.list2 <- row_order(HM2)  #Extract clusters (output is a list)
lapply(rcl.list2, function(x) length(x)) #check/confirm size gene clusters



library(magrittr) # needed to load the pipe function '%%'


clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(LOC = rownames(df_ZS[rcl.list[[i]],]), #rownames(tf.log)[clusterlist[[i]]] if cluster somehow only has one gene. https://www.biostars.org/p/465304/
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

colnames(ref_n)


clu_df
clu_df2 <- left_join(clu_df,ref_all,by="LOC")
count <- sum(!is.na(clu_df2$gene.family))
print(count)


clu_df2 <- left_join(ref_pi,clu_df, by="LOC") #%>% left_join(y=ref_pi[,1:4],by=c("LOC"))  %>% left_join(ref_n,by=c("LOC","annotation","MapMan.annot","gene.family")) 


write.csv(clu_df2,file="clusterextracts.csv",row.names=FALSE)



clu_list <- split(clu_df, clu_df$Cluster) #close

clu_list <- lapply(unique(clu_df$Cluster), function(cluster) {
  locations <- clu_df$LOC[clu_df$Cluster == cluster]
  names(locations) <- names(clu_df$Cluster)
  return(locations)                                         
}) #better. the cluster next inside function() could have been named anything. So this code is essentially saying, apply this function to each unique value from the 'Cluster' column of the data frame_df, put each LOC that matches to that unique value from the "Cluster" value into a elemenT, which is ordered/identified by the index
names(clu_list) <- unique(clu_df$Cluster)

names(clu_list) <- names(rcl.list)


#ok, now ready for GO terms.
library(clusterProfiler)

#build the c#build the c#build the comparison
all_GO_comp <- compareCluster(geneCluster = clu_list,
                              
                              fun = "enrichGO",
                              
                              OrgDb = "org.Csativa.eg.db",
                              
                              ont="BP",
                              
                              pAdjustMethod = 'fdr',
                              
                              keyType = 'GID', # this was used when creating org.Csativa.eg.db
                              
                              qvalueCutoff  = 0.01,
                              
                              minGSSize = 10 # min gene number per GO term, will determine which modules might fall out, soemthing to play with
                              
)


# problem is that the qvalues dont scale nicely in the dotplot, so first

# use cutoff, e.g. logcut =10 and generate the -logs


logcut = 20

  

logs = as.data.frame(all_GO_comp@compareClusterResult$qvalue)

colnames(logs) = 'qvalue'



logs = logs %>%
  
  mutate(nlogq = case_when(
    
    -log(qvalue) > logcut ~ logcut, # everything above cutoff
    
    -log(qvalue) <= logcut ~ -log(qvalue)))

#now add those logs back to object to use.
all_GO_comp@compareClusterResult$nlogq = logs$nlogq

options(enrichplot.colours = c("#408ea4","#c5242a"))

dotplot(all_GO_comp,
        
        showCategory=20, # top20 for every module
        
        color='nlogq',
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        font.size =  10,
        
        x = 'Cluster',
        title='All Go Modules') +
  
  scale_color_viridis_c(option = "viridis",direction=-1)

dev.off()
ggsave('all_GO_modules_top20.pdf', width = 25, height = 25, unit = 'cm')


all_GO_summary = as.data.frame(all_GO_comp)

write.csv(all_GO_summary, 'all_GO_summary.csv')

view(clu_df2)
dev.off()

# first need to prepare a GoSemSim object first:
GOdata = godata(
  
  OrgDb = "org.Csativa.eg.db",
  
  keytype = "GID",
  
  ont = 'BP',
  
  computeIC = TRUE,
  
  processTCSS = TRUE,
  
  cutoff = NULL
  
)

# now use the above to simplify

# need to specify clusterProfiler::simplify because of some other function with the same name :(

# without semData = GOdata it only works using measure = 'Wang' and semdata = NULL


all_GO_comp_simple = clusterProfiler::simplify(

  all_GO_comp,

  cutoff = 0.7, # lower number = fewer GOs

  by = "qvalue",

  select_fun = min,

  measure = "Rel",

  semData = GOdata)



dotplot(all_GO_comp_simple,
        
        showCategory=10, # top10 for every module
        
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        color="nlogq",
        
        font.size = 14,
        
        x = 'Cluster',
        title= 'Simplified') +
  
  #scale_color_viridis_c(option = "cividis")+
  
  #scale_color_viridis_c(option = "inferno")
  
  scale_color_viridis_c(option = "viridis")


ggsave('Simplified_GO.pdf', width = 35, height = 35, unit = 'cm')

all_GO_comp_simple_drop = dropGO(all_GO_comp_simple, level = c(1,2,3))


dotplot(all_GO_comp_simple_drop,
        
        showCategory=20, # top20 for every module
        
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        color='nlogq',
        
        font.size = 5,
        
        x = 'Cluster',
        title= 'Simplified and dropping high level GO terms') +
  
  scale_color_viridis_c(option = "viridis")

ggsave('Simpliedanddroppinghighlevels_GO.pdf', width = 25, height = 35, unit = 'cm')

write.csv(all_GO_comp_simple_drop, 'all_GO_comp_simple_drop.csv')











#senescing plots
library(org.Csativa.eg.db) 
senescing_genes <- read.csv("OL0PGENES.csv")
str(senescing_genes)
DEGSTPM <- read.csv("DEGs LFC.csv")
str(DEGSTPM)
colnames(DEGSTPM)
senescing_genes_df <- left_join(senescing_genes,DEGSTPM[,c(1,2,3,4,5,6,16)],by="LOC")
colnames(senescing_genes_df)
senescing_genes_df$Regulation <- ifelse(senescing_genes_df$logFC_OLP0 >0,"Upregulated", "Downregulated") #conditional, assign the row as Upregulated if >0 and Downregulated if <0 in a new column called Regulation
senescing_genes_df$Regulation <- ifelse(senescing_genes_df$logFC_OLP0 > 0, "Upregulated", "Downregulated")
senescing_genes_df_GO <- senescing_genes_df[,c(1,14)]
colnames(senescing_genes_df_GO)
senescing_genes_df_GO$Regulation <- factor(senescing_genes_df_GO$Regulation, levels=c("Upregulated","Downregulated"))
levels(senescing_genes_df_GO$Regulation)
#now pivot wider into list for GO term plots
wide_senescing_genes_df_GO <- pivot_wider(data=senescing_genes_df_GO,names_from=Regulation,values_from = LOC,values_fn = list(LOC = list))
gene_list <- as.list(wide_senescing_genes_df_GO)
gene_list <- lapply(gene_list, function(x) unlist(x, use.names = FALSE))
gene_list <- gene_list[c("Upregulated", "Downregulated")]

gc()
dev.off()
pdf("Sesnescing_genes.pdf",height=14)
set.seed(124)
class(df_senescing)
hf=Heatmap(df_senescing,
           col =colfun2,
           use_raster=FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           cluster_rows=TRUE,
           cluster_row_slices = FALSE,
           #clustering_method_columns="ward.D",
           clustering_method_rows = "ward.D",
           row_dend_reorder = TRUE,
           row_km = 5,
           row_km_repeats = 100,
           #column_km=9,
           #column_km_repeats=100,
           row_gap = unit(2, "mm"),
           row_names_gp = gpar(fontsize = 6),
           #column_split = factor(col_order$Grouping),
           #cluster_column_slices = FALSE
           column_gap=unit(2,"mm"),
           column_title="OL:0P Genes",
           column_title_gp=gpar(fontsize=15,fontface="bold")
)
HM=draw(hf)
dev.off()
#testing to see if I get the same clusters from a tpm heatmap
pdf("Sesnescing_genes2.pdf",height=7)
set.seed(124)
hb=Heatmap(df_senescingjustOL0P_log2tpm,
           col =colfun,
           use_raster=FALSE,
           cluster_columns=FALSE,
           show_row_names = FALSE,
           #cluster_rows=TRUE,
           cluster_row_slices = FALSE,
           #clustering_method_columns="ward.D",
           #clustering_method_rows = "ward.D",
           #row_dend_reorder = TRUE,
           show_row_dend = FALSE,
           row_km = 5,
           row_km_repeats = 100,
           #column_km=9,
           #column_km_repeats=100,
           row_gap = unit(2, "mm"),
           row_names_gp = gpar(fontsize = 6),
           #column_split = factor(col_order$Grouping),
           #cluster_column_slices = FALSE
           column_gap=unit(2,"mm"),
           column_title="OL:0P Genes(TPM)",
           column_title_gp=gpar(fontsize=15,fontface="bold")
)
HM=draw(hb)
dev.off()

r.dend <- row_dend(HM)  #If needed, extract row dendrogram
rcl.list <- row_order(HM)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x)) #check/confirm size gene clusters


library(magrittr) # needed to load the pipe function '%%'


clu_df <- lapply(names(rcl.list), function(i){
  out <- data.frame(LOC = rownames(df_senescing[rcl.list[[i]],]),
                    Cluster = paste0("cluster", i),
                    stringsAsFactors = FALSE)
  return(out)
}) %>%  #pipe (forward) the output 'out' to the function rbind to create 'clu_df'
  do.call(rbind, .)

clu_df
clu_df2 <- left_join(clu_df,ref_all,by="LOC")
count <- sum(!is.na(clu_df2$gene.family))
print(count)
write.csv(clu_df2,file="OL0PGENES.csv",row.names=FALSE)

clu_list <- lapply(unique(clu_df$Cluster), function(cluster) {
  locations <- clu_df$LOC[clu_df$Cluster == cluster]
  names(locations) <- names(clu_df$Cluster)
  return(locations)
}) #better. the cluster next inside function() could have been named anything. So this code is essentially saying, apply this function to each unique value from the 'Cluster' column of the data frame_df, put each LOC that matches to that unique value from the "Cluster" value into a elemenT, which is ordered/identified by the index
names(clu_list) <- unique(clu_df$Cluster) #includes clusterprefix and number

names(clu_list) <- names(rcl.list) #includes only cluster number

#no clustering
clu_list5 <- df3_with_na
clu_list5 <- clu_list5 %>%
  mutate(pseudocluster = "pseudocluster")

colnames(clu_list5)
clu_list5 <- clu_list5[,c(24,1)]
colnames(clu_list5)
clu_list5 <- clu_list5[,2]
char_list <- list(pseudocluster=clu_list5)

all_GO_comp <- compareCluster(geneCluster = gene_list,
                              
                              fun = "enrichGO",
                              
                              OrgDb = "org.Csativa.eg.db",
                              
                              ont="BP",
                              
                              pAdjustMethod = 'fdr',
                              
                              keyType = 'GID', # this was used when creating org.Csativa.eg.db
                              
                              qvalueCutoff  = 0.01,
                              
                              minGSSize = 10 # min gene number per GO term, will determine which modules might fall out, soemthing to play with
                              
)

logcut = 20



logs = as.data.frame(all_GO_comp@compareClusterResult$qvalue)

colnames(logs) = 'qvalue'



logs = logs %>%
  
  mutate(nlogq = case_when(
    
    -log(qvalue) > logcut ~ logcut, # everything above cutoff
    
    -log(qvalue) <= logcut ~ -log(qvalue)))

#now add those logs back to object to use.
all_GO_comp@compareClusterResult$nlogq = logs$nlogq

options(enrichplot.colours = c("#408ea4","#c5242a"))

dotplot(all_GO_comp,
        
        showCategory=40, # top20 for every module
        
        color='nlogq',
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        font.size = 7,
        
        x = 'Cluster',
        title='All Go Modules') +
  
  scale_color_viridis_c(option = "viridis",direction=-1)
ggsave('OL0P_allgenes_top20.pdf', width = 25, height = 25, unit = 'cm')

# first need to prepare a GoSemSim object first:
GOdata = godata(
  
  OrgDb = org.Csativa.eg.db,
  
  keytype = "GID",
  
  ont = 'BP',
  
  computeIC = TRUE,
  
  processTCSS = TRUE,
  
  cutoff = NULL
  
)
all_GO_comp_simple = clusterProfiler::simplify(
  
  all_GO_comp,
  
  cutoff = 1, # lower number = fewer GOs
  
  by = "qvalue",
  
  select_fun = min,
  
  measure = "Rel",
  
  semData = GOdata)



dotplot(all_GO_comp_simple,
        
        showCategory=15, # top10 for every module
        
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        color="nlogq",
        
        font.size = 7,
        
        x = 'Cluster',
        title= 'Simplified: OL0P') +
  
  #scale_color_viridis_c(option = "cividis")+
  
  #scale_color_viridis_c(option = "inferno")
  
  scale_color_viridis_c(option = "viridis")

ggsave('OL0P_downvsup_nodrop.pdf', width = 35, height = 30, unit = 'cm')


all_GO_comp_simple_drop = dropGO(all_GO_comp_simple, level = c(1,2,3))


dotplot(all_GO_comp_simple_drop,
        
        showCategory=20, # top20 for every module
        
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        color='nlogq',
        
        font.size = 8,
        
        x = 'Cluster',
        title= 'Simplified and dropping high level GO terms:OL0P') +
  
  scale_color_viridis_c(option = "viridis")
ggsave('OL0P_simplified_drop_highlevel.pdf', width = 25, height = 25, unit = 'cm')

gene_GO_mapping <- all_GO_comp_simple@compareClusterResult

# View the gene-to-GO mapping
View(all_GO_comp_simple_drop)

gene_GO_mapping <- all_GO_comp_simple_drop@compareClusterResult

df_expanded <- gene_GO_mapping %>%
  separate_rows(geneID, sep = "/") %>%
  arrange(geneID)
df_expanded <-dplyr:: rename(df_expanded,"LOC"="geneID")
colnames(senescing_genes_df)
senescing_genes_df_GO$
library("dplyr")
df_expanded_annot <- dplyr:: left_join(df_expanded,senescing_genes_df_GO,by="LOC", relationship = "many-to-many")

df_expanded_annot <- left_join(df_expanded_annot,senescing_genes_df[,c(1,9)],by="LOC")

write.csv(df_expanded_annot,file="OL0P_expandedannot.csv",row.names = FALSE)


#%>% left_join(senescing_genes_df[,c(1,3)],by="LOC")

View(df_expanded_annot)
senescing_genes_df_GO %>%
  count(LOC) %>%
  filter(n > 1)
write.csv(df_expanded_annot,file="df_expanded_annot.csv",row.names = FALSE)

cnetplot(all_GO_comp_simple,
         showCategory = 10,  # Show top 10 categories
         circular = FALSE,   # Set to TRUE for a circular layout
         colorEdge = TRUE,
         max.overlaps=5)   # Edge color represents gene-category relationships
data(geneList)

#tauscoring####
#START HERE FOR HEATMAP PLOTS IF USING THE TPM10 FILTER.
dataHM_base <- read.csv("data_heatmap_base_TPM10.csv") 
dataHM_ZS <- read.csv("data_heatmap_zscored.csv")
colnames(dataHM_base)
tau0P <- dataHM_base[,c(1,2,3,4,8,12,16,20)]
test <- tau0P
tau025P <- dataHM_base[,c(1,2,3,5,9,13,17,21)]
tau1P <- dataHM_base[,c(1,2,3,6,10,14,18,22)]
tau2P <- dataHM_base[,c(1,2,3,7,11,15,19,23)]
tau0P[, 4:ncol(tau0P)] <- apply(tau0P[, 4:ncol(tau0P)], 2, ts)
df1 <- DEG
write.csv(tau0P,file="tautest.csv",row.names=FALSE)
tau_score <- function(x) {
  N <- length(x)
  tau <- sum(1 - x) / (N - 1)
  return(tau)
}

ts <- function(x) {
  t<-sum(1-x/max(x))/(length(x)-1)
  return(tau)
}
t<-sum(1-x/max(x))/(length(x)-1)



tau_scores <- apply(tau0P,1, tau_score)
test[, 4:ncol(test)] <- apply(test[, 4:ncol(test)], 1,calculate_t )

 df1[, 4:ncol(df1)] <- apply(df1[, 4:ncol(df1)], 2, function(x) log2(x))  #log transformations
df1[, 4:ncol(df1)] <- apply(df1[, 4:ncol(df1)], 2, function(x) log2(x + 1))  #log transformations
                 
devtools::install_github("AllenInstitute/scrattch.hicat")
colnames(test)
library(scrattch.hicat)
print(test)
colnames(test)
tester <- test %>% 
  tibble:: column_to_rownames(var = "LOC")
tester <- as.matrix(tester)
tautest <- calc_tau(tester, byRow = TRUE)
colnames(tautest)
tautest2 <- as.data.frame(tautest)
tautest2$LOC <- rownames(tautest2)
colnames(tautest2)

# Example dataframe
# Example data
data <- matrix(c(10, 5, 8, 12, 6, 7,
                 3, 9, 11, 4, 2, 10
                 8, 6, 4, 2, 3, 7,
                 4, 9, 6, 5, 8, 1), 
               ncol = 6, byrow = TRUE)

colnames(data) <- paste0("V", 1:6)  # Assign column names V1, V2, ..., V6
rownames(data) <- paste0("Row", 1:nrow(data))  # Assign row names Row1, Row2, ..., RowN

# Convert matrix to data frame for clarity
df <- as.data.frame(data)

calculate_t_row <- function(row) {
  t <- sum(1 - row / max(row)) / (length(row) - 1)
  t <- round(t,2)
  return(t)
}

tau0P
tau025P
tau1P
tau2P
# Apply calculate_t_row function across rows
test$tau0P <- apply(tau0P[,4:ncol(tau0P)], 1, calculate_t_row)
test$tau025P <- apply(tau025P[,4:ncol(tau025P)], 1, calculate_t_row)
test$tau1P <- apply(tau1P[,4:ncol(tau1P)], 1, calculate_t_row)
test$tau2P <- apply(tau2P[,4:ncol(tau2P)], 1, calculate_t_row)
colnames(test)

test <- test[,c(1,2,3,9,10,11,12)]
write.csv(test,file="tau.csv",row.names=FALSE)
tausumm <- read.csv("tau.csv")
toplot <- tausumm[,1:4]
colnames(tausumm)
view(toplot)


melted_tau= melt(tausumm,id.vars=c("LOC","Feature","Annotation"))

breaks=c(0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.80,0.85,0.90,0.95,1.0)

ggplot(melted_tau, aes(x=value))+
  geom_histogram(bins=20,color="black",fill="#69b3a2")+
  facet_wrap(~variable,ncol=2,scales="free")+
  xlim(0,NA)+
  ylim(0,500)+
  theme_bw()+
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        axis.text.x = element_text(size = 14, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 14, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 14, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=16),
        legend.title=element_text(size=18,vjust=0.5)
  )+
  ggtitle("Bin size=0.05 | Tau score based on Same treatment but across organs")

geom_bar() +
  scale_x_binned(n.breaks = 10)
ggplot(tausumm,aes(x=tau0P))+
  geom_histogram(bins=20,color="black",fill="#69b3a2")+
  ggtitle("Bin size = 0.05 | 0P") +
  theme_bw() +
  theme(
    plot.title = element_text(size=15)
  )

ggplot(tausumm,aes(x=tau0P))+
  geom_bar() +
  scale_x_binned(n.breaks = 10)
  ggtitle("Bin size = 0.05 | 0P") +
  theme_bw()
  +geom_text(stat='count', aes(label=..count..), vjust=-1)
  
ggplot(melted_tau,aes(x=value)) +
  geom_histogram( breaks = breaks, color = "black",fill="#69b3a2") + 
  scale_x_continuous(breaks = c(0,0.25,0.50,0.75,1.0))+
  geom_text(stat='bin', aes(label=after_stat(count)),breaks=breaks, vjust=-1.5)+
  ylim(0,500)+
  facet_wrap(~variable,ncol=2,scales="free")+
  theme_bw()+
  theme(axis.title.x = element_text(color="black", size=12, face="bold"),
        axis.title.y = element_text(color="black", size=12, face="bold"),
        axis.text.x = element_text(size = 14, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 14, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 14, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.text = element_text(size=16),
        legend.title=element_text(size=18,vjust=0.5)
  )+
  ggtitle("Bin size=0.05 | Tau score based on Same treatment but across organs")
`dev.off()
ggsave(file="taufacets.pdf",units="cm",height=25, width=35)


##Heatmaps for trait correlations
#Starting trait heatmap to show patterns in one mass file.

#first load in relevant DF yes?
#checking row.names match
WGCNA5 <- read.csv("Traits_WGCNA2.csv")
WGCNA6 <- WGCNA5[-9,]
View(WGCNA6$Combined)
rownames(WGCNAfiltered)
rownames(WGCNA6) <- NULL 

WGCNA6 <- WGCNA6 %>% 
  tibble:: column_to_rownames(var = "Combined")

collectGarbage()

table(rownames(WGCNAfiltered) == rownames(WGCNA6))


View(WGCNA5)

# Apply z-scoring to each column of the dataframe

# Add rownames as the first column of the data

WGCNA7 <- as.data.frame(apply(WGCNA5[,2:ncol(WGCNA5)], 2, function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)))

rownames(WGCNA7) <- WGCNA5[,1]

#
# Apply log transformation (log1p handles log(0) by returning -Inf)
WGCNA8 <- as.data.frame(apply(WGCNA5[, 2:ncol(WGCNA5)], 2, function(x) log2(x+1)))

# Set rownames using the first column of WGCNA5
rownames(WGCNA8) <- WGCNA5[, 1]


#make heatmap

WGCNA7mat <- data.matrix(WGCNA7)
WGCNA8mat <- data.matrix(WGCNA8)


colfun2 = circlize::colorRamp2(c(-1.5,0,1.5), #gives range of the scale
                               c("#408ea4","white","#c5242a")) #for zscoring


range_values <- range(WGCNA8mat, na.rm = TRUE)
print(range_values)
colfun3 <- circlize::colorRamp2(c(range_values[1], range_values[2]), 
                                c("#408ea4", "#c5242a"))

print(rownames(WGCNA7mat))
col_order=data.frame(ID=c("LR:P0","LR:P0.25","LR:P1","LR:P2",
                          "OLP:0","OLP:0.25","OLP:1","OLP:2",
                          "MLP:0","MLP:0.25","MLP:1","MLP:2",
                          "YSP:0","YSP:0.25","YSP:1","YSP:2",
                          "TFP:0","TFP:0.25","TFP:1","TFP:2"),
                     Grouping=c(rep("LR",4),rep("OL",4),rep("ML",4),rep("YS",4),rep("TF",4))
)
col_order$Grouping=factor(col_order$Grouping,levels=c("LR","OL","ML","YS","TF"))
levels(col_order$Grouping)



roworder= data.frame(ID=c("TPM.LR_P0",   "TPM.LR_P025", "TPM.LR_P1",   "TPM.LR_P2", "TPM.OL_P0","TPM.OL_P025", "TPM.OL_P1",   "TPM.OL_P2",   "TPM.ML_P0",   "TPM.ML_P025", "TPM.ML_P1",   "TPM.ML_P2",     "TPM.YS_P0", "TPM.YS_P025", "TPM.YS_P1",  "TPM.YS_P2","TPM.TF_P0","TPM.TF_P025", "TPM.TF_P1","TPM.TF_P2"),
                     Grouping=c(rep("LR",4),rep("OL",4),rep("ML",4),rep("YS",4),rep("TF",4))
)



pdf("TRAITHM2.PDF",height=8)
b=Heatmap(t(WGCNA7mat),
          col =colfun2,
          #use_raster=FALSE,
          cluster_columns=FALSE,
          # cluster_rows=TRUE,
          show_row_names = TRUE,
          #cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
          cluster_rows=TRUE,
          #clustering_method_columns="ward.D",
          #clustering_method_rows = "ward.D",
          row_dend_reorder = FALSE,
          # row_km = 12,
          #column_km=9,
          #column_km_repeats=100,
          # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
          row_gap = unit(2, "mm"),
          #row_split=factor(row_order$Grouping), #row_split=paste(row_order$Grouping)
          row_title_gp = gpar(font = 2,fontsize=6),
          row_names_gp = gpar(fontsize = 6),
          column_split = factor(roworder$Grouping),
          cluster_column_slices = FALSE,
          column_gap=unit(2,"mm"),
          column_title="Pi Regulatory genes", 
          #column_title_side = "bottom",
          column_title_gp=gpar(fontsize=15,fontface="bold")
)
draw(b)
dev.off()
ide
View(datExpr_selected)
view(datExpr_selectedZS)


#Plotting HM for RM
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/")
redmodule <- read.csv("Red Module/Redmodulegenes3.csv")
redmodulegrouping <- redmodule[,c(24,10,13),drop=FALSE]
redmodulegrouping <- redmodule[,c(36,16,19),drop=FALSE] #if useing redmodulegenes2
colnames(redmodulegrouping)
redmodulegrouping <- redmodule[,c(36,16,19),drop=FALSE] #if useing redmodulegenes2
redmodulegrouping$Grouping <- factor(redmodulegrouping$Grouping)
levels(redmodulegrouping$Grouping)
redmodule <- dplyr::rename(redmodule, LOC = shared.name) 
redmodule <- redmodule[,"LOC",drop=FALSE]
setwd("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
dataHM_ZS <- read.csv("data_heatmap_zscored.csv_TPM10.csv")
dataFC <- read.csv("filtALL_no_0P_cond.csv")
colnames(dataFC)
dataFC <- dataFC[, !grepl("FDR", colnames(dataFC))]
colnames(dataFC)
dataFC <- dataFC[,c(1,2,3,15,16,17,12,13,14,10,11,12,7,8,9,4,5,6)]
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/")
write.csv(dataFC,"./Red Module/dataFC2.csv",row.names=FALSE)
#now go into the actual file and add columns for the 1 Pi treatment groups manuall
#load back in
dataFC <- read.csv("./Red Module/dataFC2.csv")
colnames(dataFC)=colnames(dataHM_ZS)
write.csv(dataFC,file="dataFC3.csv",row.names=FALSE)
dataFC <- read.csv("dataFC3.csv",header=TRUE)
getwd()
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/Red Module/")
#OK, now let's filter for genes we care about.
redmodulegenes <- left_join(redmodule,dataHM_ZS, by="LOC")
redmodulegenes <- redmodulegenes[,-c(2,3)] %>% 
  tibble:: column_to_rownames(var = "LOC") 

redmodulegenes <- data.matrix(redmodulegenes)

dataFCgenes <- left_join(redmodule,dataFC, by="LOC")
dataFCgenes <- dataFCgenes[,-c(2,3)] %>% 
  tibble:: column_to_rownames(var = "LOC") %>% round(1)

dataFCgenes <- data.matrix(dataFCgenes)


colfun2 = circlize::colorRamp2(c(-1.5,0,1.5), #gives range of the scale
                               c("#408ea4","white","#c5242a")) #for zscoring


#draw first HM
Heatmap(redmodulegenes,
                 col =colfun2,
                 #use_raster=FALSE,
                 cluster_columns=FALSE,
                 #cluster_rows=TRUE,
                 show_row_names = TRUE,
                 cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
                 cluster_rows=TRUE,
                 #clustering_method_columns="ward.D",
                 #clustering_method_rows = "ward.D",
                 row_dend_reorder = FALSE,
                 # row_km = 12,
                 #column_km=9,
                 #column_km_repeats=100,
                 # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
                 #row_gap = unit(2, "mm"),
                 row_split=factor(redmodulegrouping$Grouping), #row_split=paste(row_order$Grouping)
                 row_title_gp = gpar(font = 2,fontsize=5),
                 row_names_gp = gpar(fontsize = 5),
                 column_split = factor(col_order$Grouping),
                 cluster_column_slices = FALSE,
                 column_gap=unit(2,"mm"),
                 column_title="Pi Regulatory genes", 
                 #column_title_side = "bottom",
                 column_title_gp=gpar(fontsize=15,fontface="bold")
)


#ok, lets think about making this HM more informative.

redmodulegrouping <- redmodule[,c(24,10,13),drop=FALSE]
redmodulegrouping$Grouping <- factor(redmodulegrouping$Grouping)
levels(redmodulegrouping$Grouping)
#draw first HM

pdf("./HeatmapRedModule_geneID.pdf",height=12,width=6)

Heatmap(redmodulegenes,
        col =colfun2,
        row_labels = redmodulegrouping$gene.family.custom,  # Replace row names with gene names
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Ensure dataFCgenes has valid values
          if (!is.na(dataFCgenes[i, j])) {
            grid.text(as.character(dataFCgenes[i, j]), x, y, 
                      gp = gpar(fontsize = 5,fontface="bold", col = "black"))
          }
        },
        use_raster=TRUE,
        cluster_columns=FALSE,
        #cluster_rows=TRUE,
        show_row_names = TRUE,
        cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
        cluster_rows=TRUE,
        #clustering_method_columns="ward.D",
        #clustering_method_rows = "ward.D",
        row_dend_reorder = FALSE,
        # row_km = 12,
        #column_km=9,
        #column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        #row_gap = unit(2, "mm"),
        row_split=factor(redmodulegrouping$Grouping), #row_split=paste(row_order$Grouping)
        row_title_gp = gpar(font = 2,fontsize=5),
        row_names_gp = gpar(fontsize = 5,fontface="italic"),
        column_split = factor(col_order$Grouping),
        cluster_column_slices = FALSE,
        column_gap=unit(2,"mm"),
        column_title="Red Module", 
        #column_title_side = "bottom",
        column_title_gp=gpar(fontsize=15,fontface="bold"),
        column_names_gp=gpar(fontsize=8)
)

dev.off()
getwd()


#SL list
getwd()
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA")
ref_sl <- read.csv("sl_genes.csv")
colnames(ref_sl)
View(ref_sl)
mlist <- read.csv("mastersheetWGCNA.csv")
colnames(mlist)
ref_slmlist <- left_join(ref_sl,mlist[,c(1,3)],by="LOC")
View(ref_slmlist)

ref_sl2 <- dplyr::rename(ref_sl, Grouping = Role) 
ref_sl2$Grouping
ref_sl2$Grouping <- factor(ref_sl2$Grouping)
levels(ref_sl2$Grouping)
ref_sl2$Grouping <- factor(ref_sl2$Grouping,levels=c("SL precursor","SL synthesis","SL degradation","Sequestration","Perception","Corepressor","TF"))
levels(ref_sl2$Grouping)
ref_sl2 <- ref_sl2[-c(14,15,19,20,21),]

#now load in stuff for LFC
setwd("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
DEGtable_mod <- read.csv("DEG_TPM_all_no0POL.CSV",stringsAsFactors = TRUE,header=TRUE) #5341 DEG


ref_sl2genes <- left_join(ref_sl2,dataHM_ZS, by="LOC")

ref_sl2genes <- ref_sl2genes[,-c(2:6)]%>% 
  tibble:: column_to_rownames(var = "LOC") 
ref_sl2genes <- ref_sl2genes%>%  filter_at(vars(1:ncol(.)), any_vars(!is.na(.)))
ref_sl2genes <- data.matrix(ref_sl2genes)


dataFCgenessl2 <- left_join(ref_sl2,dataFC, by="LOC")
colnames(dataFCgenessl2)
dataFCgenessl2 <- dataFCgenessl2[,-c(2:6)] %>% 
  tibble:: column_to_rownames(var = "LOC") %>% round(1)
dataFCgenessl2 <- dataFCgenessl2%>%  filter_at(vars(1:ncol(.)), any_vars(!is.na(.)))
dataFCgenessl2 <- data.matrix(dataFCgenessl2)


getwd()
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA")
pdf("./unsigned_v6/Greenyellow Module/SL_hm3_LOC.pdf",height=6,width=8)

Heatmap(ref_sl2genes,
        col =colfun2,
        row_labels = ref_sl2$LOC,  # Replace row names with gene names
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Ensure dataFCgenes has valid values
          if (!is.na(dataFCgenessl2[i, j])) {
            grid.text(as.character(dataFCgenessl2[i, j]), x, y, 
                      gp = gpar(fontsize = 8,fontface="bold", col = "black"))
          }
        },
        use_raster=TRUE,
        cluster_columns=FALSE,
        #cluster_rows=TRUE,
        show_row_names = TRUE,
        cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
        #cluster_rows=TRUE,
        #clustering_method_columns="ward.D",
        #clustering_method_rows = "ward.D",
        row_dend_reorder = FALSE,
        # row_km = 12,
        #column_km=9,
        #column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        #row_gap = unit(2, "mm"),
        row_split=factor(ref_sl2$Grouping), #row_split=paste(row_order$Grouping)
        row_title_gp = gpar(font = 2,fontsize=8),
        row_names_gp = gpar(fontsize = 8,fontface="italic"),
        column_split = factor(col_order$Grouping),
        cluster_column_slices = FALSE,
        column_gap=unit(2,"mm"),
        column_title="Red Module", 
        #column_title_side = "bottom",
        column_title_gp=gpar(fontsize=15,fontface="bold"),
        column_names_gp=gpar(fontsize=8)
)

dev.off()



#Plotting HM for GY
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/")
redmodule <- read.csv("Red Module/Redmodulegenes2.csv")
redmodule <- read.csv("Greenyellow Module/Greenmodulegenes6.csv")
colnames(redmodule)
redmodulegrouping <- redmodule[,c(23,10,12),drop=FALSE]
redmodulegrouping <- redmodule[,c(36,16,19),drop=FALSE] #if useing redmodulegenes2 #shared.name,gene.family.custom,Grouping
redmodulegrouping <- redmodule[,c(26,12,14),drop=FALSE] #for greenmodules5
redmodulegrouping$Grouping <- factor(redmodulegrouping$Grouping)
levels(redmodulegrouping$Grouping)
redmodule <- dplyr::rename(redmodule, LOC = shared.name) 
redmodule <- redmodule[,"LOC",drop=FALSE]
setwd("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/")
dataHM_ZS <- read.csv("data_heatmap_zscored.csv_TPM10.csv")
dataFC <- read.csv("filtALL_no_0P_cond.csv")
colnames(dataFC)
dataFC <- dataFC[, !grepl("FDR", colnames(dataFC))]
colnames(dataFC)
dataFC <- dataFC[,c(1,2,3,15,16,17,12,13,14,10,11,12,7,8,9,4,5,6)]
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/")
write.csv(dataFC,"./Red Module/dataFC2.csv",row.names=FALSE)
#now go into the actual file and add columns for the 1 Pi treatment groups manuall
#load back in
dataFC <- read.csv("./Red Module/dataFC2.csv")
colnames(dataFC)=colnames(dataHM_ZS)
write.csv(dataFC,file="dataFC3.csv",row.names=FALSE)
#START HERE FOR dataFC file since already processed.
dataFC <- read.csv("dataFC3.csv",header=TRUE)
getwd()
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/unsigned_v6/Red Module/")
#OK, now let's filter for genes we care about.
redmodulegenes <- left_join(redmodule,dataHM_ZS, by="LOC")
redmodulegenes <- redmodulegenes[,-c(2,3)] %>% 
  tibble:: column_to_rownames(var = "LOC") 

redmodulegenes <- data.matrix(redmodulegenes)

dataFCgenes <- left_join(redmodule,dataFC, by="LOC")
dataFCgenes <- dataFCgenes[,-c(2,3)] %>% 
  tibble:: column_to_rownames(var = "LOC") %>% round(1)

dataFCgenes <- data.matrix(dataFCgenes)


colfun2 = circlize::colorRamp2(c(-2,0,2), #gives range of the scale
                               c("#408ea4","white","#c5242a")) #for zscoring

pdf("./HeatmapGYModule_LOC4.pdf",height=12,width=6)
#draw first HM
Heatmap(redmodulegenes,
        col =colfun2,
        #use_raster=FALSE,
        cluster_columns=FALSE,
        #cluster_rows=TRUE,
        show_row_names = TRUE,
        cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
        cluster_rows=TRUE,
        #clustering_method_columns="ward.D",
        #clustering_method_rows = "ward.D",
        row_dend_reorder = FALSE,
        # row_km = 12,
        #column_km=9,
        #column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        #row_gap = unit(2, "mm"),
        row_split=factor(redmodulegrouping$Grouping), #row_split=paste(row_order$Grouping)
        row_title_gp = gpar(font = 2,fontsize=5),
        row_names_gp = gpar(fontsize = 5),
        column_split = factor(col_order$Grouping),
        cluster_column_slices = FALSE,
        column_gap=unit(2,"mm"),
        column_title="Pi Regulatory genes", 
        #column_title_side = "bottom",
        column_title_gp=gpar(fontsize=15,fontface="bold")
)
draw()
dev.off()

#ok, lets think about making this HM more informative.
#draw HM to include LFC

pdf("./HeatmapGYModule_geneID3.pdf",height=12,width=6)

Heatmap(redmodulegenes,
        col =colfun2,
        row_labels = redmodulegrouping$gene.family.custom,  # Replace row names with gene names
        cell_fun = function(j, i, x, y, width, height, fill) {
          # Ensure dataFCgenes has valid values
          if (!is.na(dataFCgenes[i, j])) {
            grid.text(as.character(dataFCgenes[i, j]), x, y, 
                      gp = gpar(fontsize = 5,fontface="bold", col = "black"))
          }
        },
        use_raster=TRUE,
        cluster_columns=FALSE,
        #cluster_rows=TRUE,
        show_row_names = TRUE,
        cluster_row_slices = FALSE, #IMPORTANT TO HAVE IT OFF IF YOU WANT TO REODER PER YOUR LIST.
        cluster_rows=TRUE,
        #clustering_method_columns="ward.D",
        #clustering_method_rows = "ward.D",
        row_dend_reorder = FALSE,
        # row_km = 12,
        #column_km=9,
        #column_km_repeats=100,
        # row_km_repeats = 100,#divide into clusters by kmeans, repeat 100 use consensus
        #row_gap = unit(2, "mm"),
        row_split=factor(redmodulegrouping$Grouping), #row_split=paste(row_order$Grouping)
        row_title_gp = gpar(font = 2,fontsize=5),
        row_names_gp = gpar(fontsize = 5,fontface="italic"),
        column_split = factor(col_order$Grouping),
        cluster_column_slices = FALSE,
        column_gap=unit(2,"mm"),
        column_title="GY Module", 
        #column_title_side = "bottom",
        column_title_gp=gpar(fontsize=15,fontface="bold"),
        column_names_gp=gpar(fontsize=8)
)

dev.off()
getwd()

