BiocManager::install("VennDetail")
BiocManager::install("hicVennDiagram")
library(devtools)
install_github("jokergoo/ComplexHeatmap")
install
install.packages("tidyverse")
install.packages("UpSetR")
library(sleuth)
library(tidyverse)
library(dplyr)
library(gridExtra)
library(grid)
library(UpSetR)

setwd("C:/Users/BWeeYang/Documents/Kallisto Output HP01")
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/RNAseq/HP0001/Kallisto Output HP01")
setwd("D:/Kallisto Output HP01/")


####################################################################################################################
####################################################################################################################

#annotation file  filtering out/transforming/keeping elements that will be useful for the analysis. want to group down to gene level, while retaining information in the features and annotation descriptions of what the transcript variants were.
annotraw=read.csv("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/GCFfeature table.csv",header=T, stringsAsFactors =FALSE)
annotraw <- annotraw[,c(15,1,14)]
colnames(annotraw)
annotraw <- annotraw %>% filter(str_detect(X..feature,"RNA"))
collapse <- annotraw %>% group_by(symbol)  %>%
  dplyr::summarise(X..feature = paste(unique(X..feature), collapse = ', '), #make sure u r using the dplyr version of summarise. plyr has a different logic.
            name = paste(unique(name), collapse = ', '))

collapse <- dplyr::rename(collapse,"LOC"="symbol")
write.csv(collapse,"collapse.csv",row.names=FALSE)
view(collapse)

####################################################################################################################
####################################################################################################################

#loading in base file from tximport + filtered annotation file.
annotbase=read.csv("C:/Users/BWeeYang/Dropbox/Benjamin/Experiments/P Series/HP_0001_E1 Han NW/Data/RNAseq/kallisto-out/CS10_tximport.csv",header=T, stringsAsFactors =FALSE) #using the same target_id to gene_id file in tximport script. have to remember to ask oli how he got that.

collapse=read.csv("collapse.csv",header=T, stringsAsFactors =FALSE)


annot <- left_join(annotbase,collapse,by="LOC") # adding annotations to the base file. there will be repeats because there are multiple transcripts variants mapping to the same gene.
#annot=annot %>% filter(str_detect(X..feature,"gene"))
#annot <- annot %>% filter(str_detect(X..feature,"mRNA"))
#annot <- annot %>% filter(str_detect(product_accession,"X"))
#annot <- rename(annot,"target_id"="product_accession")
annot <- rename(annot,"Short_description"="name")
annot <- rename(annot,"feature"="X..feature")
view(annot)
colnames(annot)
#annot <- annot[,c(11,15,1,14)]
write.csv(annot,"annot.csv",row.names = FALSE)

colnames(annot)

####################################################################################################################
####################################################################################################################
#START HERE
annot_header=c("target_id","feature","Short_description") #this has to match whatever your output min results table column headers. target_id is used because it collapses the target_id and LOC column headers into target_id, matching transcript IDs to gene IDs and summarzing at gene level. So your gene IDs will be retained in the target_id column. I'm not sure why but including the LOC argument here messes things up. have to remember to document what it messes up.

annot <- read.csv("annot.csv",header=TRUE,stringsAsFactors = FALSE) 
colnames(annot)

#load sample definition file, three columns:
#A samplebas
#B condition
#C path to each file for the above, needs to be in the above base_dir
#IMPORTANT: the cwarnigsondition gets sorted alphabetically and then the first is the 'control' treatment
#therefore make sure that it is first in order, e.g. name this condition aaahighP
base_dir <- "C:/Users/BWeeYang/Documents/Kallisto Output HP01"
base_dir

s2c=read.csv("sleuthTFsampleinfo.csv",header=T, stringsAsFactors =FALSE)
s2c

     
#now everything is defined to actually start quantification

#Now the "sleuth object" can be constructed. This requires four commands that 
#(1) load the kallisto processed data into the object 
#(2) estimate parameters for the sleuth response error measurement (full) model 
#(3) estimate parameters for the sleuth reduced model, and 
#(4) perform differential analysis (testing). 

#first prep data, load processed data into the sleuth object 'so'
#this filters for lowly expressed genes etc. this filter setting can be adjusted 
#by adding parameters to the function
#aggregation is important if kallisto has been run on isoforms, but you want to have gene level results
#as kallisto has been run on representative transcripts, not necessary here


# define the filter, min_reads: how many reads in a sample, min_prop: in what proportion of the sample
# proportion is important, as there might be a condition where a gene is not expressed
#no. of minimal reads, better to have it early on. min_prop-> out of total no. reps, what proportion 


basic_filter <- function (row, min_reads = 5, min_prop = 0.47) 
{
  mean(row >= min_reads) >= min_prop
}

#custom filter if using 2 factors. 
design_filter <- function(design, row, min_reads=5, min_prop = 0.47){
  sum(apply(design, 2, function(x){
    y <- as.factor(x);
    return(max(tapply(row, y, function(f){sum(f >= min_reads)})/
                 tapply(row, y, length)) == 1 
           || basic_filter(row, min_reads, min_prop)
    )
  })) > 0}

#(1) load processed kallisto data into sleuth object 'so',  #need to ask oli why not incorporating full_model=design argument
so47 <- sleuth_prep(s2c, 
                  filter_fun = basic_filter, #filter_fun = function(x)(design_filter(so47$design,x)), if using different filter https://www.biostars.org/p/443205/
                  ~ condition, # could add other treatments, e.g condition + genotype
                  target_mapping = annot,
                  extra_bootstrap_summary = TRUE,
                  aggregation_column = "LOC", # defines a gene level analysis, leave out for isoforms
                  gene_mode=TRUE, #needed to add this https://github.com/pachterlab/sleuth/issues/164
                  transform_fun_counts = function(x) log2(x + 0.5))



#(2) fit the full model, model gene expression to condition(s) as defined above
so <- sleuth_fit(so47)


#to check the models:
models(so)
#this gives the table showing what the models(s) look(s) like, it basically shows the conditions
#except the 'control' everything is calculated againts - this is the (intercept)
#check that this is what you wanted


#3fit the reduced model, where gene expression is not dependent on any factor
so <- sleuth_fit(so, ~1, 'reduced')

models(so) 

#this performs the test to identify transcripts with a significantly better fit with the “full” model
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)

tests(so)




#remove control need oli's help with this bit.
conditions = unique(s2c$condition)
view(conditions)
cond2test = conditions[-1]
cond_number = length(cond2test)
view(cond_number)

for (condition in cond2test)
{
  current=paste("condition",condition,sep="")
  so <- sleuth_wt(so, which_beta = current)
  results_table <- sleuth_results(so, current, test_type = 'wd')
  
  #this rest of this for loop only summarises and merges data tables
  results_table <- results_table[order(results_table$target_id),]
  write.csv(results_table, file=paste0(condition,".csv", sep=""),row.names = FALSE)
  
  #first create a data.frames with the right row number and annot in the first loop
  if(condition == conditions[2]) # i.e. first diffex condition
  {
    sum_table = results_table %>% dplyr::select(all_of(annot_header))
    sum_table = sum_table[order(sum_table$target_id),]
  } 
  
  small_table = results_table %>% select(target_id, qval, b)
  small_table = small_table[order(small_table$target_id),]
  header <- c("LOC", paste0("FDR_", condition,sep=""), paste0("logFC_", condition,sep=""),row.names = FALSE)
  colnames(small_table) <-header
  write.csv( small_table, file=paste0(condition,"_small.csv"))
  
  
  #create summary for all conditions
  sum_table = sum_table[order(sum_table$target_id),]
  sum_table <- data.frame(sum_table, small_table[,-1])
  
  
}

sum_table = sum_table[order(sum_table$target_id),]
write.csv( sum_table, file=paste("eCol_summaryTF.csv"),row.names = FALSE) #row.names=F, quote=F) ,not sure what this does 
####################################################################################################################
####################################################################################################################



#Wald_test for DE
#not necessary as loop above already does it but here it is anyway. so <- sleuth_wt(so, which_beta = current)


#once done, can choose to clear environment and reload in your table outputs of choice for further cleanup and filtering.
#but first lets filter for p<0.05 and then logfc>1.5 and logfc<-1.5. Probably best to do it across the rows so don't lose information
colsumLR=read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/eCol_summaryLR.csv",stringsAsFactors = F,header=T)
colsumOL=read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/eCol_summaryOL.csv",stringsAsFactors = F,header=T)
colsumML=read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/eCol_summaryML.csv",stringsAsFactors = F,header=T)
colsumYS=read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/eCol_summaryYS.csv",stringsAsFactors = F,header=T)
colsumTF=read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/eCol_summaryTF.csv",stringsAsFactors = F,header=T)

colnames(colsumLR)
view(colsumLR)
filtcolsumLR <- colsumLR %>% filter(if_any(c(5:10),complete.cases))
filtcolsumLR <- colsumLR %>% filter(if_any(c(5:10), negate = TRUE, ~complete.cases(.)))
colnames(colsumLR)


#ok this approach sucks. Alternative approach.Per organ type, load in smalltables first then combine individual dataframes into a list. "D:/Sleuth Results/ if working from TD
annot_genelevel=read.csv(file ="annot_genelevel.csv",stringsAsFactors = F, header= T)
head(annot_genelevel)
#start here if looping
OL0P <- read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/bbOLP0_small.csv",stringsAsFactors = F,header=T)
OL0P <-OL0P %>% filter(if_any(c(3:4),complete.cases))
filtOL0P <- OL0P %>% filter(FDR_bbOLP0< 0.05, logFC_bbOLP0> 1|logFC_bbOLP0< -1)
filtOL0P <- left_join(filtOL0P,annot_genelevel,by="LOC") # you better check at this satge that no.observations don't change or it means your joing table failed miserably dude
colnames(filtOL0P)
filtOL0P <- select(filtOL0P,LOC,name,FDR_bbOLP0,logFC_bbOLP0)
write.csv(filtOL0P,file="filt_OL0P.csv",row.names=FALSE)


OL025P <- read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/ccOLP025_small.csv",stringsAsFactors = F,header=T)
OL025P <-OL025P %>% filter(if_any(c(3:4),complete.cases))
colnames(OL025P)
filtOL025P <- OL025P %>% filter(FDR_ccOLP025< 0.05, logFC_ccOLP025> 1|logFC_ccOLP025< -1)
filtOL025P <- left_join(filtOL025P,annot_genelevel,by="LOC") # you better check at this satge that no.observations don't change or it means your joing table failed miserably dude
colnames(filtOL025P)
filtOL025P <- select(filtOL025P,LOC,name,FDR_ccOLP025,logFC_ccOLP025)
write.csv(filtOL025P,file="filt_OL025P.csv",row.names=FALSE)

OL02P <- read.csv(file = "C:/Users/BWeeYang/Documents/Sleuth Results/ddOLP2_small.csv",stringsAsFactors = F,header=T)
OL02P <-OL02P %>% filter(if_any(c(3:4),complete.cases))
colnames(OL02P)
filtOL2P <- OL02P %>% filter(FDR_ddOLP2< 0.05, logFC_ddOLP2> 1|logFC_ddOLP2< -1)
filtOL2P <- left_join(filtOL2P,annot_genelevel,by="LOC") # you better check at this satge that no.observations don't change or it means your joing table failed miserably dude
colnames(filtOL2P)
filtOL2P <- select(filtOL2P,LOC,name,FDR_ddOLP2,logFC_ddOLP2)
write.csv(filtOL2P,file="filt_OL2P.csv",row.names=FALSE)


#making summarized frames.

OLall <- full_join(filtOL0P, filtOL025P, by = c("LOC","name")) %>%
  full_join(filtOL2P, by = c("LOC","name"))

filtML0P <- read.csv("filt_ML0P.csv",header = T)
filtML025P <- read.csv("filt_ML025P.csv",header = T)
filtML2P <- read.csv("filt_ML2P.csv",header = T)
MLall <- full_join(filtML0P, filtML025P, by = c("LOC","name"))# %>%
  full_join(filtML2P, by = c("LOC","name"))

filtYS0P <- read.csv("filt_YS0P.csv",header = T)
filtYS025P <- read.csv("filt_YS025P.csv",header = T)
filtYS2P <- read.csv("filt_YS2P.csv")
YSall <- full_join(filtYS0P, filtYS025P, by = c("LOC","name")) #%>%
  full_join(filtYS2P, by = c("LOC","name"))

filtTF0P <- read.csv("filt_TF0P.csv",header = T)
filtTF025P <- read.csv("filt_TF025P.csv",header = T)
filtTF2P <- read.csv("filt_TF2P.csv",header = T)
TFall <- full_join(filtTF0P, filtTF025P, by = c("LOC","name")) # %>%
  full_join(filtTF2P, by = c("LOC","name"))


filtLR0P <- read.csv("filt_LR0P.csv",header = T)
filtLR025P <- read.csv("filt_LR025P.csv",header = T)
filtLR2P <- read.csv("filt_LR2P.csv",header = T)
LRall <- full_join(filtLR0P, filtLR025P, by = c("LOC","name")) %>%
  full_join(filtLR2P, by = c("LOC","name"))

#custom summarized frame.
OLall_no_0P <- full_join(filtOL025P,filtOL2P, by = c("LOC","name",))
#conditional merge to include only genes that have been DE in 025 and 2P, without the specific OL0P genes
OLall_cond_0P <- left_join(OLall_no_0P,filtOL0P,by = c("LOC","name"))
write.csv(OLall_cond_0P,file="filtOLall_noDEGOL0P.csv",row.names=FALSE)

filtALL_no_0P <- full_join(TFall,YSall, by=c("LOC","name")) %>%
  full_join(MLall,by=c("LOC","name"))%>%
  full_join(OLall_cond_0P,by=c("LOC","name"))%>%
  full_join(LRall,by=c("LOC","name"))
colnames(extra)
rejig <- select(extra,1,2,27,21,22,23,24,25,26,19,20,15,16,17,18,22,15,11,12,13,14,7,8,9,10,3,4,5,6) #note to self that I manually added in the empty 2P rows for relevant organs e.g. ML YS and TF which did not have ANY DEG for 2P, using the original filt_all as base.
colnames(rejig)
write.csv(rejig,"filtALL_no_0P_cond.csv",row.names=FALSE)
#official summarized frame.
filtALL <- full_join(TFall,YSall, by=c("LOC","name")) %>%
  full_join(MLall,by=c("LOC","name"))%>%
  full_join(OLall,by=c("LOC","name"))%>%
  full_join(LRall,by=c("LOC","name"))

colnames(filtALL)
head(annot_genelevel)

#If want to include the features back in
extra <- left_join(filtALL_no_0P,annot_genelevel,by=c('LOC',"name"))
colnames(extra)
head(extra)            
colnames(filtALL)
write.csv(extra,"filt_ALL.csv",row.names=FALSE)
filtMasterRef

#ok dipshit, let's do this. Want to get a count of the DEGs ONLY expressed in old leaves, 0P
ALL <- read.csv("DEG_TPM_all_atleast10TPM.csv")
OLALL <- read.csv("filt_ALL.csv")
AllFC <- left_join(ALL,OLALL,by=c('LOC'))
AllFC <- AllFC[,-(2:25)]
#View(AllFC)
colnames(AllFC)
#Run the next two steps only if you are dealing with just the OL tissue type dataframe
# Assuming your data frame is named df
AllFC <- AllFC[, !grepl("FDR", colnames(AllFC))]
colnames(AllFC)
write.csv(AllFC,"filt_ALL_TPM10.csv",row.names = FALSE) #10327 genes that are DE and TPM>10 
# Ru
AllFC <- AllFC[!(is.na(AllFC$logFC_bbOLP0) & is.na(AllFC$logFC_bbOLP0)), ]
colnames(AllFC)


#Run from here if dealing with all organs dataframe that has FDR columns as well
# Filter for genes differentially expressed ONLY in respective organ types
#create column name lsit
filtered_df <- AllFC[
  !is.na(AllFC$logFC_ccMLP025) & AllFC$logFC_ccMLP025 != 0 & #prefilters df so that condition of interest will only keep rows with DE
    rowSums(!is.na(AllFC[, grepl("logFC_", colnames(AllFC))]) & AllFC[, grepl("logFC_", colnames(AllFC))] != 0) == 1, #All other logFC_ columns must be NA or 0, because we have already done the prefilter df step
]


# Count of genes with logFC_bbOLP0 > 0
count_greater <- sum(filtered_df$logFC_ccMLP025 > 0)

# Count of genes with logFC_bbOLP0 < 0
count_less <- sum(filtered_df$logFC_ccMLP025 < 0)

# Print the results
cat("Count of genes with  > 0:", count_greater, "\n")
cat("Count of genes with  < 0:", count_less, "\n")
colnames(filtered_df)

write.csv(filtered_df,file="filt_OL0P_only_TPM10.csv",row.names=FALSE)


###automates the above
# Create a list to store the filtered dataframes
filtered_list <- list()

# Loop over each logFC column (excluding the first column which is 'LOC')
logFC_columns <- grep("^logFC_", colnames(AllFC), value = TRUE)

# Loop through each logFC column and apply the filtering criteria
for (col in logFC_columns) {
  # Filter the dataset based on the logFC column
  filtered_df <- AllFC[
    !is.na(AllFC[[col]]) & AllFC[[col]] != 0 &
      rowSums(!is.na(AllFC[, grepl("logFC_", colnames(AllFC))]) & AllFC[, grepl("logFC_", colnames(AllFC))] != 0) == 1, 
  ]
  
  # Store the filtered data frame in the list with the column name as the list name
  filtered_list[[col]] <- filtered_df
}

# Print the list of filtered data frames
names(filtered_list)  # Display the names of the data frames in the list

##test
# Loop through each data frame in the filtered_list
for (df_name in names(filtered_list)) {
  
  # Get the filtered data frame from the list
  filtered_df <- filtered_list[[df_name]]
  
  # Get the respective column name for logFC (it's the same as the df_name)
  logFC_column <- df_name
  
  # Count the number of genes with logFC > 0
  count_greater <- sum(filtered_df[[logFC_column]] > 0, na.rm = TRUE)
  
  # Count the number of genes with logFC < 0
  count_less <- sum(filtered_df[[logFC_column]] < 0, na.rm = TRUE)
  
  # Print the results for the current data frame and column
  cat("For", df_name, ":\n")
  cat("Count of genes with", logFC_column, "> 0:", count_greater, "\n")
  cat("Count of genes with", logFC_column, "< 0:", count_less, "\n")
  cat("\n")
}


#getting counts for all DEG and also TPM>10. OL0p included. So doesn't care that it's organ:treatment specific, provides overview without going into organ:treatment specificity i.e. is there a subset of genes solely expressed in an organ:treatment combo
filt_ALL_no_0P <- read.csv("filt_ALL_TPM10.csv")
filt_ALL_no_0P <- filt_ALL_no_0P[, !grepl("FDR", colnames(filt_ALL_no_0P))]
colnames(filt_ALL_no_0P)

###manual way of getting counts
count_greater <- sum(filt_ALL_no_0P$logFC_ddOLP2  > 0,na.rm=TRUE)
count_less <- sum(filt_ALL_no_0P$logFC_ddOLP2  < 0,na.rm=TRUE)
# Print the results
cat("Count of genes with  > 0:", count_greater, "\n")
cat("Count of genes with  < 0:", count_less, "\n")
colnames(filt_ALL_no_0P)


###Automates the above
# List of columns with 'logFC' that you want to analyze
logFC_columns <- grep("^logFC", colnames(filt_ALL_no_0P), value = TRUE)

# Loop over each logFC column to calculate the counts for values > 0 and < 0
count_results <- sapply(logFC_columns, function(col) {
  count_greater <- sum(filt_ALL_no_0P[[col]] > 0, na.rm = TRUE)
  count_less <- sum(filt_ALL_no_0P[[col]] < 0, na.rm = TRUE)
  c(Count_Greater = count_greater, Count_Less = count_less)
})

# Print the results
for (col in logFC_columns) {
  cat(paste("For column", col, ":\n"))
  cat("Count of genes with logFC > 0:", count_results["Count_Greater", col], "\n")
  cat("Count of genes with logFC < 0:", count_results["Count_Less", col], "\n\n")
}

####################################################################################################################
####################################################################################################################
#getting DEG counts #piping from  earlier. these dataframes omit the FDR adjusted pvalues
filtOLall <- full_join(filtOL0P, filtOL025P, by = c("LOC","name")) %>%
  full_join(filtOL2P, by = c("LOC","name"))
filtOLall <- filtOLall[,-c(3,5,7)]
view(filtOLall)
write.csv(filtOLall,file='filt_OLALL.csv',row.names=FALSE)



#get counts into single dataframe
df_list <- list(filtTFall,filtYSall,filtMLall,filtOLall,filtLRall)

#cycle through list for DEG counts

counts_list <- lapply(df_list, function(df) {
  df %>%
    summarise(across(.cols = (ncol(df) - 2):ncol(df), .fns = ~sum(!is.na(.)))) %>%
    pivot_longer(cols = everything(), names_to = "treatment", values_to = "DEG")
})

print(counts_list)
t <- bind_rows(counts_list,.id="organ")
view(t)
# Convert the list of counts data frames into a single data frame
counts_df <- bind_rows(counts_list, .id = "organ")
counts_df <- counts_df %>%
  mutate(organ = recode(organ, 
                             "1" = "TF", 
                             "2" = "YS",
                             "3" = "ML",
                        "4" = "OL",
                        "5"="LR"))

print(counts_df)
# Define the repeating sequence
repeating_seq <- rep(c(0, 0.25, 2), length.out = nrow(counts_df))
print(repeating_seq)

# Replace the values in the 'column' column with the repeating sequence
counts_df <- counts_df %>%
  mutate(Treatment = repeating_seq) %>% select(1,4,3)
print(counts_df)

counts_df$organ <- factor(counts_df$organ, levels = c("TF", "YS", "ML", "OL", "LR"))
counts_df$Treatment <- as.factor(counts_df$Treatment)
levels(counts_df$organ)
levels(counts_df$Treatment)

write.csv(counts_df,file="counts_df.csv",row.names=FALSE)
counts_df <- read.csv("counts_df.csv",stringsAsFactors = TRUE,header = TRUE)
counts_df$Treatment <- as.factor(counts_df$Treatment)
scientific_colors <- c("0" = "#2CA000",  # Blue
                       "0.25" = "#FFF000", # Orange
                       "2" = "#BF0000")   # Green

ggplot(counts_df,aes(x=Treatment,y=DEG,fill=Treatment))+
  geom_bar(stat="identity")+
  facet_wrap(~organ,nrow=1)+
  geom_text(aes(label = DEG), vjust = -0.01,size=4.)+
  theme_bw()+
  theme(axis.title.x = element_text(color="black", size=15, face="bold"),
                     axis.title.y = element_text(color="black", size=15, face="bold"),
                     axis.text.x = element_text(size = 15, angle =0, vjust =0.5, face = "bold"),
                     axis.text.y = element_text(size = 15, vjust =0.5, face = "bold"),
                     strip.text.x = element_text(size = 15, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
                     panel.grid.major.x = element_blank(),
                     panel.grid.minor.y = element_blank()+
  theme(legend.key.size = unit(15, 'cm'), #change legend key size
                legend.key.height = unit(15, 'cm'), #change legend key height
                legend.key.width = unit(15,'cm'), #change legend key width
                legend.title = element_text(size=14), #change legend title font size
                legend.text = element_text(size=10)) #change legend text font size
  )
ggsave("DEGbarplots.pdf",height=20,width=20,units = "cm")


#DEGs using custom list built from previous section workflow to show down/upregulated genes
UpDown <- read.csv("DEG_forbarplots.csv")
UpDown2 <- read.csv("DEG_forbarplots_TPM10.csv")
print(UpDown)

# Melt the data for visualization
# Melt the data for visualization
df_long <- UpDown2 %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Count") %>%
  mutate(IsSpecial = !is.na(IsolateUp) | !is.na(IsolateDown)) # Highlight special rows
df_long$Treatment <- as.factor(df_long$Treatment)
df_long$Organ <- factor(df_long$Organ, levels=c("LR","OL","ML","YS","TF"))
levels(df_long$Organ)
df_long
# Plot
a=ggplot(df_long, aes(x = Treatment, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge",alpha=0.65) +
  geom_point(data = df_long %>% filter(IsSpecial),
             aes(x = Treatment, y = ifelse(Direction == "Up", IsolateUp, IsolateDown),color=Direction),
             #color ="white" , 
             size = 3, shape = 4) +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5) + # Add gene count labels
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) + # Flip the colors
  labs(title = "Upregulated and Downregulated Genes by Organ and Treatment",
       x = "Organ:Treatment",
       y = "Gene Count",
       fill = "Regulation") +
  facet_wrap(~ Organ,nrow=1) + # Add facet_wrap for organ type
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(size = 15, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 15, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 15, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()+
          theme(legend.key.size = unit(15, 'cm'), #change legend key size
                legend.key.height = unit(15, 'cm'), #change legend key height
                legend.key.width = unit(15,'cm'), #change legend key width
                legend.title = element_text(size=14), #change legend title font size
                legend.text = element_text(size=10))+ #change legend text font size
  theme(axis.text.x = element_text(angle = 45, hjust = 1)))
a

#alternative with barplots
a=ggplot(df_long, aes(x = Treatment, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = "dodge",alpha=0.65) +
  geom_bar(data = df_long,
           aes(x = Treatment, y = ifelse(Direction == "Up", IsolateUp, IsolateDown)),
           #color ="white" , 
           stat = "identity", position = "dodge") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.8), vjust = -0.5) + # Add gene count labels
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) + # Flip the colors
  labs(title = "Upregulated and Downregulated Genes by Organ and Treatment",
       x = "Organ:Treatment",
       y = "Gene Count",
       fill = "Regulation") +
  facet_wrap(~ Organ,nrow=1) + # Add facet_wrap for organ type
  theme_bw() +
  theme(axis.title.x = element_text(color="black", size=15, face="bold"),
        axis.title.y = element_text(color="black", size=15, face="bold"),
        axis.text.x = element_text(size = 15, angle =0, vjust =0.5, face = "bold"),
        axis.text.y = element_text(size = 15, vjust =0.5, face = "bold"),
        strip.text.x = element_text(size = 15, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()+
          theme(legend.key.size = unit(15, 'cm'), #change legend key size
                legend.key.height = unit(15, 'cm'), #change legend key height
                legend.key.width = unit(15,'cm'), #change legend key width
                legend.title = element_text(size=14), #change legend title font size
                legend.text = element_text(size=10))+ #change legend text font size
          theme(axis.text.x = element_text(angle = 45, hjust = 1)))
df_long

#DEG barplot with up down bars
a=ggplot(df_long, aes(x = Treatment, y = ifelse(Direction == "Up", Count, -Count), fill = Direction)) + #If Direction is "Up", the y value will be set to Count (positive value). If Direction is not "Up" (presumably "Down"), the y value will be set to -Count (negative value).
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) + 
  geom_bar(data = df_long,
           aes(x = Treatment, y = ifelse(Direction == "Up", IsolateUp, -IsolateDown)),
           #color ="white" , 
           stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_dodge(width = 0.1), vjust = ifelse(df_long$Direction == "Up", -0.6, 1.2)) + # Adjust text position for up and down
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) + # Flip the colors
  labs(title = "Upregulated and Downregulated Genes by Organ and Treatment",
       x = "Organ:Treatment",
       y = "Gene Count",
       fill = "Regulation") +
  facet_wrap(~ Organ, nrow = 1)+  # Add facet_wrap for organ type
  ylim(-2000,2000)+
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = 15, face = "bold"),
        axis.title.y = element_text(color = "black", size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 15, vjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold", margin = margin(0.1, 0, 0.1, 0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))

a
ggsave(a,file="UpDownDEGS_regular_ylim.pdf")
#DEG barplot with up down bars, flipped coords
b=ggplot(df_long, aes(y = Treatment, x = ifelse(Direction == "Up", Count, -Count), fill = Direction)) + #If Direction is "Up", the y value will be set to Count (positive value). If Direction is not "Up" (presumably "Down"), the y value will be set to -Count (negative value).
  geom_bar(stat = "identity", position = "stack", alpha = 0.65) + 
  geom_bar(data = df_long,
           aes(y = Treatment, x = ifelse(Direction == "Up", IsolateUp, -IsolateDown)),
           #color ="white" , 
           stat = "identity", position = "stack") +
  geom_text(aes(label = Count), position = position_dodge(width = -0.6), vjust = ifelse(df_long$Direction == "Up", -0.15, 0.9)) + # Adjust text position for up and down
  scale_fill_manual(values = c("Down" = "#408ea4", "Up" = "#c5242a")) + # Flip the colors
  labs(title = "Upregulated and Downregulated Genes by Organ and Treatment",
       y = "Organ:Treatment",
       x = "Gene Count",
       fill = "Regulation") +
  facet_wrap(~ Organ, ncol = 1)+  # Add facet_wrap for organ type
  #ylim(-2000,2000)+
  theme_bw() +
  theme(axis.title.x = element_text(color = "black", size = 15, face = "bold"),
        axis.title.y = element_text(color = "black", size = 15, face = "bold"),
        axis.text.x = element_text(size = 15, angle = 45, vjust = 0.5, face = "bold"),
        axis.text.y = element_text(size = 15, vjust = 0.5, face = "bold"),
        strip.text.x = element_text(size = 15, face = "bold", margin = margin(0.1, 0, 0.1, 0, "cm")),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.key.size = unit(1, 'cm'),
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))
b
ggsave(b,file="UpDownDEGS_flipped coords.pdf",width=25,units="cm")

# Man I'm talented at R. Just kidding I'm a master at making GPT generate boilerplate code for me teehee


####################################################################################################################
####################################################################################################################
#ok time to make new dataframes for upsetplots
LRupset <- filtLRall
MLupset <- filtMLall
OLupset <- filtOLall
YSupset <- filtYSall
TFupset <- filtTFall
view(OLupset)


# Iterate over Column3 and Column4-  repeat unique dataframe from XXupset series 

for (col_name in c("logFC_bbTFP0", "logFC_ccTFP025","logFC_ddTFP2")) {
  # Replace values in the current column with corresponding values from Column1 only if they are not NA
  not_na_indices <- !is.na(TFupset[[col_name]])
  TFupset[[col_name]][not_na_indices] <- TFupset$LOC[not_na_indices]
  view(TFupset)
}
#lol u didnt run the code again so your upsetall is now just raw logFCs u idiot. but then again shouldn't need them as LOCs/binary as the next steps convert that for u...maybe for venndiagram? xD


#OK now have individual frames..
UpsetALL <- full_join(TFupset, YSupset, by=c("LOC","name"))%>% 
  full_join(MLupset, by=c("LOC","name"))%>% 
  full_join(OLupset, by=c("LOC","name"))%>% 
  full_join(LRupset, by=c("LOC","name"))
view(UpsetALL)
colnames(UpsetALL)


write.csv(UpsetALL, "Upset_table_merged_all.csv",row.names=FALSE) #this is mastersheet going forward. difference between this and upsetAll is that this one only has logFC values.
UpsetALL <- read.csv("Upset_table_merged_all.csv",stringsAsFactors = FALSE,header = TRUE)
#
DEGS2 <- read.csv("DEG_TPM_all_no_0P_atleast10TPM.csv",stringsAsFactors = FALSE,header=TRUE)
DEGS2merge <- DEGS2[,1:2]
colnames(DEGS2merge)
Upsetfiltered <- left_join(DEGS2merge,UpsetALL,by="LOC")
view(Upsetfiltered)
#now prepDEGS2merge#now preparing for upset plotting.
test2=UpsetALL[,-2]
colnames(test2)
test2=Upsetfiltered[,-c(2:3)]

for (i in 2:16) {  #change range to match.
  # Replace values with 1 if they match the corresponding value in Column1, else replace with 0
  test2[[i]] <- ifelse(test2[[i]] == test2$LOC,1,0)
  view(test2)
}
test2 <- replace(test2, is.na(test2), 0) # replace NA values with 0
view(test2)


write.csv(test2, "Upset_table_IPMB_binary.csv",row.names=FALSE)#this one is for the filtered list for IPMB
write.csv(test2, "Upset_table_merged_all_binary.csv",row.names=FALSE)
test2=read.csv("Upset_table_IPMB_binary.csv",stringsAsFactors = FALSE,header=TRUE)
view(testing)
view(test3)
test3=test2[,-1]
colnames(test3)
names <- c("TF:P0","TF:P025","TF:P2","YS:P0","YS:P025","YS:P2","ML:P0","ML:P025","ML:P2","OL:P0","OL:P025","OL:P2","LR:P0","LR:P025","LR:P2")
test3=test3%>%rename_with(~names,everything())
test4=select(test3,13,14,15,10,11,12,7,8,9,4,5,6,1,2,3)
colnames(test4)
colheader=colnames(test)
colheader=colheader[-1]
print(colheader)


#plotting upsets
#first, across per treatments, across different organ types.See gene overlap/consistency across treatments
#first filter for per organ.so..

colnames(test3)

#
tf2plot <- test3[,-c(10,11,12)] #tester plot
upset(test3,
      order.by = "freq",
      nintersects=13,
      mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ",
      nsets=7,
      empty.intersections = "off"
)


tf2plot <- test3[,c(1,4,7,10,13)] #0P
colnames(tf2plot)
P0_organs=upset(tf2plot,
      order.by = "freq",
      nintersects=13,
      mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
      )
P0_organs
colnames(test3)
tf2plot <- test3[,c(1,4,7,13)] #0P,without OL
colnames(tf2plot)
P0_organs_wOL=upset(tf2plot,
        order.by = "freq",
        nintersects=9,
        mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
)
#grid.text("Distribution of DE Proteins", x = 0.65, y = 0.95,gp = gpar(fontsize = 10)) doesn't come out the way I want it to
dev.off()
P0_organs
P0_organs_wOL

colnames(test3)

tf2plot <- test3[,c(2,5,8,11,14)] #025P
colnames(tf2plot)
P025_organs=upset(tf2plot,
                  order.by = "freq",
                  nintersects=6,
                  mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
)


colnames(test3)
tf2plot <- test3[,c(2,5,8,14)] #025P,without OL
colnames(tf2plot)
P025_organs_wOL=upset(tf2plot,
                order.by = "freq",
                empty.intersections = "off",
                nintersects=6,
                mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
)

P025_organs
P025_organs_wOL

colnames(test3)




tf2plot <- test3[,c(3,6,9,12,15)] #2P
colnames(tf2plot)
P2_organs=upset(tf2plot,
                  order.by = "freq",
                  nintersects=13,
                  mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
)

colnames(test3)
tf2plot <- test3[,c(3,6,9,12,15)] #2P,without OL
colnames(tf2plot)
P2_organs_wOL=upset(tf2plot,
                      order.by = "freq",
                      nintersects=9,
                      mainbar.y.label = "Organ Intersections", sets.x.label = "DEGs Per Organ"
)

P2_organs
P2_organs_wOL


colnames(test3)

tf2plot <- test3[,c(1,2,3)] #0P
colnames(tf2plot)

organs_TF=upset((subset(test3[,c(1:3)])),
                order.by = "freq",
                empty.intersections = "on",
                nintersects=13,
)
organs_TF

organs_YS=upset((subset(test3[,c(4:6)])),
                order.by = "freq",
                nintersects=13,
)
organs_YS

organs_ML=upset((subset(test3[,c(7:9)])),
                order.by = "freq",
                nintersects=13,
)
organs_ML


organs_OL=upset((subset(test3[,c(10:12)])),
                order.by = "freq",
                nintersects=3,
)
organs_OL

organs_LR=upset((subset(test3[,c(13:15)])),
                order.by = "freq",
                nintersects=6,
)
organs_LR



####################################################################################################################
####################################################################################################################
#Extracting intersect Lists
testing=read.csv("Upset_table_merged_all_binary.csv",stringsAsFactors = FALSE,header=TRUE)
view(testing)
colnames(testing)
alpha=testing[,c(1,14,15,16)] # need column 1 as the identifier
colnames(alpha)


data_with_intersection <- alpha %>%
  tidyr::unite(col = "intersection", -c("LOC"), sep = "")

colnames(data_with_intersection)
head(data_with_intersection)


#display upset plot as intergers chr
intquery <- data_with_intersection %>%
  group_by(intersection) %>%
  dplyr::summarise(n = n()) %>%
  arrange(desc(n))
intquery

#extract out into a list. the binary numbers are your unique identifiers for the set of interest.
query <- data_with_intersection %>%
  group_by(intersection) %>%
  dplyr::summarise(list = list(LOC)) %>%
  mutate(list = setNames(list, intersection)) %>%
  pull(list)

intquery

#extract the specific interaction identifier
intersect_of_interest <- query$"111"
intersect_of_interest_DF <- setNames(data.frame(intersect_of_interest), "LOC")#save as df
colnames(intersect_of_interest_DF)
intersect_of_interest_DF <- left_join(intersect_of_interest_DF,annot_genelevel,by="LOC")


query[["01110"]]

####################################################################################################################
####################################################################################################################
#u moron, u jinxed urself. luckily u also saved a table with LOCs already converted! haha..haha..hahahahahahahahadsfgfasdf.
venntable <- read.csv("Upset_table_merged_all.csv")
head(venntables)



names <- c("TF_P0","TF_P025","TF_P2","YS_P0","YS_P025","YS_P2","ML_P0","ML_P025","ML_P2","OL_P0","OL_P025","OL_P2","LR_P0","LR_P025","LR_P2")
venntables= venntable%>% dplyr::rename_with(~names,.cols=c(3:17))
venntables=venntables[,-c(1,2)]
colnames(venntables)
list <- venntables[,c(1,2,3)]

print(as.data.frame(t(lapply(UpsetALL[, 3:17], function(col) sum(!is.na(col))))))#sanity check
sapply(list, function(x) sum(!is.na(x))) #sanity checks make sure x))) #sanity checks make sure 

vennlist <- lapply(list, function(x) na.omit(as.character(x))) #function(x) represents the function being applied here is an anonymous function defined using function(x), where x represents each element of venntables. lapply applies the fucntion each element of the list and returns a list of the results. In this case the function(x is applying na.omit and as.character

library(ggVennDiagram)

Uno <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("P0 without OL")
Dos <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("PO")
Tres <-ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("YS")
Quattro <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("ML")
Pento <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("OL")
Sixto <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("LR")
Hepto <- ggVennDiagram(vennlist) + scale_fill_gradient(low="grey90",high = "red")+ggtitle("TF")
Uno+ggtitle("P0 without OL")
Dos+ggtitle("P0")
Tres
Quattro
Pento
Sixto
Hepto
grid=grid.arrange(Uno,Dos,Tres,Quattro,Pento,Sixto,Hepto,ncol=3,nrow=3)

ggsave(grid,file="venn_sanity_checks.pdf",width=70,height=60,units="cm")
dev.off()
####################################################################################################################
####################################################################################################################

# Column to check for duplicates
checkingdupes=read.csv(file = "Upset_table_merged_all.csv",header=T, stringsAsFactors = F)
#or
checkingdupes=view(test2) 
column_to_check <- checkingdupes$LOC

# Check if any entries in the column are repeats


any_repeats <- any(duplicated(column_to_check))


# Output the result
if (any_repeats) {
  cat("There are repeated entries in the column.")
} else {
  cat("There are no repeated entries in the column.")
}

#
if (any_repeats) {
  cat("There are repeated entries in the column.\n")
  
  # Get the indices of duplicated rows
  duplicated_indices <- which(duplicated(column_to_check))
  
  # Display the duplicated rows
  duplicated_rows <- checkingdupes[duplicated_indices, ]
  print(duplicated_rows)
} else {
  cat("There are no repeated entries in the column.\n")
}

####################################################################################################################
####################################################################################################################



#RUBBISH

testing <- ifelse(is.na(filtLRall$logFC_bbLRP0),filtLRall$LOC, filtLRall$logFC_ddLRP0)
view(testing)
data$Column3 <- ifelse(is.na(data$Column3), data$Column1, data$Column3)
data$Column4 <- ifelse(is.na(data$Column4), data$Column1, data$Column4)

ggsave("filteredbarplots.png",height=29, width= 35, units='cm')

ggplot(counts_df, aes(x = condition, y = DEG,fill=condition)) +
  geom_bar(stat = "identity") +
  labs(x = "Treatments (mM P)", y = "DEGs(p<0.05 & -1>FDR>1)", title = "LR DEGs") +
  theme_minimal() +
  scale_fill_manual(values = scientific_colors,labels = c("logFC_bbLRP0" = "0", "logFC_ccLRP025" = "0.25", "logFC_ddLRP2" = "2"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  geom_text(aes(label = DEG), vjust = -0.5)+
  scale_x_discrete(labels = c("logFC_bbLRP0" = "0", "logFC_ccLRP025" = "0.25", "logFC_ddLRP2" = "2"))#+
#scale_fill_discrete(labels = c("logFC_bbLRP0" = "0", "logFC_ccLRP025" = "0.25", "logFC_ddLRP2" = "2"))

ggsave("test.pdf",width=39,height=50,units="cm")

png(paste0(Sys.Date(),"upset.png"), units="cm", width=8, height=8, res=300)
dev.off()
#making a list
a<- filtLR0P[,2] 
b <- filtLR025P[,2]
c <- filtLR2P[,2]
view(a)
view(c)
#LRList=list(filtLR0P,filtLR025P,filtLR2P)
LRList=list(a,b,c)

upset(fromList(LRList),order.by="frea")

plot(LRList, type = "upset")

#filtering for most significant genes



sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)

sleuth_significant <- dplyr::filter(results_table, qval <= 0.05)
sleuth_significant <- sleuth_significant[order(sleuth_significant$qval),]

head(sleuth_significant, 20)

a=plot_bootstrap(so, "XM_030624680.1", units = "est_counts", color_by = "condition")+ ggtitle("SQD2")
b=plot_bootstrap(so, "XM_030655136.1", units = "est_counts", color_by = "condition")+ ggtitle("SPX2")
c=grid.arrange(a,b,ncol=2,nrow=1)
a
dev.off()

sleuth_live(so)
dev.off()


#expression heatmap
oe <- sleuth_wt(so, 
                which_beta = 'conditionbbbLRP0')

# output results oe only necesary if you have more than 2 condition in your analysis aside from the control!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

sleuth_results_oe <- sleuth_results(oe, 
                                    test = 'conditionbbbLRP0', 
                                    show_all = TRUE)


sig_transcripts <- sleuth_results_oe %>% 
  filter(qval < 0.05)

plot_transcript_heatmap(so,transcripts=sleuth_significant$target_id[1:5],units="tpm")
print(heatmap)


ggplot2::ggsave("heatmap.pdf", width = 30, height = 25, units = "cm")


pdf(file = "PCA_plot.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) 


plot_pca(so, color_by = 'condition',text_labels=TRUE)
plot_pc_variance(so)
plot_sample_heatmap(so)
plot_group_density(so, 
                   use_filtered = TRUE, 
                   units = "est_counts",
                   trans = "log", 
                   grouping = "condition")

sig_transcripts <- sleuth_results_so %>% 
  filter(qval < 0.05)

plot_transcript_heatmap(so, 
                        transcripts = sleuth_significant$target_id[1:25],units="tpm")

ggsave("heatmap.pdf", width = 30, height = 25, units = "cm")
pdf(file = "heatmap.pdf",   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 8) 



dev.off()


quit()
