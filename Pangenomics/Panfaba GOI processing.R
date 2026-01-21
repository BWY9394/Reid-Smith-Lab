#load packages.
#Core packages
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggforce)
library(ggrepel)
library(ggsci)
library(multcompView)
library(emmeans)
library(multcomp)
library(gridExtra)
library(patchwork)
library(stringr)
library(Hmisc)
library(dplyr)
library(plyr)


#start here
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/Reid Lab")

#load in dataframe for cleaning
Roy <- read.csv("./Panfaba work/Roy_Raw_Nod_Genes.csv",header=TRUE)
View(Roy)

#Great
expanded_Roy <- Roy %>%
  # Separate the 'value' column by comma, creating new rows for each
  separate_rows(MT4, sep = ",") %>%
  # Optional: Trim any leading/trailing whitespace that might result from the split
  mutate(value = trimws(MT4)) 


#Now we merge datasets
MT4_T5 <- read.csv("./Panfaba work/MT4-MT5 Custom query.csv",header=TRUE)
View(MT4_T5)
MT4_T5 <- MT4_T5[,c(1,7)]
colnames(MT4_T5)

MT4_T5 <- MT4_T5 %>%
  dplyr::rename(MT4 = objectB,
                MT5 = X.objectA) %>%
  dplyr::distinct(MT4, .keep_all = TRUE)



Mergetest <- left_join(expanded_Roy,MT4_T5, by= "MT4") #235
