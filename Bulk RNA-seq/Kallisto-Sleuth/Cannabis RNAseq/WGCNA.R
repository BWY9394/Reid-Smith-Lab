

#Load in packages
library(WGCNA)
library(flashClust)
library(splitstackshape)
library(ggplot2)
library(corrplot)
library(reshape2)
library(tidyverse)


# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = 4)


#WD
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/")

###Preprocessing


#load in TPM file with raw means only
TPMmeans=read.csv("TPM_mean.csv")

TPMmeans$LRtotal=rowSums(TPMmeans[,2:5])
TPMmeans$MLtotal=rowSums(TPMmeans[,6:9])
TPMmeans$OLtotal=rowSums(TPMmeans[,10:13])
TPMmeans$TFtotal=rowSums(TPMmeans[,14:17])
TPMmeans$YStotal=rowSums(TPMmeans[,18:21])
colnames(TPMmeans)

DEGTPMmeansfilt <- TPMmeans[,c(1,22,23,24,25,26)] #29283 genes
TPMfilter10TPM <- DEGTPMmeansfilt[rowSums(DEGTPMmeansfilt[, -1] > 10) > 0, ] #18,665 genes NOT DEG yet but filtered to keep genes only >10 in at least 1 organ type

TPM_All_10TPM <- left_join(TPMfilter10TPM, TPMmeans, by="LOC")
TPM_All_10TPM <- TPM_All_10TPM[,-c(22:26)]
write.csv(TPM_All_10TPM, file="TPMfilter10TPM.csv",row.names=FALSE)





###Start here for pipeline

TPM_All_10TPM <- read.csv("TPMfilter10TPM.csv",stringsAsFactors = TRUE)
TPM_All_10TPM_transpose <- as.data.frame(t(TPM_All_10TPM[-c(1:1)]))
                                     

#extract gene names to import to new matrix/dataframe. Pretty handy function, have to bear in mind next time.
names(TPM_All_10TPM_transpose) = TPM_All_10TPM$LOC
names(TPM_All_10TPM_transpose)

#extract treatment IDs to import to new matrix/dataframe
rownames(TPM_All_10TPM_transpose) = names(TPM_All_10TPM)[-c(1:1)]
rownames(TPM_All_10TPM_transpose)

#TPM_All_10TPM_transpose<- apply(TPM_All_10TPM_transpose, 2, function(x) log2(x + 1))


##check to see if genes are all good
gsg = goodSamplesGenes(TPM_All_10TPM_transpose, verbose = 3)
gsg$goodSamples
gsg$allOK


#cluster samples to get an idea how samples are related etc.
sampleTree = flashClust(dist(TPM_All_10TPM_transpose), method = "average")


#plot tree, this can be then used to split samples by cutting at defined height
par(cex = 0.6); #controls size of text elements in plot
par(mar = c(0,4,2,0)) # Adjusts the margins of the plot, bottom marging, left margin, top marging, right margin
plot(sampleTree, main = "Sample clustering to detect outliers",
     sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut, change hight to see where tree will be cut
abline(h = 56000, col = "red")

dev.copy2pdf(file = "./WGCNA_sampleclustering_tree.pdf")


clust = cutreeStatic(sampleTree, cutHeight = 56000, minSize = 3)
table(clust)
print(clust) #9th entry is 0, as in this is cluster 0 and it has one treatment in it, which coincides with the row order.


keepSamples = (clust==1) #assigning clust 1 with 19 entries to be kept for downstream analysis, all samples except for OL0P
print(keepSamples)
WGCNAfiltered= TPM_All_10TPM_transpose[keepSamples, ] #creating a df which holds the samples, excluding the outlier samples, in the same transposed formatt

str(WGCNAfiltered)
adjustedrownames=c("TPM.LR_P0","TPM.LR_P025","TPM.LR_P1","TPM.LR_P2","TPM.OL_P025","TPM.OL_P1", "TPM.OL_P2","TPM.ML_P0","TPM.ML_P025","TPM.ML_P1", "TPM.ML_P2","TPM.YS_P0","TPM.YS_P025","TPM.YS_P1", "TPM.YS_P2", "TPM.TF_P0", "TPM.TF_P025","TPM.TF_P1", "TPM.TF_P2") #adjusting order
WGCNAfiltered <- WGCNAfiltered[adjustedrownames,,drop=FALSE]



#load in trait data
WGCNA <- read.csv("Traits_WGCNA.csv",stringsAsFactors = TRUE)
colnames(WGCNA$Organ)
levels(WGCNA$Organ)
rownames(WGCNAfiltered)
WGCNA$Organ <- factor(WGCNA$Organ,levels=c("Roots","Mature Leaf","Old Leaf","Top Flower","Young Stem"))
WGCNA2 <- WGCNA[,4:6]
colnames(WGCNA2)
WGCNA3 <- pivot_wider(WGCNA2,names_from=variable,values_from =mean )
View(WGCNA3)
rownames(WGCNAfiltered)
rownames(WGCNA3)
dim(WGCNA3)
WGCNA4 <- WGCNA3 %>%
   filter(!str_detect(Combined, "0.5"))

write.csv(WGCNA4,file="Traits_WGCNA2.csv",row.names=FALSE)


#checking row.names match
WGCNA5 <- read.csv("Traits_WGCNA2.csv")
WGCNA6 <- WGCNA5[-9,]
print(WGCNA6$Combined)
rownames(WGCNAfiltered)
rownames(WGCNA6) <- NULL 

WGCNA6 <- WGCNA6 %>% 
   tibble:: column_to_rownames(var = "Combined")

collectGarbage()
WGCNA6 <- WGCNA6[adjustedrownames,,drop=FALSE]


table(rownames(WGCNAfiltered) == rownames(WGCNA6))


#Curated gene list for later
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

# Check for duplicates and print the repeating entries
repeated_entries <- unique(ref_all$LOC[duplicated(ref_all$LOC)])
cat("Repeating entries: ", ifelse(length(repeated_entries) > 0, paste(repeated_entries, collapse = ", "), "None"), "\n")

#delete repeaters
ref_all <- ref_all[-c(122,111,186,228,229,328,335),]

# Check for duplicates and print the repeating entries
repeated_entries <- unique(ref_all$LOC[duplicated(ref_all$LOC)])
cat("Repeating entries: ", ifelse(length(repeated_entries) > 0, paste(repeated_entries, collapse = ", "), "None"), "\n")

write.csv(ref_all, "ref_all_removedrepeats.csv",row.names = FALSE)
view(ref_all)
#load in already adj32usted DF  
ref_all <- read.csv("ref_all_removedrepeats.csv")





#now all input data is establiushed, use output dir

dir.create("WGCNA") #not needed if rerunning 
setwd("C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/")
getwd()


# next we generate the dendrogram and relate it to traits, also outliers are identified
# sample network based on squared Euclidean distance note that we transpose the data
# this means you create distance for genotypes and NOT genes (that is how the table is organised)
# A = adjacency(t(datExprRoot), type = "distance")

A = adjacency(t(WGCNAfiltered), type = "signed hybrid")
k = as.numeric(apply(A, 2, sum)) - 1
# standardized connectivity, 'normalises' Root genotypes to each other 
Z.k = scale(k)
thresholdZ.k = -2.5  # often -2.5, vary to test, means 5 SD from other genotypes
# the color vector indicates outlyingness (red)
outlierColor = ifelse(Z.k < thresholdZ.k, "red", "black")
# calculate the cluster tree using flahsClust or hclust
sampleTree = hclust(as.dist(1 - A), method = "average")

# Convert traits to a color representation: where red indicates high values
traitColors = data.frame(numbers2colors(WGCNA6, signed = FALSE))
dimnames(traitColors)[[2]] = paste(names(WGCNA6), sep = "")
datColors = data.frame(outlierC = outlierColor, traitColors)
# Plot the sample dendrogram and the colors underneath. 
#The plot can indicate if a trait is simple, i.e. one genotype is mostly responsible for variation of one trait
# or if the trait is complex, a more even distribution aross genotypes
plotDendroAndColors(sampleTree, 
                    groupLabels = names(datColors),
                    colors = datColors,
                    cex.colorLabels = 0.8, 
                    cex.dendroLabels = 0.9, 
                    cex.rowText = 0.2,
                    main = "Sample dendrogram and trait heatmap")


dev.copy2pdf(file = "sample_trait_heatmap.pdf")

#now start with expression networks
# Choose a set of soft thresholding powers, important parameter for Root downstream analyses
powers = c(1:20)
sft = pickSoftThreshold(WGCNAfiltered,
                        powerVector = powers,
                        networkType = "signed hybrid",
                        corFnc = cor,
                        corOptions=list(use="p"),
                        #corOptions=list(maxPOutliers=0.1),
                        verbose = 5)
sft$fitIndices


sizeGrWindow(12, 9)
par(mfrow = c(1, 2))

# SFT index as a function of different powers
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)",
     ylab = "SFT, signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red")
# this line corresponds to using an R^2 cut-off of h
abline(h = 0.85, col = "red")

#create connectivity plot
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n",
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")
#pick a threshold that is at saturation, here sft = 10
dev.copy2pdf(file = "sft_optimisation.pdf")


##now we are ready to detect modules
#first we define most important parameters:
netType = "signed hybrid";
#netType = "signed";
#netType = "unsigned";
clustering = "pearson";
#clustering = "bicor";
pow = 13;#picked 10, but could go with 11 or 12 as well.
minModSize = 200;
maxBlokSize = 21000;
mergingThresh = 0.25;#need to ask Oli what to base this on
split = 4; #sesnitivity of module detection

getwd()
#calculate pairwise correlations between genes and creates a WGCN (Weighted gene co-expression network), creating modules assigned to numbers, MEs(module eigengenes), the latter of which are the first principal component of the expression data for all genes in a given module
#
net = blockwiseModules(WGCNAfiltered, 
                       corType = clustering,
                       maxPOutliers = 0.1, #for bicor
                       maxBlockSize = maxBlokSize,
                       networkType = netType,
                       power = pow,
                       minModuleSize = minModSize,
                       deepSplit = split,
                       mergeCutHeight = mergingThresh,
                       numericLabels = TRUE, #assigns modules into numbers rather than just colors
                       saveTOMs = TRUE,
                       pamRespectsDendro = FALSE,
                       saveTOMFileBase = "COG007_TOM",
                       verbose=5,
                       nthreads=4
      
)


plot(net$dendrograms[[1]], main = "Module Dendrogram", xlab = "", sub = "")
#load(file = "C:/Users/BWeeYang/OneDrive - LA TROBE UNIVERSITY/Desktop/WGCNA/COG007_TOM-block.1.RData") # this will load back in as "TOM"
#dimtom=as.matrix(TOM)
#dim(net)
#net2 <- TOM
#net2=as.matrix(TOM)
#dim(TOM_matrix)  # Get dimensions of the matrix
#net2[1:5, 1:5]
#ls(net2)
#show number of modules created
#able(net$colors)
#heatmap(net2, symm = TRUE, col = topo.colors(100)) #DO NOT PLOT THE HEATMAP OR R WILL CRASH LOL


#extract modules and calculate eigenvalues
net$colors

moduleLabelsAutomatic = net$colors
moduleColorsAutomatic = labels2colors(moduleLabelsAutomatic) #this assigns the colors according to the module number

#extract MEs
#The goal of calculating the module eigengene is to reduce the gene expression data from each module (which contains many genes) to a single value that represents the overall expression profile of that module for each sample.

MEsAutomatic = net$MEs #creates a dataframe where each column represents a module eigengene. 
dim(MEsAutomatic) #gives dimension of data frame (no. rows and no. columns respectively)
moduleNumber = dim(MEsAutomatic)[2] #
blocknumber = 1 

#defining parameters for traits, one by one
#GCNA_nitrate = as.data.frame(WGCNA6$nitrate)
#ames(WGCNA_nitrate) = "WGCNA_nitrate"
#GS.WGCNA_nitrate = as.numeric(cor(WGCNAfiltered, WGCNA_nitrate, use = "p"))
#GS.WGCNA_nitrate.Color = numbers2colors(GS.WGCNA_nitrate, signed = T)

#make list for the colors
traits <- colnames(WGCNA6)
trait_colors <- list()

for (trait in traits) {
   # Extract trait data
   trait_data <- as.data.frame(WGCNA6[[trait]])
   #colnames(trait_data) <- paste0("WGCNA_", trait) #only necessary if you need to have a genotype difference. but maybe doesn't make sense to input this here. why not at end of collected dataframe for colors or group labels?
   
   # Calculate correlation
   GS_trait <- as.numeric(cor(WGCNAfiltered, trait_data, use = "p"))
   
   # Generate color representation
   trait_colors[[trait]] <- numbers2colors(GS_trait, signed = TRUE)
}

# Combine all colors into a data frame

datColors <- data.frame(
   moduleColorsAutomatic,
   trait_colors
)[net$blockGenes[[blocknumber]], ] #uset this if no comparisons being made i.e. genotypes


#datColors2 <- data.frame(
   #moduleColorsAutomatic,
   #setNames(trait_colors, paste0("GS.COG07_", traits, ".Color"))
#)[net$blockGenes[[blocknumber]], ]

# Generate groupLabels dynamically
groupLabels <- c("Module colors", traits) #if no comparisons
#groupLabels2 <- c("Module colors", paste0("GS.COG007_", traits))


# Plot dendrogram with colors #possible the next step here is to look at subsetting specific comparisons only. But should this be done earlier?
plotDendroAndColors(net$dendrograms[[blocknumber]],
                    colors = datColors,
                    groupLabels = groupLabels,
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = FALSE,
                    main = paste("mergingThresh:", mergingThresh,
                                 ", cluster method:", clustering,
                                 ", Net type:", netType, "\n",
                                 "power: ", pow,
                                 ", min module size: ", minModSize,
                                 ", deepSplit: ", split,
                                 ", number of modules: ", moduleNumber),
                    guideHang = 0.05)

dev.copy2pdf(file = "modules_trait_heatmap_13modules.pdf")


moduleColorsWGCNA = moduleColorsAutomatic
nGenes = ncol(WGCNAfiltered)
nSamples = nrow(WGCNAfiltered)

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(WGCNAfiltered, moduleColorsWGCNA)$eigengenes
MEsWGCNA = orderMEs(MEs0) 
modTraitCor = cor(MEsWGCNA, WGCNA6, use = "p")
modTraitP = corPvalueStudent(modTraitCor, nSamples)

write.csv(modTraitCor, "module_trait_correlations.csv")
write.csv(modTraitP, "module_trait_pvalues.csv")

MEList = moduleEigengenes(WGCNAfiltered, colors = moduleColorsWGCNA)
MEs = MEList$eigengenes

traits <- colnames(WGCNA6)
trait_data_list <- list()  # Initialize an empty list to store the data frames


#make list for just the raw mean data
for (trait in traits) {
   # Extract trait data and store in the list with the trait name as the key
   trait_data_list[[trait]] <- as.data.frame(WGCNA6[[trait]])
}

#make list for p values
trait_correlations <- list()
for (trait in colnames(WGCNA6)) {
   # Extract numeric trait data
   trait_data <- as.numeric(WGCNA6[[trait]])
   
   # Calculate correlation
   GS_trait <- as.numeric(cor(WGCNAfiltered, trait_data, use = "p"))
   
   # Store correlation values directly in the list
   trait_correlations[[trait]] <- GS_trait
}
MET <- orderMEs(cbind(MEs, do.call(cbind, trait_data_list)))


MET=orderMEs(cbind(MEs))


plotEigengeneNetworks(MET,"",marDendro=c(0,4,1,2),
                      marHeatmap=c(3,4,1,2),
                      cex.lab=0.8,
                      xLabelsAngle=90,
                      signed = TRUE,
                      plotAdjacency = T)

dev.copy2pdf(file = "module_correlation.pdf")
#101224 stopping here

#Since we have a moderately large number of modules and traits, a 
# suitable graphical representation will help in reading the table. We 
# color code each association by the correlation value: Will display 
# correlations and their p-values

#display correlation and P value in same matrix
#textMatrix = paste(signif(modTraitCor, 2),"\n(", signif(modTraitP, 1), ")", sep = "")

#only show P values, colours already give correlation data
textMatrix = paste(signif(modTraitP, 1))

dim(textMatrix) = dim(modTraitCor)

#
# Create the textMatrix with both correlation values and p-values
textMatrix <- paste0(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")")

# Ensure the textMatrix has the same dimensions as the correlation matrix
dim(textMatrix) <- dim(modTraitCor)

# Convert to matrix with the same dimensions as modTraitCor2
textMatrix <- matrix(textMatrix, nrow = nrow(modTraitCor), ncol = ncol(modTraitCor))

# Display the correlation values in modTraitCor within a heatmap plot
sizeGrWindow(12, 15)
par(cex = 0.6);
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = modTraitCor, xLabels = names(WGCNA6),
               yLabels = names(MEsWGCNA),
               ySymbols = names(MEsWGCNA),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = TRUE,
               cex.text = 0.9,
               zlim = c(-1,1),
               main = paste("Module - trait relationships"))


dev.copy2pdf(file = "modules_trait_heatmap_clean_slighedit.pdf")

# calculate the module membership values (aka. module eigengene based # connectivity kME): 
datKME = signedKME(WGCNAfiltered, MEsWGCNA)

#Intramodular analysis: identifying genes with high GS and MM
colorOfColumn = substring(names(datKME), 4)

#Gene significance versus module eigengene plots
selectModules = colorOfColumn[1:14]


#plot verbosescatterplot
# Open the PNG device with desired dimensions
#pdf("GSVskME_sulfate2.pdf", width = 2400, height = 1800)  # Set width and height as needed
#par(mfrow = c(4, length(selectModules)/4))
sizeGrWindow(24, 20)
par(mfrow = c(4, 4)) 


#Creates scatter plots that visualize relationship between module membership(kME) and gene significance (GS), or between module membership and trait correlations
for (module in selectModules)
{
   column = match(module, colorOfColumn)
   restModule = moduleColorsWGCNA == module
   #do plot
   verboseScatterplot(datKME[restModule, column],
                      trait_correlations[["Organic_P"]][restModule], #using the correlations generated into trait_correaltions
                      xlab = paste("Module Membership ", module, "module"),
                      ylab = "GS.Organic P", #G.S=Gene significance
                      main = paste("kME.", module, "vs. GS","\n"), #"\n" adds a line break after the title finishes
                      abline = TRUE,
                      col = module,
                      cex = 0.55, #dot size,
                      cex.axis = 1.75, #axis labels
                      cex.lab = 1.75, #axis legend
                      cex.main = 1.75) #title text size
}
   
wgc
dev.copy2pdf(file = "GSVskME_traits_OrganicP.pdf", width = 26, height = 20)
print(traits)
dev.off()  # Close the graphics device after saving
#try other traits and other modules, how do you make R do this ?

dev.copy2pdf(file = "GSVskME_sulfate2.pdf",width=24,height=18)

#plot MEs for modules
column = length(colorOfColumn)
selectModules = sort(colorOfColumn[1:column])
sizeGrWindow(25, 15)
#par(mfrow = c(3, length(selectModules)/4))
par(mfrow = c(5, 3), mar=c(1,3,1,1))

moduleNumbers=as.matrix(table(moduleColorsWGCNA))

i=0
for (whichmodule in selectModules)
{ 
   ME=MEsWGCNA[, paste("ME",whichmodule,sep="")]
   i=i+1
   #set margins of grphics window
   #create a heatmap whose columns correspond to the 'arrays'
   # and whose rows correspond to genes
   barplot(ME, 
           col=whichmodule,
           main=paste(whichmodule,"(",moduleNumbers[i,1],")"),
           cex.main=1,
           ylab="eigengene expression",
           xlab=paste(whichmodule),
           ylim=c(-1,1),
           cex.axis = 1
   )
}

dev.copy2pdf(file = "MEs_WGCNA.pdf",height=15, width=25)

#merging of tables for module membership
annot <- read.csv("annot_genelevel.csv")
genes=names(WGCNAfiltered)
View(genes)
Gene_annot = match(genes, annot$LOC) #match returns numberical index (posiiton) of the first occurence of each element from genes vector in the comparison vector (annot$LOc)
View(Gene_annot)
datGS.Traits = data.frame(cor(WGCNAfiltered, WGCNA6, use = "p"))
names(datGS.Traits) = paste("cor", names(datGS.Traits), sep = ".")
datOutput = data.frame(
                       LOC = names(WGCNAfiltered),
                       name = annot[Gene_annot, "name"],
                       moduleColorsWGCNA,
                       datKME
)
test =data.frame( moduleColorsWGCNA,
                    datKME)
test$LOC <- rownames(datKME)



#save output
datOutputmerged=left_join(datOutput, ref_all, by="LOC")
write.csv(datOutputmerged,file="WGCNA_color_modules_output.csv",row.names=FALSE)


#use venn diagrams to visualize spread of genes over colors?








#plot MEs for modules
column = length(colorOfColumn)
selectModules = sort(colorOfColumn[1:column])
sizeGrWindow(11.7, 8.3)
#par(mfrow = c(3, length(selectModules)/4))
par(cex = 1);
par(mfrow = c(5, 3), mar=c(1,3,1,1))

moduleNumbers=as.matrix(table(moduleColorsWGCNA))

i=0
for (whichmodule in selectModules)
{ 
   ME=MEsWGCNA[, paste("ME",whichmodule,sep="")]
   i=i+1
   #set margins of grphics window
   #create a heatmap whose columns correspond to the 'arrays'
   # and whose rows correspond to genes
   barplot(ME, 
           col=whichmodule,
           main=paste(whichmodule,"(",moduleNumbers[i,1],")"),
           cex.main=1,
           ylab="eigengene expression",
           xlab=paste(whichmodule),
           ylim=c(-1,1),
           cex.axis = 0.5
   )
}

dev.copy2pdf(file = "MEs_WGCNA.pdf")
#test for normality
hist(WGCNA6$phosphate)
shapiro.test(WGCNA6$phosphate) # >0.05 is normal
qqnorm(WGCNA6$phosphate)
qqline(WGCNA6$phosphate, col = "red") #data follows line = roughly normal
WGCNA6_scaled <- data.frame(scale(WGCNA6))
hist(WGCNA6_scaled$phosphate)
shapiro.test(WGCNA6_scaled$phosphate) # >0.05 is normal
qqnorm(WGCNA6_scaled$phosphate)
qqline(WGCNA6_scaled$phosphate, col = "red")
#generate a cytoscape file, this will only work for small modules as files gets large, a lot of edges
#this generates an object with the two files for cytoscape, the edgefiel for connections and the node file for annotations
selected = c("black")

datExpr_selected = WGCNAfiltered[,which(moduleColorsWGCNA %in% selected)]

probes = names(WGCNAfiltered)
TOM_selected = TOMsimilarityFromExpr(WGCNAfiltered, 
                                     power = 18,
                                     networkType = "signed hybrid",
                                     TOMType="signed")

dimnames(TOM_selected) = list(probes, probes)
cyt = exportNetworkToCytoscape(TOM_selected,
                               edgeFile = paste("CytoscapeInput-edges-", paste(module, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(module, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.5, #the higher the fewer genes
                               nodeNames = probes,
                               altNodeNames = probes)

nodes = cyt$nodeData
edges = cyt$edgeData
write.csv(edges, "COG007_edges.csv")
write.csv(nodes, "COG007_nodes.csv")
par(mfrow = c(1,2),mar=c(4,4,3,2)+.3);
cex1 = 0.9


####GENE Ontology FOR EACH SPEARATE MODULE? or perhaps lets just do the lot.
#load in appropriate packagesLOC115694687
library("ggforce")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library(org.Csativa.eg.db) 
library("clusterProfiler")
library(GOSemSim)
library(enrichplot)
library(viridis)

#let's build the df required.
testDF <- read.csv("unsigned_v2/WGCNA_color_modules_output_unsignedv2.csv")
testDF <- read.csv("signed_h_lowersft/signed_h_lowersft/WGCNA_color_modules_output.csv")
testDF <- read.csv("signed_h_highsft_split25/WGCNA_color_modules_output_signedh_split25.csv")
testDF <- read.csv("unsigned_v3/WGCNA_color_modules_output_unsigned_v3.csv")
testDF <- read.csv("./unsigned_v4/WGCNA_color_modules_output_unsignedv4.csv")
testDF <- read.csv("./unsigned_v5/WGCNA_color_modules_output_unsignedv5.csv")
testDF <- read.csv("./unsigned_v6/WGCNA_color_modules_output_unsignedv6.csv")
testDF <- read.csv("./unsigned_v7/WGCNA_color_modules_output_unsignedv7.csv")
colnames(testDF)
testDF <- testDF[,c(1,3)] #get rid of anything besides gene names + associated module color
testDF <- testDF[,c(2,4)]
colnames(testDF)
testDF$moduleColorsWGCNA2 <- factor(testDF$moduleColorsWGCNA2)



#code for isolating specific gene clusters
testDF2 <- subset(testDF, moduleColorsWGCNA2== "red")


all_GO_comp <- compareCluster(geneCluster = gene_list,
                              
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



# now add those logs to the object to use

all_GO_comp@compareClusterResult$nlogq = logs$nlogq


dotplot(all_GO_comp,
        
        showCategory=20, # top20 for every module
        
        includeAll=FALSE, # this is important to keep the ones that are overlapping
        
        color = 'nlogq',
        
        font.size = 6,
        
        x = 'Cluster') +
   
   scale_color_viridis_c(option = "viridis")



ggsave('all_GO_modules_top20.pdf', width = 20, height = 25, unit = 'cm')




GOdata = godata(
   
   OrgDb = org.Csativa.eg.db,
   
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
   
   cutoff = 0.3, # lower number = fewer GOs
   
   by = "qvalue",
   
   select_fun = min,
   
   measure = "Rel",
   
   semData = GOdata)



dotplot(all_GO_comp_simple,
        
        showCategory=10, # top10 for every module
        
        includeAll=TRUE, # this is important to keep the ones that are overlapping
        
        color = 'nlogq',
        
        font.size = 8,
        
        x = 'Cluster') +
   
   #scale_color_viridis_c(option = "cividis")+
   
   #scale_color_viridis_c(option = "inferno")+
   
   scale_color_viridis_c(option = "viridis")+
   ggtitle("unsigned_v6")

ggsave('./unsigned_v6/unsigned_v6.pdf', width = 18, height = 25, unit = 'cm')

# Extract gene-to-GO mapping
gene_GO_mapping <- all_GO_comp_simple@compareClusterResult

# View the gene-to-GO mapping
View(gene_GO_mapping)
df_expanded <- gene_GO_mapping %>%
   separate_rows(geneID, sep = "/") %>%
   arrange(geneID)
df_expanded <-dplyr:: rename(df_expanded,"LOC"="geneID")
df_expanded_annot <- left_join(df_expanded,annot,by="LOC")
View(df_expanded_annot)

write.csv(df_expanded_annot,file="df_expanded_annot.csv",row.names = FALSE)

cnetplot(all_GO_comp_simple,
         showCategory = 10,  # Show top 10 categories
         circular = FALSE,   # Set to TRUE for a circular layout
         colorEdge = TRUE,
         max.overlaps=5)   # Edge color represents gene-category relationships
data(geneList)

library(org.Csativa.eg.db)

#learning to build your own annotation db for organism of choice
library(biomaRt)


#find plant marts from bioMart
listEnsemblGenomes()
#we want plants_mart, so creater vector for it and load from biomart
ensembl_plants <- useEnsemblGenomes(biomart = "plants_mart")

#search function
searchDatasets(ensembl_plants, pattern = "Cannabis")


# Connect to Ensembl
ensembl_cs10<- useEnsemblGenomes(biomart = "plants_mart",
                                 dataset = "csfemale_eg_gene") # Change species if needed
# View all available attributes
attributes <- listAttributes(ensembl_cs10)
View(attributes)  # Check for gene name-related attributes

# View all available filters
filters <- listFilters(ensembl_cs10)
View(filters)  # Check for gene ID-related filters
# Convert LOC IDs to gene names
genes <- c("LOC115713503", "LOC115709945")  # Replace with your list

result <- getBM(attributes = c("entrezgene_id", "external_gene_name"), 
                filters = "entrezgene_id", 
                values = genes, 
                mart = ensembl_cs10)
library(biomaRt)

# Get a small sample of gene names
test_result <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                     filters = "external_gene_name",
                     values="SPX2",
                     mart = ensembl_cs10 
                     )  # Limit to 10 results

print(test_result)

print(result)
keys <- keytypes(org.Csativa.eg.db)
bitr(genes, fromType="GID", toType=keys, OrgDb=org.Csativa.eg.db)
bitr(genes, fromType="ENTREZID", toType=c("ALIAS","EVIDENCE","GID","GO","GOALL","SYMBOL","ONTOLOGYALL","PMID","REFSEQ"), OrgDb=org.Cs.eg.db)
KEYS
library(clusterProfiler)
library(org.Csativa.eg.db)
print(result)
x <- c("GPX3",  "GLRX",   "LBP",   "CRYAB", "DEFB1", "HCLS1",   "SOD2",   "HSPA2",
       "ORM1",  "IGFBP1", "PTHLH", "GPC3",  "IGFBP3","TOB1",    "MITF",   "NDRG1",
       "NR1H4", "FGFR3",  "PVR",   "IL6",   "PTPRM", "ERBB2",   "NID2",   "LAMB1",
       "COMP",  "PLS3",   "MCAM",  "SPP1",  "LAMC1", "COL4A2",  "COL4A1", "MYOC",
       "ANXA4", "TFPI2",  "CST6",  "SLPI",  "TIMP2", "CPM",     "GGT1",   "NNMT",
       "MAL",   "EEF1A2", "HGD",   "TCN2",  "CDA",   "PCCA",    "CRYM",   "PDXK",
       "STC1",  "WARS",  "HMOX1", "FXYD2", "RBP4",   "SLC6A12", "KDELR3", "ITM2B")
eg = bitr(x, fromType="G", toType="ENTREZID", OrgDb="org.Hs.eg.db")
print(eg)
keytypes(org.Cs.eg.db)
keytypes(org.Csativa.eg.db)
library(AnnotationHub)
hub <- AnnotationHub()
query(hub,  c("cannabis","orgdb") )
yes
org.Cs.eg.db <- hub[["AH114845"]]
columns(org.Cs.eg.db)
k <- keys(org.Cs.eg.db)
length(k)

AnnotationDbi::saveDb(org.Cs.eg.db, file = "org.Cs.eg.sqlite")
org.Cs.eg.db <- AnnotationDbi::loadDb("org.Cs.eg.sqlite")

#Checking for Ren matches

library(ggstatsplot)
install.packages("(ggstatsplot)")


mlist <- read.csv("mastersheetWGCNA.csv")
rengenes <- read.csv("unsigned_v6/Ren Genes.csv")


View(mlist)
checkrenmodules <- left_join(rengenes, mlist[,c(1,2,3,23,18,20)],by="LOC")
View(checkrenmodules)
checkrenmodules2 <- checkrenmodules %>%
   mutate(Annotation=ifelse(Annotation !="notDE", "DE", Annotation)) # ifelse(condition, value_if_TRUE, value_if_FALSE) is an R function for conditional replacement.Annotation != "notDE": Checks if each value in Annotation is not "notDE".
checkrenmodules3 <- checkrenmodules2 %>%
   mutate(Annotation=ifelse(is.na(Annotation), "TPM<10", Annotation))
view(checkrenmodules3)
rednetwork <- read.csv("./unsigned_v6/Red Module/Redmodulegenes3.csv")
rednetwork <-dplyr:: rename(rednetwork,"LOC"="name")
colnames(rednetwork)
greenyellow <- read.csv("./unsigned_v6/Greenyellow Module/Greenmodulegenes4.csv")
colnames(greenyellow)
greenyellow <-dplyr:: rename(greenyellow,"LOC"="name")
checkrenmodules5 <- left_join(checkrenmodules3, rednetwork[,c(23,19,16,33,15)],by="LOC")
checkrenmodules6 <- left_join(checkrenmodules5, greenyellow[,c(19,14,12,23,10)],by="LOC")

write.csv(checkrenmodules6,file="unsigned_v6/Ren Genes_R_mod.csv",row.names=FALSE)
checkrenmodules3$Selection <- factor(checkrenmodules3$Selection)
checkrenmodules3$Annotation <- factor(checkrenmodules3$Annotation)
checkrenmodules3$moduleColorsWGCNA2 <- factor(checkrenmodules3$moduleColorsWGCNA2)
str(checkrenmodules3)
write.csv(checkrenmodules3, file="unsigned_v6/Ren genes ploting dataframe.csv",row.names = FALSE)
checkrenmodules3 <- read.csv("unsigned_v6/Ren genes ploting dataframe.csv")
checkrenmodules4 <- checkrenmodules3[,c(1,8,9)]
view(checkrenmodules3)

#Fisher's
ContingencyTable_DE <- subset(checkrenmodules4, Annotation == "DE"  )
test <- fisher.test(ContingencyTable_DE[,c(1:2)])
test <- fisher.test(ContingencyTable_DE)

checkrenmodules3
checkrenmodules4
#get counts
library(dplyr)
checkrenmodules3 <- read.csv("unsigned_v6/Ren Genes_R_mod.csv")
data_summary <- checkrenmodules3 %>%
   count(Selection, moduleColorsWGCNA2, Annotation) #count() is a shortcut for group_by() + summarise(n = n()) It groups the data by the specified columns (Selection, moduleColorsWGCNA2, and Annotation). Then, it counts how many times each combination occurs and creates a new column called n (which contains the counts).

view(data_summary)
str(data_summary)
write.csv(data_summary, file="unsigned_v6/Ren genes ploting dataframe.csv",row.names = FALSE)


data_summary_wide <- tidyr::pivot_wider(data_summary,names_from = moduleColorsWGCNA2, values_from = n, values_fill = 0)
View(data_summary_wide)
write.csv(data_summary_wide,file="data_summary_wide.csv",row.names = FALSE)
data_summary_wide <- read.csv("data_summary_wide.csv")

#now subsetting for DE or notDE 
subsetforplot <- subset(data_summary_wide, Annotation == "DE")
subsetforplot <- subsetforplot[,-2] %>%
   tibble:: column_to_rownames(var="Selection")


str(subsetforplot)

mosaicplot(subsetforplot,
           main = "Mosaic plot",
           color = TRUE
)

#chi+fisher's
chisq.test(subsetforplot)$expected
tes <- fisher.test(subsetforplot, simulate.p.value = TRUE, B = 1e6)
tes$p.value

#extract p-value to data frame to save if needed
df <- data.frame(p_value = tes[["p.value"]])

####
#now for whole counts i.e. not separating by selection criteria
checkrenmodules4

data_summary2 <- checkrenmodules3 %>%
   count( moduleColorsWGCNA2, Annotation) #count() is a shortcut for group_by() + summarise(n = n()) It groups the data by the specified columns (Selection, moduleColorsWGCNA2, and Annotation). Then, it counts how many times each combination occurs and creates a new column called n (which contains the counts).
write.csv(data_summary_wide2,file="countsforRengenesaccordingtoWGCNAmodules.csv",row.names=FALSE)
data_summary_wide2 <- tidyr::pivot_wider(data_summary2,names_from = moduleColorsWGCNA2, values_from = n, values_fill = 0)
data_summary_wide2



#now subsetting for DE or notDE 
subsetforplot2 <- data_summary_wide2[-3,]
subsetforplot2 <- tibble:: column_to_rownames(subsetforplot2,var="Annotation")
subsetforplot2

str(subsetforplot2)

mosaicplot(subsetforplot2,
           main = "Mosaic plot",
           color = TRUE
)

#chi+fisher's
chisq.test(subsetforplot2)$expected
tes <- fisher.test(subsetforplot2, simulate.p.value = TRUE, B = 1e6)
tes$p.value

#extract p-value to data frame to save if needed
df <- data.frame(p_value = tes[["p.value"]])



#ok now doing on a per-module basis for 
modulecounts <- c(901,3573,2365,1641,280,63,496,789,490,1473,4905,1689)
rencounts <- c(35,191,148,75,7,1,33,48,24,34,202,74)
modulecolors <- c("black",	"blue",	"brown",	"green",	"greenyellow",	"grey",	"magenta",	"pink",	"purple",	"red",	"turquoise",	"yellow",)
# Create the dataframe
df <- data.frame(
   Module = c("WGCNA", "yellow"),
   Hemp = c(18665, 1689),
   Ren = c(992, 71)
)
df <- tibble:: column_to_rownames(df,var="Module")

# Print the dataframe
print(df)

chisq.test(df)$expected
df_test <- fisher.test(df)
df_test$p.value




#statplots

View(checkrenmodules4)
subsetforplot2 <- subset(checkrenmodules3, Annotation == c("DE"))#for selection x module comparison
subsetforplot2 <- subset(checkrenmodules3, Annotation %in% c("DE","notDE"))#for DE/nonDE x module comparison
subsetforplot2 <- subsetforplot2[,-3] #for selection x module comparison
subsetforplot2 <- subsetforplot2[,-1] #for DE/nonDE x module comparison
colnames(subsetforplot2)
View(subsetforplot2)
library("ggstatsplot")
e=ggbarstats(
   subsetforplot2, Annotation, moduleColorsWGCNA2,
   results.subtitle = FALSE,
   subtitle = paste0(
      "DE and nonDEgenes together: Fisher's exact test", ", p-value = ",
      ifelse(tes$p.value < 0.001, "< 0.001", round(tes$p.value, 3))
   )
)
e
ggsave(e,file="./unsigned_v6/DEnonDE_Distribution_Rengenes_Modules.pdf",height=25,width=35,units="cm")


a= ggbarstats(
   subsetforplot2, Selection, moduleColorsWGCNA2,
   results.subtitle = FALSE,
   subtitle = paste0(
      "DE genes: Fisher's exact test", ", p-value = ",
      ifelse(tes$p.value < 0.001, "< 0.001", round(tes$p.value, 3))
   )
)

library(gridExtra)

c=grid.arrange(a, b,e, d, layout_matrix = rbind(c(1, 2,3), c(4, 4)))

ggsave(c,file="Rengenes_fishers_collated2.pdf", ,height=25,width=45,units="cm")
#normal barplots

ggplot(subset(checkrenmodules3, Selection == "hemp_v_basal" ), aes(x=moduleColorsWGCNA2,y=count("Annotation"), grouping=moduleColorsWGCNA2)+
       geom_bar()
       )

d=ggplot(data_summary, aes(x = moduleColorsWGCNA2, y = n, fill = Selection, grouping=Selection,)) +
   geom_bar(stat = "identity", position = "dodge") +
   theme_minimal() +
   labs(
      x = "Module Color",
      y = "Count",
      fill = "Annotation Category",
      title = "Annotation Counts per Module Color"
   ) +
   geom_text(aes(label = n), position = position_dodge(width = 0.9), vjust = -0.5, size = 3)+
   facet_wrap(~Annotation,ncol=3)+
   theme_bw()+
   theme(axis.title.x = element_text(color="black", size=10, face="bold"),
         axis.title.y = element_text(color="black", size=10, face="bold"),
         axis.text.x = element_text(size = 10, angle =45, vjust =0.5, face = "bold"),
         axis.text.y = element_text(size = 10, vjust =0.5, face = "bold"),
         strip.text.x = element_text(size = 10, face="bold", margin = margin(0.1,0,0.1,0, "cm")),
         #panel.grid.major.x = element_blank(),
         panel.grid.minor.y = element_blank(),
         panel.grid.major.x=element_blank()
   )

ggsave(file="Rengenes_barplot2.pdf", height=25,width=35,units="cm")
