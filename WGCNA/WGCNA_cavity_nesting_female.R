if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)

#BiocManager::install("impute")
library(impute)

#BiocManager::install("GO.db")
library(GO.db)

#BiocManager::install("preprocessCore")
library(preprocessCore)

#install.packages("WGCNA")
library("WGCNA")

#install.packages("readr")
library(readr)

#install.packages("ggforce")
library(ggforce)

#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

#install.packages("dplyr")
library(dplyr)

if (!("magrittr" %in% installed.packages())) {
  install.packages("magrittr")
}
library(magrittr)

if (!("limma" %in% installed.packages())) {
  BiocManager::install("limma")
}
library(limma)

if (!("ggplot2" %in% installed.packages())) {
  install.packages("ggplot2")
}
library(ggplot2)

#setwd("~/Downloads/WGCNA - new cavity nesting project/")


setwd("~/Dropbox/Rosvall_Postdoc/Cavity Nesting/Cavity transcriptomics/WGCNA/WGCNA CN Signed Hybrid/Females")

############################################################################
#   Read in count data
############################################################################
countdata<-as.matrix(read.csv("CN_raw_counts_females.csv"))

# Save protein ID, gene name and Gene description in separate variable
geneInfo <- countdata[,1:3]

# remove first three columns so that our matrix is only numbers
ncol(countdata)
countdata <- countdata[,4:64]
head(countdata)
rnames <- geneInfo[,2]
cnames <- colnames(countdata)

#reformat the matrix so that it reads the counts as numeric
rawCounts<-matrix(as.numeric(unlist(countdata)),nrow=nrow(countdata))
rownames(rawCounts) <- rnames
colnames(rawCounts) <- cnames


############################################################################
#   Read in meta data
############################################################################

# Read in the meta data
metadata <- read.csv(file = "CN_WGCNA_traits_females.csv", header=T)
all.equal(cnames, metadata[,1])
rownames(metadata) <- metadata[,1]


############################################################################
#   QC outlier detection
############################################################################

gsg <- goodSamplesGenes(t(rawCounts))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# remove genes that are detectd as outliers (a total of 625 genes)
data1 <- rawCounts[gsg$goodGenes == TRUE,]
#there were no bad genes or bad samples!

# detect outlier samples - hierarchical clustering - method 1
htree1 <- hclust(dist(t(data1)), method = "average")
plot(htree1)

#Think asbout removing ETF6 and HSF3

# pca - method 2
pca <- prcomp(t(data1))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

#After looking at the PCA I would probably remove HSF3

### NOTE: If there are batch effects observed, correct for them before moving ahead

# exclude outlier samples
samples.to.be.excluded <- c('ETF6', 'HSF3')
data.subset <- data1[,!(colnames(data1) %in% samples.to.be.excluded)]

#Also need to subset the metadata so everything matches up
sub.metadata <- metadata[!(rownames(metadata) %in% samples.to.be.excluded),]

############################################################################
#   Create DESeq data object and generate normalized counts
############################################################################

# Now we have to take our raw counts and normalize them
dds_datExpr <- DESeqDataSetFromMatrix(countData = data.subset, colData = sub.metadata, design =~ Nest.Type)

#remove genes with low counts
#Since there are 12 samples per species we will do 90% of 12
dds90 <- dds_datExpr[rowSums(counts(dds_datExpr) >= 15) >= 11,]
nrow(dds90) # 10609 genes

#Run DESeq
des_datExpr <- DESeq(dds90)

vsd<-varianceStabilizingTransformation(des_datExpr)

# Transpose data
normalizedCounts <- assay(vsd) %>%
  t()

############################################################################
#   Run soft power threshold and blockwise modules
############################################################################

# To identify which genes are co-expressed similarly across samples, WGCNA first creates a weighted network to define which genes are most similar in expression.
# The measure of similarity it uses is based on a correlation matrix, but it requires the definition of a threshold value, which in turn depends on a "power" parameter that defines the exponent used when transforming the correlation values. 
# The choice of power parameter will affect the number of modules identified. 
# We can use this pickSoftThreshold function to identify a good choice for our power parameter. 
# This code might take a little while to run (10 minutes or so)
sft <- pickSoftThreshold(normalizedCounts, dataIsExpr = T, networkType = "signed hybrid")

#?pickSoftThreshold()

# Calculate the signed R^2, a measure of the model fit. We will plot this value to figure out what our "power" threshold should be. 
sft_df <- data.frame(sft$fitIndices) %>%
  dplyr::mutate(model_fit = -sign(slope) * SFT.R.sq)

# Now we'll plot the "model fit" (i.e., the signed R^2) vs. the power threshold to decide which power threshold gives us the best model fit. 
ggplot(sft_df, aes(x = Power, y = model_fit, label = Power)) + # This line passes our data to the program ggplot
  geom_point() + geom_text(nudge_y = 0.1) + # This tells it that we want it to be a dotplot, and that we want each point to have a label right next to it
  geom_hline(yintercept = 0.80, col = "red") + # This plots a horizontal red line on the same axes, according to the formula y = 0.8 (which is the signed R^2 threshold that the program recommends)
  ylim(c(min(sft_df$model_fit), 1.05)) + # This sets the bounds of our y axis
  xlab("Soft Threshold (power)") +  # These next three lines set the labels for our axes, and makes the formatting of the graph nicer
  ylab("Scale Free Topology Model Fit, signed R^2") + 
  ggtitle("Scale Independence") + theme_classic()

sft$powerEstimate

bwnet8 <- blockwiseModules(normalizedCounts, # This line passes our "normalized_counts" dataset that we created on lines 191-192 to the blockwiseModules function
                           maxBlockSize = 5000, # We will split the data into chunks of 5000 genes because our computer can't do all 20,000+ at once
                           TOMType = "signed", # Type of matrix we have
                           power = 8, # I chose a power parameter of 14, based on the plot we generated in the previous step. Do you understand why?
                           numericLabels = T, # Modules will be numbered
                           randomSeed = 1234) # Since this process is a random search, we won't necessarily get the same results every time. If we set a seed though, these results will be reproducible. 



############################################################################
#   Module Eigengenes
############################################################################


module_eigengenes8 <- bwnet8$MEs

metadata$Nest.Type <- relevel(as.factor(metadata$Nest.Type), ref= "Open")

adjacency8 = adjacency(normalizedCounts, power = 8, type = "unsigned") #Calculating the adjacency matrix

TOM8 = TOMsimilarity(adjacency8) 

dissTOM8 = 1-TOM8 ##Calculating the dissimilarity

# Call the hierarchical clustering function
geneTree8 = hclust(as.dist(dissTOM8), method = "average")

plot(geneTree8, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)


minModuleSize = 30;
# Module identification using dynamic tree cut:
dynamicMods8 = cutreeDynamic(dendro = geneTree8, distM = dissTOM8,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize);
table(dynamicMods8)

dynamicColors8 = labels2colors(dynamicMods8)

# Plot the dendrogram and colors underneath
sizeGrWindow(8,6)
plotDendroAndColors(geneTree8, dynamicColors8, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")


#Clustering the modules
#Calculating eigengenes
MEList8 = moduleEigengenes(normalizedCounts, colors = dynamicColors8)
MEs8 = MEList8$eigengenes

#Calculating the module dissimilarity eigengenes
MEDiss8 = 1-cor(MEs8)

#Clustering the eigengenes modules
METree8 = hclust(as.dist(MEDiss8), method = "average")

#Plotting the result
plot(METree8, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
#Grouping the clusters from a cut-off
MEDissThres = 0.25
#Plotting a cut-off line
abline(h=MEDissThres, col = "red")


#Merge modules that are siilar to one another
merge8 = mergeCloseModules(normalizedCounts, dynamicColors8, cutHeight = MEDissThres, verbose = 3)

# The merged module colors
mergedColors8 = merge8$colors

# Eigengenes of the new merged modules:
mergedMEs8 = merge8$newMEs

write.csv(mergedMEs8, file="eigengenes_female.csv")

plotDendroAndColors(geneTree8, cbind(dynamicColors8, mergedColors8),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Rename to moduleColors
moduleColors8 = mergedColors8

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels8 = match(moduleColors8, colorOrder)-1;

MEs8 = mergedMEs8

# Save module colors and labels for use in subsequent parts
save(MEs8, moduleLabels8, moduleColors8, geneTree8, file = "Cavity_Module_Data_Female.Rdata")


############################################################################
#   Module trait correlations
############################################################################

#Relating modules to characteristics and identifying important genes
#Defining the number of genes and samples
nGenes1 = ncol(normalizedCounts)
nSamples1 = nrow(normalizedCounts)

#get module eigengenes of new merged modules
MEs0.8 = moduleEigengenes(normalizedCounts, moduleColors8)$eigengenes
MEs.8 = orderMEs(MEs0.8)

#code check
dim(MEs.8)
dim(sub.metadata)
rownames(MEs.8)
rownames(metadata) <- metadata$Individual

# Change all of the metadata data traits into numeric vectors and save as new df
metadata2 <- data.frame(sub.metadata[,8:9])
metadata2$Species.Common <- as.numeric(as.factor(sub.metadata[,2]))
metadata2$Species.Scientific <- as.numeric(as.factor(sub.metadata[,3]))
metadata2$Nest.Type.Specific <- as.numeric(as.factor(sub.metadata[,4]))
metadata2$Nest.Type <- as.numeric(as.factor(sub.metadata[,5]))
#metadata2$Sex <- as.numeric(as.factor(metadata[,6]))
metadata2$Family <- as.numeric(as.factor(sub.metadata[,7]))

# correlate modules with traits and calculate pvalues
moduleTraitCor8 = bicor(MEs.8, metadata2, use = "p",robustY = FALSE,maxPOutliers =0.1);
moduleTraitPvalue8 = corPvalueStudent(moduleTraitCor8, nSamples1)

#Displaying correlations and its p-values
textMatrix8 =  paste(signif(moduleTraitCor8, 2), "\n(",
                     signif(moduleTraitPvalue8, 1), ")", sep = "")

dim(textMatrix8) = dim(moduleTraitCor8)

par(mar = c(6, 8.5, 3, 3))
par(mar = c(1, 3, 2, 2))

dev.off()
#Displaying the correlation values in a heatmap plot

labeledHeatmap(Matrix = moduleTraitCor8,
               xLabels = names(metadata2),
               yLabels = names(MEs.8),
               ySymbols = names(MEs.8),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix8,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()


geneModuleMembership = as.data.frame(cor(normalizedCounts, MEs.8, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples1))

geneInfo = data.frame(Gene.name = colnames(normalizedCounts),
                      moduleColor = moduleColors8, geneModuleMembership)

write.csv(geneInfo, file="geneInfo.female.Pow8.csv")



############################################################################
#   Exporting networkt to cytoscape
############################################################################

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(normalizedCounts, power = 8); 
# Is this 'normalizedCounts' the equivalent of datExpr? Need to add 

# Read in annotation file
annot = read.csv("CN_annGeneNames.csv")

# Select modules
modules = c("sienna4")

#Select module probes 
probes=colnames(normalizedCounts) 
inModule=is.finite(match(moduleColors8,modules))
modProbes=probes[inModule]
modGenes=annot$GeneName[match(modProbes,annot$ZebraFinchProteinID)]

#Select the corresponding Topological Overlap 
modTOM=TOM[inModule,inModule]
dimnames(modTOM)=list(modProbes,modProbes)

# Export the netLwork into edge and node list files Cytoscape can read

cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors8[inModule]);

chooseTopHubInEachModule(normalizedCounts,colorh = "sienna4") 
