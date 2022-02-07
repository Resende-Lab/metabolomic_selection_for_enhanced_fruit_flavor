# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#
# This script was adapted from Langfelder and Horvath
# https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/
#

###
# Metabolic Network Analysis
###

library(WGCNA)

###
# Tomato Data
###


species = "tomato"
type = "unsigned"
ntraits = 5
resdir = "results/fig1/"
  
set.seed(1234)
if (species == "tomato"){
  tom = read.csv("data/input/tom_imputed_scaled.csv",check.names = F)
} else {
  tom = read.csv("data/input/bb_imputed_scaled.csv",check.names = F)
}

rownames(tom) = paste(c(1:nrow(tom)), tom[,1])
tom = tom[, -(1:(ntraits+1))]

toclust = dist(scale(tom), method = "euclidean", diag = TRUE, upper = TRUE)
toclust = toclust/max(toclust)

sampleTree = hclust(toclust, method = "average");
clust = cutreeStatic(sampleTree, cutHeight = 0.7, minSize = 10)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples = (clust == 1)
tom = tom[keepSamples, ]

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE)

# Choose a set of soft-thresholding powers
powers = c(c(1:50))

# Call the network topology analysis function
sft = pickSoftThreshold(tom, RsquaredCut=0.6,powerVector = powers, verbose = 5,networkType=type)
sft$powerEstimate

# Plot the results:
pdf(paste0(resdir, sprintf("%s_%s_Power.pdf", species, type)), width = 9, height = 5);
sizeGrWindow(9, 5)
par(mfrow = c(1, 2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2",type = "n",
     main = paste("Scale independence"));
text(sft$fitIndices[, 1], -sign(sft$fitIndices[,3]) * sft$fitIndices[, 2],
     labels = powers,cex = cex1, col = "red");
abline(h = 0.7, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
dev.off()

# Constructing the gene network and identifying modules is now a simple function call:
net = blockwiseModules(tom, 
                       power = sft$powerEstimate, 
                       networkType = type,
                       TOMType = type, 
                       minModuleSize = 3, 
                       maxBlockSize = 181,
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.2,
                       numericLabels = TRUE, 
                       pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste0(resdir, sprintf("%s_TOM", species)),
                       verbose = 5)

# Save
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]


####
####
####

# Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(tom, power = sft$powerEstimate);

# Select modules
modules = unique(moduleColors);

# Select module volatiles
probes = colnames(tom)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];

# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste(resdir, species,"_",type,"CytoscapeInput-edges-", length(modules), "modules.txt", sep=""),
                               nodeFile = paste(resdir, species,"_",type,"CytoscapeInput-nodes-", length(modules), "modules.txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule])

