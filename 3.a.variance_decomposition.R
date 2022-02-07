# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

###
# Figure 3 Multivariate Analysis
###

# For ease of use, this script uses some non-free software (asreml) for a portion of the analysis. 
#   This does hinder reproducibility to a degree as many will not have a license. If this is the case for you, 
#   please try swapping out the asreml function below for the open source lme4 package and see if you get the same results
#   

# plots
library(ggplot2)
library(ggrepel)
library(corrplot)
library(gridExtra)
library(GGally)

# Manipulation
library(dplyr)
library(reshape2)
library(openxlsx)
library(readr)
library(readxl)
library(tidyverse)

# Multivariate Analysis
library(factoextra)
library(FactoMineR)
library(survminer)

# Linear Models
library(lme4)
library(lmerTest)
library(asreml)
library(AGHmatrix)


# Set results directory and seed
resdir <- "results/fig3/"
set.seed(100)

# Load data
key = read.csv("data/input/tom_metabolites_clusters_key.csv")
scaled.tom = read.csv("data/input/tom_imputed_scaled.csv", check.names = F)

# Make unique genotype ids
scaled.tom$id <- paste(1:nrow(scaled.tom), scaled.tom$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "umami", "intensity")
mets = colnames(scaled.tom)[!colnames(scaled.tom) %in% c(sensory, "id")]

# Creating a metabolite key by merging
idx.tom = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")


# non.vocs vs. VOCs
non.vocs = idx.tom %>%
  filter(GBLUP_cluster == "Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()

vocs = idx.tom %>%
  filter(GBLUP_cluster != "Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()

# Non.vocs matrix
W = scaled.tom %>%
  select(non.vocs$Metabolite) %>%
  as.matrix(.)
W_non.vocs = W%*%t(W) 

G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_non.vol=1-G/max(G, na.rm=TRUE)

# Vocs matrix
W1 = scaled.tom %>%
  select(vocs$Metabolite) %>%
  as.matrix(.)
W_vocs = W1%*%t(W1) 

G1 = as.matrix(dist(W1, method="euclidean", diag=TRUE,upper=TRUE))
G_vol=1-G1/max(G1, na.rm=TRUE)

# ASREML
for(i in 1:length(sensory)){
  #  i=1
  sensory.idx = sensory[i] 
  print(sensory.idx)
  df = data.frame(y = scaled.tom[,sensory.idx],
                  genotype = paste(scaled.tom$id,1:nrow(scaled.tom),sep="_"))
  
  df$genotype <- as.factor(df$genotype)
  
  # Names for the matrices
  rownames(G_vol) = df$genotype;colnames(G_vol) = df$genotype
  rownames(G_non.vol) = df$genotype; colnames(G_non.vol) = df$genotype
  
  library(asreml)
  fm= asreml(fixed=y~1,
             random = ~ vm(genotype, G_vol) + vm(genotype, G_non.vol), 
             na.action = na.method(y="include"), 
             #residual = ~ dsum(~units|year),
             data=df)
  
  total = sum(summary(fm)$varcomp[,1])
  cat("\n")
  print(sensory.idx)
  print(paste("VOCs --->",summary(fm)$varcomp[1,1]/total))
  print(paste("nonVOCs --->",summary(fm)$varcomp[2,1]/total))
  print(paste("Resi --->",summary(fm)$varcomp[3,1]/total))
  cat("\n\n")
}



# Per group
## Tomato

resdir <- "results/fig3/"
set.seed(100)

# Load data
key = read.csv("data/input/tom_metabolites_clusters_key.csv")
scaled.tom = read.csv("data/input/tom_imputed_scaled.csv", check.names = F)

# Make unique genotype ids
scaled.tom$id <- paste(1:nrow(scaled.tom), scaled.tom$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "umami", "intensity")
mets = colnames(scaled.tom)[!colnames(scaled.tom) %in% c(sensory, "id")]

# Creating a metabolite key by merging
idx.tom = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")


# Sugar
a1 = idx.tom %>%
  filter(GBLUP_cluster=="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite[c(2,3,6)]) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_sugar=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# Acids
a1 = idx.tom %>%
  filter(GBLUP_cluster=="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite[c(1,4,5)]) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_acids=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# AA derived
a1 = idx.tom %>%
  filter(GBLUP_cluster=="AA derived") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_AA=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# Carotenoid
a1 = idx.tom %>%
  filter(GBLUP_cluster=="Carotenoid") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_caroteinoid=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# Lipid
a1 = idx.tom %>%
  filter(GBLUP_cluster=="Lipid") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_lipid=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# Phe derived
a1 = idx.tom %>%
  filter(GBLUP_cluster=="Phe derived") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.tom %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_phe=1-G/max(G, na.rm=TRUE) + diag(0.001,209,209)

# ASREML
df.list.tom = list()
for(i in 1:length(sensory)){
  #  i=7
  sensory.idx = sensory[i] 
  print(sensory.idx)
  df = data.frame(y = scaled.tom[,sensory.idx],
                  genotype = paste(scaled.tom$id,1:nrow(scaled.tom),sep="_"))
  df$genotype <- as.factor(df$genotype)
  # Names for the matrices
  rownames(G_AA) = df$genotype;colnames(G_AA) = df$genotype
  rownames(G_sugar) = df$genotype; colnames(G_sugar) = df$genotype
  rownames(G_acids) = df$genotype; colnames(G_acids) = df$genotype
  rownames(G_caroteinoid) = df$genotype;colnames(G_caroteinoid) = df$genotype
  rownames(G_phe) = df$genotype; colnames(G_phe) = df$genotype
  rownames(G_lipid) = df$genotype; colnames(G_lipid) = df$genotype
  
  library(asreml)
  fm= asreml(fixed=y~1,
             random=~ vm(genotype,G_sugar) + vm(genotype,G_acids) + vm(genotype,G_AA) + vm(genotype,G_phe) + vm(genotype,G_caroteinoid) + vm(genotype,G_lipid), 
             na.action = na.method(y="include"), 
             #residual = ~ dsum(~units|year),
             data=df)
  fm = update(fm)
  
  total = sum(summary(fm)$varcomp[,1])
  
  cat("\n")
  print(sensory.idx)
  print(paste("sugar --->",summary(fm)$varcomp[1,1]/total))
  print(paste("acids --->",summary(fm)$varcomp[2,1]/total))
  print(paste("AA --->",summary(fm)$varcomp[3,1]/total))
  print(paste("phe --->",summary(fm)$varcomp[4,1]/total))
  print(paste("carot --->",summary(fm)$varcomp[5,1]/total))
  print(paste("lipid --->",summary(fm)$varcomp[6,1]/total))
  print(paste("Resi --->",summary(fm)$varcomp[7,1]/total))
  cat("\n\n")
  
  df.list.tom[[i]] = data.frame(sensory =  c("sugar", "acids", "AA", "phe", "carot",
                                             "lipid", "residual"),
                                trait = sensory.idx, species="tomato",
                                var.comp = (summary(fm)$varcomp[,1]/total)*100)
  
}




#############
# Blueberry #
#############

# Set results directory and seed
resdir <- "results/fig3/"
set.seed(100)


# Load data
key = read.csv("data/input/bb_metabolites_clusters_key.csv")
scaled.bb = read.csv("data/input/bb_imputed_scaled.csv", check.names = F)

# Make unique genotype ids
scaled.bb$id <- paste0("bb", 1:nrow(scaled.bb), "_", scaled.bb$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "intensity")
mets = colnames(scaled.bb)[!colnames(scaled.bb) %in% c(sensory, "id")]


# Creating a metabolite key by merging
idx.bb = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")


# non.vocs vs. VOCs
non.vocs = idx.bb %>%
  filter(GBLUP_cluster=="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()

vocs = idx.bb %>%
  filter(GBLUP_cluster!="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()

# Non.vocs matrix
W = scaled.bb %>%
  select(non.vocs$Metabolite) %>%
  as.matrix(.)

#W_non.vocs = W%*%t(W) 
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_non.vol=1-G/max(G, na.rm=TRUE)

# Vocs matrix
W1 = scaled.bb %>%
  select(vocs$Metabolite) %>%
  as.matrix(.)
#W_vocs = W1%*%t(W1) 
G1 = as.matrix(dist(W1, method="euclidean", diag=TRUE,upper=TRUE))
G_vol=1-G1/max(G1, na.rm=TRUE)

# ASREML
for(i in 1:length(sensory)){
  sensory.idx = sensory[i] 
  print(sensory.idx)
  df = data.frame(y = scaled.bb[,sensory.idx],
                  genotype = scaled.bb$id)

  df$genotype <- as.factor(df$genotype)
  
  # Names for the matrices
  rownames(G_vol) = df$genotype;colnames(G_vol) = df$genotype
  rownames(G_non.vol) = df$genotype; colnames(G_non.vol) = df$genotype
  
  library(asreml)
  fm = asreml(fixed=y~1,
             random=~ vm(genotype, G_vol) + vm(genotype, G_non.vol), 
             na.action = na.method(y = "include"), 
             #residual = ~ dsum(~units|year),
             data = df)
  
  total = sum(summary(fm)$varcomp[,1])
  
  print(sensory.idx)
  cat("\n")
  print(paste("VOCs --->",summary(fm)$varcomp[1,1]/total))
  print(paste("nonVOCs --->",summary(fm)$varcomp[2,1]/total))
  print(paste("Resi --->",summary(fm)$varcomp[3,1]/total))
  cat("\n\n")
}



## Blueberry

# Load data
key = read.csv("data/input/bb_metabolites_clusters_key.csv")
scaled.bb = read.csv("data/input/bb_imputed_scaled.csv", check.names = F)

# Make unique genotype ids
scaled.bb$id <- paste0("bb", 1:nrow(scaled.bb), "_", scaled.bb$id)

# Gathering names of traits and metabolites
sensory = c("liking","sweetness", "sour", "intensity")
mets = colnames(scaled.bb)[!colnames(scaled.bb) %in% c(sensory, "id")]


# Creating a metabolite key by merging
idx.bb = data.frame(Metabolite = mets) %>%
  merge(., key, by="Metabolite")


table(idx.bb$GBLUP_cluster)

# Sugar
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite[c(2,3,4,5)]) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_sugar=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# Acids
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Acid/Sugar") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite[c(1,6)]) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_acids=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# AA derived
a1 = idx.bb %>%
  filter(GBLUP_cluster=="AA derived") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_AA=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# Carotenoid/Terpenoid
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Carotenoid/Terpenoid") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_caroteinoid=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# Lipid
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Lipid") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_lipid=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# Phe derived
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Phe derived") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_phe=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)

# Ester
a1 = idx.bb %>%
  filter(GBLUP_cluster=="Ester") %>%
  select(Metabolite) %>%
  droplevels()
W = scaled.bb %>%
  select(a1$Metabolite) %>%
  as.matrix(.)
G = as.matrix(dist(W, method="euclidean", diag=TRUE,upper=TRUE))
G_ester=1-G/max(G, na.rm=TRUE) + diag(0.001,244,244)


# ASREML
df.list = list()
for(i in 1:length(sensory)){
  
  sensory.idx = sensory[i] 
  print(sensory.idx)
  df = data.frame(y = scaled.bb[,sensory.idx],
                  genotype = scaled.bb$id)
  df$genotype <- as.factor(df$genotype)
  
  # Names for the matrices
  rownames(G_AA) = df$genotype;colnames(G_AA) = df$genotype
  rownames(G_sugar) = df$genotype; colnames(G_sugar) = df$genotype
  rownames(G_acids) = df$genotype; colnames(G_acids) = df$genotype
  rownames(G_caroteinoid) = df$genotype;colnames(G_caroteinoid) = df$genotype
  rownames(G_phe) = df$genotype; colnames(G_phe) = df$genotype
  rownames(G_lipid) = df$genotype; colnames(G_lipid) = df$genotype
  rownames(G_ester) = df$genotype; colnames(G_ester) = df$genotype
  
  
  fm= asreml(fixed=y~1,
             random=~ vm(genotype,G_sugar) + vm(genotype,G_acids) + vm(genotype,G_AA) + vm(genotype,G_phe) + vm(genotype,G_caroteinoid) + vm(genotype,G_lipid) + vm(genotype,G_ester), 
             na.action = na.method(y="include"), 
             #residual = ~ dsum(~units|year),
             data=df)
  
  fm = update(fm)
  
  total = sum(summary(fm)$varcomp[,1])

  print(sensory.idx)
  cat("\n")
  print(paste("sugar --->",summary(fm)$varcomp[1,1]/total))
  print(paste("acids --->",summary(fm)$varcomp[2,1]/total))
  print(paste("AA --->",summary(fm)$varcomp[3,1]/total))
  print(paste("phe --->",summary(fm)$varcomp[4,1]/total))
  print(paste("carot --->",summary(fm)$varcomp[5,1]/total))
  print(paste("lipid --->",summary(fm)$varcomp[6,1]/total))
  print(paste("Ester --->",summary(fm)$varcomp[7,1]/total))
  print(paste("Resi --->",summary(fm)$varcomp[8,1]/total))
  cat("\n\n")
  
  
  df.list[[i]] = data.frame(sensory =  c("sugar", "acids", "AA", "phe", "carot",
                                         "lipid", "ester", "residual"),
                            trait = sensory.idx, species="blueberry",
                            var.comp = (summary(fm)$varcomp[,1]/total)*100)
  
}




a = do.call("rbind", df.list)
b = do.call("rbind", df.list.tom)
var.comp = rbind(a,b)

save(var.comp, file = paste0(resdir, "GBLUP_seed100.Rdata"))
write.csv(var.comp, file = paste0(resdir, "VarComp_seed100.csv"), row.names = F)

# Make a rough plot to verify results - we will plot results out nicely in a separate script
ggplot(data=var.comp, aes(x=trait, y=var.comp, fill=sensory)) +
  geom_bar(stat="identity", position=position_dodge())+
  facet_wrap(~species,scales = "free")+
  theme_bw() 

# Check this run against the published results
this_analysis <- read.csv(paste0(resdir, "VarComp_seed100.csv"))
published_analysis <- read.csv("data/supplemental_datasets/SD3_VarComp.csv")

summary(this_analysis$var.comp - published_analysis$var.comp)





























