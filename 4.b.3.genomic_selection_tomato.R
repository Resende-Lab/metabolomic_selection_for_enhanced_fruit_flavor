args = commandArgs(trailingOnly=TRUE)

print(args[1])
print(args[2])
print(args[3])

experiment <- as.numeric(args[1])
trait <- as.numeric(args[2])
cv <- as.numeric(args[3])

#experiment <- 1
#trait = 2
#cv = 1


library(caret)
library(BGLR)
library(randomForest)
library(rrBLUP)

library(dplyr)
library(tidyr)
library(reshape2)


# Load
load("snps.Rdata")

tom = read.csv("./data/tom_imputed_scaled_subset_aggregated_fin.csv", header = T, check.names = F, na.strings = ".")

traits <- colnames(tom)[2:6]
mets <- colnames(tom)[7:ncol(tom)]

bayes <- c("BayesA", "BayesB", "BayesC", "BRR", "BL")
mods <- c(bayes, "RandomForest", "SupportVectorMachineLinear",  
          "ElasticNet", "NeuralNet", "KernelPartialLeastSquares",
          "LinearRegression", "RKHS", "GradientBoostMachine",
          "ReleventVectorMachine", "SVMRadial", "SVMRadialSigma",
          "BayesianNeuralNet", "Xgboost")

load("snps.Rdata")
snp.data <- snp.data[rownames(snp.data) != "TS-601",]
snp.data <- snp.data[-32,]

NMOD <- length(mods)
NPHEN <- length(traits)
NMET <- length(mets)
NFOLD <- 10

# Gather Initial Seeds
set.seed(0)
seeds <- sample(1:1e4, 100)

#df <- matrix(NA, ncol = 7, nrow = 100)

seed <- seeds[experiment]
set.seed(seed)

tom <- tom[sample(1:71, 70),]
idx <- sample(rep(1:NFOLD, 7))

A <- A.mat(snp.data[as.character(tom$sequencingID), ] )
snpX <- snp.data[as.character(tom$sequencingID), ] 
snpX[is.na(snpX)] <- -1

rownames(tom) <- tom$sequencingID

# Create Result Matrices
acc <- matrix(NA, nrow = 2, ncol = 1)

rm(list = c("snp.data", "phenotypic.data"))

############
system.time({
  set.seed(paste(seed, trait, cv, sep = ""))
  trnX <- as.matrix(tom[, mets])
  trnY <- data.frame(tom[, c(trait + 1, 1)])
  trnY[idx == cv, 1] <- NA
  tstY <- as.matrix(tom[idx == cv, trait + 1])
  
  ans1 <- kin.blup(trnY, K = A, geno = "sequencingID", pheno = traits[trait])
  acc[1,1] <- cor(ans1$g[idx == cv], tstY)

  #Caret
  carX <- as.matrix(cbind(Y = tom[idx != cv, trait + 1], snpX[idx != cv, ]))
  tstX <- snpX[idx == cv, ]
  
  #Random Forest
  trnCtrl <- trainControl(method = "cv", number = 5)
  grid <- expand.grid(.mtry = c(2000, 5000, 10000))
  rf <- train(x = carX[,-1], y = (carX[,1]), method = 'rf', ntree = 1000, tuneGrid = grid, trControl = trnCtrl)
  
  res <- as.vector(predict(rf, tstX))
  acc[2,1] <- cor(res, tstY)
  

  })  
##########

results <- data.frame(experiment = experiment,
                      trait = traits[trait],
                      cv = cv,
                      model = c("gBLUP", "randomForest"), 
                      accuracy = acc[,1]
                      )


if(!file.exists("rfGenomic.txt")){
  write.table(results, "rfGenomic.txt", row.names = F, col.names = T, append = T, quote = F)
} else {
  write.table(results, "rfGenomic.txt", row.names = F, col.names = F, append = T, quote = F)
}


















