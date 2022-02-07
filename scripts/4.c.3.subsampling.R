args = commandArgs(trailingOnly=TRUE)

# Permutations
print(args[1])
print(args[2])
print(args[3])

experiment <- as.numeric(args[1])
trait <- as.numeric(args[2])
subsamp <- as.numeric(args[3])

#experiment <- 1
#trait = 2
#subsamp = 100

library(dplyr)
library(reshape2)
library(stringr)
library(foreach)
library(doParallel)

library(caret)
library(BGLR)
library(elasticnet)
library(glmnet)
library(nnet)
library(neuralnet)
library(pls)
library(randomForest)
library(kernlab)
library(brnn)
library(gbm)
library(xgboost)
library(e1071)
library(rrBLUP)

############
# Analysis #
############

#bb = read.csv("./data/bb_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")
df = read.csv("./data/input/tom_imputed_scaled_fin.csv", header = T, check.names = F, na.strings = ".")

traits <- colnames(df)[2:6]
mets <- colnames(df)[7:ncol(df)]

#trait <- traits[trait]

bayes <- c("BayesA", "BayesB", "BayesC", "BRR", "BL")
models <- c(bayes,
            "RandomForest",
            "SupportVectorMachineLinear",
            "ElasticNet",
            "NeuralNet",
            "KernelPartialLeastSquares",
            "LinearRegression",
            "RKHS",
            "GradientBoostMachine",
            "ReleventVectorMachine",
            "SVMRadial",
            "SVMRadialSigma",
            "BayesianNeuralNet",
            "gBLUP")

NPHEN <- length(traits)
#NFTRS <- length(mets)
NMOD <- 18
#NMET <- 68
NFOLD <- 10
NTEST <- 39

set.seed(0)
seeds <- sample(1:1e4, 10)


#############
# Functions #
#############

  set.seed(seeds[experiment])
  NSAMP <- nrow(df)

  idx <- sample(1:NSAMP, NTEST)
  trn <- df[-idx,]
  tst <- df[idx,]

  idx <- sample(rep(1:17, 10))
  acc <- matrix(NA, nrow = NMOD, ncol = 1)
  #  wts <- matrix(NA, nrow = NMOD, ncol = NMET)

    trnX <- as.matrix(trn[1:subsamp, mets])
    trnY <- as.matrix(trn[1:subsamp, trait + 1])

    tstX <- as.matrix(tst[, mets])
    tstY <- as.matrix(tst[, trait + 1])

    X.bglr <- rbind(trnX, tstX)
    Y.bglr <- rbind(trnY, as.matrix(rep(NA, nrow(tstY))))

    X.caret <- as.matrix(data.frame(Y = trnY, trnX, check.names = F))

    A <- A.mat(X.bglr)
    Y.gblup <- data.frame(id = row.names(trn[1:subsamp,]), pheno = trnY)

    #BAYES
    for(method in 1:5){
      ETA <- list(X = list(X = X.bglr, model = bayes[method]))
      fit <- BGLR(y = Y.bglr, ETA = ETA,
                  nIter = 30000,
                  burnIn = 10000,
                  thin = 100,
                  saveAt = paste("9.tmp/",
                                 sprintf("%03d", experiment), ".",
                                 sprintf("%03d", subsamp), ".",
                                 trait, ".", bayes[method], ".", sep = ""),
                  verbose = F)
      acc[method, ] <- cor(fit$yHat[is.na(Y.bglr)], tstY)
      #      wts[method,] <- fit$ETA$X$b
    }

    #RKHS
    G <- tcrossprod(X.bglr) / ncol(X.bglr)
    EVD <- eigen(G)
    ETA <- list(MRK = list(V = EVD$vectors, d = EVD$values, model="RKHS"))
    fit <- BGLR(y = Y.bglr, ETA = ETA, nIter = 10000, burnIn = 2000, thin = 10, verbose = FALSE)
    acc[12, ] <- cor(fit$yHat[is.na(Y.bglr)], tstY)

    #CARET MODELS
    #RANDOM FOREST
    trnCtrl <- trainControl(method = "cv", number = NFOLD)
    grid <- expand.grid(mtry = seq(10, 50, by = 2))
    rf <- train(Y ~ ., data = X.caret, method = 'rf', importance = TRUE, tuneGrid = grid, ntree = 750, trControl = trnCtrl)
    res <- as.vector(predict(rf, tstX))
    acc[6, ] <- cor(res, tstY)
    #    wts[6,] <- as.matrix(varImp(object = rf)$importance)

    #SVMLinear
    grid <- expand.grid(C = c(0.025, 0.0325, 0.05))
    svm <- train(Y ~ ., data = X.caret, method = "svmLinear", importance = TRUE, trControl = trnCtrl, tuneGrid = grid)
    res <- as.vector(predict(svm, tstX))
    acc[7, ] <- cor(res, tstY)

    #SVMRadial
    svmr <- train(Y ~ ., data = X.caret, method = "svmRadial", importance = TRUE, trControl = trnCtrl)
    res <- as.vector(predict(svmr, tstX))
    acc[15, ] <- cor(res, tstY)

    #SVMRadialSigma
    svmrs <- train(Y ~ ., data = X.caret, method = "svmRadialSigma", importance = TRUE, trControl = trnCtrl)
    res <- as.vector(predict(svmrs, tstX))
    acc[16, ] <- cor(res, tstY)

    #ENET
    grid <- expand.grid(alpha=seq(0.25, 0.75, by = 0.05),lambda = seq(0.01, 0.15, by = 0.001))
    enet <- train(Y ~ ., data = X.caret, method = "glmnet", trControl = trnCtrl, tuneGrid = grid)
    res <- predict(enet, tstX)
    acc[8, ] <- cor(res, tstY)
    #    wts[8, sub.id] <- as.matrix(varImp(object = enet)$importance)

    #NNET
    nnetGrid <-  expand.grid(size = c(5,7,10), decay = c(2.5,5,7.5))
    nnet <- train(Y ~ ., data = X.caret, method = 'nnet', tuneGrid = nnetGrid, maxit = 100000, trace = F, linout = 1)
    res <- predict(nnet, tstX)
    acc[9, ] <- cor(res, tstY)
    #    wts[9, sub.id] <- as.matrix(varImp(object = nnet)$importance)

    #BNNET
    bnnet <- train(Y ~ ., data = X.caret, method = 'brnn', maxit = 100000, trace = F, linout = 1)
    res <- predict(bnnet, tstX)
    acc[17, ] <- cor(res, tstY)

    #StochasticGradientBoostMachine
    grid <- expand.grid(.n.trees = c(500, 750, 1000),
                        .interaction.depth = c(5,7,10),
                        .shrinkage = c(0.01,0.02,0.05),
                        .n.minobsinnode = c(3,5))
    gbm <- train(Y ~ ., data = X.caret, method = "gbm", trControl = trnCtrl, tuneGrid = grid)
    res <- as.vector(predict(gbm, tstX))
    acc[13, ] <- cor(res, tstY)
    #  	wts[13, ] <- as.matrix(varImp(object = gbm)$importance)

    #ReleventVectorMachine(RadialBasisKernel)
    rvm <- train(Y ~ ., data = X.caret, method = "rvmRadial", sigma = 0.01, trControl = trnCtrl)
    res <- as.vector(predict(rvm, tstX))
    acc[14, ] <- cor(res, tstY)

    #KERNEL PARTIAL LEAST SQUARES REGRESSION
    pls <- train(Y ~ ., data = X.caret, method = 'kernelpls', trControl = trnCtrl)
    res <- predict(pls, tstX)
    acc[10, ] <- cor(res, tstY)
    #  	wts[10, ] <- as.matrix(varImp(object = pls)$importance)

    #LINEAR
    LM <- lm(Y ~ ., data = data.frame(X.caret))
    res <- as.vector(predict(LM, data.frame(tstX)))
    acc[11, ] <- cor(res, tstY)
    #    wts[11,] <- LM$coefficients[2:(NMET+1)]

    # GBLUP
    ans1 <- kin.blup(data = Y.gblup, K = A, geno = "id", pheno = "pheno")
    acc[18, ] <- cor(ans1$g[-(1:subsamp)], tstY)


  results <- data.frame(experiment = experiment,
                        trait = traits[trait],
                        subSize = subsamp,
                        model = models,
                        accuracy = acc)


if(!file.exists("resultsSub.txt")){
  write.table(results, "resultsSub.txt", row.names = F, col.names = T, append = T, quote = F)
} else {
  write.table(results, "resultsSub.txt", row.names = F, col.names = F, append = T, quote = F)
}





