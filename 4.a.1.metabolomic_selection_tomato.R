# 
# METABOLOMIC SELECTION FOR ENHANCED FRUIT FLAVOR
# Colantonio and Ferrao et al., 2022
# https://doi.org/10.1073/pnas.2115865119
#

###
# Figure 4a Metabolomic Selection Models
###

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
library(spls)

#library("qualityTools")
library(dplyr)
library(tidyr)


# We will paralellize this analysis using the foreach and doParallel packages
#  Some of the models such as XGBoost perform much better using gpu parallelization
#  This analysis is computationally intensive and would take an average laptop overnight, 
#  So during development don't run all experiments or limit the number of models tested
#  Paralellization speeds up the analysis significantly but also can crash R at times...

# Paralellization
library(foreach)
library(doParallel)

# RegisterParallel
#  Here is where we generate a cluster of R threads that can be used in the for each loops below
#  This asks for the number of threads available to the system minus 2 (to keep the system from being overwhelmed)

no_cores <- detectCores() - 2
cl <- makeCluster(no_cores)
registerDoParallel(cl)


########################
# Predictive Modeling  #
########################

#bb = read.csv("./data/bb_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")
tom = read.csv("./data/input/tom_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")

traits <- colnames(tom)[2:6]
mets <- colnames(tom)[7:74]

bayes <- c("BayesA", "BayesB", "BayesC", "BRR", "BL")
mods <- c(bayes,
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
          "Xgboost")


run <- format(Sys.time(), "tom_MS_run_%y%m%d_%I%M%p")
resdir <- "./results/fig4/a.panel/"


dir.create(paste0(resdir, run))

for(i in 1:10){
  dir.create(paste0(resdir, run, "/experiment", i))
  dir.create(paste0(resdir, run, "/experiment", i, "/9.tmp"))
  dir.create(paste0(resdir, run, "/experiment", i, "/5.accuracies"))
  dir.create(paste0(resdir, run, "/experiment", i, "/6.meanSquaredErrors"))
  dir.create(paste0(resdir, run, "/experiment", i, "/7.importanceValues"))
  dir.create(paste0(resdir, run, "/experiment", i, "/7.importanceValues/nested"))
}




####################
# Begin Experiment #
####################

NMOD <- 18
NPHEN <- 5
NMET <- 68
NFOLD <- 10

# Gather Initial Seeds
set.seed(0)
seeds <- sample(1:1e4, 10)


#####################
## Create Main Log ##
#####################

mainlog = paste0(resdir, run, "/1.mainlog.txt")
cat(strrep("#", 80), file = mainlog)
cat(paste("\n\nStarting run: Results_", run, "at", Sys.time(), "\n \n"), file = mainlog, append = T)
cat(strrep("#", 80), file = mainlog, append = T)
cat("\n", file = mainlog, append = T)


#####################
# Start Experiments #
#####################

for(experiment in 1:10){
  
  ####################  
  ## Writing to log ##   
  dir = paste0(resdir, run, "/experiment", experiment)
  tmpdir = paste0(resdir, run, "/experiment", experiment, "/9.tmp/")
  log <- paste(dir, "/1.log.txt", sep = "")
  cat(strrep("#", 80), file = log)
  cat(strrep("#", 80), file = mainlog, append = T)
  cat(paste("\nStarting experiment", experiment, "at", Sys.time(), "\n \n"), file = log, append = T)
  cat(paste("\n\nStarting experiment", experiment, "at", Sys.time(), "\n \n"), file = mainlog, append = T)
  ####################
  
  tom = read.csv("./data/input/tom_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")
  
  seed <- seeds[experiment]
  set.seed(seed)
  
  tom <- tom[sample(1:209, 200),]
  idx <- sample(rep(1:NFOLD, 20))
  
  
  ####################    
  ## Writing to log ##    
  cat(paste("The seed used for this experiment is:", seed, "\n \n"), file = log, append = T)
  cat(paste("The seed used for this experiment is:", seed, "\n \n"), file = mainlog, append = T)
  
  cat(paste("Storing data as indexed in this experiment at:", "\n", log, "\n \n"), file =  log, append = T)
  write.csv(cbind(idx, tom), paste(dir, "/2.index.csv", sep = ""), row.names = F)
  ####################
  
  
  # Create Result Matrices
  acc <- matrix(NA, nrow = NMOD, ncol = 1)
  mse <- matrix(NA, nrow = NMOD, ncol = 1)
  wts <- matrix(NA, nrow = NMOD, ncol = NMET)
  
  
  ################
  # START MODELS #
  ################
  
  ####################  
  ## Writing to log ##    
  cat(paste("Begin Cross Validation", "\n\n", sep = ""), file = log, append = T)
  cat(paste("Date and Time\tTrait\tTraitName\tTraitID\tFold\tSeed\n"), file = log, append = T)
  ####################
  
  start <- proc.time()
  
  tomresults <- foreach(trait = c(1:NPHEN), 
                        .combine = 'list', 
                        .multicombine = TRUE, 
                        .packages = c("BGLR", 
                                      "caret", 
                                      "randomForest", 
                                      "pls", 
                                      "nnet", 
                                      "glmnet", 
                                      "neuralnet", 
                                      "gbm")) %:%   
    foreach(cv = 1:NFOLD, 
            .combine = 'list', 
            .multicombine = TRUE) %dopar% {
              
              
              ## Writing to log ##    
              cat(paste(Sys.time(), "\tstart\t", traits[trait], "\t", trait, "\tfold\t", cv, "\tseed\t", paste(seed, trait, cv, sep = ""), "\n", sep = ""), file = log, append = T)
              ####################
              
              set.seed(paste(seed, trait, cv, sep = ""))
              
              trnX <- as.matrix(tom[, (2+NPHEN):(1+NPHEN+NMET)])
              trnY <- as.matrix(tom[, trait + 1])
              trnY[idx == cv] <- NA
              tstY <- as.matrix(tom[idx == cv, trait + 1])
              
              
              #BAYES
              for(method in 1:5){
                ETA <- list(X = list(X = trnX, model = bayes[method]))
                fit <- BGLR(y = trnY, ETA = ETA, nIter = 30000, burnIn = 10000, 
                            saveAt = paste0(tmpdir, "fold", cv, "_", bayes[method], "_", traits[trait],  "_"), 
                            thin = 100, verbose = F)
                acc[method,] <- cor(fit$yHat[idx == cv], tstY)
                mse[method,] <- mean((fit$yHat[idx == cv] - (tstY))^2, na.rm = TRUE)
                wts[method,] <- fit$ETA$X$b 
              }      
              
              #RKHS
              G <- tcrossprod(trnX) / ncol(trnX)
              EVD <- eigen(G)
              ETA <- list(MRK = list(V = EVD$vectors, d = EVD$values, model="RKHS"))
              fit <- BGLR(y = trnY, ETA = ETA, nIter = 10000, burnIn = 2000, 
                          saveAt = paste0(tmpdir, "fold", cv, "_RKHS_", traits[trait], "_"),
                          thin = 10, verbose = FALSE)
              acc[12,] <- cor(fit$yHat[idx == cv], tstY)
              mse[12,] <- mean((fit$yHat[idx == cv] - (tstY))^2, na.rm = TRUE)
              
              
              #CARET MODELS
              carX <- as.matrix(cbind(Y = tom[idx != cv, trait + 1], tom[idx != cv, (2+NPHEN):(1+NPHEN+NMET)]))
              tstX <- tom[idx == cv, (2+NPHEN):(1+NPHEN+NMET)]
              
              
              #RANDOM FORREST
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              grid <- expand.grid(mtry = seq(10, 50, by = 2))
              rf <- train(Y ~ ., data = carX, method = 'rf', importance = TRUE, tuneGrid = grid, ntree = 750, trControl = trnCtrl)
              res <- as.vector(predict(rf, tstX))
              acc[6,] <- cor(res, tstY)
              mse[6,] <- mean((res - tstY) ^ 2, na.rm = TRUE)
              wts[6,] <- as.matrix(varImp(object = rf)$importance)
              
              
              #SVMLinear
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              grid <- expand.grid(C = c(0.025, 0.0325, 0.05))
              svm <- train(Y ~ ., data = carX, method = "svmLinear", importance = TRUE, trControl = trnCtrl, tuneGrid = grid)
              res <- as.vector(predict(svm, tstX))
              acc[7,] <- cor(res, tstY)
              mse[7,] <- mean((res - tstY)^2, na.rm = TRUE)
              
              #SVMRadial
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              svmr <- train(Y ~ ., data = carX, method = "svmRadial", importance = TRUE, trControl = trnCtrl)
              res <- as.vector(predict(svmr, tstX))
              acc[15,] <- cor(res, tstY)
              mse[15,] <- mean((res - tstY) ^ 2, na.rm = TRUE)
              
              #SVMRadialSigma
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              svmrs <- train(Y ~ ., data = carX, method = "svmRadialSigma", importance = TRUE, trControl = trnCtrl)
              res <- as.vector(predict(svmrs, tstX))
              acc[16,] <- cor(res, tstY)
              mse[16,] <- mean((res - tstY)^2, na.rm = TRUE)
              
              #ENET
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              grid <- expand.grid(alpha=seq(0.25, 0.75, by = 0.05),lambda = seq(0.01, 0.15, by = 0.001))
              enet <- train(Y ~ ., data = carX, method = "glmnet", trControl = trnCtrl, tuneGrid = grid)
              res <- predict(enet, tstX)
              acc[8,] <- cor(res, tstY)
              mse[8,] <- mean((res - tstY)^2, na.rm = TRUE) 
              wts[8,] <- as.matrix(varImp(object = enet)$importance)
              
              #NNET
              nnetGrid <-  expand.grid(size = c(5,7,10), decay = c(2.5,5,7.5))
              nnet <- train(Y ~ ., data = carX, method = 'nnet', tuneGrid = nnetGrid, maxit = 100000, trace = F, linout = 1)
              res <- predict(nnet, tstX)
              acc[9,] <- cor(res, tstY)
              mse[9,] <- mean((res - tstY)^2, na.rm = TRUE) 
              wts[9,] <- as.matrix(varImp(object = nnet)$importance)
              
              #BNNET
              bnnet <- train(Y ~ ., data = carX, method = 'brnn', maxit = 100000, trace = F, linout = 1)
              res <- predict(bnnet, tstX)
              acc[17,] <- cor(res, tstY)
              mse[17,] <- mean((res - tstY)^2, na.rm = TRUE) 
              
              #StochasticGradientBoostMachine
              grid <- expand.grid(.n.trees = c(500, 750, 1000), .interaction.depth = c(5,7,10), .shrinkage = c(0.01,0.02,0.05), .n.minobsinnode = c(3,5))
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              gbm <- train(Y ~ ., data = carX, method = "gbm", trControl = trnCtrl, tuneGrid = grid)
              res <- as.vector(predict(gbm, tstX))
              acc[13,] <- cor(res, tstY)
              mse[13,] <- mean((res - tstY)^2, na.rm = TRUE) 
              wts[13,] <- as.matrix(varImp(object = gbm)$importance)
              
              #ReleventVectorMachine(RadialBasisKernel)
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              rvm <- train(Y ~ ., data = carX, method = "rvmRadial", sigma = 0.01, trControl = trnCtrl)
              res <- as.vector(predict(rvm, tstX))			
              acc[14,] <- cor(res, tstY)
              mse[14,] <- mean((res - tstY)^2, na.rm = TRUE) 
              
              #KERNEL PARTIAL LEAST SQUARES REGRESSION
              trnCtrl <- trainControl(method = "cv", number = NFOLD)
              pls <- train(Y ~ ., data = carX, method = 'kernelpls', trControl = trnCtrl)
              res <- predict(pls, tstX)
              acc[10, ] <- cor(res, tstY)
              mse[10, ] <- mean((res - tstY)^2, na.rm = TRUE) 
              wts[10, ] <- as.matrix(varImp(object = pls)$importance)
              
              #xgboost    
              #      grid <- expand.grid(nrounds = 1000, lambda = 1, alpha = 0, eta = 1)
              #      trnCtrl <- trainControl(method = "cv", number = 10)
              #      xbs <- train(Y ~ ., data = carX, method = "xgbLinear", trControl = trnCtrl, tuneGrid = grid)
              #      res <- as.vector(predict(xbs, tstX))
              #      acc[18,] <- cor(res, tstY)
              
              #LINEAR
              LM <- lm(Y ~ ., data = data.frame(carX))
              res <- as.vector(predict(LM, data.frame(tstX)))
              acc[11,] <- cor(res, tstY)
              mse[11,] <- mean((res - tstY)^2, na.rm = TRUE)
              wts[11,] <- LM$coefficients[2:(NMET+1)]
              
              
              ####################        
              ## Writing to log ##    
              cat(paste(Sys.time(), "\tend\t", traits[trait], "\t", trait, "\tfold\t", cv, "\tseed\t", paste(seed, trait, cv, sep = ""), "\n", sep = ""), file = log, append = T)
              ####################
              
              list(acc, mse, wts)
            }
  

  
  
  ##########################
  ## Recording Experiment ##  
  ##########################  
  
  ####################  
  ## Writing to log ##   
  time <- proc.time() - start
  cat("\n\n\n", file = log, append = T)
  cat(paste(Sys.time(), "Finishing Experiment\t", experiment, "\n"), file = log, append = T)
  cat(paste("Model Training elapsed:\t", time[3]/60 , " Minutes\n \n"), file = log, append = T)
  cat(paste(names(time), sep = "\t"), file = log, append = T)
  cat("\n", file = log, append = T)
  cat(paste(time, sep = "\t"), file = log, append = T)
  cat("\n\n\n", file = log, append = T)
  ####################
  
  ############################
  ## Calculating Accuracies ##  
  ############################  
  
  accs <- rep(list(matrix(NA, nrow = NMOD, ncol = NFOLD)), NPHEN)
  mses <- rep(list(matrix(NA, nrow = NMOD, ncol = NFOLD)), NPHEN)
  
  wtss <- rep(list(matrix(NA, nrow = NMET, ncol = NFOLD)), NMOD)
  wtss <- rep(list(wtss), NPHEN)
  
  rs <- matrix(NA, nrow = NMOD, ncol = NPHEN)
  rsmse <- matrix(NA, nrow = NMOD, ncol = NPHEN)
  rswts <- rep(list(matrix(NA, nrow = NMET, ncol = NMOD)), NPHEN)
  
  
  for(i in 1:NPHEN){
    for(ii in 1:NFOLD){
      accs[[i]][, ii] <- tomresults[[i]][[ii]][[1]][1:18]
      mses[[i]][, ii] <- tomresults[[i]][[ii]][[2]][1:18]
      
      for(iii in 1:NMOD){
        wtss[[i]][[iii]][, ii] <- tomresults[[i]][[ii]][[3]][iii, ]
      }
    }
  }
  
  
  for(i in 1:NPHEN){
    rs[,i] <- rowMeans(accs[[i]])
    rsmse[,i] <- rowMeans(mses[[i]])
    
    
    colnames(accs[[i]]) <- paste0("fold", 1:NFOLD)
    rownames(accs[[i]]) <- mods
    write.csv(accs[[i]], file = paste0(dir, "/5.accuracies/accuracy_",  traits[i], ".csv"))
    
    colnames(mses[[i]]) <- paste0("fold", 1:NFOLD)
    rownames(mses[[i]]) <- mods
    write.csv(mses[[i]], file = paste0(dir, "/6.meanSquaredErrors/mse_", traits[i], ".csv"))
    
    
    for(ii in 1:NMOD){
      rswts[[i]][, ii] <- rowMeans(wtss[[i]][[ii]])
      
      colnames(wtss[[i]][[ii]]) <- paste0("fold", 1:NFOLD)
      rownames(wtss[[i]][[ii]]) <- mets
      write.csv(wtss[[i]][[ii]], file = paste0(dir, "/7.importanceValues/nested/", mods[ii], "_", traits[i], ".csv"))
    }
    
    colnames(rswts[[i]]) <- mods
    rownames(rswts[[i]]) <- mets
    write.csv(rswts[[i]], file = paste0(dir, "/7.importanceValues/weights_", traits[i],".csv"))
  }
  
  
  colnames(rs) <- colnames(tom)[2:(NPHEN+1)]
  rownames(rs) <- mods
  
  res <- cbind(rowMeans(rs))
  res <- cbind(res[order(res, decreasing = T),])
  
  colnames(rsmse) <- colnames(tom)[2:(NPHEN + 1)]
  rownames(rsmse) <- mods
  
  
  ####################  
  ## Writing to log ##    
  cat(paste("Writing Final Accuracies to file\n", paste(dir, "/FinalAccuracies.csv", sep = ""), "\n \n"), file = log, append = T)
  write.csv(rs, file = paste(dir, "/3.FinalAccuracies.csv", sep = ""))
  
  sink(log, append = T)
  print(rs)
  sink()
  
  cat(paste("\n\n\nWriting Final Mean Squared Errors to file\n", paste(dir, "/FinalMeanSquaredErrors.csv", sep = ""), "\n \n"), file = log, append = T)
  write.csv(rsmse, file = paste(dir, "/4.FinalMeanSquaredErrors.csv", sep = ""))
  
  sink(log, append = T)
  print(rsmse)
  sink()
  
  cat(paste("\n\nFinal Rankings for experiment ", experiment, "are as follows:\n"), file = log, append = T)
  
  sink(log, append = T)
  print(res)
  sink()
  
  cat(paste("\n\n", Sys.time(), "Finished Experiment\t", experiment, "\n \n"), file = log, append = T)
  cat(strrep("#", 80), file = log, append = T)
  cat("\n\n\n\n", file = log, append = T)
  ####################
  
  #########################
  ## Writing to Main log ##
  cat(paste("Final Accuracies for experiment \t", experiment, "are as follows:", "\n \n"), file = mainlog, append = T)
  
  sink(mainlog, append = T)
  print(rs)
  sink()
  
  cat(paste("\n\nFinal Rankings for experiment ", experiment, "are as follows:\n"), file = mainlog, append = T)
  
  sink(mainlog, append = T)
  print(res)
  sink()
  
  cat(paste("\n\n", Sys.time(), "Finished Experiment\t", experiment, "\n \n"), file = mainlog, append = T)
  cat(strrep("#", 80), file = mainlog, append = T)
  cat("\n\n\n\n", file = mainlog, append = T)
  #######################  
  
  
  #######################
  ## End of Experiment ##  
  #######################
  
  # Delete Temp Files
  # delete a directory -- must add recursive = TRUE
  unlink(tmpdir, recursive = TRUE)
  rm(list = ls()[!ls() %in% c("traits", "mods", "mets", "resdir", "bayes", "run", "NMOD", "NMET", "NPHEN", "NFOLD", "seeds", "cl", "mainlog")])
}




results <- list()

for(i in 1:10){
  dir = paste0(resdir, "200408_0638PM", "/experiment", i)
  results[[i]] <- read.csv(file = paste(dir, "/3.FinalAccuracies.csv", sep = ""))  
}

#results[[1]]$experiment <- 1
res <- results[[1]]

for(i in 1:10){
  #  results[[i]]$experiment <- i
  res <- cbind(res, results[[i]][,2:8])
}

rmns <- rowMeans(res[,2:dim(res)[2]])

#idk?? not ideal
fin <- as.matrix(cbind(as.character(res$X), as.numeric(rmns)))
fin <- fin[1:17,]
fin <- as.matrix(fin)
fin[,2] <- as.numeric(as.character(fin[,2]))

fin <- fin[order(fin[,2], decreasing = T),]

fin <- as.data.frame(fin, rownames = T)






# #############
# CUDA XGBOOST 04/11/20
# ############


####################
# Begin Experiment #
####################

NMOD <- 18
NPHEN <- 7
NMET <- 68
NFOLD <- 10

# Gather Initial Seeds
set.seed(0)
seeds <- sample(1:1e4, 10)


#####################
## Create Main Log ##
#####################

#  mainlog = paste0(resdir, run, "/1.mainlog.txt", sep = "")
#  cat(strrep("#", 80), file = mainlog)
#  cat(paste("\n\nStarting run: Results_", run, "at", Sys.time(), "\n \n"), file = mainlog, append = T)
#  cat(strrep("#", 80), file = mainlog, append = T)
#  cat("\n", file = mainlog, append = T)


#####################
# Start Experiments #
#####################

# Some of these paths need updating for github

for(experiment in 1:10){
  
  ####################  
  ## Writing to log ##   
  dir = paste0(resdir, run, "/experiment", experiment)
  #  log <- paste(dir, "/1.log.txt", sep = "")
  #  cat(strrep("#", 80), file = log)
  #  cat(strrep("#", 80), file = mainlog, append = T)
  #  cat(paste("\nStarting experiment", experiment, "at", Sys.time(), "\n \n"), file = log, append = T)
  #  cat(paste("\n\nStarting experiment", experiment, "at", Sys.time(), "\n \n"), file = mainlog, append = T)
  ####################
  
  tom = read.csv("./data/tom_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")
  
  seed <- seeds[experiment]
  set.seed(seed)
  
  tom <- tom[sample(1:209, 200),]
  idx <- sample(rep(1:NFOLD, 20))
  
  
  ####################    
  ## Writing to log ##    
  #  cat(paste("The seed used for this experiment is:", seed, "\n \n"), file = log, append = T)
  #  cat(paste("The seed used for this experiment is:", seed, "\n \n"), file = mainlog, append = T)
  
  #  cat(paste("Storing data as indexed in this experiment at:", "\n", log, "\n \n"), file =  log, append = T)
  #  write.csv(cbind(idx, tom), paste(dir, "/2.index.csv", sep = ""), row.names = F)
  ####################
  
  
  # Create Result Matrices
  acc <- matrix(NA, nrow = NMOD, ncol = 1)
  mse <- matrix(NA, nrow = NMOD, ncol = 1)
  wts <- matrix(NA, nrow = NMOD, ncol = NMET)
  
  
  ################
  # START MODELS #
  ################
  
  ####################  
  ## Writing to log ##    
  #  cat(paste("Begin Cross Validation", "\n\n", sep = ""), file = log, append = T)
  #  cat(paste("Date and Time\tTrait\tTraitName\tTraitID\tFold\tSeed\n"), file = log, append = T)
  ####################
  
  start <- proc.time()
  
  
  # Create orthogonal array
  tdo <- taguchiDesign("L27_3")
  
  # Parameters
  eta <- c(0.005, 0.02, 0.1)
  dep <- c(1, 2, 3)
  ss=c(0.4, 0.6, 0.8)
  gamma=c(0.01, 0.1, 1)
  mcw=c(1, 2, 5)
  lambda = c(0.005, 0.01, 0.02)
  alpha = c(0.02, 0.1, 0.2)
  nr = c(200, 300, 400)
  partrees = c(2, 3, 5)
  
  par <- tdo@design[,c(1:9)]
  res2 <- matrix(NA, nrow = 8, ncol = dim(par)[1])
  res <- matrix(NA, nrow = 8*nr, ncol = dim(par)[1])
  varimp <- rep(list(rep(list(), 10)), 8)
  acc <- matrix(NA, nrow = 8, ncol = 10)
  cors <- matrix(nrow = 10, ncol = 8)
  start <- proc.time()
  
  # Perform sampling, assess optimal parameters, train final model, and compute accuracies
  for(trait in 1:1){  
    for(cv in 1:10){
      for(i in 1:(dim(par)[1])){
        
        set.seed(paste(seed, trait, cv, sep = ""))
        
        trnX <- as.matrix(tom[, (2+NPHEN):(1+NPHEN+NMET)])
        trnY <- as.matrix(tom[, trait + 1])
        trnY[idx == cv] <- NA
        tstY <- as.matrix(tom[idx == cv, trait + 1])
        
        
        fit <- xgb.cv(data = as.matrix(tom[idx != cv, (2+NPHEN):(1+NPHEN+NMET)]), 
                      label = as.matrix(tom[idx != cv, trait + 1]), 
                      nfold = 10,
                      eta = eta[par[i,1]],
                      eval_metric = 'rmse',
                      max_depth = dep[par[i,2]],
                      subsample = ss[par[i,3]],
                      gamma = gamma[par[i,4]],
                      min_child_weight = mcw[par[i,5]],
                      nrounds = nr[par[i,6]],
                      lambda = lambda[par[i,7]],
                      alpha = alpha[par[i,8]],
                      num_parallel_tree = partrees[par[i,9]],
                      nthread = 4,
                      #                       updater = "grow_gpu_hist",
                      tree_method = 'gpu_hist', verbose = F)
        
        res2[trait, i] <- fit$evaluation_log$test_rmse_mean[nr[par[i,6]]]      
        rm(fit)
        gc()
        
      }
      
      df <- as.data.frame(tdo)
      opt <- vector()
      
      for(i in 1:9){
        opt[i] <- which.min(tapply(res2[(trait),], df[,i+3], mean))
      }
      
      fit <- xgboost(data = as.matrix(tom[idx != cv, 9:76]), 
                     label = as.matrix(tom[idx != cv, trait + 1]), 
                     eta = eta[opt[1]],
                     eval_metric = 'rmse',
                     max_depth = dep[opt[2]],
                     subsample = ss[opt[3]],
                     gamma = gamma[opt[4]],
                     min_child_weight = mcw[opt[5]],
                     nrounds = nr[opt[6]],
                     lambda = lambda[opt[7]],
                     alpha = alpha[opt[8]],
                     num_parallel_tree = partrees[opt[9]],
                     nthread = 4,
                     #                   updater = "grow_gpu_hist",                   
                     tree_method = 'gpu_hist', verbose = F)
      
      pred <- predict(fit, as.matrix(tom[idx == cv, 9:76]), reshape =T)
      acc[trait, cv] <- cor(pred, tom[idx == cv, trait + 1])
      varimp[[trait]][[cv]] <- xgb.importance(colnames(tom)[9:dim(tom)[2]], model = fit)
    }
  }
  
  proc.time() - start

  
  
  acc <- matrix(NA, nrow = 7, ncol = 10)
  accs <- rep(list(acc), 10)
  
  varimp <- rep(list(rep(list(), 10)), 7)
  varimps <- rep(list(varimp), 10)
  mainlog <- "./data/xgboost_tomato_results/1.mainlog.txt"
  # Perform sampling, assess optimal parameters, train final model, and compute accuracies
  
  for(experiment in 1:10){
    cat(strrep("#", 80), file = mainlog, append = T)
    cat(paste("\n\nStarting experiment", experiment, "at", Sys.time(), "\n \n"), file = mainlog, append = T)
    
    tom = read.csv("./data/input/tom_imputed_scaled.csv", header = T, check.names = F, na.strings = ".")
    
    seed <- seeds[experiment]
    set.seed(seed)
    
    tom <- tom[sample(1:209, 200),]
    idx <- sample(rep(1:NFOLD, 20))
    
    start <- proc.time()
    
    
    for(trait in 1:7){  
      for(cv in 1:10){
        
        set.seed(paste(seed, trait, cv, sep = ""))
        
        fit <- xgboost(data = as.matrix(tom[idx != cv, (2+NPHEN):(1+NPHEN+NMET)]), 
                       label = as.matrix(tom[idx != cv, trait + 1]), 
                       nrounds = 1000, 
                       eta = 0.1,
                       max_depth = 2,
                       subsample = 0.5,
                       gamma = 0.01,
                       min_child_weight = 5,
                       lambda = 0.0075,
                       alpha = 0.2,
                       num_parallel_tree = 20,
                       tree_method = 'gpu_hist', verbose = F)
        
        pred <- predict(fit, as.matrix(tom[idx == cv, (2+NPHEN):(1+NPHEN+NMET)]), reshape =T)
        accs[[experiment]][trait, cv] <- cor(pred, tom[idx == cv, trait + 1])
        varimps[[experiment]][[trait]][[cv]] <- xgb.importance(colnames(tom)[(2+NPHEN):(1+NPHEN+NMET)], model = fit)
        
        rm(fit)
        gc()
      }
    }
    
    proc.time() - start
    
    
    #accs[[experiment]] <- acc[1:7,]
    
    #colnames(accs[[experiment]]) <- folds
    rownames(accs[[experiment]]) <- colnames(tom)[2:8]
    #########################
    ## Writing to Main log ##
    
    cat(paste("Final Accuracies for experiment \t", experiment, "are as follows:", "\n \n"), file = mainlog, append = T)
    
    sink(mainlog, append = T)
    print(accs[[experiment]][order(rowMeans(accs[[experiment]]), decreasing = T), order(colMeans(accs[[experiment]]), decreasing = T)])
    sink()
    
    cat(paste("\n\nFinal Rankings for experiment ", experiment, "are as follows:\n"), file = mainlog, append = T)
    
    sink(mainlog, append = T)
    print(mean(rowMeans(accs[[experiment]])))
    sink()
    
    cat(paste("\n\n", Sys.time(), "Finished Experiment\t", experiment, "\n \n"), file = mainlog, append = T)
    cat(paste("Elapsed time:\n"), file = mainlog, append = T)
    cat((proc.time() - start)[3] / 60, file = mainlog, append = T)
    cat(" minutes \n\n", file = mainlog, append = T)
    
    cat(strrep("#", 80), file = mainlog, append = T)
    cat("\n\n\n\n", file = mainlog, append = T)
    #######################  
    
  }
}


####
# Processing results presented in the paper
####

######
#  04/12/20
#######
  
  for(i in 1:10){
    write.csv(accs[[i]], file = paste("./data/xgboost_tomato_results/exp", i, "Accuracies.csv", sep = ""))
  }
  
  results <- list()
  
  for(i in 1:10){
    dir = paste0(resdir, "200408_0638PM", "/experiment", i, sep = "")
    results[[i]] <- read.csv(file = paste(dir, "/3.FinalAccuracies.csv", sep = ""))  
    results[[i]][18, 2:8] <- rowMeans(accs[[i]])
  }
  
  #results[[1]]$experiment <- 1
  res <- results[[1]]
  
  for(i in 1:10){
    #  results[[i]]$experiment <- i
    res <- cbind(res, results[[i]][,2:8])
  }
  
  rmns <- rowMeans(res[,2:dim(res)[2]])
  
  fin <- as.matrix(cbind(as.character(res$X), as.numeric(rmns)))
  fin <- as.matrix(fin)
  fin[,2] <- as.numeric(as.character(fin[,2]))
  fin <- fin[order(fin[,2], decreasing = T),]
  fin <- as.data.frame(fin, rownames = T)
  
  
  results[[1]]$experiment <- 1
  res <- results[[1]]
  
  for(i in 1:10){
    results[[i]]$experiment <- i
    res <- rbind(res, results[[i]])
  }
  
  cmns <- colMeans(res[,2:dim(res)[2]])
  
  fin <- as.matrix(cbind(as.character(res$X), as.numeric(rmns)))
  fin <- as.matrix(fin)
  fin[,2] <- as.numeric(as.character(fin[,2]))
  fin <- fin[order(fin[,2], decreasing = T),]
  fin <- as.data.frame(fin, rownames = T)
  
  
  avgs <- data.frame(matrix(NA, ncol = 8, nrow = 18))
  
  for(i in 1:18){
    res <- results[[1]][i,]
    
    for(ii in 1:10){
      #  results[[i]]$experiment <- i
      res <- rbind(res, results[[ii]][i,])
    }
    
    avgs[i, 2:8] <- colMeans(res[,2:8])
  }
  
  
  colnames(avgs) <- colnames(results[[1]])
  avgs$X <- mods  
  
  fin2 <- data.frame(cbind(avgs$X, rowMeans(avgs[,2:8])))
  fin2 <- fin2[order(fin2[,2], decreasing = T),]
  fin2
  
  avgs <- avgs[order(rowMeans(avgs[, 2:8]), decreasing = T), c(1, 1 + order(colMeans(avgs[, 2:8]), decreasing = T))]
  avgs[19,2:8] <- colMeans(avgs[2:8])
  avgs[,9] <- rowMeans(avgs[2:8])
  
  avgs[19,1] <- "TraitAverages"
  colnames(avgs)[c(1,9)] <- c("Models", "ModelAverages")
  
  write.csv(avgs, "./data/xgboost_tomato_results/FinalAccuracies")
  #avgs <- avgs[c(19, 1:18), c(1, 9, 2:8)]

  

































