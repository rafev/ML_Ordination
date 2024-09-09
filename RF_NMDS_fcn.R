RF_ordination = function(Meth.array, limma.out, pheno_file, phenotypes, lfc){
  
  X = t(Meth.array[unique(limma.out[limma.out$adj.P.Val <= 0.05 & abs(limma.out$logFC) >= lfc,]$Row.names), ])
  
  Ytmp = pheno_file[pheno_file$phenotype %in% phenotypes, c("phenotype", "SampleID")]
  Ytmp$phenotype = paste0("HC", Ytmp$phenotype)
  Ytmp$phenotype = as.factor(as.character(Ytmp$phenotype))
  
  XX = X[row.names(X) %in% Ytmp$SampleID, ]
  
  XX = merge(XX, Ytmp, by.x = "row.names", by.y = "SampleID")
  rownames(XX) = XX$Row.names
  XX = XX[,2:ncol(XX)]
  XX = relocate(XX, phenotype)
  
  #rf = randomForest::randomForest(x = XX, y = YY, proximity = TRUE, ntree = 1000, importance = TRUE)
  
  set.seed(1234)
  inTrain = createDataPartition(
    y = XX[,"phenotype"],
    p = .75,
    list = FALSE
  )
  Train = XX[inTrain,]
  Valid = XX[-inTrain,]
  
  folds = 10
  repeats = 3
  cctrl1 = trainControl(method = "repeatedcv",
                        number = folds,     # number of folds
                        repeats = repeats,  # number of repeats
                        classProbs = T,     # assigns class probability
                        summaryFunction = ifelse(length(phenotypes) == 2, defaultSummary, multiClassSummary),
                        sampling = "up", # method of subsampling within Train's sampling routine
                        verboseIter = T)
  rf = train(x = Train[,2:ncol(Train)],
             y = Train[[1]],
             method = "rf", trControl = cctrl1,
             metric = ifelse(length(phenotypes) == 2, "Accuracy", "MeanBalancedAccuracy"),
             family = if(length(phenotypes) > 2) {"multinomial"},
             nTree = 1000, importance = TRUE, proximity = TRUE)
  rf.predict = predict(rf, newdata = Valid, type = "raw")
  varImp_plt = randomForest::varImpPlot(rf$finalModel)
  
  # -- Calculate CV performance
  if (length(phenotypes > 2)) {
    
    auc_res = list()
    for (i in phenotypes) {
      Valid.tmp = gsub("HC", "", Valid[[1]])
      auc_res[[i]] = pROC::roc((ifelse(Valid.tmp == i, 1, 2)),
                               as.numeric(rf.predict))$auc
    }
    auc_res = mean(unlist(auc_res))
    
  } else {
    auc_res = pROC::roc(as.numeric(Valid), as.numeric(rf.predict))$auc
  }
  
  # -- Calculate proximity matrix on whole
  proximity = predict(object = rf$finalModel, newdata = XX, type = "prob", proximity = TRUE)$proximity
  rf.mds = stats::cmdscale(1 - proximity, eig = TRUE, k = 3)
  rf.mds = data.frame(rf.mds$points)
  rf.mds$class = XX[[1]]
  rf.mds$col = ifelse(rf.mds$class == "HC1","#42B540FF",
                      ifelse(rf.mds$class == "HC2", "red",
                             ifelse(rf.mds$class == "HC3", "#00486BFF",
                                    ifelse(rf.mds$class == "HC4", "#8B4500", "purple"))))
  
  # -- Calculate proximity matrix on validation
  proximity = predict(object = rf$finalModel, newdata = Valid, type = "prob", proximity = TRUE)$proximity
  rf.mds_valid = stats::cmdscale(1 - proximity, eig = TRUE, k = 3)
  rf.mds_valid = data.frame(rf.mds_valid$points)
  rf.mds_valid$class = Valid[[1]]
  rf.mds_valid$col = ifelse(rf.mds_valid$class == "HC1","#42B540FF",
                            ifelse(rf.mds_valid$class == "HC2", "red",
                                   ifelse(rf.mds_valid$class == "HC3", "#00486BFF",
                                          ifelse(rf.mds_valid$class == "HC4", "#8B4500", "purple"))))
  
  
  Output = list("varImpPlt" = varImp_plt,
                "auc" = auc_res,
                "MDS" = rf.mds,
                "Validation" = rf.mds_valid)
  return(Output)
  
}