# Read idat files
readIdatsObj = function(resource_name){
  rdata_path <- paste0("./rdata/", resource_name, ".rdata")
  idats <- readRDS(rdata_path)
  return(idats)
}

# Read saved model
readModel = function(model_name){
  model_path <- paste0("./models/", model_name, ".rmodel")
  model <- readRDS(model_path)
  return(model)
}

# Save model
saveModel = function(model){
  model_name <- format(Sys.time(), "%Y%m%d%H%M%S")
  model_path <- paste0("./models/", model_name, ".rmodel")
  saveRDS(result, file = model_path)
}

# Split test and train data
trainTestSplit = function(X, Y, test_size_rate){
  data_size <- length(Y)
  test_size <- as.integer(data_size*test_size_rate)

  test_Y <- Y[1:test_size]
  train_Y <- Y[(test_size+1):data_size]
  test_X <- X[1:test_size,]
  train_X <- X[(test_size+1):data_size,]

  return(list(train_Y, test_Y, train_X, test_X))
}

rocAUC = function(predicted, y){
  pred <- prediction(predicted, y)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  return(auc)
}

vec.include = function(v, query){
  which_v <- c()
  for(q in query){
    which_v <- append(which_v, which(v == q))
  }
  return(which_v)
}

saveCoef = function(fit_model, cpg_ids){
  coef = coef(fit_model, s = fit_model$lambda.min)
  cpg_ids <- cpg_ids[which(coef[-1] != 0)]
  cpg_ids <- c("intercept",cpg_ids)
  cpg <- t(rbind(cpg_ids, coef@x))
  write.csv(cpg, "./data/coef/coef.csv")
}

closestGene = function(){
  cpg_data <- read.csv(paste0("./data/tsd_cpgs.csv"))[-1, ]
  
  cpg_id <- c()
  gene <- c()
  tss <- c()
  distance <- c()
  for(ncpg in 1:length(cpg_data[,1])){
    cpg <- cpg_data[ncpg, ]
    res <- findClosestGene(as.character(cpg$chr), cpg$position, genome="hg19", position="txStart")
    cpg_id <- append(cpg_id, as.character(cpg$cpg.id))
    gene <- append(gene, as.character(res$geneName[1]))
    tss <- append(tss, res$txStart[1])
    distance <- append(distance ,res$Distance[1])
  }
  cpg.data.frame <- data.frame(CpG.ID=cpg_id, Gene.Name=gene, TSS=tss, Distance=distance)
  write.csv(cpg.data.frame, "./data/closest_gene_from_cpg.csv")
}