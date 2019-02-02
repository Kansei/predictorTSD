library(methylumi)
library(lumi)
library(gplots)
library(ROCR)
library(boot)
source("./lib/asmn/norm_factors.R")
source("./lib/asmn/normalize_asmn.R")

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

# Convert sleep status to 0 or 1
convertSleepStatus2Num = function(attributes_data){
  sorted_data <- attributes_data[order(attributes_data$Assay.Name),]
  converted_data <- as.numeric(sorted_data[1] == "total acute sleep deprivation")

  return(converted_data)
}

# Convert Beta-value to M-value
convertBeta2M = function(beta){
  m <- beta2m(beta)
  return(m)
}

# All sample mean normalize
normalize = function(idats){
  # Create normalization factors.
  normfactors <- norm_factors(methylumidata = idats, type = "methylumi", controldata = NULL)
  # Perform all sample mean normalization.
  normdata <- normalize_asmn(normfactors = normfactors, methylumidata = idats, type = "methylumi", rawdata = NULL)
  return(normdata)
}

# Except nan, inf
exceptMissingValue = function(X){
  cpg_num <- length(X[,1])
  nan <- which(is.nan(X)) %% cpg_num
  inf <- which(is.infinite(X)) %% cpg_num
  missing_value <- sort(unique(c(nan,inf)))
  if(missing_value[1] == 0){
    missing_value[1] <- cpg_num
  }
  excepted_X <- X[-missing_value, ]
  return(excepted_X)
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

preprocessingIdats = function(idats){
  # Perform all sample mean normalization
  norm_idats <- normalize(idats)
  # Pick out Beta-value from idats
  beta_value <- norm_idats@assayData[["betas"]]
  # Convert beta-value to M-value
  m_value <- convertBeta2M(beta_value)
  # Except missing value
  excepted_m_value <- exceptMissingValue(m_value)
  # Transpose matrix to reformat (barcode_name x cpg_sites)
  X <- t(excepted_m_value)
  return(X)
}

rocAUC = function(predicted, y){ 
  pred <- prediction(predicted, y)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  return(auc)
}
