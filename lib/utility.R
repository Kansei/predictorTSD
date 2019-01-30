suppressPackageStartupMessages(library('methylumi'))
suppressPackageStartupMessages(library('lumi'))
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

# All sample normalize
normalize = function(idats){
  normfactors <- norm_factors(controldata = NULL, subjects = NULL, methylumidata = idats, type = "methylumi")
  normdata <- normalize_asmn(normfactors = normfactors, rawdata = NULL, featuredata = NULL, methylumidata = idats, type = "methylumi")
  return(normdata)
}

# Except nan, inf
exceptMissingValue = function(X){
  cpg_num <- length(X[,1])
  nan <- which(is.nan(X)) %% cpg_num
  inf <- which(is.infinite(X)) %% cpg_num
  missing_value <- sort(unique(c(nan,inf)))
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
