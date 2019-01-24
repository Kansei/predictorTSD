suppressPackageStartupMessages(library('methylumi'))
suppressPackageStartupMessages(library('lumi'))
source("./lib/asmn/norm_factors.R")
source("./lib/asmn/normalize_asmn.R")

# Read idat files
readIdats = function(resource_name){
  resource_path <- paste("./data/", resource_name, sep = "")
  idatPath <- paste(resource_path, "/idat", sep = "")
  barcodes <- scan(paste(resource_path, "/barcodes.txt", sep = ""), what = character(), sep = "\n")
  idats <- methylumIDAT(barcodes = barcodes, idatPath=idatPath)

  return(idats)
}

# Convert sleep status to 0 or 1
convertSleepStatus2Num = function(attributes_data){
  sorted_data <- attributes_data[order(attributes_data$Assay.Name),]
  converted_data <- as.numeric(sorted_data[1] == "total acute sleep deprivation")

  return(converted_data)
}

# convrt Beta-value to M-value
convertBeta2M = function(beta){
  m <- beta2m(beta)
  return(m)
}

# all sample normalize
normalize = function(idats){
  normfactors <- norm_factors(controldata = NULL, subjects = NULL, methylumidata = idats, type = "methylumi")
  normdata <- normalize_asmn(normfactors = normfactors, rawdata = NULL, featuredata = NULL, methylumidata = idats, type = "methylumi")
  return(normdata)
}

# except nan, inf
exceptMissingValue = function(X){
  cpg_num <- length(X[,1])
  nan <- which(is.nan(X)) %% 485577
  inf <- which(is.infinite(X)) %% 485577
  missing_value <- sort(unique(c(nan,inf)))
  excepted_X <- X[-missing_value, ]
  return(excepted_X)
}

# split test and train data
trainTestSplit = functoin(X, Y, test_size_rate){
  data_size <- length(Y)
  test_size <- as.integer(data_size*test_size_rate)
  train_size <- data_size - test_size
  
  train_X <- X[1:train_size]
  test_X <- X[train_size+1:data_size]
  train_Y <- Y[1:train_size]
  test_Y <- Y[ttest_size+1:data_size]
  
  return(list(train_X, test_X, train_Y, test_Y))
}
