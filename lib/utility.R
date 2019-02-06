library(methylumi)
library(lumi)
library(gplots)
library(ROCR)
library(boot)
library(HDCI)
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
  except_cpg <- sort(unique(c(nan,inf)))
  if(except_cpg[1] == 0){
    except_cpg[1] <- cpg_num
  }
  print("number of except CpG site")
  print(length(except_cpg))
  excepted_X <- X[-except_cpg, ]
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

boLasso = function(x, y, B){
  result <- bootLasso(x, y, B = B, intercept = TRUE)
  selected_coef_position <- which(result$Beta != 0)
  if(length(selected_coef_position) == 0){
    stop("failed train\n")
  }

  print("CpG site:")
  print(colnames(x[, selected_coef_position]))
  print("Coef:")
  print(result$Beta[selected_coef_position])
  print("Interval:")
  print(result$interval[1, selected_coef_position])

  return(selected_coef_position)
}

tsdSpecificCpG = function(all_cpg_ids){
  tsd_specific_cpg_ids <- scan("./data/DifferentiallyMethylatedProbes.txt", what = character(), sep = "\n")
  which_tsd_cpgs <- c()
  for(cpg_id in tsd_specific_cpg_ids){
    which_tsd_cpgs <- append(which_tsd_cpgs, which(all_cpg_ids == cpg_id))
  }
  return(which_tsd_cpgs)
}

allCpG = function(){
  return(TRUE)
}

cpgSelection.function = function(function_name){
  switch(function_name,
         "boLasso" = boLasso,
         "tsdSpecificCpG" = tsdSpecificCpG,
         "allCpG"= allCpG,
         NULL
         )
}

prediction.lassoBinomal = function(){
  predicted <- predict.glmnet(learned_model, s="lambda.min", newx=test_X)
  return(predicted)
}

prediction.lassoBinomal.cv = function(learned_model, test_X){
  predicted <- predict.cv.glmnet(learned_model, s="lambda.min", newx=test_X)
  return(predicted)
}

prediction.logistic = function(learned_model, test_X){
  predicted <- predict.glm(learned_model, newdata = data.frame(test_X) ,type = "responise")
  return(predicted)
}

prediction.function = function(function_name){
  switch(function_name,
    "lassoBinomal" = prediction.lassoBinomal,
    "lassoBinomal.cv" = prediction.lassoBinomal.cv,
    "logistic" = prediction.logistic,
     NULL
  )
}

rocAUC = function(predicted, y){
  pred <- prediction(predicted, y)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  return(auc)
}
