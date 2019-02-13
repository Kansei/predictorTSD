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
  na <- which(is.na(X)) %% cpg_num
  nan <- which(is.nan(X)) %% cpg_num
  inf <- which(is.infinite(X)) %% cpg_num
  
  if((length(na) == 0)&&(length(nan) == 0)&&(length(inf) == 0)){
    print("except CpG is nothing")
    return(X)
  }
  
  except_cpg <- sort(unique(c(na,nan,inf)))
  if(except_cpg[1] == 0){
    except_cpg[1] <- cpg_num
  }
  print("number of except CpG site")
  print(length(except_cpg))
  excepted_X <- X[-except_cpg, ]
  return(excepted_X)
}

standardize = function(data){
  return(scale(data))
}

preprocessingIdats = function(idats, s = T, m = T, n = T){
  # Perform all sample mean normalization
  if(n){
    idats <- normalize(idats)
  }
  # Pick out Beta-value from idats
  beta.matrix <- idats@assayData[["betas"]]
  
  # Convert beta-value to M-value
  if(m){
    methyl.matrix <- convertBeta2M(beta.matrix)  
  }else{
    methyl.matrix <- beta.matrix 
  }
  # Except missing value
  methyl.matrix <- exceptMissingValue(methyl.matrix)
  
  # Transpose matrix to reformat (barcode_name x cpg_sites)
  if(s){
    X <- standardize(t(methyl.matrix))
  }else{
    X <- t(methyl.matrix)
  }
  return(X)
}
