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

normalize = function(){
  
}
