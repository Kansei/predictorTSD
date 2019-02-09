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
