library(glmnet)

model.lassoBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}

model.lassoBinomial = function(x, y){
  result <- glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}


