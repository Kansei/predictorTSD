library(glmnet)
library(HDCI)

model.lassoBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}

model.lassoBinomial = function(x, y){
  result <- glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}

model.bootLasso = function(x, y, B){
  result <- bootLasso(x, y, B = B, intercept = TRUE)
  return(result)
}

model.function = function(model_name){
  switch(model_name,
         "lassoBinomal.cv" = model.lassoBinomial.cv,
         "lassoBinomal" = model.lassoBinomial,
         "bootLasso" = model.bootLasso,
         NULL
         )
}
