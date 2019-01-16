library(glmnet)

model.lassoBinomial = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=1)
  return(result)
}

