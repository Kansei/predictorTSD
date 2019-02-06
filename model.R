library(glmnet)
library(glm2)

model.lassoBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=1, standardize = TRUE, nfolds=10)
  coef = coef(result, s = result$lambda.min)
  print(colnames(x)[which(coef != 0)])
  return(result)
}

model.lassoBinomial = function(x, y){
  result <- glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}

model.logistic = function(x, y){
  data <- data.frame(y, x)
  glm.output <- glm(y ~ .,  family = binomial(link = "logit"), data = data)
  glm.cv.output <- cv.glm(data = data, glm.output, K=5)
  print("Glm.cv error, adjusted error")
  print(glm.cv.output$delta)
  return(glm.output)
}

model.function = function(model_name){
  switch(model_name,
         "lassoBinomal.cv" = model.lassoBinomial.cv,
         "lassoBinomal" = model.lassoBinomial,
         "logistic" = model.logistic,
         NULL
         )
}
