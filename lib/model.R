model.lassoBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=1, standardize = TRUE, nfolds=length(y), type.measure="deviance")
  print("lambda:")
  print(result$lambda.min)
  coef = coef(result, s = result$lambda.min)
  print(colnames(x)[which(coef != 0)])
  return(result)
}

model.lassoBinomial = function(x, y){
  result <- glmnet(x, y, family="binomial", alpha=1, standardize = TRUE)
  return(result)
}

model.ridgeBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=0, standardize = TRUE, nfolds=length(y), type.measure="deviance")
  print("lambda:")
  print(result$lambda.min)
  coef = coef(result, s = result$lambda.min)
  print(colnames(x)[which(coef[-1] != 0)])
  return(result)
}

model.elasticNetBinomial.cv = function(x, y){
  result <- cv.glmnet(x, y, family="binomial", alpha=0.5, standardize = TRUE, nfolds=length(y), type.measure="deviance", grouped=F)
  print("lambda:")
  print(result$lambda.min)
  coef <- coef(result, s = result$lambda.min)
  cpg_id <- c("intercept", colnames(x)[which(coef[-1] != 0)])
  coef <- coef[coef != 0]
  write.csv(data.frame(cpg_id, coef), "./data/coef.csv")
  
  return(result)
}

model.adaptiveLassoBinomial.cv = function(x, y){
  fit.ridge.cv <- cv.glmnet(x, y, family="binomial", alpha=0, standardize = TRUE, nfolds=length(y), type.measure="deviance")
  best_ridge_coef <- as.numeric(coef(fit.ridge.cv, s = fit.ridge.cv$lambda.min))[-1]

  result <- cv.glmnet(x, y, family="binomial", alpha=1, nfolds=10, penalty.factor = 1 / abs(best_ridge_coef))
  print("lambda:")
  print(result$lambda.min)
  coef = coef(result, s = result$lambda.min)
  print(colnames(x)[which(coef[-1] != 0)])
  
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
         "lassoBinomial.cv" = model.lassoBinomial.cv,
         "lassoBinomial" = model.lassoBinomial,
         "ridgeBinomial.cv" = model.ridgeBinomial.cv,
         "elasticNetBinomial.cv" = model.elasticNetBinomial.cv,
         "adaptiveLassoBinomial.cv" = model.adaptiveLassoBinomial.cv,
         "logistic" = model.logistic,
         NULL
         )
}
