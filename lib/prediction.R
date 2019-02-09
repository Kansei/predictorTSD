prediction.glmnet = function(fit_model, test_X){
  predicted <- predict.glmnet(fit_model, newx=test_X)
  return(predicted)
}

prediction.cv.glmnet = function(fit_model, test_X){
  predicted <- predict.cv.glmnet(fit_model, s="lambda.min", newx=test_X)
  return(predicted)
}

prediction.logistic = function(fit_model, test_X){
  predicted <- predict.glm(fit_model, newdata = data.frame(test_X) ,type = "response")
  return(predicted)
}

prediction.function = function(function_name){
  switch(function_name,
    "lassoBinomal" = prediction.glmnet,
    "lassoBinomal.cv" = prediction.cv.glmnet,
    "ridgeBinomial.cv" = prediction.cv.glmnet,
    "elasticNetBinomial.cv" = prediction.cv.glmnet,
    "adaptiveLassoBinomial.cv" = prediction.cv.glmnet,
    "logistic" = prediction.logistic,
     NULL
  )
}
