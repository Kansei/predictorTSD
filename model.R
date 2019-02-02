library(glmnet)
library(HDCI)
library(glm2)

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
  selected_coef_position <- which(result$Beta != 0)
  if(length(selected_coef_position) == 0){
    stop("failed train\n")
  }
  
  print("CpG site:")
  print(colnames(X[, selected_coef_position]))
  print("Coef:")
  print(result$Beta[selected_coef_position])
  print("p-value:")
  print(result$interval[1, selected_coef_position])
  data <- data.frame(train_Y, train_X[, selected_coef_position])
  
  glm.output <- glm(train_Y ~ .,  family = binomial(link = "logit"), data = data)
  glm.cv.output <- cv.glm(data = data, glm.output, K=10)
  print("Glm.cv error, adjusted error")
  print(glm.cv.output$delta)
  return(glm.output)
}

model.function = function(model_name){
  switch(model_name,
         "lassoBinomal.cv" = model.lassoBinomial.cv,
         "lassoBinomal" = model.lassoBinomial,
         "bootLasso" = model.bootLasso,
         NULL
         )
}
