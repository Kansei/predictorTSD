source("./lib/utility.R")
library(glmnet)

resource_name = 'E-MTAB-4664'

idats <- readIdats(resource_name)

# To reformat matrix (barcode_name x cpg_sites)
beta_value <- t(idats@assayData[["betas"]])
m_value <- convertBeta2M(beta_value)
X <- m_value

attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))

Y <- convertSleepStatus2Num(attributes_data[5:6])

lasso.model.cv <- cv.glmnet(x = X, y = Y, family = "binomial", alpha = 1)

