source("./lib/utility.R")
source("./model.R")

resource_name = 'E-MTAB-4664'
idats <- readIdats(resource_name)

# To reformat matrix (barcode_name x cpg_sites) <- normalize(idats)
beta_value <- norm_idats@assayData[["betas"]]
m_value <- convertBeta2M(beta_value)
excepted_m_value <- exceptMissingValue(m_value)
X <- t(excepted_m_value)

attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))
Y <- convertSleepStatus2Num(attributes_data[5:6])

train_test <- trainTestSplit(X, Y, 0.3)
train_X = train_test[1]
test_X = train_test[2]
train_Y = train_test[3]
test_Y = train_test[4]

result <- model.lassoBinomial.cv(X, Y)

paramater <- coef(result, s="lambda.min")

paramater[which(paramater != 0),1]

prd <- predict(result, newx = X, type = "class", s = "lambda.min")