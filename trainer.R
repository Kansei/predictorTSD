source("./header.R")

resource_name = 'E-MTAB-4664'
# Read idats object
idats <- readIdatsObj(resource_name)
# Data preprocessing
X <- preprocessingIdats(idats)

# Read sample's attribute csv file
attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))
# Convert sleep status to 0/1
Y <- convertSleepStatus2Num(attributes_data[5:6])

# Split test and train data
test_size_rate = 0.3
data_size <- length(Y)
test_size <- as.integer(data_size*test_size_rate)
test_Y <- Y[1:test_size]
train_Y <- Y[(test_size+1):data_size]
test_X <- X[1:test_size,]
train_X <- X[(test_size+1):data_size,]

function_name = "tsdSpecificCpG"
which_cpg <- cpgSelection.function(function_name)(colnames(train_X))
# boLasso(train_X, train_Y, 10)
# tsdSpecificCpG(colnames(train_X))
# allCpG()

# Choose model
# logistic, lassoBinomal.cv, lassoBinomal, ridgeBinomial.cv, elasticNetBinomial.cv, adaptiveLassoBinomial.cv
model_name = "logistic"
fit_model <- model.function(model_name)(train_X[, which_cpg], train_Y)

predicted <- prediction.function(model_name)(fit_model, test_X[, which_cpg])
print(predicted)

auc <- rocAUC(predicted, test_Y)
print("AUC")
print(auc)

# saveModel(result)
