source("./lib/utility.R")
source("./model.R")

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
cpgSelection <- cpgSelection.function(function_name)
which_cpg <- cpgSelection(colnames(train_X))
# boLasso(train_X, train_Y, 10)
# tsdSpecificCpG(colnames(train_X))
# allCpG()

# Choose model
# "logistic, "lassoBinomal.cv", "lassoBinomal"
model_name = "logistic"
model <- model.function(model_name)
result <- model(train_X[,which_cpg], train_Y)

pred.glm <- predict.glm(result, newdata = data.frame(test_X[, which_cpg]) ,type = "response")
auc <- rocAUC(pred.glm, test_Y)
print("AUC")
print(auc)

# saveModel(result)
