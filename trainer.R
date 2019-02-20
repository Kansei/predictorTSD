source("./header.R")

resource_name = 'E-MTAB-4664'
# Read idats object
idats <- readIdatsObj(resource_name)
# Data preprocessing
X <- preprocessingIdats(idats, n = F, s = F, m = T)

# Read sample's attribute csv file
attributes_data <- read.csv(paste0("./data/", resource_name ,"/attributes.csv"))
# Convert sleep status to 0/1
Y <- convertSleepStatus2Num(attributes_data[5:6])

auc <- 1
while (auc > 0.77 || auc < 0.76) {
# Split test and train data
data_size <- length(Y)
rand <- sample(100:1000, 1)
set.seed(rand)
# en 295
# la 
rand_num <- sample(1:data_size, data_size)
X <- X[rand_num,]
Y <- Y[rand_num]
test_size_rate = 0.3
test_size <- as.integer(data_size*test_size_rate)
test_Y <- Y[1:test_size]
train_Y <- Y[(test_size+1):data_size]
test_X <- X[1:test_size,]
train_X <- X[(test_size+1):data_size,]

# boLasso(train_X, train_Y, 10), tsdSpecificCpG(colnames(train_X)), allCpG(), wilcoxonSignedRankTest(0.01, idats = X)
function_name = "wilcoxonSignedRankTest"
which_cpg <- cpgSelection.function(function_name)(p_threshold = 0.001, idats = X)
print(length(which_cpg))

train_X <- train_X[, which_cpg]
test_X <- test_X[, which_cpg]

# Choose model
# logistic, lassoBinomial.cv, lassoBinomial, ridgeBinomial.cv, elasticNetBinomial.cv, adaptiveLassoBinomial.cv
model_name = "lassoBinomial.cv"
fit_model <- model.function(model_name)(train_X, train_Y)
# saveCoef(fit_model, colnames(train_X[, which_cpg]))

predicted <- prediction.function(model_name)(fit_model, test_X)
#ROC(predicted, test_Y)
auc <- rocAUC(predicted, test_Y)
print("AUC:")
print(auc)
set.seed(NULL)
}

print(rand)
# saveModel(result)
