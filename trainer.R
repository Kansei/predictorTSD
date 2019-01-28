source("./lib/utility.R")
source("./model.R")

resource_name = 'E-MTAB-4664'
idats <- readIdatsObj(resource_name)

# To reformat matrix (barcode_name x cpg_sites)
norm_idats <- normalize(idats)
beta_value <- norm_idats@assayData[["betas"]]
m_value <- convertBeta2M(beta_value)
excepted_m_value <- exceptMissingValue(m_value)
X <- t(excepted_m_value)

attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))
Y <- convertSleepStatus2Num(attributes_data[5:6])

# split test and train data
test_size_rate <- 0.3
data_size <- length(Y)
test_size <- as.integer(data_size*test_size_rate)
test_Y <- Y[1:test_size]
train_Y <- Y[(test_size+1):data_size]
test_X <- X[1:test_size,]
train_X <- X[(test_size+1):data_size,]

result <- model.lassoBinomial.cv(train_X, train_Y)

paramater <- coef(result, s="lambda.min")
# paramater[which(paramater != 0), 1]

model_name <- paste(format(Sys.time(), "%Y%H%M%S"), ".rmodel", sep = "")
model_path <- paste("./models/", model_name, sep = "")
saveRDS(result, file = model_path)

