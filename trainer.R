source("./lib/utility.R")
source("./model.R")

resource_name = 'E-MTAB-4664'
idats <- readIdats(resource_name)

# To reformat matrix (barcode_name x cpg_sites)
norm_idats <- normalize(idats)
beta_value <- norm_idats@assayData[["betas"]]
m_value <- convertBeta2M(beta_value)
excepted_m_value <- exceptMissingValue(m_value)
X <- t(excepted_m_value)

attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))
Y <- convertSleepStatus2Num(attributes_data[5:6])

result <- model.lassoBinomial.cv(X, Y)

coef(result, s="lambda.min")