source("./lib/utility.R")
source("./model.R")

resource_name = 'E-MTAB-4664'

idats <- readIdats(resource_name)

# To reformat matrix (barcode_name x cpg_sites)
beta_value <- t(idats@assayData[["betas"]])
m_value <- convertBeta2M(beta_value)i
X <- normalize(m_value)

attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))

Y <- convertSleepStatus2Num(attributes_data[5:6])

result <- model.lassoBinomial(X, Y)

