suppressPackageStartupMessages(library('methylumi'))
# Read idat files

resource_name = "E-MTAB-4664"

resource_path <- paste("./data/", resource_name, sep = "")
idatPath <- paste(resource_path, "/idat", sep = "")
barcodes <- scan(paste(resource_path, "/barcodes.txt", sep = ""), what = character(), sep = "\n")
idats <- methylumIDAT(barcodes = barcodes, idatPath=idatPath)

rdata_name <- paste(resource_name, ".rdata", sep = "")
rdata_path <- paste("./rdata/", rdata_name, sep = "")

saveRDS(idats, file = rdata_path)

