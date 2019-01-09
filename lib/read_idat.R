suppressPackageStartupMessages(require('methylumi'))

idatPath <- "../data/E-MTAB-4664"

barcodes <- scan("../data/E-MTAB-4664/barcodes.txt", what = character(), sep = "\n")

idats <- methylumIDAT(barcodes = barcodes, idatPath=idatPath)
