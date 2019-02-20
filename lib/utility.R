# Read idat files
readIdatsObj = function(resource_name){
  rdata_path <- paste0("./rdata/", resource_name, ".rdata")
  idats <- readRDS(rdata_path)
  return(idats)
}

# Read saved model
readModel = function(model_name){
  model_path <- paste0("./models/", model_name, ".rmodel")
  model <- readRDS(model_path)
  return(model)
}

# Save model
saveModel = function(model){
  model_name <- format(Sys.time(), "%Y%m%d%H%M%S")
  model_path <- paste0("./models/", model_name, ".rmodel")
  saveRDS(result, file = model_path)
}

# Split test and train data
trainTestSplit = function(X, Y, test_size_rate){
  data_size <- length(Y)
  test_size <- as.integer(data_size*test_size_rate)

  test_Y <- Y[1:test_size]
  train_Y <- Y[(test_size+1):data_size]
  test_X <- X[1:test_size,]
  train_X <- X[(test_size+1):data_size,]

  return(list(train_Y, test_Y, train_X, test_X))
}

rocAUC = function(predicted, y){
  pred <- prediction(predicted, y)
  perf <- performance(pred, "tpr", "fpr")
  plot(perf)
  auc.tmp <- performance(pred,"auc")
  auc <- as.numeric(auc.tmp@y.values)
  return(auc)
}

vec.include = function(v, query){
  which_v <- c()
  for(q in query){
    which_v <- append(which_v, which(v == q))
  }
  return(which_v)
}

saveCoef = function(fit_model, cpg_ids){
  coef = coef(fit_model, s = fit_model$lambda.min)
  cpg_ids <- cpg_ids[which(coef[-1] != 0)]
  cpg_ids <- c("intercept",cpg_ids)
  cpg <- t(rbind(cpg_ids, coef@x))
  write.csv(cpg, "./data/coef/coef.csv")
}

closestGene = function(){
  cpg_data <- read.csv(paste0("./data/tsd_cpgs.csv"))[-1, ]
  
  cpg_id <- c()
  gene <- c()
  tss <- c()
  distance <- c()
  for(ncpg in 1:length(cpg_data[,1])){
    cpg <- cpg_data[ncpg, ]
    res <- findClosestGene(as.character(cpg$chr), cpg$position, genome="hg19", position="txStart")
    cpg_id <- append(cpg_id, as.character(cpg$cpg.id))
    gene <- append(gene, as.character(res$geneName[1]))
    tss <- append(tss, res$txStart[1])
    distance <- append(distance ,res$Distance[1])
  }
  cpg.data.frame <- data.frame(CpG.ID=cpg_id, Gene.Name=gene, TSS=tss, Distance=distance)
  write.csv(cpg.data.frame, "./data/closest_gene_from_cpg.csv")
}

ROC = function(score, actual, add=FALSE,
               col="black", col.area="", type="l", pch=16) {
  o = order(score, decreasing=TRUE)
  fp = tp = fp_prev = tp_prev = 0
  nF = sum(actual == FALSE)
  nT = sum(actual == TRUE)
  score_prev = -Inf
  ber_min = Inf
  area = 0
  rx = ry = numeric(length(o))
  n = 0
  for (i in seq_along(o)) {
    j = o[i]
    if (score[j] != score_prev) {
      area = area + (fp - fp_prev) * (tp + tp_prev) / 2
      n = n + 1
      rx[n] = fp/nF
      ry[n] = tp/nT
      ber = (fp/nF + 1 - tp/nT)/2
      if (ber < ber_min) {
        ber_min = ber
        th = score_prev
        rx_best = fp/nF
        ry_best = tp/nT
      }
      score_prev = score[j]
      fp_prev = fp
      tp_prev = tp
    }
    if (actual[j] == TRUE) {
      tp = tp + 1
    } else {
      fp = fp + 1
    }
  }
  
  area = area + (fp - fp_prev) * (tp + tp_prev) / 2
  n = n + 1
  rx[n] = fp/nF  # = 1
  ry[n] = tp/nT  # = 1
  if (!add) {
    plot(NULL, xlim=c(0,1), ylim=c(0,1), asp=1,
         xlab="False Positive", ylab="True Positive",
         xaxs="i", yaxs="i")
    abline(h=(1:9)/10, v=(1:9)/10, col=gray(0.9))
    abline(0, 1, col=gray(0.4))
    abline(h=0:1, v=0:1)
  }
  t = (rx_best + ry_best)/2
  abline(ry_best-rx_best, 1, col=gray(0.8))
  lines(c(rx_best, t), c(ry_best, t), col=gray(0.8))
  if (col.area != "") {
    polygon(c(rx[1:n],1), c(ry[1:n],0), col=col.area)
  }
  lines(rx[1:n], ry[1:n], type=type, lwd=2, col=col, pch=pch, xpd=TRUE)
  cat("AUC =", area/(nF*nT), "th =", th, "\n")
  cat("BER =", (rx_best + (1-ry_best))/2,
      "OR =", (ry_best/(1-ry_best))/(rx_best/(1-rx_best)), "\n")
  print(table(score >= th, actual, dnn=c("Predicted","Actual")))
  invisible(list(rx=rx, ry=ry, AUC=area/(nF*nT), th=th))
}