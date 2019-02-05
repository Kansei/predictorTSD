source("./lib/utility.R")

resource_name = 'E-MTAB-4664'
idats <- preprocessingIdats(readIdatsObj(resource_name))
all_cpg_ids <- colnames(idats)
all_assay_names <- rownames(idats)

tsd_specific_cpg_ids <- scan("./data/DifferentiallyMethylatedProbes.txt", what = character(), sep = "\n")

### CpGの部位とIDの行列を作成
which_tsd_cpgs <- c()
for(cpg_id in tsd_specific_cpg_ids){
  which_tsd_cpgs <- append(which_tsd_cpgs, which(all_cpg_ids == cpg_id))
}
tsd_specific_cpgs <- t(rbind(tsd_specific_cpg_ids, which_tsd_cpgs))

## サンプルをsleepとTSDで分ける
attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))[4:6]
which_tsd_sample <- which(attributes_data$environmental.stress == "total acute sleep deprivation")

tsd_samples <- attributes_data[which_tsd_sample, ]
tsd_samples <- tsd_samples[order(tsd_samples$individual),]
sleep_samples <- attributes_data[-which_tsd_sample, ]
sleep_samples <- sleep_samples[order(sleep_samples$individual),]

## cpgでfor文を回す
for(ncpg in 1:nrow(tsd_specific_cpgs)){
  cpg_id <- as.character(tsd_specific_cpgs[ncpg, 1])
  which_cpg <- as.numeric(tsd_specific_cpgs[ncpg, 2])

  max_m <- max(idats[, which_cpg])
  min_m <- min(idats[, which_cpg])

  png(filename=paste0("./data/plot/methylation_bw_sleep_and_TSD/", cpg_id, ".png"))
#
  plot(c(0,1), c(min_m ,max_m),
    type="n",
    ylim=c(min_m, max_m),
    xaxp=c(-1,1,2),
    xlab="TSD or Sleep", ylab="M-value",
    main = cpg_id
  )
  for(nsample in 1:nrow(tsd_samples)){
   tsd_sample <- tsd_samples$Assay.Name[nsample]
   sleep_sample <- sleep_samples$Assay.Name[nsample]

   tsd_m <- idats[which(all_assay_names == tsd_sample), which_cpg]
   sleep_m <- idats[which(all_assay_names == sleep_sample), which_cpg]
   points(c(0,1), c(sleep_m, tsd_m), type="o")
  }
  dev.off()
}
