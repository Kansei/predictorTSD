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

except_tsd_samples <- c("9297962119_R04C01", "9297962119_R03C01", "9297962042_R04C01", "9297962042_R03C01")
which_except_tsd_sample <- c()
for(sample in except_tsd_samples){
  which_except_tsd_sample <- append(which_except_tsd_sample, which(tsd_samples$Assay.Name == sample))
}
tsd_samples <- tsd_samples[-which_except_tsd_sample, ]

except_sleep_samples <- c("9297962119_R02C02", "9297962119_R02C01", "9297962042_R02C02", "9297962042_R02C01")
which_except_sleep_sample <- c()
for(sample in except_sleep_samples){
  which_except_sleep_sample <- append(which_except_sleep_sample, which(sleep_samples$Assay.Name == sample))
}
sleep_samples <- sleep_samples[-which_except_sleep_sample, ]

## cpgでfor文を回す
for(ncpg in 1:nrow(tsd_specific_cpgs)){
  cpg_id <- as.character(tsd_specific_cpgs[ncpg, 1])
  which_cpg <- as.numeric(tsd_specific_cpgs[ncpg, 2])

  max_m <- max(idats[, which_cpg])
  min_m <- min(idats[, which_cpg])

  png(filename=paste0("./data/plot/methylation_bw_sleep_and_TSD/", cpg_id, ".png"))
  par(mai=c(.8,.8,.8,1.2))
  plot(c(0,1), c(min_m, max_m),
    type="n",
    ylim=c(min_m, max_m),
    xaxp=c(-1,1,2),
    xlab="TSD or Sleep", ylab="M-value",
    main=cpg_id
  )

  labels <- c()
  sample_range <- 1:nrow(tsd_samples)

  for(nsample in sample_range){
    tsd_sample_name <- tsd_samples$Assay.Name[nsample]
    sleep_sample_name <- sleep_samples$Assay.Name[nsample]
    labels <- append(labels, tsd_samples$individual[nsample])

    tsd_m <- idats[which(all_assay_names == tsd_sample_name), which_cpg]
    sleep_m <- idats[which(all_assay_names == sleep_sample_name), which_cpg]

    difference_m <- sleep_m - tsd_m
    if(difference_m > 0){
      color <- 1
    } else {
      color <- 2
    }

    points(c(0,1), c(sleep_m, tsd_m), type="o", pch=nsample, col=color)
  }

  par(xpd=T)
  legend(par()$usr[2], par()$usr[4], legend = labels, title = "Sample Id", col = 1, pch = sample_range, bg = "transparent")
  dev.off()
}
