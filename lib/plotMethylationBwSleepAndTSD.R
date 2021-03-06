source("./lib/utility.R")

plotMvalue = function(args){
  resource_name = 'E-MTAB-4664'

  idats <- preprocessingIdats(readIdatsObj(resource_name))
  all_cpg_ids <- colnames(idats)
  all_assay_names <- rownames(idats)

  if(args[1] == "rand"){
    dir <- "rand/"
    rand_cpgs <- sample(1:length(all_cpg_ids), as.numeric(args[2]))
    plot_cpg_ids <- c()
    for(cpg in rand_cpgs){
      plot_cpg_ids <- append(plot_cpg_ids, all_cpg_ids[cpg])
    }
  } else if(args[1] == "tsd-specific"){
    dir <- "beta/standardize/"
    plot_cpg_ids <- scan("./data/DifferentiallyMethylatedProbes.txt", what = character(), sep = "\n")
  } else {
    dir <- "custom/"
    plot_cpg_ids <- args
  }

  which_tsd_cpgs <- vec.include(all_cpg_ids, plot_cpg_ids)

  plot_cpgs <- t(rbind(plot_cpg_ids, which_tsd_cpgs))

  attributes_data <- read.csv(paste("./data/", resource_name ,"/attributes.csv", sep = ""))[4:6]
  which_tsd_sample <- which(attributes_data$environmental.stress == "total acute sleep deprivation")

  tsd_samples <- attributes_data[which_tsd_sample, ]
  tsd_samples <- tsd_samples[order(tsd_samples$individual),]
  sleep_samples <- attributes_data[-which_tsd_sample, ]
  sleep_samples <- sleep_samples[order(sleep_samples$individual),]

  except_tsd_samples <- c("9297962119_R04C01", "9297962119_R03C01", "9297962042_R04C01", "9297962042_R03C01")
  except_sleep_samples <- c("9297962119_R02C02", "9297962119_R02C01", "9297962042_R02C02", "9297962042_R02C01")

  tsd_samples <- tsd_samples[-vec.include(tsd_samples$Assay.Name, except_tsd_samples), ]
  sleep_samples <- sleep_samples[-vec.include(sleep_samples$Assay.Name, except_sleep_samples), ]

  for(ncpg in 1:nrow(plot_cpgs)){
    cpg_id <- as.character(plot_cpgs[ncpg, 1])
    which_cpg <- as.numeric(plot_cpgs[ncpg, 2])

    max_m <- max(idats[, which_cpg])
    min_m <- min(idats[, which_cpg])

    png(filename=paste0("./data/plot/methylation_bw_sleep_and_TSD/", dir, cpg_id, ".png"))
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
}
