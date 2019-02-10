library(exactRankTests)

wilcoxonSignedRankTest = function(p_threshold = 0.05){
  resource_name = 'E-MTAB-4664'

  idats <- preprocessingIdats(readIdatsObj(resource_name))
  cpg_ids <- colnames(idats)
  assay_names <- rownames(idats)

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

  which_tsd <- vec.include(assay_names, tsd_samples$Assay.Name)
  which_sleep <- vec.include(assay_names, sleep_samples$Assay.Name)

  which_cpg <- c()
  for(ncpg in 1:length(cpg_ids)){
    result <- wilcox.exact(x=idats[which_tsd, ncpg], y=idats[which_sleep, ncpg],paired=T)
    if(result$p.value <= p_threshold){
      which_cpg <- append(which_cpg, ncpg)
    }
  }

  return(t(rbind(cpg_ids[which_cpg], which_cpg)))
}
