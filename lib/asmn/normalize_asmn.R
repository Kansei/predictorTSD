# @keywords internal
# @export
multiplybygreen <-
  function(onerow, normfactors){
    green.factors <- normfactors[[2]]
    out <- onerow*green.factors
    return(out) 
  }

# @keywords internal
# @export
multiplybyred <-
  function(onerow, normfactors){
    red.factors <- normfactors[[1]]
    out <- onerow*red.factors
    return(out) 
  }


#' Perform all sample mean normalization.
#' 
#' This function normalizes either raw data output from BeadStudio or data of type \code{MethyLumiSet}.
#' 
#' @param normfactors The output from \code{\link{norm_factors}}
#' @param rawdata A \code{data.frame}, or \code{matrix} containing the methylation values the user wishes to normalize. Must contain column names that specify the detection signal type (A or B) and have the number of rows equal the number of CpG sites to be analyzed.
#' @param featuredata A \code{data.frame} containing information for each row of \code{rawdata} specifying the ID of the CpG sites, the Infinium design type (I or II), and the color channel (red or green). This determines which normalization factors are used. 
#' @param methylumidata A \code{MethyLumiSet} data object containing the methylated/unmethylated sites as well as color channel identifiers accessed using \code{fData}.  See the vignette for this package for an example using this type of data. 
#' @param type One of either "raw" (default) or "methylumi" indicating the type of data supplied by the user. 
#' @export
#' @import methylumi stats Biobase
#' @examples
#' ids <- seq(from = 10, to = 59)
#' n <- length(ids)
#' reds.s <- data.frame(matrix(rep(round(abs(rnorm(n, 4, 1))), 30), nrow=30))
#' names(reds.s) <- paste("X", ids, ".Signal_Red", sep = "")
#' greens.s <- data.frame(matrix(rep(round(abs(rnorm(n, 5, 2))), 30), nrow = 30))
#' names(greens.s) <- paste("X", ids, ".Signal_Grn", sep = "")
#' dat <- data.frame(reds.s, greens.s)
#' dat <- dat[,order(names(dat))]
#' indices <- sample(1:30, 30, replace = FALSE)
#' TargetID.s <- rep(NA, 30)
#' TargetID.s[indices[1:9]] <- "NORM_A"
#' TargetID.s[indices[10:18]] <- "NORM_T"
#' TargetID.s[indices[19:25]] <- "NORM_C"
#' TargetID.s[indices[25:30]] <- "NORM_G"
#' controldata.s <- data.frame(TargetID = TargetID.s, dat)
#' normfactors <- norm_factors(controldata=controldata.s, subjects=NULL)
#' ncpg <- 100
#' IlmnIDs.s <- paste("cg00", seq(1:ncpg), sep = "") 
#' Infinium_Design_type.s <- sample(c("I", "II"), size=ncpg, replace = TRUE)
#' Color_Channel.s <- vector()
#' Color_Channel.s[Infinium_Design_type.s == "I"] <- sample(c("Red", "Grn"), size = length(Color_Channel.s[Infinium_Design_type.s == "I"]), replace = TRUE)
#' featuredata.s <- data.frame(IlmnIDs = IlmnIDs.s, Infinium_Design_Type = Infinium_Design_type.s, Color_Channel = Color_Channel.s)
#' signalA.s <- data.frame(matrix(rep(round(abs(rnorm(n, 4000, 2000))), ncpg), nrow=ncpg))
#' names(signalA.s) <- paste("X", ids, ".Signal_A", sep = "")
#' signalB.s <- data.frame(matrix(rep(round(abs(rnorm(n, 5000, 2000))), ncpg), nrow = ncpg))
#' names(signalB.s) <- paste("X", ids, ".Signal_B", sep = "")
#' dat <- data.frame(signalA.s, signalB.s)
#' mydata.s <- dat[,order(names(dat))]
#' 
#' newbetas <- normalize_asmn(normfactors=normfactors, rawdata=mydata.s, featuredata=featuredata.s, methylumidata = NULL, type = "raw")

normalize_asmn <- function(normfactors, rawdata, featuredata, methylumidata = NULL, type = "raw"){
     if(is.null(rawdata) & is.null(methylumidata)) stop("Must supply some control data")
  if(is.null(rawdata) ==FALSE & type %in% "methylumi") stop("If using raw control data, type must be 'raw'")
  if(is.null(methylumidata) ==FALSE & type %in% "raq") stop("If using MethyLumiM data, type must be 'methylumi'")
  if(is.null(normfactors)) stop("Must supply normalization factors")
    if(type %in% "raw"){
        signalA <- grep("Signal_A", colnames(rawdata))
        signalB <- grep("Signal_B", colnames(rawdata))
        signals <- colnames(rawdata)[sort(c(signalA, signalB))]
        ids <- as.character(featuredata$IlmnID)
        bs2 <- grep("II", featuredata$Infinium_Design_Type)
        gr <- grep("Grn", featuredata$Color_Channel)
        rd <- grep("Red", featuredata$Color_Channel)
        ids.reorder <- c(ids[bs2], ids[rd], ids[gr])
        cat("Subsetting data", "\n")
        bs2sigA <- as.matrix(rawdata[bs2, signalA])
        bs2sigB <- rawdata[bs2, signalB]
        rdsigA <- rawdata[rd, signalA]
        rdsigB <- rawdata[rd, signalB]
        grsigA <- rawdata[gr, signalA]
        grsigB <- rawdata[gr, signalB]
        cat("Multiplying factors", "\n")
        bs2sigA.norm <- apply(bs2sigA, 1, multiplybyred, normfactors = normfactors)
        bs2sigB.norm <- apply(bs2sigB, 1, multiplybygreen, normfactor = normfactors)
        rdsigA.norm <- apply(rdsigA, 1, multiplybyred, normfactors = normfactors)
        rdsigB.norm <- apply(rdsigB, 1, multiplybyred, normfactors = normfactors)
        grsigA.norm <- apply(grsigA, 1, multiplybygreen, normfactors = normfactors)
        grsigB.norm <- apply(grsigB, 1, multiplybygreen, normfactors = normfactors)
        sigAs <- cbind(bs2sigA.norm, rdsigA.norm, grsigA.norm)
        sigBs <- cbind(bs2sigB.norm, rdsigB.norm, grsigB.norm)
        cat("Normalizing", "\n")
        sigAs.n <- apply(sigAs, c(1, 2), function(z) {
            max(0, z)
        })
        sigBs.n <- apply(sigBs, c(1, 2), function(z) {
            max(0, z)
        })
        norm.betas <- data.frame(t(sigBs.n/(sigAs.n + sigBs.n + 100)))
        rownames(norm.betas) <- ids.reorder
        cat("Arranging data", "\n")
        norm.betas <- norm.betas[order(ids.reorder),]
        return(norm.betas)
    }
    if(type %in% "methylumi"){
      Red.factor <- normfactors$Red
      Grn.factor <- normfactors$Green
      Grn <- list( M1=methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Grn'), ],
                   U1=unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Grn'), ],
                   M2=methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Both'), ] )
      Grn <- lapply(Grn, function(y) sweep(y, 2, FUN="*", Grn.factor))
      methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Grn'), ] <- Grn$M1
      unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Grn'), ] <- Grn$U1
      methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Both'), ] <- Grn$M2
      
      Red <- list( M1=methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Red'), ],
                   U1=unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Red'), ],
                   U2=unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Both'), ] )
      Red <- lapply(Red, function(y) sweep(y, 2, FUN="*", Red.factor))
      methylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Red'), ] <- Red$M1
      unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Red'), ] <- Red$U1
      unmethylated(methylumidata)[ which(fData(methylumidata)$COLOR_CHANNEL=='Both'), ] <- Red$U2
      betas(methylumidata) <- methylated(methylumidata) / total.intensity(methylumidata)
      return(methylumidata)
    }
}
