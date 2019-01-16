#' Create normalization factors.
#' 
#' Using either raw data from BeadStudio or \code{MethyLumiSet} data, this function creates normalization factors using either all subjects' means or a subset. 
#' 
#' @param controldata A \code{data.frame} containing the red and green color channels for all control samples in raw format output from BeadStudio. May also contain detection p-values for each subject. Must contain a column called TargetID, which contains the CpG site identifiers. Either this or methylumidata must be supplied. 
#' @param subjects Optional. User can specify the index or names of a given subject or range of subjects to use in the creation of the normalization factors as oppposed to using all the samples.
#' @param methylumidata A \code{MethyLumiSet} object containing control data (which can be accessed using the \code{normctls} function). See the vignette for this package for an example using this type of data. 
#' @param type One of either "raw" (default) or "methylumi" indicating the type of data supplied by the user. 
#' @export
#' @import methylumi stats Biobase
#' @return A list of length two containing the normalization factors for each subject in each color channel. 
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
#'
#' normfactors <- norm_factors(controldata=controldata.s, subjects=NULL, type = "raw")
#' normfactors <- norm_factors(controldata=controldata.s, subjects=1, type = "raw")

norm_factors <- function(controldata, subjects = NULL, methylumidata = NULL, type = "raw"){
  if(is.null(controldata) & is.null(methylumidata)) stop("Must supply some control data")
  if(is.null(controldata) ==FALSE & type %in% "methylumi") stop("If using raw control data, type must be 'raw'")
  if(is.null(methylumidata) ==FALSE & type %in% "raw") stop("If using MethyLumiSet data, type must be 'methylumi'")
  if(type %in% "methylumi"){
  		if (is.null(subjects)){
  			controls <- normctls(methylumidata)
  			green_means <- mean(controls$Cy3)
			  red_means <- mean(controls$Cy5)
  		} else if(is.null(subjects) %in% FALSE) {
  			controls <- normctls(methylumidata)
  			green_means <- mean(controls$Cy3[subjects, drop = FALSE])
			  red_means <- mean(controls$Cy5[subjects, drop = FALSE])
  	}	
   green_ind_means <- colMeans(controls$Cy3)
   red_ind_means <- colMeans(controls$Cy5)
   green_factors <- green_ind_means/green_means
   red_factors <- red_ind_means/red_means
   normfactors_means <- list(red_factors, green_factors)
   names(normfactors_means) <- c("Red", "Green")
   return(normfactors_means)
  }
  else if(type %in% "raw"){
  	if (is.null(subjects)){
  		green_columns <- grep("Grn", names(controldata))
  		red_columns <- grep("Red", names(controldata))
  		red_ind_means <- apply(controldata[controldata$TargetID %in% c("NORM_A", "NORM_T"), red_columns], 2, mean)
  		green_ind_means <- apply(controldata[controldata$TargetID %in% c("NORM_G", "NORM_C"), green_columns], 2, mean)
  		red_means <- mean(as.matrix(controldata[controldata$TargetID %in% c("NORM_A", "NORM_T"), red_columns]))
 		green_means <- mean(as.matrix(controldata[controldata$TargetID %in%c("NORM_G", "NORM_C"), green_columns]))
  }
  else if(is.null(subjects) %in% FALSE){
  		green_columns <- grep("Grn", names(controldata))
		red_columns <- grep("Red", names(controldata))
  		green_columns_s <- green_columns[subjects, drop = FALSE]
  		red_columns_s <- red_columns[subjects, drop = FALSE]
  		red_ind_means <- apply(controldata[controldata$TargetID %in% c("NORM_A", "NORM_T"), red_columns], 2, mean)
  		green_ind_means <- apply(controldata[controldata$TargetID %in% c("NORM_G", "NORM_C"), green_columns], 2, mean)
 		red_means <- mean(as.matrix(controldata[controldata$TargetID %in% c("NORM_A", "NORM_T"), red_columns_s]))
 		green_means <- mean(as.matrix(controldata[controldata$TargetID %in% c("NORM_G", "NORM_C"), green_columns_s]))
 		}
  green_factors <- green_ind_means/green_means
  red_factors <- red_ind_means/red_means
  normfactors_means <- list(red_factors, green_factors)
  names(normfactors_means) <- c("Red", "Green")
  return(normfactors_means)
 }
}