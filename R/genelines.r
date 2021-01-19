#' @title Genelines
#'
#' @description Calculate changes between treatments
#'
#' @usage genelines(obj, condition, comparisons = NULL)
#'
#' @param obj A data frame or matrix of normalised gene counts; genes are rows, samples are columns. Row names should correspond to Gene IDs.
#' @param condition A list of the treatment group for each column in obj.
#' @param comparisons A list of desired comparisons, in the format \emph{list("1-2", "2-3", "1-3")}. This would be equivalent to treatment 2 - treatment 1, treatment 3 - treatment 2, and treatment 3 - treatment 1. If NULL (the default), comparisions are performed between successive treatment levels in \emph{condition}. Note that treatment names must NOT contain dashes (-), as this is used to separate treatment pairs.
#'
#' @return A list of two data frames: the first (means) contains the mean value of each treatment group for the scaled and centred counts of each gene, and the second (slopes) contains the differences between pairs of treatment levels.
#'
#' @details Note that although features are referred to as genes for convenience, this function is best suited for non-gene data. If you have RNAseq data you are advised to use \code{\link{genelinesLog2FC}}.
#' Treatment levels are compared in sequential pairs, i.e. treatment 2 - treatment 1, treatment 3 - treatment 2, unless otherwise specified by *comparisons*. This is equivalent to the slope of a line drawn between the treatment means.
#'
#' @seealso \code{\link{graph.genelines}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0), ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' genelines(genes, treats)
#'
#' @export
genelines<-function(obj, condition, comparisons = NULL){
  condition<-as.factor(condition)
  #Initialise outputs
  if (!is.null(comparisons)) {
    changes <- list(
      means = matrix(
        NA,
        ncol = nlevels(condition),
        nrow = nrow(obj),
        dimnames = list(row.names(obj), levels(condition))
      ),
      slopes = matrix(
        NA,
        ncol = length(comparisons),
        nrow = nrow(obj),
        dimnames = list(row.names(obj), rep(NA, length(comparisons)))
      )
    )
  }
  else {
    changes <- list(
      means = matrix(
        NA,
        ncol = nlevels(condition),
        nrow = nrow(obj),
        dimnames = list(row.names(obj), levels(condition))
      ),
      slopes = matrix(
        NA,
        ncol = nlevels(condition) - 1,
        nrow = nrow(obj),
        dimnames = list(row.names(obj), rep(NA, nlevels(condition) - 1))
      )
    )
  }
  #Scale and centre genes
  obj<-scale(t(obj))
  #Calculate mean for each treatment
  for (i in 1:ncol(obj)){
    agg<-aggregate(obj[,i], list(condition), mean)
    changes$means[i,]<-agg$x
  }
  #Calculate change between treatments
  if (!is.null(comparisons)) {
    complevels <- vector()
    for (j in 1:length(comparisons)) {
      c1 <- gsub("(.+)-(.+)", "\\1", comparisons[j])
      c2 <- gsub("(.+)-(.+)", "\\2", comparisons[j])
      changes$slopes[, j] <- changes$means[, c2] - changes$means[,
                                                                 c1]
      colnames(changes$slopes)[j] <- paste("C", j, "_",
                                           j + 1, sep = "")
      complevels <- c(complevels, c1, c2)
    }
  }
  else{
    for (j in 1:nlevels(condition) - 1) {
      colnames(changes$slopes)[j] <- paste("C", j, "_", j + 1, sep = "")
      changes$slopes[, j] <- changes$means[, j + 1] - changes$means[, j]
    }
  }
  changes$means<-as.data.frame(changes$means)
  changes$slopes<-as.data.frame(changes$slopes)
  if (exists("complevels")) {
    changes$means <- changes$means[, names(changes$means) %in%
                                     complevels]
  }
  return(changes)
}
