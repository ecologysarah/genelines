#' @title Pattern Changes
#'
#' @description Find genes with a particular expression pattern
#'
#' @usage change.pattern(changes, pattern, tolerance)
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param pattern A named list of the desired values for each comparison. Names should correspond to the column names in \emph{changes$slopes}.
#' @param tolerance A single number corresponding to how much above/below each value in \emph{pattern} a gene may be and still be returned.
#'
#' @return A list of two data frames: the first contains the mean value of each treatment group for the scaled and centred counts of each gene, and the second contains the differences between pairs of treatment levels.
#'
#' @details The differences between pairs of treatment levels is equivalent to the slope of a line drawn between the treatment means. \emph{change.pattern} allows the user to define a set of numerical relationships between treatments and find which genes follow this pattern. The amount of deviation allowed from the specified pattern is determined by \emph{tolerance}.Note that if \emph{changes was produced by  \code{\link{genelinesLog2FC}}}, the slopes correspond to log2 fold changes.
#'
#' @seealso \code{\link{change.direction}}
#' @seealso \code{\link{change.similar}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0),
#' ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' changes<-genelines(genes, treats)
#' change.pattern(changes, pattern = list(C1_2 = 0.5), tolerance = 0.2)
#'
#' @export
change.pattern<-function(changes, pattern, tolerance){
  if(length(pattern) != ncol(changes$slopes)){
    stop("Error: number of patterns must equal the number of comparisons.")
  }
  for (i in 1:length(pattern)){
    slopeno<-rownames(changes$slopes)[changes$slopes[,names(pattern[i])] >= (pattern[[i]]-tolerance) & changes$slopes[,names(pattern[i])] <= (pattern[[i]]+tolerance)]
  changes$means<-changes$means[rownames(changes$means) %in% slopeno,]
  changes$slopes<-changes$slopes[rownames(changes$slopes) %in% slopeno,]
  }
  return(changes)
}
