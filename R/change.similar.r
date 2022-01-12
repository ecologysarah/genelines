#' @title Similarity Changes
#'
#' @description Find genes with an expression pattern similar to a gene of interest
#'
#' @usage change.similar(changes, gene, tolerance)
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param gene A character vector (in quotes) corresponding to one of the gene names in \emph{rownames(changes$slopes)}.
#' @param tolerance A single number corresponding to how much above/below the value for \emph{gene} another gene may be and still be returned.
#'
#' @return A list of two data frames: the first contains the mean value of each treatment group for the scaled and centred counts of each gene, and the second contains the differences between pairs of treatment levels.
#'
#' @details The differences between pairs of treatment levels is equivalent to the slope of a line drawn between the treatment means. \emph{change.similar} allows the user to choose a gene of interest and find which genes follow the same expression pattern. The amount of deviation allowed from the target gene is determined by \emph{tolerance}.
#'
#' @seealso \code{\link{change.pattern}}
#' @seealso \code{\link{change.direction}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0),
#' ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' changes<-genelines(genes, treats)
#' change.similar(changes, gene = "g12", tolerance = 0.2)
#'
#' @export
change.similar<-function(changes, gene, tolerance){
  for (i in 1:ncol(changes$slopes)){
    target<-changes$slopes[rownames(changes$slopes)==gene,i]
    slopeno<-rownames(changes$slopes)[changes$slopes[,i] >= target-tolerance & changes$slopes[,i] <= target+tolerance]
  changes$means<-changes$means[rownames(changes$means) %in% slopeno,]
  changes$slopes<-changes$slopes[rownames(changes$slopes) %in% slopeno,]
  }
  return(changes)

}
