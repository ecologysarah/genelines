#' @title Direction Changes
#'
#' @description Find and plot genes with all possible expression patterns for two comparisons
#'
#' @usage combined.directions(changes, tolerance, label.plots=NULL)
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param tolerance A single number corresponding to how much above/below each value in \code{\link{change.direction}} a gene may be and still be returned.
#' @param label.plots A single keyword for position (as supplied to \code{\link[graphics]{legend}}) for labelling the plots with the number of genes in each group. Defaults to NULL.
#'
#' @return A list of eight, containing the genes that follow the patterns flat-up,   flat-down, up-up, up-flat, up-down,  down-up, down-flat, down-down. Within each item are two data frames: the first contains the mean value of each treatment group for each gene, and the second contains the differences between pairs of treatment levels.
#'
#' @details The differences between pairs of treatment levels is equivalent to the slope of a line drawn between the treatment means. \code{\link{change.direction}} allows the user to define a single set of relationships between treatments and find which genes follow this pattern. The amount of deviation allowed from the specified pattern is determined by \emph{tolerance}.
#'
#' @seealso \code{\link{genelines}}
#' @seealso \code{\link{genelinesLog2FC}}
#' @seealso \code{\link{change.direction}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=900, mean=100, sd=10), 0),
#' ncol = 9, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treat1",3), rep("treat2",3)))
#' changes<-genelines(genes, treats)
#' combined.directions(changes, tolerance = 0.2, label.plots="center")
#' @import graphics
#'
#'@export
combined.directions<-function(changes, tolerance, label.plots=NULL){
  if(ncol(changes$slopes) != 2){
    stop("all.directions is currently implemented for two comparisons only.")
  }
  genesets<-list(
  FlatUp=change.direction(changes, list(C1_2 = "flat", C2_3 = "up"), tolerance),
  FlatDown=change.direction(changes, list(C1_2 = "flat", C2_3 = "down"), tolerance),
  UpUp=change.direction(changes, list(C1_2 = "up", C2_3 = "up"), tolerance),
  UpFlat=change.direction(changes, list(C1_2 = "up", C2_3 = "flat"), tolerance),
  UpDown=change.direction(changes, list(C1_2 = "up", C2_3 = "down"), tolerance),
  DownUp=change.direction(changes, list(C1_2 = "down", C2_3 = "up"), tolerance),
  DownFlat=change.direction(changes, list(C1_2 = "down", C2_3 = "flat"), tolerance),
  DownDown=change.direction(changes, list(C1_2 = "down", C2_3 = "down"), tolerance)
  )

  present<-sapply(genesets, function(x)nrow(x$means)>0)
  #plot the genelines for each group with > 0 genes
  par(mfrow = c(3,ceiling(sum(present)/3)))
  for(i in names(genesets)[present]){
    graph.genelines(genesets[[i]])
    mtext(i, side=3)
    if (!is.null(label.plots)){
      legend(label.plots, legend = paste("n =", nrow(genesets[[i]]$means)), bty="n")
    }
  }
  par(mfrow = c(1, 1))
  return(genesets)
}
