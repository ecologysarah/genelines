#' @title graph.genelines
#'
#' @description Plotting function for genelines
#'
#' @usage graph.genelines(changes, ...)
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param ... Additional parameters to be passed on to \code{\link{lines}}
#' @return A plot showing the change in mean value of each treatment group for the scaled and centred counts of each gene.
#'
#' @seealso \code{\link{genelines}}
#' @seealso \code{\link{genelinesLog2FC}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0),
#' ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' changes<-genelines(genes, treats)
#' graph.genelines(changes)
#' @import graphics
#'
#'@export
graph.genelines<-function(changes, ...){
  plot(1, type="n", xlab="", ylab="Relative expression", ylim=c(min(changes$means), max(changes$means)), xlim=c(1, ncol(changes$means)), xaxt = "n", xaxs="i", yaxs="i")
  axis(1, at=seq(1,ncol(changes$means),1), labels=names(changes$means))
  for (i in 1:nrow(changes$means)){
    lines(x=seq(1,ncol(changes$means),1), y=as.numeric(changes$means[i,]), ...)
  }
}
