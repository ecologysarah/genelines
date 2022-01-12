#' @title Direction Changes
#'
#' @description Find genes with a given direction of expression change
#'
#' @usage change.direction(changes, directions, mode="tol", tolerance=0.1)
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param directions A named list of the desired directions for each comparison. Names should correspond to the column names in \emph{changes$slopes}. Can contain the values "up", "down" and "flat".
#' @param mode How to decide whether a gene should be assigned to a particular direction. Defaults to "tol", whereby the parameter \emph{tolerance} is used to classify directions. Alternatively, "DE" counts a gene as "flat" only if it was not differentially expressed. Option "DE" will only work if \emph{changes} was produced by \code{\link{genelinesLOG2FC}} with \emph{indicate.DE=TRUE}.
#' @param tolerance A single number corresponding to how much above/below each value in \emph{directions} a gene may be and still be returned: e.g. a tolerance of 0.1 means that genes with slopes -0.1-0.1 are regarded as "flat", > 0.1 as "up" and < -0.1 as "down". Only used if \emph{mode = "tol"}. Defaults to 0.1.
#'
#' @return A list of two data frames: the first contains the mean value of each treatment group for the scaled and centred counts of each gene, and the second contains the differences between pairs of treatment levels.
#'
#' @details The differences between pairs of treatment levels is equivalent to the slope of a line drawn between the treatment means. \emph{change.direction} allows the user to define a set of relationships between treatments and find which genes follow these relationships. The amount of deviation allowed is determined by \emph{tolerance}, i.e. a tolerance of 0.2 means that a slope of > 0.2 will count as "up", a slope of < -0.2 will count as "down", and a slope between -0.2 to 0.2 (inclusive) will count as "flat".
#'
#' @seealso \code{\link{change.pattern}}
#' @seealso \code{\link{change.similar}}
#' @seealso \code{\link{combined.directions}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0),
#' ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' changes<-genelines(genes, treats)
#' change.direction(changes, directions = list(C1_2 = "up"), tolerance = 0.2)
#' change.direction(changes, directions = list(C1_2 = "down"), tolerance = 0.2)
#' change.direction(changes, directions = list(C1_2 = "flat"), tolerance = 0.2)
#' @import data.table
#'
#' @export
change.direction<-function(changes, directions, mode = "tol", tolerance=0.1){
  if(!mode %in% c("tol", "DE")){
    stop("Error: mode must be either 'tol' or 'DE'.")
  }
  if(length(directions) != ncol(changes$slopes)){
    stop("Error: number of directions must equal the number of comparisons.")
  }

  if (mode=="tol"){
    for (i in 1:length(directions)){
      if (directions[[i]] == "up") {slopeno<-rownames(changes$slopes)[`>`(changes$slopes[,names(directions[i])], 0+tolerance)]}
      else if (directions[[i]] == "down") {slopeno<-rownames(changes$slopes)[`<`(changes$slopes[,names(directions[i])], 0-tolerance)]}
      else if (directions[[i]] == "flat") {slopeno<-rownames(changes$slopes)[between(changes$slopes[,names(directions[i])], 0-tolerance, 0+tolerance)]}
      else {stop("Error: direction must be a list whose values are either up, down, or flat.")}
      changes$means<-changes$means[rownames(changes$means) %in% slopeno,]
      changes$slopes<-changes$slopes[rownames(changes$slopes) %in% slopeno,]
  }
}
  else {
    if(is.null(changes$DE)){
      stop("Error: object 'changes' does not have DE information.")
    }
    for (i in 1:length(directions)){
      if (directions[[i]] == "up") {
        slopeno<-rownames(changes$slopes)[changes$DE[,names(directions[i])] == TRUE]
        slopeno<-slopeno[changes$slopes[match(slopeno, rownames(changes$slopes)),names(directions[i])] > 0]}
      else if (directions[[i]] == "down") {
        slopeno<-rownames(changes$slopes)[changes$DE[,names(directions[i])] == TRUE]
        slopeno<-slopeno[changes$slopes[match(slopeno, rownames(changes$slopes)),names(directions[i])] < 0]}
      else if (directions[[i]] == "flat") {slopeno<-rownames(changes$slopes)[changes$DE[,names(directions[i])] == FALSE]}
      else {stop("Error: direction must be a list whose values are either up, down, or flat.")}
      changes$means<-changes$means[rownames(changes$means) %in% slopeno,]
      changes$slopes<-changes$slopes[rownames(changes$slopes) %in% slopeno,]
  }
}
  return(changes)
}
