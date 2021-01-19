#' @title geneclust
#'
#' @description Perform cluster analysis on genelines output.
#'
#' @usage geneclust(changes, k, cor.method = "pearson", clust.method = "complete")
#'
#' @param changes A list as returned by \code{\link{genelines}} or \code{\link{genelinesLog2FC}}.
#' @param k The desired number of clusters.
#' @param cor.method The correlation method for the distance matrix, as passed to \code{\link[stats]{cor}}.Defaults to "pearson".
#' @param clust.method The clustering method, as passed to \code{\link[stats]{hclust}}.Defaults to "complete".
#'
#' @return A list of cluster membership for each gene and a plot of each cluster.
#'
#' @seealso \code{\link{genelines}}
#' @seealso \code{\link{genelinesLog2FC}}
#'
#' @examples
#' genes<-matrix(round(rnorm(n=600, mean=100, sd=10), 0), ncol = 6, dimnames = list(paste0("g", seq(1:100))))
#' treats<-as.factor(c(rep("control",3), rep("treated",3)))
#' changes<-genelines(genes, treats)
#' geneclust(changes, k=3)
#'
#' @export
geneclust<-function(changes, k, cor.method = "pearson", clust.method = "complete"){
  #Hierachical clustering
  d<-as.dist(1-cor(t(changes$means),  method = cor.method))
  clust<-hclust(d, method = clust.method)
  #Cut the tree
  groups<-cutree(clust, k=k)

  #list which genes are in each group
  #graph the genelines for each group
  par(mfrow = c(3, ceiling(k/3)))
  for(h in 1:k){
    genegroup<-names(groups)[groups==h]
    subs<-changes$means[rownames(changes$means) %in% genegroup,]
    plot(1, type="n", xlab="", ylab="", ylim=c(min(changes$means), max(changes$means)), xlim=c(1, ncol(changes$means)), xaxt = "n", xaxs="i", yaxs="i")
    axis(1, at=seq(1,ncol(changes$means),1), labels=names(changes$means))
    for (i in 1:nrow(subs)){
      lines(as.numeric(subs[i,]))
    }
  }
    par(mfrow = c(1, 1))
    return(groups)
}
