#' @title Genelines for log2 fold change data
#'
#' @description Calculate expression changes between treatments
#'
#' @usage genelinesLog2FC(data, condition, sig.level=0.05, FC, degs="inclusive")
#'
#' @param data A list of DESeq results objects, in the order you wish to make the comparisons.
#' @param condition A vector of the treatment group levels, in the same order they appear in data. Note that treatment names must NOT contain dashes (-), as this is used to separate treatment pairs.
#' @param sig.level The significance cutoff used to decide a differentially expressed gene. Defaults to 0.05; can be removed by setting to NULL.
#' @param FC The fold change cutoff used to decide a differentially expressed gene. Defaults to 0 (no cutoff).
#' @param degs Whether to return genes differentially expressed in ANY comparison ("inclusive") or ALL comparisons ("exclusive"). Defaults to inclusive.
#'
#' @return A list of two data frames: the first (means) contains the mean value of each treatment group relative to the first treatment group, and the second (slopes) contains the log2 fold changes between pairs of treatment levels.
#'
#' @details Treatment levels are compared in the order supplied by *data*, i.e. treatment 2 - treatment 1, treatment 3 - treatment 2. The 'means' are for compatibility with other *genelines* functions only. All genes are assigned a value of 0 for treatment 1. Genes are only included in the results if they meet the *sig.level* and *FC* criteria.
#'
#' @seealso \code{\link{graph.genelines}}
#'
#' @examples
#' genes<-list(
#' deData1=data.frame(log2FoldChange=round(rnorm(n=100, mean=100, sd=10), 0)),
#' deData2=data.frame(log2FoldChange=round(rnorm(n=100, mean=100, sd=10), 0)))
#' rownames(genes$deData1)<-paste0("g", seq(1:100))
#' rownames(genes$deData2)<-paste0("g", seq(1:100))
#' genelinesLog2FC(genes, c("control", "treat1", "treat2"))
#'
#' @export
genelinesLog2FC<-function(data, condition, sig.level=0.05, FC=0, degs="inclusive"){
  #Create function to identify DE genes
  indicateDEGS<-function(deData){
    DE<-!(deData$log2FoldChange < FC & deData$log2FoldChange > -FC)
    if(!is.null(sig.level)){
      DE[deData$padj > sig.level]<-FALSE
    }
    DE[is.na(DE)]<-FALSE
    return(DE)
  }
  #Mark DE genes in each dataset
  i<-1
  while (i <= length(data)) {
    data[[i]]$DE<-indicateDEGS(data[[i]])
    #data[[i]]<-data[[i]][data[[i]]$DE==TRUE,]
    i<-i+1
  }
  marked<-cbind(sapply(data, function(x)x$DE))
  #Combine datasets
  log2fc<-data.frame(V1=data[[1]]$log2FoldChange, row.names = rownames(data[[1]]))
  i<-2
  while (i <= length(data)) {
    log2fc[,i]<-data[[i]]$log2FoldChange[match(rownames(log2fc), rownames(data[[i]]))]
    i = i + 1
  }
  #Remove non-DE genes
  if(degs=="inclusive"){
    log2fc<-log2fc[rowSums(marked)!=0,]
  }
  else if(degs=="exclusive"){
    log2fc<-log2fc[rowSums(marked)==length(data),]
  }
  else {stop("Error: 'degs' must be either 'exclusive' or 'inclusive'.")}
  log2fc<-log2fc[!is.na(rowSums(log2fc)),]
  if(nrow(log2fc)==0){
    stop("Error: No genes are differentially expressed. Try changing the filtering criteria.")
  }
  #Initialise outputs
  changes <- list(
  means = matrix(
    NA,
    ncol = length(condition),
    nrow = nrow(log2fc),
    dimnames = list(row.names(log2fc), condition)
  ),
  slopes = as.matrix(
    log2fc,
    dimnames = list(row.names(log2fc), condition[-1])
  )
)
  #Calculate mean for each treatment
  changes$means[,1]<-0
  for (i in 1:ncol(changes$slopes)){
    changes$means[,i+1]<-changes$means[,i]+changes$slopes[,i] #1?
  }
  i<-1
  while (i <= ncol(changes$slopes)){
  colnames(changes$slopes)[i]<-paste0("C", i, "_", i+1)
  i<-i+1
}
  changes$means<-as.data.frame(changes$means)
  changes$slopes<-as.data.frame(changes$slopes)
  return(changes)
}
