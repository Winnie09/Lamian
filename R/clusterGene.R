#' Cluster genes based on their temporal patterns  of gene expression (constant time test) or the temporal patterns of gene expression group difference (sample covariate test).
#' 
#' This function is used to cluster genes based on their temporal patterns  of gene expression (constant time test) or the temporal patterns of gene expression group difference (sample covariate test).
#' 
#' @export
#' @import ggplot2 RColorBrewer splines gridExtra viridis
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param testobj object returned from lamian.test(). 
#' @param gene a character vector of gene names. It can be of length 1 or > 1.
#' @param type One of c('Time', 'Variable'). 
#' @param k.auto logical. If FALSE (default), users need to specify k as the number of clusters. If TRUE, k will be automatically determined by elbow's method.
#' @param k a numeric number. The number of clusters. Only useful when k.auto = FALSE.
#' @param method the clustering method, either "kmeans" for k-means clustering or "hierarchical" for hierarchical clustering.
#' @param scale.difference logical. If FALSE, then do not standarize the group difference, but only scale by the maximum value of the group difference absolute values. If TRUE, then standardize the group difference before doing the clustering. 
#' @examples
#' data(mantestobj)
#' clu = clusterGene(testobj = mantestobj, k = 4)
clusterGene <- function(testobj, gene, type = 'variable', k.auto = FALSE,  k=5, method = 'kmeans', scale.difference = F){
  if (toupper(type) == 'TIME'){
    clu <- cluster_gene(testobj = testobj, 
                        gene = gene,
                        k = k,
                        type = type,
                        method = method)
  } else if (toupper(type) == 'VARIABLE'){
    if ('DDGType' %in% names(testobj)){
      DDGType <- testobj$DDGType
    } else {
      DDGType <- getDDGType(testobj)
    }
    if (k > 2){
      clu <- cluster_gene(testobj, gene = names(DDGType)[!DDGType %in% c('nonDDG', 'meanSig')], type = 'variable', k=k-2, scale.difference = scale.difference, method = method, k.auto = k.auto)
    } else {
      clu = NULL
    }
      
    design = testobj$design
    cellanno = testobj$cellanno
    meandiff <- sapply(c(0,1), function(i){
      as <- rownames(design[design[,2]==i, ])
      if ('expr' %in% names(testobj)){
        rowMeans(testobj$expr[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
      } else {
        rowMeans(testobj$expr.ori[names(DDGType)[DDGType == 'meanSig'], cellanno[cellanno[,2] %in% as,1]])
      }
      
    })
    large0 <- rownames(meandiff)[meandiff[,1] >= meandiff[,2]]
    large1 <- rownames(meandiff)[meandiff[,1] < meandiff[,2]]
    
    clu2 <- rep(max(clu)+1, length(large0))
    names(clu2) <- large0
    clu3 <- rep(max(clu)+2, length(large1))
    names(clu3) <- large1
    clu = c(clu, clu2, clu3)  
  }
  clu <- clu[gene]
  return(clu)
}
