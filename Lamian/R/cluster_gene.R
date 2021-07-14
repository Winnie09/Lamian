#' Support function for external clusterGene(). Cluster genes based on their temporal patterns  of gene expression (constant time test) or the temporal patterns of gene expression group difference (sample covariate test).
#'
#' This function is used to support the external function clusterGene(). cluster genes based on their temporal patterns  of gene expression (constant time test) or the temporal patterns of gene expression group difference (sample covariate test).
#'
#' @import ggplot2 RColorBrewer splines gridExtra viridis
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param testobj object returned from lamian.test().
#' @param gene a character vector of gene names. It can be of length 1 or > 1.
#' @param type One of c('Time', 'Variable').
#' @param k.auto logical. If FALSE (default), users need to specify k as the number of clusters. If TRUE, k will be automatically determined by elbow's method.
#' @param k a numeric number. The number of clusters. Only useful when k.auto = FALSE.
#' @param method the clustering method, either "kmeans" for k-means clustering or "hierarchical" for hierarchical clustering.
#' @param scale.difference logical. If FALSE, then do not standarize the group difference, but only scale by the maximum value of the group difference absolute values. If TRUE, then standardize the group difference before doing the clustering.
cluster_gene <- function(testobj,
                         gene,
                         k,
                         k.auto = FALSE,
                         type = 'time',
                         method = 'kmeans',
                         scale.difference = FALSE) {
  if (type == 'time') {
    if ('populationFit' %in% names(testobj)) {
      fit <- testobj$populationFit
    } else {
      fit <- getPopulationFit(testobj, gene = gene, type = type)
    }
  } else if (type == 'variable') {
    if ('covariateGroupDiff' %in% names(testobj)) {
      fit <- testobj$covariateGroupDiff
    } else{
      fit <- getCovariateGroupDiff(testobj = testobj, gene = gene)
    }
  }
  if (scale.difference) {
    mat.scale <- scalematrix(fit[gene, , drop = FALSE])
  } else {
    max <- apply(abs(fit[gene, , drop = FALSE]), 1, max)
    mat.scale <- fit[gene, , drop = FALSE] / max
  }
  
  if (method == 'kmeans') {
    # set.seed(12345)
    if (k.auto) {
      clu <- mykmeans(mat.scale, maxclunum = 20)$cluster
    } else {
      clu <- kmeans(mat.scale, k, iter.max = 1000)$cluster
    }
  } else if (method == 'hierarchical') {
    clu <- cutree(hclust(dist(mat.scale)), k = k)
  }
  # order clusters by genes' max expr position
  v <- sapply(unique(clu), function(i) {
    tmp <- fit[clu == i, ]
    mean(apply(tmp, 1, which.max))
  })
  names(v) <- unique(clu)
  trans <- cbind(as.numeric(names(v)), rank(v))
  n <- names(clu)
  clu <- trans[match(clu, trans[, 1]), 2]
  names(clu) <- n
  return(clu)
}
