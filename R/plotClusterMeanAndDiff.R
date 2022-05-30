#' Plot the cluster mean of population-level fitting values and group difference.
#'
#' This function is used for plotting the cluster mean of population-level fitting values and group difference.
#'
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @import ggplot2 RColorBrewer gridExtra reshape2 grDevices
#' @export
#' @param testobj output object from lamian_test(). It is a list.
#' @param cluster users can defined the clusters of the cells. It is a numeric vector whose names are cell names. By default it is retrieved from testobj object.
#' @param free.scale1 logical. If TRUE, free scale the y-axis of the cluster mean plot.
#' @param free.scale2 logical. If TRUE, free scale the y-axis of the group difference plot.
#' @examples
#' data(mantestobj)
#' plotClusterMeanAndDiff(testobj = mantestobj)

plotClusterMeanAndDiff <- function(testobj, 
                                   cluster = testobj[['cluster']],
                                   free.scale1 = TRUE,
                                   free.scale2 = FALSE){
  ## only works for Covariate Test. 
  if ('populationFit' %in% names(testobj)){
    fit <- testobj$populationFit
  } else {
    print("The object testobj should contain populationFit!")
  }
  ## give populationFit cell names as colnames
  colnames(testobj$populationFit[[1]]) <- colnames(testobj$populationFit[[1]]) <- colnames(testobj$expr) 
  a1 <- ifelse(free.scale1, 'free_y', 'fixed') 
  a2 <- ifelse(free.scale2, 'free_y', 'fixed') 
  int <- intersect(rownames(fit[[1]]), names(cluster))
  clu <- cluster[int]
  pd <- lapply(1:length(fit), function(i){
    mat <- fit[[i]][int, ]
    tmp <- sapply(sort(unique(clu)), function(i){
      colMeans(mat[clu == i, , drop = FALSE])
    })
    tmp <- melt(tmp)
    tmp[,1] = testobj[['pseudotime']][as.character(tmp[,1])]
    colnames(tmp) <- c('pseudotime', 'cluster', 'populationFitClusterMean')
    tmp <- data.frame(tmp, type = names(fit)[i])
  })
  pd <- do.call(rbind, pd)
  pd$cluster <- factor(pd$cluster)
  p1 <- ggplot(data = pd) + geom_line(aes(x = pseudotime, y = populationFitClusterMean, color = type), size = 1)+
    theme_classic() + 
    facet_wrap(~cluster, nrow = length(unique(pd$cluster)), scales = a1)+
    theme(legend.position = 'none')  +
    ylab('Model-fitted expression patterns for groups')+
    xlab('Pseudotime')
  if (length(unique(pd$type)) < 8){
    p1 <- p1 + scale_color_brewer(palette = 'Dark2')
  } else {
    p1 <- p1 + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$type))))
  }
  
  if ('covariateGroupDiff' %in% names(testobj)){
    fit <- testobj$covariateGroupDiff
  } else {
    fit <- getCovariateGroupDiff(testobj = testobj, gene = int)
  }
  fit <- fit[names(clu), , drop=FALSE]
  colnames(fit) <- seq(1, ncol(fit))
  
  
  tmp <- sapply(sort(unique(clu)), function(i){
    m <- colMeans(fit[clu == i, , drop = FALSE])
  })
  colnames(tmp) <- sort(unique(clu))
  pd2 <- melt(tmp)
  
  colnames(pd2) <- c('pseudotime', 'cluster', 'covariateGroupDiff')
  pd2$cluster <- factor(pd2$cluster)
  p2<- ggplot(data = pd2) + geom_line(aes(x = pseudotime, y = covariateGroupDiff), size = 1)+
    theme_classic() +
    facet_wrap(~cluster, nrow = length(unique(pd2$cluster)), scales = a2)+
    xlab('Pseudotime') + ylab('Group difference')
  if (length(unique(pd$cluster)) < 8){
    p2 <- p2 + scale_color_brewer(palette = 'Set1')
  } else {
    p2 <- p2 + scale_color_manual(values = colorRampPalette(brewer.pal(8,'Dark2'))(length(unique(pd$cluster))))
  }
  grid.arrange(p1,p2,ncol=2)
}

