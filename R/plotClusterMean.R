#' Plot the cluster mean of population-level fitting values.
#'
#' This function is used for plotting the cluster mean of population-level fitting values.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @import ggplot2 RColorBrewer gridExtra viridis reshape2
#' @export
#' @return a plot
#' @param testobj output object from lamian.test(). It is a list.
#' @param cluster users can defined the clusters of the cells. It is a numeric vector whose names are cell names. By default it is retrieved from testobj object.
#' @param type one of c('Time', 'Variable').
#' @param facet logical. If TRUE, facet each cluster.
#' @examples
#' data(mantestobj)
#' plotClusterMean(mantestobj, type = 'variable')

plotClusterMean <- function(testobj,
                            cluster,
                            type = 'time',
                            facet = FALSE) {
  ## type: "time" (default) or "variable.
  ## facet: plot each cluster indiviidually. Lnly works when type == 'time'. When type == 'variable', the plot will be facet anyway for each cluster.
  if ('populationFit' %in% names(testobj)) {
    fit <- testobj$populationFit
  } else {
    print("The object testobj should contain populationFit!")
  }
  if ('cluster' %in% names(testobj))
    cluster <- testobj[['cluster']]
  if (toupper(type) == 'VARIABLE') {
    int <- intersect(rownames(fit[[1]]), names(cluster))
    clu <- cluster[int]
    pd <- lapply(seq_len(length(fit)), function(i) {
      mat <- fit[[i]][int,]
      tmp <- sapply(sort(unique(clu)), function(i) {
        colMeans(mat[clu == i, , drop = FALSE])
      })
      tmp <- melt(tmp)
      colnames(tmp) <-
        c('pseudotime', 'cluster', 'populationFitClusterMean')
      tmp <- data.frame(tmp, type = names(fit[i]))
    })
    pd <- do.call(rbind, pd)
    pd$cluster <- factor(pd$cluster)
    p <-
      ggplot(data = pd) + geom_line(aes(x = pd[,1], y = pd[,3], color = type),
                                    size = 1) +
      theme_classic() +
      theme(axis.text.x = element_blank()) +
      facet_wrap( ~cluster, scale = 'free') + 
      xlab('Pseudotime') + 
      ylab('Expression')
    if (length(unique(pd$type)) < 8) {
      p <- p + scale_color_brewer(palette = 'Dark2')
    } else {
      p <-
        p + scale_color_manual(values = grDevices::colorRampPalette(brewer.pal(8, 'Dark2'))(length(unique(pd$type))))
    }
  } else if (toupper(type) == 'TIME') {
    int <- intersect(rownames(fit), names(cluster))
    clu <- cluster[int]
    mat <- fit[int,]
    tmp <- sapply(sort(unique(clu)), function(i) {
      colMeans(mat[clu == i, , drop = FALSE])
    })
    pd <- melt(tmp)
    colnames(pd) <-
      c('pseudotime', 'cluster', 'populationFitClusterMean')
    pd$pseudotime <- as.numeric(pd$pseudotime)
    pd$cluster <- factor(pd$cluster)
    p <-
      ggplot(
        data = pd,
        aes(
          x = pseudotime,
          y = populationFitClusterMean,
          group = cluster,
          color = cluster
        )
      ) +
      geom_smooth(size = 1) +
      theme_classic() +
      theme(axis.text.x = element_blank())
    if (facet) {
      p <- p + facet_wrap( ~ cluster)
    }
    if (length(unique(pd$cluster)) < 8) {
      p <- p + scale_color_brewer(palette = 'Set1')
    } else {
      p <-
        p + scale_color_manual(values = grDevices::colorRampPalette(brewer.pal(9, 'Set1'))(length(unique(pd$cluster))))
    }
  }
  
  print(p)
}
