#' Plot the cluster mean or of group differnce for a set of genes or plot the group difference for individual genes.
#'
#' This function is used for plotting the cluster mean of group difference or plotting the group difference for individual genes.
#'
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @import ggplot2 RColorBrewer reshape2
#' @export
#' @param testobj output object from lamian_test(). It is a list.
#' @param gene a vector of gene names.
#' @param cluster users can defined the clusters of the cells. It is a numeric vector whose names are cell names. By default it is retrieved from testobj object.
#' @param each logical. If TRUE, the plot for each gene will be generated.
#' @param sep a string in the gene names that needs to replaced with blank.
#' @param reverse logical. If FALSE (default), group difference is calculated by group 1 - group 0. If TRUE, group difference is calculated by group 0 - group 1.
#' @param facet.grid logical. If FALSE (default), seperate clusters by facet_wrap, otherwise by facet_grid.
#' @param free.scale logical. If TRUE (default), the y-axis is on free scale.
#' @param axis.text.blank logical. If TRUE, leave axis text as blank.
#' @examples
#' data(mantestobj)
#' plotClusterDiff(testobj = mantestobj)
plotClusterDiff <- function(testobj,
                            gene = names(testobj$cluster),
                            cluster = testobj[['cluster']],
                            each = FALSE,
                            sep = '',
                            reverse = FALSE,
                            facet.grid = FALSE,
                            free.scale = TRUE,
                            axis.text.blank = FALSE) {
  a <- ifelse(free.scale, 'free', 'fixed')
  if ('covariateGroupDiff' %in% names(testobj)) {
    fit <- testobj$covariateGroupDiff[gene, , drop = FALSE]
  } else {
    fit <-
      getCovariateGroupDiff(testobj = testobj,
                            gene = gene,
                            reverse = reverse)
  }
  colnames(fit) <- seq(1, ncol(fit))
  
  if (each) {
    pd <- melt(fit)
    colnames(pd) <- c('gene', 'pseudotime', 'covariateGroupDiff')
    pd[, 1] <- sub(sep, '', pd[, 1])
    pd$gene <-
      factor(as.character(pd$gene), levels = sub(sep, '', gene))
    if (facet.grid) {
      p <-
        ggplot(data = pd) + geom_line(aes(x = pd[,2], y = pd[,3])) +
        theme_classic() +
        facet_grid( ~ gene, scales = a)
    } else {
      p <-
        ggplot(data = pd) + geom_line(aes(x = pd[,2], y = pd[,3])) +
        theme_classic() +
        facet_wrap( ~ gene, scales = a)
    }
    
  }  else {
    clu = cluster[gene]
    tmp <- sapply(sort(unique(clu)), function(i) {
      m <- colMeans(fit[clu == i, , drop = FALSE])
    })
    colnames(tmp) <- sort(unique(clu))
    pd <- melt(tmp)
    
    colnames(pd) <- c('pseudotime', 'cluster', 'covariateGroupDiff')
    pd$cluster <- factor(pd$cluster)
    p <-
      ggplot(data = pd) + geom_line(aes(x = pd[,1], y = pd[,3], color = pd[,2])) +
      theme_classic() + scale_x_continuous(breaks = c(min(pd$pseudotime), max(pd$pseudotime)))
    if (length(unique(pd$cluster)) < 8) {
      p <- p + scale_color_brewer(palette = 'Dark2')
    } else {
      p <-
        p + scale_color_manual(values = grDevices::colorRampPalette(brewer.pal(8, 'Dark2'))(length(unique(pd$cluster))))
    }
  }
  if (axis.text.blank) {
    p <-
      p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p <-
      p + scale_x_continuous(breaks = c(min(pd$pseudotime), max(pd$pseudotime)))
  }
  print(p)
}
