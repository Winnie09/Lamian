#' Plot enriched GO terms for each of the clusters.
#'
#' This function is used to plot enriched GO terms for each of the clusters.
#'
#' @import ggplot2 RColorBrewer splines gridExtra viridis reshape2
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param  goRes output from GOEnrich. A list of GO enrichment result. Length same as number of clusters.
#' @param  n number of top GO terms shown for each cluster
#' @param  sortByFDR if TRUE (default), sort GO terms by FDR first, and then negative FC. If FALSE,  soft by negative FC and then FDR.
#' @param fdr.cutoff FDR cutoff for choosing "statistically significant" GO terms. Will be combined with those choosing by fc.cutoff.
#' @param fc.cutff fold change cutoff for choosing "statistically significant" GO terms. Will be combined with those choosing by fdr.cutoff.
#' @examples
#' data(mantestobj)
#' plotGOEnrich(goRes = mantestobj$goRes, n = 2, fdr.cutoff = 1, fc.cutoff = 2)


plotGOEnrich <-
  function(goRes,
           n = 5,
           sortByFDR = TRUE,
           fdr.cutoff = 0.05,
           fc.cutoff = 2) {
    d <- lapply(names(goRes), function(i) {
      cbind(Cluster = i, goRes[[i]])
    })
    d <- do.call(rbind, d)
    if (sortByFDR) {
      d <- d[order(d$FDR, -d$FC), ]
    } else {
      d <- d[order(-d$FC, d$FDR), ]
    }
    cd <- do.call(rbind, sapply(sort(unique(d$Cluster)), function(i) {
      tmp <- d[d$Cluster == i, ]
      tmp <-
        tmp[tmp$FDR < fdr.cutoff & tmp$FC > fc.cutoff, , drop = FALSE]
      if (nrow(tmp) > 0) {
        tmp[seq_len(min(n, nrow(tmp))), ]
      } else {
        NULL
      }
    }, simplify = FALSE))
    ut <- unique(cd$Term)
    d <- d[d$Term %in% ut, c('Cluster', 'Term', 'FDR', 'FC')]
    d <- d[complete.cases(d),]
    d <- d[d$FDR < fdr.cutoff & d$FC > fc.cutoff, , drop = FALSE]
    
    dmat <- dcast(d, Term ~ Cluster)
    rownames(dmat) <- dmat[, 1]
    dmat <- as.matrix(dmat[, -1, drop = FALSE])
    dmat <- is.na(dmat)
    
    v <- sapply(seq_len(nrow(dmat)), function(i)
      which.min(dmat[i, ]))
    names(v) <- rownames(dmat)
    
    pd <- melt(dmat)
    colnames(pd) <- c('Term', 'Cluster', 'enc')
    pd$Cluster = as.character(pd$Cluster)
    pd$Term <-
      factor(as.character(pd$Term), levels = names(v[rev(order(v))]))
    pd$enc <- ifelse(pd$enc, 'Non-significant', 'Significant')
    pd$Cluster <- factor(pd$Cluster, levels = unique(pd$Cluster)[order(as.numeric(unique(pd$Cluster)))])
    p <-
      ggplot(pd, aes(x = pd[,2], y = pd[,1], fill = pd[,3])) + geom_tile() + theme_classic() + scale_fill_manual(values =
                                                                                                               c('darkblue', 'orange')) +
      theme(legend.position = 'right') +
      scale_y_discrete(position = "right") +
      xlab('Cluster') + ylab('GO Terms')
    return(p)
  }
