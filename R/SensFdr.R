#' Calculate Sensitivity, real FDR, and reported FDR
#'
#' This function is used to calculate Sensitivity, real FDR, and reported FDR.
#'
#' @import ggplot2 RColorBrewer splines gridExtra viridis
#' @return a datafrome. Columns are c('Sensitivity','Real_FDR','Reported_FDR').
#' @export
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param TruePositive a character vector of true positive genes.
#' @param Statistics: a dataframe or matrix, should contain a column of fdr , for example c('FStat','P.Value','adj.P.Val'). Row names are gene names.
#' @examples
#' stat = data.frame(fdr = abs(rnorm(10)), fc = rnorm(10))
#' rownames(stat) = letters[seq(1,10)]
#' a = SensFdr(TruePositive = letters[seq(1,5)], Statistics = stat)
SensFdr <- function(TruePositive, Statistics) {
  ## find the column of fdr
  col.na <- which(colSums(is.na(Statistics)) == nrow(Statistics))
  if (length(col.na) > 0)
    Statistics <- Statistics[, -col.na, drop = FALSE]
  fdrchar <-
    intersect(
      colnames(Statistics),
      c(
        'adj.P.Val',
        'adj.pvalue',
        'fdr',
        'FDR',
        'Fdr',
        'adj.p',
        'adj.P',
        'adj.Pval',
        'fdr.overall'
      )
    )
  fdrcol <- which(colnames(Statistics) == fdrchar)
  ## if not ordered by significance, then rank by significance
  
  if (sum((diff(Statistics[, fdrcol]) < 0) + 0) > 0) {
    Statistics <-
      Statistics[order(Statistics[, fdrcol]), , drop = FALSE]
  }
  Order <- rownames(Statistics)
  ## calculate sensitivity, realfdr, reported fdr
  perf <- t(sapply(seq(1, length(Order)), function(i) {
    num <- sum(Order[seq(1, i)] %in% TruePositive)
    c(num / length(TruePositive), (i - num) / i, Statistics[i, fdrcol])
  }))
  ## reorder real fdr so that it is monotonically increasing
  if (nrow(perf) > 1) {
    for (i in (nrow(perf)):2) {
      if (perf[i - 1, 2] > perf[i, 2])
        perf[i - 1, 2] <- perf[i, 2]
    }
  }
  colnames(perf) <- c('Sensitivity', 'Real_FDR', 'Reported_FDR')
  rbind(c(0, 0, 0), perf)
}

