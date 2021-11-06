#' Perform k-means clustering
#'
#' Perform k-means clustering for the cells (or genes) on any cell (or gene) by low-dimensional reprensentation coordinates matrix.
#'
#' @import parallel
#' @return k-means clustering result
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param matrix cell by principal component (pc) matrix for cell clustering.
#' @param number.cluster the number of clusters in cell clustering that will be used in trajectory inference. If NA (default),the number of clusters will be determined automatically by elbow's method.
#' @param maxclunum the maximum number of clusters in the elbew's method.
#' @param seed seed applied before kmeans() to address it reproducibility issue.
#' @examples
#' data(mandata)
#' a = mykmeans(matrix = mandata$expr, number.cluster = 2)
mykmeans <-
  function(matrix,
           number.cluster = NA,
           maxclunum = 30,
           seed = 12345) {
    library(parallel)
    if (is.na(number.cluster)) {
      rss <- mclapply(seq_len(maxclunum), function(clunum) {
        ## set.seed(seed)
        tmp <- kmeans(matrix, clunum, iter.max = 1000)
        tmp$betweenss / tmp$totss
      }, mc.cores = 30)
      rss <- unlist(rss)
      # number.cluster <- which(diff(rss) < 1e-2)[1]
      x <- 2:maxclunum
      number.cluster <-
        x[which.min(sapply(seq_len(length(x)), function(i) {
          x2 <- pmax(0, x - i)
          sum(lm(rss[-1] ~ x + x2)$residuals ^ 2)  ## check this
        }))]
    }
    ## set.seed(seed)
    clu <- kmeans(matrix, number.cluster)
    return(clu)
  }
