#' Perform branch proportion test.
#'
#' This function is used to perform branch proportion test.
#'
#' @import stats
#' @return  a vector of p.values.
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param data a dataframe or matrix where the each row is a branch, each column is a sample (patient). Values are the branch proportion. This data can be the second element of the output (a list) of evaluate_uncertainty().
#' @param design: a data frame. Number of rows should be the same as the number of unique samples. Each row is a sample. First column is the sample names. Second column is the covariate realization valuels for each of the samples that needs to be tested on.
#' @examples
#' a <- branchPropTest(data = matrix(rnorm(24), nrow =3), design = data.frame(sample = paste0('BM', seq(1,8)), sex = c(rep(1,4), rep(0,4))))

branchPropTest <- function(data, design) {
  id1 = which(design[, 2] == 0)
  id2 = which(design[, 2] == 1)
  apply(data, 1, function(i)
    t.test(i[id1], i[-c(id2)])$p.value)
}

