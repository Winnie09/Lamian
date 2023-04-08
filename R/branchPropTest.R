#' Perform branch proportion test.
#'
#' This function is used to perform branch proportion test.
#'
#' @import stats
#' @return  a vector of p.values.
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param data If "method = 't.test'", then the data is a dataframe or matrix where the each row is a branch, each column is a sample (patient). Values are the branch proportion. This data can be the second element of the output (a list) of evaluate_uncertainty(). If "method == 'multinom'", then the data is a dataframe where each row is a cell, each of the three columns are cell names, branch names, and sample names.
#' @param design: a data frame. Number of rows should be the same as the number of unique samples. Each row is a sample. First column is the sample names. Second column is the covariate realization valuels for each of the samples that needs to be tested on.
#' @examples
#' a <- branchPropTest(data = matrix(rnorm(24), nrow =3), design = data.frame(sample = paste0('BM', seq(1,8)), sex = c(rep(1,4), rep(0,4))), method = 't.test')

branchPropTest <- function(data, design, method = 't.test') {
  if (method == 't.test'){
    id1 = which(design[, 2] == 0)
    id2 = which(design[, 2] == 1)
    res <- apply(data, 1, function(i)
      t.test(i[id1], i[id2])$p.value)
    
  }
  else if (method == 'multinom'){
    colnames(data) <- c('cell', 'branch', 'sample')
    data$covariate <- design[match(data[,3], rownames(design)), 2]
    test <- nnet::multinom(branch ~ covariate, data = data)
    z <- summary(test)$coefficients/summary(test)$standard.errors
    p <- (1 - pnorm(abs(z), 0, 1)) * 2
    res <- p
  }
  return(res)
}

