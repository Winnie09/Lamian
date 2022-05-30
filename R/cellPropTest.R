#' Perform cell proportion test.
#'
#' This function is used to perform cell proportion test.
#'
#' @import stats
#' @return  a list. A lamian_test output object except there is one more element of statistics in the list.
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param ncores number of cores. if ncores > 1, then implement in parallel mode.
#' @param test.type the type of test. test.type = 'Time', use TCD test. test.type = 'Variable', use XCD test.
#' @param testvar a numeric number indicating the column in the design matrix that needs to be tested while controlling for other columns (not intercept). Default is 2. testvar = 2 means the second column in the design needs to be tested.
#' @examples
#' data(mandata)
#' a <- cellPropTest(cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, ncores = 1)

cellPropTest <- function(cellanno,
                         pseudotime,
                         design = NULL,
                         ncores = detectCores(),
                         test.type = 'variable',
                         testvar = 2) {
  ptw <-
    cut(pseudotime, seq(min(pseudotime), max(pseudotime), length.out = 100), include.lowest = T)
  ptdat <-
    table(ptw, cellanno[match(names(pseudotime), cellanno[, 1]), 2])
  ptdat <-
    t(t(ptdat) / colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples.
  ptdat <- as.data.frame(ptdat)
  colnames(ptdat) <- c('pt', 's', 'prop')
  ptdat[, 1] <- match(ptdat[, 1], levels(ptw))
  
  ptdat$cell <- paste0('cell', 1:nrow(ptdat))
  ptexpr <- t(ptdat[, c('prop', 'prop'), drop = F])
  colnames(ptexpr) <- ptdat$cell
  
  ptpt <- ptdat$pt
  names(ptpt) <- ptdat$cell
  ptcellanno <-
    data.frame(
      cell = ptdat$cell,
      sample = ptdat$s,
      stringsAsFactors = F
    )
  res <-
    lamian_test(
      expr = ptexpr,
      cellanno = ptcellanno,
      pseudotime = ptpt,
      design = design,
      testvar = testvar,
      ncores = ncores,
      test.type = test.type,
      test.method = 'permutation',
      ncores.fit = 1,
      fix.all.zero = F
    )
  
  res$statistics <- res$statistics[1, 2:3]
  return(res)
}
