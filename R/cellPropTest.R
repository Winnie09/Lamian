#' Perform cell proportion test.
#'
#' This function is used to perform cell proportion test.
#'
#' @import stats
#' @return  a list. A lamian.test output object except there is one more element of statistics in the list.
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param ncores number of cores. if ncores > 1, then implement in parallel mode.
#' @param permuiter the number of permutation. Value will be passed to the "permuiter" of lamian.test().
#' @param EMmaxiter the number of EM iterations. Value will be passed to the "EMmaxiter" of lamian.test().
#' @examples
#' data(mandata)
#' a <- cellPropTest(cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, ncores = 1, permuiter = 2, EMmaxiter = 3)

cellPropTest <-
  function(cellanno,
           pseudotime,
           design = NULL,
           permuiter = 100,
           EMmaxiter = 100,
           ncores = detectCores()) {
    ptw <-
      cut(pseudotime, seq(min(pseudotime), max(pseudotime), length.out = 100), include.lowest = TRUE)
    ptdat <-
      table(ptw, cellanno[match(names(pseudotime), cellanno[, 1]), 2])
    ptdat <-
      t(t(ptdat) / colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples.
    ptdat <- as.data.frame(ptdat)
    colnames(ptdat) <- c('pt', 's', 'prop')
    ptdat[, 1] <- match(ptdat[, 1], levels(ptw))
    
    ptdat$cell <- paste0('cell', seq_len(nrow(ptdat)))
    ptexpr <- t(ptdat[, c('prop', 'prop'), drop = FALSE])
    colnames(ptexpr) <- ptdat$cell
    
    ptpt <- ptdat$pt
    names(ptpt) <- ptdat$cell
    res <-
      lamian.test(
        expr = ptexpr,
        cellanno = data.frame(cell = ptdat$cell, sample = ptdat$s),
        pseudotime = ptpt,
        design = design,
        ncores = ncores,
        test.type = 'Variable',
        demean = FALSE,
        test.method = 'permutation',
        ncores.fit = 1,
        fix.all.zero = FALSE,
        permuiter = permuiter,
        EMmaxiter = EMmaxiter
      )
    res$statistics <- res$statistics[1, 2:3]
    return(res)
  }
