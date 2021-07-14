#' Obtain the covariate group difference of the genes in sample covariate test.
#'
#' This function is used to obtain the covariate group difference of the genes in sample covariate test.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a matrix of gene by cell (ordered by pseudotime). Values are covariate group difference.
#' @export
#' @param testobj output object of the function tespt().
#' @param gene a character vector of genes.
#' @param output a gene by pseudotime matrix. Entries are group difference w.r.t the variable, i.e., the unit-covariate incremental difference.
#' @param reverse logitcal. FALSE (default) when group 1 - group 0.  TRUE when group 0 - group 1.
#' @examples
#' data(mantestobj)
#' a = getCovariateGroupDiff(testobj = mantestobj, gene = rownames(expr)[seq_len(2)])

getCovariateGroupDiff <- function(testobj,
                                  gene,
                                  reverse = FALSE) {
  knotnum = testobj$knotnum[gene]
  pseudotime = seq(1, max(testobj$pseudotime))
  beta <- lapply(testobj$parameter[gene], function(i) {
    if (reverse) {
      -i$beta
    } else {
      i$beta
    }
  })
  names(beta) <- gene
  philist <-
    lapply(min(knotnum):max(knotnum), function(num.knot) {
      if (num.knot == 0) {
        # phi <- cbind(1, bs(pseudotime))
        phi <- bs(pseudotime, intercept = TRUE)
      } else {
        knots = seq(min(pseudotime), max(pseudotime), length.out = num.knot + 2)[2:(num.knot +
                                                                                      1)]
        # phi <- cbind(1, bs(pseudotime, knots = knots))
        phi <- bs(pseudotime, knots = knots, intercept = TRUE)
      }
    })
  names(philist) <- as.character(min(knotnum):max(knotnum))
  
  fit <- lapply(gene, function(i) {
    id <- (seq_len(length(beta[[i]]) / 2)) * 2
    a <- philist[[as.character(knotnum[i])]] %*% beta[[i]][id]
  })
  fit <- do.call(cbind, fit)
  colnames(fit) <- gene
  return(t(fit))
}
