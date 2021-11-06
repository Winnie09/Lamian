#' Obtain the population-level fitting of the genes.
#'
#' This function is used to obtain the population-level fitting pattern along pseudotime of the genes, using the estimated parameters from the lamian model.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @import splines
#' @export
#' @return a matrix, population fitting of the genes (rows) by cells (ordered by pseudotime)
#' @param testobj output object of the function tespt().
#' @param gene a vector of gene names
#' @param type one of c('Time', 'Variable). Case insensitive.
#' @examples
#' data(mantestobj)
#' a <- getPopulationFit(testobj = mantestobj, gene = rownames(mantestobj$populationFit[[1]])[seq(1,3)], type = 'variable')

getPopulationFit <- function(testobj,
                             gene = NULL,
                             type = 'time') {
  ## if type = 'time', then return population fit (a vector for a gene; or a gene by num.cell matrix) for constant test (test on time)
  ## if type = 'variable', then return population fit for all levels of that character (a matrix, columns are population fit for each level in the variabel).
  ## gene: a vector of gene names.
  design = testobj$design
  pseudotime = testobj$pseudotime
  knotnum = testobj$knotnum
  pseudotime = pseudotime[order(pseudotime)]
  type <- toupper(type)
  if (sum(design[, 1]) != nrow(design)) {
    print(
      "The first column of design matrix should be all 1s (intercept)! Using the first column as the variable column ..."
    )
    design = cbind(intercept = 1, design)
  }
  colnames(design)[1] <- 'intercept'
  if (is.null(gene))
    gene <- rownames(testobj$statistics)
  if (type == 'TIME') {
    design = design[, 1, drop = FALSE]
  } else {
    variable = colnames(design)[2]
    design <- unique(design[, c('intercept', variable)])
    rownames(design) <-
      paste0(variable, '_', unique(design[, variable]))
  }
  
  fitlist <- lapply(gene, function(g) {
    # beta <- lapply(testobj$parameter[g], function(i) {
    #   i$beta
    # })
    # names(beta) <- g
    beta <- testobj$parameter[[g]]$beta
    x <- sapply(row.names(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE]) ###
    }, simplify = FALSE)
    
    pt <- seq(1, max(pseudotime))
    if (knotnum[g] == 0) {
      # phi <- cbind(1, bs(pt))
      phi <- bs(pt, intercept = TRUE)
    } else {
      knots = seq(min(pt), max(pt), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      # phi <- cbind(1, bs(pt, knots = knots))
      phi <- bs(pt, knots = knots, intercept = TRUE)
    }
    if (!exists('variable')) {
      if (ncol(phi) == ncol(x[[1]])) {
        fit <- t(phi %*% t(x[[1]]) %*% beta)[1, ]
      } else {
        fit <- t(phi %*% x[[1]] %*% beta)[1, ]
      }
    } else {
      fit <- lapply(x, function(i) {
        if (ncol(phi) == nrow(i)) {
          phi %*% i %*% beta
        } else {
          phi %*% t(i) %*% beta
        }
      })
      names(fit) <- names(x)
    }
    return(fit)
  })
  
  names(fitlist) <- gene
  if (type == 'VARIABLE') {
    fitres <- lapply(names(fitlist[[1]]), function(i) {
      tmp <- t(sapply(fitlist, function(j) {
        j[[i]]
      }))
    })
    names(fitres) <- names(fitlist[[1]])
  } else if (type == 'TIME') {
    fitres <- t(do.call(cbind, fitlist))
    if (ncol(fitres) == length(testobj$pseudotime)) {
      colnames(fitres) <- names(testobj$pseudotime)
    }
  }
  return(fitres)
}
