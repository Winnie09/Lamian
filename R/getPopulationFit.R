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
#' @param downsample logical. If TRUE, the cells will be downsampled, default number is num.timepoint = 1e3; if FALSE, the cells will not be downsampled.
#' @param num.timepoint the number of time points used to fit the curve. Default is the minimum of 1e3 and max(pseudotime). This argument will reduce the running time of both this function and the downstream plot function if there are too many cells. 
#' @examples
#' data(mantestobj)
#' a <- getPopulationFit(testobj = mantestobj, gene = rownames(mantestobj$populationFit[[1]])[seq(1,3)], type = 'time')

getPopulationFit <- function(testobj,
                             gene = NULL,
                             type = 'time',
                             downsample = FALSE,
                             num.timepoint = 1e3){
  type <- toupper(type)
  if (!'testvar' %in% names(testobj)) {
    testvar <- testobj$testvar <- 2
  } else {
    testvar = testobj$testvar
  }
    
  if (type == 'VARIABLE'){
    design = testobj$design[, c(1, testobj$testvar)] ## design for multi
  } else {
    design = testobj$design
  }
  knotnum = testobj$knotnum
  pseudotime = testobj$pseudotime
  pseudotime = pseudotime[order(pseudotime)]
  if (downsample){
    if (num.timepoint > length(pseudotime)){
      print('Number of existing cells exceeds target downsampling number. Ignoring downsampling...')
      pt = pseudotime
    } else {
      print('Downsampling the cells...')
      pt <- round(seq(1, length(pseudotime), length.out = min(num.timepoint, length(pseudotime)))) ## downsample the cells; select equidistant pseudotime time points  
    }
  } else {
    pt = pseudotime
  }
    
  
  if (sum(design[, 1]) != nrow(design)){
    print("The first column of design matrix should be all 1s (intercept)! We are inserting an all-1 column as the first column ...")
    design = cbind(intercept = 1, design)
    colnames(design)[1] <- 'intercept'
  }

  if (is.null(gene)) gene <- rownames(testobj$statistics)
  
  if (type == 'TIME') {
    design = design[, 1, drop = FALSE]
  } else {
    variable = colnames(design)[2]
    design <- unique(design[, c(colnames(design)[1], variable)])
    rownames(design) <- paste0(variable, '_', unique(design[, variable]))
  }
  
  fitlist <- lapply(gene, function(g){
    tmp = matrix(testobj$parameter[[g]]$beta, ncol = knotnum[g]+4)
    if (type == 'VARIABLE'){
      beta = as.vector(tmp[c(1,testvar), ]) ### subset the beta values of the intercept and the test covariate for multi
    } else {
      beta = as.vector(tmp[1, ]) ### subset the beta values of the intercept and the test covariate for multi
    }
      
    x <- sapply(row.names(design), function(i) {
      kronecker(diag(knotnum[g] + 4), design[i, , drop = FALSE]) ###
    }, simplify = FALSE)
    
    
    if (knotnum[g] == 0) {
      # phi <- cbind(1, bs(pt))
      phi <- bs(pt, intercept = TRUE)
    } else {
      knots = seq(min(pt), max(pt), length.out = knotnum[g] + 2)[2:(knotnum[g] + 1)]
      # phi <- cbind(1, bs(pt, knots = knots))
      phi <- bs(pt,knots = knots, intercept = TRUE)
    }
    
    if (type == 'VARIABLE') {
      fit <- lapply(x, function(i) {
        if (ncol(phi) == nrow(i)){
          phi %*% i %*% beta
        } else {
          phi %*% t(i) %*% beta
        }
      })
    } else {
      i = x[[1]]
      if (ncol(phi) == nrow(i)){
        fit <- phi %*% i %*% beta
      } else {
        fit <- phi %*% t(i) %*% beta
      }
    }
    return(fit)
  })
  
  names(fitlist) <- gene
  if (type == 'VARIABLE'){
    fitres <- lapply(names(fitlist[[1]]), function(i){
      tmp <- t(sapply(fitlist, function(j){
        j[[i]]
      }))
      colnames(tmp) = names(pseudotime)[pt]
    })  
    names(fitres) <- names(fitlist[[1]])
  } else if (type == 'TIME'){
    fitres <- t(do.call(cbind, fitlist))
    rownames(fitres) <- gene
    if (ncol(testobj$expr) == ncol(fitres)) colnames(fitres) <- colnames(testobj$expr)
    colnames(fitres) = names(pseudotime)[pt]  
  }
  
  
  return(fitres)
}
