#' Perform model fitting in one permutation setting.
#'
#' This function is used to perform model fitting in one permutation setting.
#'
#' @author Wenpin Hou <wh2526@cumc.columbia.edu>
#' @return a list of two elements. The first is for full model, while the second is for null model. Each element is a list containing data and parameter estimates.
#' @import stats Matrix splines parallel
#' @param iter a numeric number indicating which iteration it is in a permutation test. If iter = 1, it is  fitting the original data (not permuted).
#' @param diffType one of overall, meanDiff, or trendDiff.
#' @param gene a vector of gene names. It is used when users want to externally specify the genes (or specify the order or genes).
#' @param testvar a numeric number indicating the column in the design matrix that needs to be tested while controlling for other columns (not intercept). Its value when calls the function will be used. e.g. testvar = 2 means the second column in the design needs to be tested.
#' @param test.type the type of test. One of c('Time', 'Variable). Case insensitive. 
#' @param testvar a numeric number indicating the column in the design matrix that needs to be tested while controlling for other columns (not intercept). Default is 2. testvar = 2 means the second column in the design needs to be tested.
#' @param expr gene by cell expression matrix. The expression values should have been library-size-normalized and log-transformed. They can either be imputed or non-imputed.
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param ncores the number of cores to be used. If ncores > 1, it will be implemented in parallel mode.
#' @param verbose.output logical. If TRUE, print intermediate information.

fitfunc <- function(iter, diffType = 'overall', gene = rownames(expr), test.type = 'Time', testvar=testvar, maxknotallowed = 10, EMmaxiter=100, EMitercutoff=0.05, ncores=1, expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design, verbose.output = FALSE) {
  expr <- expr[gene, ,drop=FALSE]
  if (verbose.output) print(paste0('iter ', iter, '\n'))
  if (toupper(test.type)=='TIME') {
    if (iter == 1){
      fitres.full <- fitpt(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design[,1,drop=FALSE], testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=ncores, model=-1)
      fitres.null <- fitpt.m0(expr=expr, cellanno=cellanno, pseudotime=pseudotime, design=design[,1,drop=FALSE], EMmaxiter=EMmaxiter, maxknotallowed = maxknotallowed, EMitercutoff=EMitercutoff)
      return(list(fitres.full = fitres.full, fitres.null = fitres.null))
    } else {
      perpsn <- lapply(rownames(design), function(s){
        tmpid <- cellanno[cellanno[,2] == s, 1]  # subset cells
        tmppsn <- pseudotime[names(pseudotime) %in% tmpid] # subset time
        names(tmppsn) <- sample(names(tmppsn)) # permute time
        tmppsn
      })
      names(perpsn) <- NULL
      perpsn <- unlist(perpsn)
      perpsn <- perpsn[names(pseudotime)]
      sampcell <- as.vector(unlist(lapply(unique(cellanno[,2]), function(p){
        sample(which(cellanno[,2] == p), replace = T)
      }))) ### bootstrap
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      perpsn <- perpsn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(perpsn) <- paste0('cell_',1:length(perpsn))
      tryCatch(fitres.full <- fitpt(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design,testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=1, model = -1), warning = function(w){}, error = function(e) {})
      tryCatch(fitres.null <- fitpt.m0(expr=perexpr, cellanno=percellanno, pseudotime=perpsn, design=design, EMmaxiter=EMmaxiter, maxknotallowed = maxknotallowed, EMitercutoff=EMitercutoff), warning = function(w){}, error = function(e) {})
      if (exists('fitres.full') & exists('fitres.null')) {
        if (verbose.output) print(paste0('iter ', iter, ' success!'))
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else {
        if (verbose.output) print(paste0('iter ', iter, ' try again!'))
        return(NULL)
      }
    }
  } else if (toupper(test.type)=='VARIABLE'){
    if (verbose.output) print('testing Variable ...')
    if (diffType == 'overall'){
      mod.full = 3
      mod.null = 1
    } else if (diffType == 'meanDiff'){
      mod.full = 2
      mod.null = 1
    } else if (diffType == 'trendDiff'){
      mod.full = 3
      mod.null = 2
    }
    if (iter == 1){
      fitres.full <- fitpt(expr, cellanno, pseudotime, design, testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=1, model = mod.full) 
      fitres.null <- fitpt(expr, cellanno, pseudotime, design, testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=1, model = mod.null, knotnum = fitres.full[[2]])
      if (exists('fitres.full') & exists('fitres.null')) {
        if (verbose.output) print(paste0('iter ', iter, ' success!'))
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else {
        if (verbose.output) print(paste0('iter ', iter, ' try again!'))
        return(NULL)
      }
    } else {
      dn <- paste0(design[,testvar],collapse = '_')
      perdn <- dn
      # make sure perdesign is full rank
      while (TRUE) {
        perid <- sample(1:nrow(design))
        perdesign <- design
        perdesign[, testvar] <- design[perid, testvar]
        rnk <- as.integer(Matrix::rankMatrix(perdesign))
        perdn_new <- paste0(perdesign[, testvar], collapse = '_')
        if (perdn_new != dn && rnk == ncol(perdesign)) {
          perdn <- perdn_new
          break 
        }
      }
      row.names(perdesign) <- row.names(design)
      sampcell <- sample(1:ncol(expr),replace=T) ## boostrap cells
      perexpr <- expr[,sampcell,drop=F]
      percellanno <- cellanno[sampcell,,drop=F]
      psn <- pseudotime
      psn <- psn[sampcell]
      colnames(perexpr) <- percellanno[,1] <- names(psn) <- paste0('cell_',1:length(psn))
      
      
      fitres.full <- fitpt(perexpr, percellanno, psn, perdesign,testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=ncores, model = mod.full)
      fitres.null <- fitpt(perexpr, percellanno, psn, perdesign,testvar=testvar, maxknotallowed = maxknotallowed, EMmaxiter=EMmaxiter, EMitercutoff=EMitercutoff, ncores=ncores, model = mod.null, knotnum = fitres.full[[2]])
      if (exists('fitres.full') & exists('fitres.null')) {
        if (verbose.output) print(paste0('iter ', iter, ' success!'))
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else {
        if (verbose.output) print(paste0('iter ', iter, ' try again!'))
        return(NULL)
      }
    }
  }
}

