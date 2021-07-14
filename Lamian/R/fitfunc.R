#' Perform model fitting in one permutation setting.
#'
#' This function is used to perform model fitting in one permutation setting.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list of two elements. The first is for full model, while the second is for null model. Each element is a list containing data and parameter estimates.
#' @export
#' @import stats Matrix splines parallel
#' @importFrom stats cov2cor
#' @importFrom stats decompose
#' @importFrom stats toeplitz
#' @importFrom stats update
#' @importFrom stats spectrum
#' @param expr gene by cell expression matrix. The expression values should have been library-size-normalized and log-transformed. They can either be imputed or non-imputed.
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param verbose logical. If TRUE, print intermediate information.
#' @param iter a numeric number indicating which iteration it is in a permutation test. If iter = 1, it is  fitting the original data (not permuted).
#' @param diffType one of overall, meanDiff, or trendDiff.
#' @param gene a vector of gene names. It is used when users want to externally specify the genes (or specify the order or genes).
#' @param test.type the type of test. One of c('Time', 'Variable). Case insensitive.
#' @param ncores the number of cores to be used. If ncores > 1, it will be implemented in parallel mode.
#' @examples
#' data(mandata)
#' a = fitfunc(iter = 1, diffType = 'overall', gene = rownames(mandata$expr), test.type = 'Time', EMmaxiter=5, EMitercutoff=10, verbose=FALSE, ncores=1, expr= mandata$expr, cellanno= mandata$cellanno, pseudotime=mandata$pseudotime, design=mandata$design)


fitfunc <-
  function(iter,
           diffType = 'overall',
           gene = rownames(expr),
           test.type = 'Time',
           EMmaxiter = 100,
           EMitercutoff = 0.1,
           verbose = FALSE,
           ncores = 1,
           expr = expr,
           cellanno = cellanno,
           pseudotime = pseudotime,
           design = design) {
    ## this function serves the function lamian.test().
    ## return(list(fitres.full = fitres.full, fitres.null = fitres.null))
    ## ncores = 1 or otherwise meaningless since the upper function is running in parallel
    expr <- expr[gene, , drop = FALSE]
    if (toupper(test.type) == 'TIME') {
      if (iter == 1) {
        fitres.full <-
          fitpt(
            expr = expr,
            cellanno = cellanno,
            pseudotime = pseudotime,
            design = design[, 1, drop = FALSE],
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose,
            ncores = ncores,
            model = 1
          )
        fitres.null <-
          fitpt.m0(
            expr = expr,
            cellanno = cellanno,
            pseudotime = pseudotime,
            design = design[, 1, drop = FALSE],
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose
          )
        return(list(fitres.full = fitres.full, fitres.null = fitres.null))
      } else {
        perpsn <- lapply(rownames(design), function(s) {
          tmpid <- cellanno[cellanno[, 2] == s, 1]  # subset cells
          tmppsn <-
            pseudotime[names(pseudotime) %in% tmpid] # subset time
          names(tmppsn) <- sample(names(tmppsn)) # permute time
          tmppsn
        })
        names(perpsn) <- NULL
        perpsn <- unlist(perpsn)
        perpsn <- perpsn[colnames(expr)]
        sampcell <-
          as.vector(unlist(lapply(unique(cellanno[, 2]), function(p) {
            sample(which(cellanno[, 2] == p), replace = TRUE)
          }))) ### bootstrap
        perexpr <- expr[, sampcell, drop = FALSE]
        percellanno <- cellanno[sampcell, , drop = FALSE]
        perpsn <- perpsn[sampcell]
        colnames(perexpr) <-
          percellanno[, 1] <-
          names(perpsn) <- paste0('cell_', seq_len(length(perpsn)))
        tryCatch(
          fitres.full <-
            fitpt(
              expr = perexpr,
              cellanno = percellanno,
              pseudotime = perpsn,
              design = design,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose = verbose,
              ncores = 1,
              model = 1
            ),
          warning = function(w) {
          },
          error = function(e) {
          }
        )
        tryCatch(
          fitres.null <-
            fitpt.m0(
              expr = perexpr,
              cellanno = percellanno,
              pseudotime = perpsn,
              design = design,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose = verbose
            ),
          warning = function(w) {
          },
          error = function(e) {
          }
        )
        if (exists('fitres.full') & exists('fitres.null')) {
          return(list(fitres.full = fitres.full, fitres.null = fitres.null))
        } else {
          return(NULL)
        }
      }
    } else if (toupper(test.type) == 'VARIABLE') {
      if (diffType == 'overall') {
        mod.full = 3
        mod.null = 1
      } else if (diffType == 'meanDiff') {
        mod.full = 2
        mod.null = 1
      } else if (diffType == 'trendDiff') {
        mod.full = 3
        mod.null = 2
      }
      if (iter == 1) {
        fitres.full <-
          fitpt(
            expr,
            cellanno,
            pseudotime,
            design,
            maxknotallowed = 10,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose,
            ncores = 1,
            model = mod.full
          )
        fitres.null <-
          fitpt(
            expr,
            cellanno,
            pseudotime,
            design,
            maxknotallowed = 10,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose,
            ncores = 1,
            model = mod.null,
            knotnum = fitres.full[[2]]
          )
        if (exists('fitres.full') & exists('fitres.null')) {
          return(list(fitres.full = fitres.full, fitres.null = fitres.null))
        } else {
          return(NULL)
        }
      } else {
        dn <- paste0(as.vector(design), collapse = '_')
        perdn <- dn
        while (perdn == dn) {
          perid <- sample(seq_len(nrow(design)))
          perdesign <- design[perid, , drop = FALSE]
          perdn <- paste0(as.vector(perdesign), collapse = '_')
        }
        row.names(perdesign) <- row.names(design)
        sampcell <-
          sample(seq_len(ncol(expr)), replace = TRUE) ## boostrap cells
        perexpr <- expr[, sampcell, drop = FALSE]
        percellanno <- cellanno[sampcell, , drop = FALSE]
        psn <- pseudotime[colnames(expr)]
        psn <- psn[sampcell]
        colnames(perexpr) <-
          percellanno[, 1] <-
          names(psn) <- paste0('cell_', seq_len(length(psn)))
        
        
        fitres.full <-
          fitpt(
            perexpr,
            percellanno,
            psn,
            perdesign,
            maxknotallowed = 10,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose,
            ncores = ncores,
            model = mod.full
          )
        fitres.null <-
          fitpt(
            perexpr,
            percellanno,
            psn,
            perdesign,
            maxknotallowed = 10,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose = verbose,
            ncores = ncores,
            model = mod.null,
            knotnum = fitres.full[[2]]
          )
        if (exists('fitres.full') & exists('fitres.null')) {
          return(list(fitres.full = fitres.full, fitres.null = fitres.null))
        } else {
          return(NULL)
        }
      }
    }
  }
