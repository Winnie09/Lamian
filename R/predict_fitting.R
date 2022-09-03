#' Predict population fitting.
#'
#' This function will predict the population fitting for genes using the output from lamian_test().
#'
#' @param testObj the output object from lamian_test().
#' @param gene a vector of genes that need to do the prediction.
#' @param test.type One of c('Time', 'Variable').
#' @return a gene by cell (or pseudotime) expression matrix
#' @author Wenpin Hou <whou10@jhu.edu>
predict_fitting <-
  function(testObj,
           gene = NULL,
           test.type = 'time') {
    ## make the cells order according to pseudotime order
    ## fitting
    if ('expr.ori' %in% names(testObj)){
      expr <- testObj$expr.ori
    } else {
      expr <- testObj$expr
    }
      
    knotnum = testObj$knotnum[gene]
    design = testObj$design
    cellanno = testObj$cellanno
    pseudotime = testObj$pseudotime[colnames(expr)]
    if (is.null(gene))
      gene <- rownames(expr)
    philist <- lapply(sort(unique(knotnum)), function(num.knot) {
      if (num.knot == 0) {
        # phi <- cbind(1,bs(pseudotime))
        phi <- bs(pseudotime, intercept = TRUE)
      } else {
        knots = seq(min(pseudotime), max(pseudotime), length.out = num.knot + 2)[2:(num.knot +
                                                                                      1)]
        # phi <- cbind(1,bs(pseudotime,knots = knots))
        phi <- bs(pseudotime, knots = knots, intercept = TRUE)
      }
    })
    
    names(philist) <- as.character(sort(unique(knotnum)))
    as <- row.names(design)
    sname <-
      sapply(as, function(i)
        cellanno[cellanno[, 2] == i, 1], simplify = FALSE)
    
    pred <- lapply(unique(knotnum), function(num.knot) {
      genesub <- names(knotnum)[knotnum == num.knot]
      B <- t(sapply(genesub, function(g) {
        testObj$parameter[[g]]$beta
      }))
      
      omega <- t(sapply(genesub, function(g) {
        testObj$parameter[[g]]$omega
      }))
      
      phi <- philist[[as.character(num.knot)]]
      phi <-
        sapply(as, function(ss)
          phi[sname[[ss]], ], simplify = FALSE)
      
      if (toupper(test.type)  == 'TIME') {
        xs <- sapply(row.names(design), function(i) {
          kronecker(diag(num.knot + 4), design[i, 1, drop = FALSE])
        }, simplify = FALSE)
      } else if (toupper(test.type) == 'VARIABLE') {
        xs <- sapply(row.names(design), function(i) {
          kronecker(diag(num.knot + 4), design[i,])
        }, simplify = FALSE)
      }
      phiphi <- sapply(as, function(s) {
        t(phi[[s]]) %*% phi[[s]]
      }, simplify = FALSE)
      phiX <- sapply(as, function(s) {
        phi[[s]] %*% t(xs[[s]])
      }, simplify = FALSE)
      
      predtmp <- sapply(as, function(s) {
        sexpr <- expr[genesub, , drop = FALSE]
        sexpr_phibx <-
          sexpr[genesub, cellanno[, 2] == s, drop = FALSE] - B[genesub, ] %*% t(phiX[[s]])
        
        nb <- num.knot + 4
        oinv <- sapply(genesub, function(g) {
          chol2inv(chol(matrix(omega[g, , drop = FALSE], nrow = nb)))
        }, simplify = FALSE)
        
        Jchol <- sapply(genesub, function(g) {
          chol(phiphi[[s]] + oinv[[g]])
        }, simplify = FALSE)
        
        Jsolve <- sapply(genesub, function(g) {
          chol2inv(Jchol[[g]])
        })
        K <- tcrossprod(t(phi[[s]]), sexpr_phibx)
        JK <-
          rowsum((Jsolve * K[rep(seq_len(nb), nb), , drop = FALSE]), rep(seq_len(nb), each =
                                                                           nb)) ## u's poterior mean
        t(phi[[s]] %*% JK)
      }, simplify = FALSE)
      predtmp <- do.call(cbind, predtmp)
    })
    pred <- do.call(rbind, pred)
    pred <- pred[gene, colnames(expr), drop = FALSE]
    if ('populationFit' %in% names(testObj)){
      populationFit = testObj$populationFit
    } else{
      populationFit <- getPopulationFit(testObj, gene, type = testObj$test.type)
    }
      
    if (toupper(test.type) == 'TIME') {
      return(pred + populationFit[gene, , drop = FALSE])
    } else {
      l <- lapply(populationFit, function(i) {
        pred + i[gene, pseudotime , drop = FALSE]
      })
      return(l)
    }
    
  }
