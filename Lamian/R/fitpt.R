#' Perform the model fitting.
#'
#' This function is used to perform the model fitting.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list containing parameters and number of knots that obtained from the model for genes.
#' @export
#' @import stats Matrix splines matrixcalc parallel
#' @importFrom matrixcalc %s%
#' @param expr gene by cell expression matrix. The expression values should have been library-size-normalized and log-transformed. They can either be imputed or non-imputed.
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param maxknotallowed a numeric number. Max number of knots applied in the B-spline fitting.
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param verbose logical. If TRUE, print intermediate information.
#' @param ncores number of cores used when performing the spline fitting for genes.
#' @param model a numeric number (1,2,or 3) indicating which model from the nesting models will be used for fitting.  Model 0 is implemented by the function fitpt.m0().
#' @param knotnum If NULL (default), this function will automatically select the optimal number of knots for each gene. If specified by a numeric vector whose names are gene names, this function will used the specified number of knots for those genes. This argument is used to speed up the fitting if the number of knots has been known.
#' @examples
#' data(mandata)
#' a = fitpt(expr = mandata$expr, cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, maxknotallowed=5, EMmaxiter=10, EMitercutoff=10, verbose=FALSE, ncores=1, model = 1)
fitpt <-
  function(expr,
           cellanno,
           pseudotime,
           design,
           maxknotallowed = 10,
           EMmaxiter = 100,
           EMitercutoff = 0.1,
           verbose = FALSE,
           ncores = 1,
           model = 3,
           knotnum = NULL) {
    pseudotime <- pseudotime[colnames(expr)]
    cellanno <- cellanno[match(colnames(expr), cellanno[, 1]), ]
    sname <-
      sapply(row.names(design), function(i)
        cellanno[cellanno[, 2] == i, 1], simplify = FALSE)
    design = as.matrix(design)
    
    philist <- lapply(0:maxknotallowed, function(num.knot) {
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
    names(philist) <- as.character(0:maxknotallowed)
    
    maxknot <- 0
    testpos <- 1
    
    while (testpos & maxknot < maxknotallowed) {
      maxknot <- maxknot + 1
      phi <- philist[[as.character(maxknot)]]
      testpos <- mean(sapply(names(sname), function(ss) {
        matrixcalc::is.positive.definite(crossprod(phi[sname[[ss]], , drop = FALSE]))
      })) == 1
    }
    maxknot <- maxknot - 1
    sexpr <-
      sapply(names(sname), function(ss)
        expr[, sname[[ss]], drop = FALSE], simplify = FALSE)
    
    ## automatically select knotnum for each gene
    if (is.null(knotnum)) {
      bicfunc <- function(num.knot) {
        phi <- philist[[as.character(num.knot)]]
        ll <- sapply(names(sname), function(ss) {
          phiss <- phi[sname[[ss]], , drop = FALSE]
          dif2 <-
            sexpr[[ss]] - sexpr[[ss]] %*% (phiss %*% base::chol2inv(base::chol(crossprod(phiss)))) %*% t(phiss)
          dif2 <- rowSums(dif2 * dif2)
          s2 <- dif2 / (length(sname[[ss]]) - ncol(phi))
          log(2 * pi * s2) * nrow(phiss) + dif2 / s2
        })
        if (is.vector(ll)) {
          sum(ll, na.rm = TRUE) + log(nrow(phi)) * ((ncol(phi) + 1) * sum(!is.na(ll)))
        } else {
          rowSums(ll, na.rm = TRUE) + log(nrow(phi)) * ((ncol(phi) + 1) * rowSums(!is.na(ll)))
        }
      }
      
      if (ncores != 1) {
        bic <- mclapply(0:maxknot, bicfunc, mc.cores = ncores)
        bic <- do.call(cbind, bic)
      } else {
        bic <- sapply(0:maxknot, bicfunc)
      }
      
      rm('sexpr')
      if (is.vector(bic)) {
        knotnum <- c(0:maxknot)[which.min(bic)]
      } else {
        knotnum <- c(0:maxknot)[apply(bic, 1, which.min)]
      }
      
      names(knotnum) <- rownames(expr)
    } else {
      knotnum <- knotnum[rownames(expr)]
    }
    
    sfit <- function(num.knot) {
      #print(paste0('num.knot = ', num.knot))
      gid <- names(which(knotnum == num.knot))
      #print(gid)
      sexpr <-
        expr[gid, , drop = FALSE] ## !!! double check, should be list len =S
      phi <- philist[[as.character(num.knot)]]
      phicrossprod <- apply(phi, 1, tcrossprod)
      phicrossprod <-
        sapply(names(sname), function(ss)
          phicrossprod[, sname[[ss]]], simplify = FALSE)
      phi <-
        sapply(names(sname), function(ss)
          phi[sname[[ss]], ], simplify = FALSE)
      
      if (model == 0) {
        xs <- sapply(row.names(design), function(i) {
          as.matrix(1, nrow = 1, ncol = 1)
        }, simplify = FALSE)
        phi <- sapply(phi, function(i) {
          i[, 1, drop = FALSE]
        }, simplify = FALSE)
      } else if (model == 1) {
        xs <- sapply(row.names(design), function(i) {
          kronecker(diag(num.knot + 4), design[i, 1, drop = FALSE])
        }, simplify = FALSE)
      } else if (model == 2) {
        # xs <- sapply(row.names(design), function(i) {
        #   tmp <- kronecker(diag(num.knot + 4), design[i, ])
        #   tmp <- tmp[-seq(4, nrow(tmp), 2), ]
        # }, simplify = FALSE)
        xs <- sapply(row.names(design), function(i) {
          ## change X
          tmp <- kronecker(diag(num.knot + 4), c(1, 1))
          tmp <- tmp[-seq(4, nrow(tmp), 2),]
          tmp[1, ] <- design[i, 2]
          tmp
        }, simplify = FALSE)
      } else if (model == 3) {
        xs <- sapply(row.names(design), function(i) {
          kronecker(diag(num.knot + 4), design[i,])
        }, simplify = FALSE)
      }
      
      as <- names(phi)
      nb <- ncol(phi[[1]])
      cn <- sapply(as, function(s)
        nrow(phi[[s]]))
      
      phiphi <- sapply(as, function(s) {
        t(phi[[s]]) %*% phi[[s]]
      }, simplify = FALSE)
      phiX <- sapply(as, function(s) {
        phi[[s]] %*% t(xs[[s]])
      }, simplify = FALSE)
      
      phiXTphiX <- sapply(as, function(s) {
        t(phiX[[s]]) %*% (phiX[[s]])
      }, simplify = FALSE)
      
      ## initialize
      B1 <- Reduce('+', phiXTphiX)
      B2 <-
        Reduce('+', sapply(as, function(s)
          sexpr[, cellanno[, 2] == s, drop = FALSE] %*% phiX[[s]], simplify = FALSE))
      B <- B2 %*% solve(B1)
      
      indfit <- sapply(as, function(s) {
        sexpr[, cellanno[, 2] == s, drop = FALSE] %*% (phi[[s]] %*% base::chol2inv(base::chol(crossprod(phi[[s]]))))
      }, simplify = FALSE)
      
      s2 <- matrix(sapply(as, function(s) {
        tmp <-
          sexpr[, cellanno[, 2] == s, drop = FALSE] - indfit[[s]] %*% t(phi[[s]])
        n <- ncol(tmp)
        (rowMeans(tmp * tmp) - (rowMeans(tmp)) ^ 2) * (n) / (n - 1)
      }),
      nrow = length(gid),
      dimnames = list(gid, as))
      alpha <- 1 / apply(s2, 1, var) * rowMeans(s2) ^ 2 + 2
      eta <- (alpha - 1) * rowMeans(s2)
      
      diffindfit <- sapply(as, function(s) {
        indfit[[s]] - B %*% xs[[s]]
      }, simplify = FALSE)
      
      omega <- t(sapply(rownames(sexpr), function(g) {
        m <- sapply(diffindfit, function(i)
          i[g, , drop = FALSE])
        m <- m - rowMeans(m)
        tcrossprod(m) / (ncol(m) - 1) + diag(nrow(m)) * 0.01
      }))
      
      iter <- 0
      gidr <- rownames(sexpr)
      all <-
        matrix(
          -Inf,
          nrow = nrow(sexpr),
          ncol = 1,
          dimnames = list(rownames = gidr)
        )
      #etalist <- alphalist <- omegalist <- Nlist <- Jslist <- list()
      while (iter < EMmaxiter && length(gidr) > 0) {
        oinv <- sapply(gidr, function(g) {
          base::chol2inv(base::chol(matrix(omega[g, ], nrow = nb)))
        }, simplify = FALSE)
        
        omegadet <- sapply(gidr, function(g) {
          log(det(matrix(omega[g, ], nrow = nb))) / 2
        })
        
        ll <- NULL
        for (s in as) {
          sexpr_phibx <-
            sexpr[gidr, cellanno[, 2] == s, drop = FALSE] - B[gidr, ] %*% t(phiX[[s]])
          
          L <- rowSums(sexpr_phibx * sexpr_phibx)
          
          Jchol <- sapply(gidr, function(g) {
            base::chol(phiphi[[s]] + oinv[[g]])
          }, simplify = FALSE)
          
          Jsolve <- sapply(gidr, function(g) {
            base::chol2inv(Jchol[[g]])
          })
          
          K <- tcrossprod(t(phi[[s]]), sexpr_phibx)
          
          L2eKJK <-
            2 * eta[gidr] + L - colSums(K[rep(seq_len(nb), nb), , drop = FALSE] * K[rep(seq_len(nb), each =
                                                                                          nb), , drop = FALSE] * Jsolve)
          A <- log(L2eKJK / 2) - digamma(alpha[gidr] + cn[s] / 2)
          N <- (2 * alpha[gidr] + cn[s]) / L2eKJK
          
          JK <-
            t(rowsum((Jsolve * K[rep(seq_len(nb), nb), , drop = FALSE]), rep(seq_len(nb), each =
                                                                               nb)))
          
          logdv <-
            -2 * colSums(log(sapply(Jchol, as.vector)[seq(1, nb * nb, nb + 1), , drop =
                                                        FALSE]))
          if (is.null(ll)) {
            ll <-
              alpha[gidr] * log(2 * eta[gidr]) + lgamma(cn[s] / 2 + alpha[gidr]) - cn[s] *
              log(pi) / 2 - lgamma(alpha[gidr]) - omegadet + logdv / 2 - (cn[s] / 2 +
                                                                            alpha[gidr]) * log(L2eKJK)
            B1 <- tcrossprod(N, as.vector(phiXTphiX[[s]]))
            rownames(B1) <- gidr
            B2 <-
              N * ((sexpr[gidr, cellanno[, 2] == s, drop = FALSE] - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
            omegalist <-
              t(Jsolve) + N * JK[, rep(seq_len(nb), nb)] * JK[, rep(seq_len(nb), each =
                                                                      nb)]
            sumA <- A
            sumN <- N
          } else {
            ll <-
              ll + alpha[gidr] * log(2 * eta[gidr]) + lgamma(cn[s] / 2 + alpha[gidr]) -
              cn[s] * log(pi) / 2 - lgamma(alpha[gidr]) - omegadet + logdv / 2 - (cn[s] /
                                                                                    2 + alpha[gidr]) * log(L2eKJK)
            B1 <- B1 + tcrossprod(N, as.vector(phiXTphiX[[s]]))
            B2 <-
              B2 + N * ((sexpr[gidr, cellanno[, 2] == s, drop = FALSE] - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
            omegalist <-
              omegalist + t(Jsolve) + N * JK[, rep(seq_len(nb), nb)] * JK[, rep(seq_len(nb), each =
                                                                                  nb)]
            sumA <- sumA + A
            sumN <- sumN + N
          }
        }
        
        np <- nrow(xs[[1]])
        B[gidr, ] <-
          t(sapply(gidr, function(g) {
            ## each column is a gene's all betas
            base::chol2inv(chol(matrix(B1[g, ], nrow = np))) %*% B2[g, ]
          }))
        
        omega[gidr, ] <- omegalist / length(as)
        
        rN <- sumN / length(as)
        rA <- sumA / length(as)
        eta[gidr] <- sapply(gidr, function(g) {
          meanN <- rN[g]
          meanA <- rA[g]
          # uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
          optim(
            eta[g],
            fn = function(eta) {
              (digamma(eta * meanN) - log(eta) + meanA) ^ 2
            },
            gr = function(eta) {
              2 * (digamma(eta * meanN) - log(eta) + meanA) * (trigamma(eta * meanN) *
                                                                 meanN - 1 / eta)
            },
            lower = 1e-10,
            method = 'L-BFGS-B'
          )$par
        })
        alpha[gidr] <- eta[gidr] * rN
        
        iter <- iter + 1
        
        llv <- all[, ncol(all)]
        llv[gidr] <- ll
        all <- cbind(all, llv)
        gidr <-
          names(which(all[, ncol(all)] - all[, ncol(all) - 1] > EMitercutoff))
        #etalist[[iter]] <- eta
        #alphalist[[iter]] <-alpha
        #omegalist[[iter]] <-omega
        #Nlist[[iter]] <- N
        #Jslist[[iter]] <- Jsolve
        rm(list = c('L', 'Jsolve', 'K'))
      }
      # print(table(apply(all,1,function(i) mean(diff(i) >= 0))))
      return(list(
        beta = B,
        alpha = alpha,
        eta = eta,
        omega = omega,
        logL = all
      ))
    }
    
    if (ncores != 1) {
      allres <- mclapply(unique(knotnum), sfit, mc.cores = ncores)
    } else {
      allres <- sapply(unique(knotnum), sfit, simplify = FALSE)
    }
    
    para <- list()
    for (i in seq_len(length(allres))) {
      for (j in row.names(allres[[i]][[1]])) {
        para[[j]] <- list(
          beta = allres[[i]][[1]][j, ],
          alpha = allres[[i]][[2]][j],
          eta = allres[[i]][[3]][j],
          omega = allres[[i]][[4]][j, ],
          ll = allres[[i]][[5]][j, ncol(allres[[i]][[5]])]
        )
      }
    }
    list(parameter = para[rownames(expr)], knotnum = knotnum[rownames(expr)])
  }
