#' Perform the model 0 (beta0 for intercept base only) fitting.
#'
#' This function is used to perform the model 0 (beta0 for intercept base only) fitting.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list containing parameters in the model 0 for genes.
#' @export
#' @import stats Matrix splines parallel matrixcalc
#' @param expr hdf5 path dir (including file name).
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param targetgene the gene names need to be fitted. if NULL, then read in h5 file "feature" of the first sample. 

#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @examples
#' data(mandata)
#' a = fitpt.m0(expr = mandata$expr, cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, EMmaxiter=5, EMitercutoff=10)

fitpt_m0_h5 <- function(expr, pseudotime, design, targetgene=NULL, boot=NULL, EMmaxiter=100, EMitercutoff=0.05) {
  # set.seed(12345)
  samp <- h5ls(expr,recursive=F)$name
  sname <- sapply(samp,function(s) as.vector(h5read(expr,paste0(s,'/barcode'))))
  if (!is.null(boot)) sname <- sapply(sname,function(i) boot[boot[,1] %in% i,2],simplify = F,USE.NAMES = T)
  gn <- h5read(expr,paste0(samp[1],'/feature'))
  if (is.null(targetgene)) targetgene <- gn
  design = as.matrix(design)
  exprreadfunc <- function(s,geneid) { ## s: the sample to read. geneid: numeric, the gene to read
    e <- h5read(expr,paste0(s,'/expr'),index=list(geneid,NULL)) ## read in the sample's gene expression on the geneid (a subset of genes)
    g <- as.vector(h5read(expr,paste0(s,'/feature'))[geneid]) ## read in the gene names
    rownames(e) <- g ##
    # if (!is.null(boot)) { ## if boot == NULL, that is, iter = 1. otherwise, redo this step: boostrap cell and rename the cells
    #   b <- as.vector(h5read(expr,paste0(s,'/barcode')))
    #   colnames(e) <- b
    #   bid <- which(boot[,1] %in% b) ## boot has all samples' cell names. bid gets the cells of this sample
    #   e <- e[,boot[bid,1],drop=F] ## get the expr
    #   colnames(e) <- boot[bid,2] ## rename the cells
    #   e[,sname[[s]],drop=F]  ## sname is the new (permuted) cell names
    # } else {
    #   e
    # }
    if (!is.null(boot)) { ## if boot == NULL, that is, iter = 1. otherwise, redo this step: boostrap cell and rename the cells
      b <- as.vector(h5read(expr,paste0(s,'/barcode')))
      bid <- which(boot[,1] %in% b) ## boot has all samples' cell names. bid gets the cells of this sample
      e <- e[,match(boot[bid,1],b)[match(sname[[s]],boot[bid,2])],drop=F] ## get the expr
      colnames(e) <- sname[[s]]
    }
    e
  }
  
  
  as <- names(sname)
  cn <- sapply(sname,length)
  
  ## initialize
  
  indfit <- sapply(as,function(s) {
    rowMeans(exprreadfunc(s,match(targetgene,gn)))
  },simplify = F)
  B <- colSums(do.call(rbind,indfit) * cn[names(indfit)])/sum(cn)
  
  s2 <- matrix(sapply(as,function(s) {
    tmp <- exprreadfunc(s,match(targetgene,gn))-indfit[[s]]
    n <- ncol(tmp)
    (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
  }),nrow=length(targetgene),dimnames=list(targetgene,as))
  alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
  eta <- (alpha-1)*rowMeans(s2)
  
  diffindfit <- sapply(as,function(s) {
    indfit[[s]] - B
  })
  s2[s2 < 0.001] <- 0.001
  omega <- rowMeans(diffindfit^2/s2)
  
  iter <- 0
  gidr <- targetgene
  oldpara <- list(beta = B, alpha = alpha, eta = eta, omega = omega)
  
  while (iter < EMmaxiter && length(gidr) > 0) {
    expr_phibx <- sapply(as,function(s) {
      exprreadfunc(s,match(gidr,gn))-B[gidr]
    },simplify = F)
    
    L <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]] * expr_phibx[[s]])
    },simplify = F)
    
    Jsolve <- matrix(sapply(as,function(s) {
      1/(cn[s] + 1/omega[gidr])
    }), nrow = length(gidr), dimnames = list(gidr, as))
    
    K <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]])
    },simplify = F)
    
    JK <- sapply(as,function(s) {
      Jsolve[,s,drop=F] * K[[s]]  ## debug here !!
    },simplify = F)
    
    L2eKJK <- sapply(as,function(s) {
      2*eta[gidr] + L[[s]] - K[[s]] * JK[[s]]
    },simplify = F)
    
    A <- matrix(sapply(as,function(s) {
      log(L2eKJK[[s]]/2)-digamma(alpha[gidr]+cn[s]/2)
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    N <- matrix(sapply(as,function(s) {
      (2*alpha[gidr]+cn[s])/L2eKJK[[s]]
    }),nrow=length(gidr),dimnames=list(gidr,as))
    
    
    B1 <- rowSums(matrix(sapply(as, function(s){
      N[,s] * cn[s]
    }), ncol = length(as), dimnames = list(gidr, as)))  ## length(gidr) * length(as) debug here !!!
    B2 <- rowSums(matrix(sapply(as, function(s){   ## debug here !!!
      N[,s,drop=F] * colSums(t(exprreadfunc(s,match(gidr,gn))) - rep(JK[[s]],each=cn[s]))
    }), ncol = length(as), dimnames = list(gidr, as)))  ## vector of length(gidr) debug here !!!
    B[gidr] <- B2/B1
    
    ## M -step:
    omega[gidr] <- rowMeans(matrix(sapply(as,function(s) {
      Jsolve[,s,drop=F] + N[,s,drop=F]*JK[[s]] * JK[[s]]  ## debug here !!!
    }), ncol=length(as), dimnames = list(gidr, as)))
    
    eta[gidr] <- sapply(gidr, function(g) {
      meanN <- mean(N[g, ])
      meanA <- mean(A[g, ])
      f <- function(eta) digamma(eta * meanN) - log(eta) + meanA
      f_lower <- f(1e-10)
      f_upper <- f(1e10)
      if (f_lower * f_upper > 0) {
        result <- tryCatch(
          uniroot(f, c(1e-20, 1e20), tol = 1e-12)$root,
          error = function(e) {
            warning(paste0("No root found for gene: ", g))
            NA
          }
        )
      } else {
        result <- uniroot(f, c(1e-10, 1e10))$root
      }
      result
    })
    alpha[gidr] <- eta[gidr] * rowMeans(N)
    
    para <- list(beta = B, alpha = alpha, eta = eta, omega = omega)
    paradiff <- sapply(names(para),function(s) {
      if (is.vector(para[[s]])) {
        abs(para[[s]]-oldpara[[s]])/abs(oldpara[[s]])
      } else {
        rowSums(abs(para[[s]]-oldpara[[s]]))/rowSums(abs(oldpara[[s]]))
      }
    })
    if (is.vector(paradiff)) paradiff <- matrix(paradiff,nrow=1,dimnames=list(gid,NULL))    
    paradiff <- apply(paradiff,1,max)
    gidr <- names(which(paradiff > EMitercutoff))
    oldpara <- para
    
    iter <- iter + 1
    
    rm(list = c('L','Jsolve', 'K'))
  }
  
  
  ll <- rowSums(matrix(sapply(as,function(s) {
    expr_phibx <- sapply(as,function(s) {
      exprreadfunc(s,match(targetgene,gn))-B[targetgene]
    },simplify = F)
    
    L <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]] * expr_phibx[[s]])
    },simplify = F)
    
    Jsolve <- matrix(sapply(as,function(s) {
      1/(cn[s] + 1/omega[targetgene])
    }), nrow = length(targetgene), dimnames = list(targetgene, as))
    
    K <- sapply(as,function(s) {
      rowSums(expr_phibx[[s]])
    },simplify = F)
    
    JK <- sapply(as,function(s) {
      Jsolve[,s,drop=F] * K[[s]]  ## debug here !!
    },simplify = F)
    
    L2eKJK <- sapply(as,function(s) {
      2*eta[targetgene] + L[[s]] - K[[s]] * JK[[s]]
    },simplify = F)
    dv <- omega[targetgene]/Jsolve[targetgene,s,drop=F]
    alpha[targetgene]*log(2*eta[targetgene])+lgamma(cn[s]/2+alpha[targetgene])-cn[s]*log(pi)/2-lgamma(alpha[targetgene])-log(dv)/2-(cn[s]/2+alpha[targetgene])*log(L2eKJK[[s]])
  }),nrow=length(targetgene),dimnames=list(targetgene,NULL)))
  
  allres  <- list(beta = B, alpha = alpha, eta = eta, omega = omega, ll = ll)
  
  para <- list()
  
  for (j in targetgene) {
    para[[j]] <- list(beta=allres[[1]][j],
                      alpha=allres[[2]][j],
                      eta=allres[[3]][j],
                      omega=allres[[4]][j],
                      ll=allres[[5]][j])
  }
  
  return(list(parameter=para))
}


