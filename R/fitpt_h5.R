#' Perform the model fitting.
#'
#' This function is used to perform the model fitting.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list containing parameters and number of knots that obtained from the model for genes.
#' @export
#' @import stats Matrix splines matrixcalc parallel rhdf5
#' @importFrom matrixcalc %s%
#' @param expr hdf5 file path and name. 
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param testvar a numeric number indicating the column in the design matrix that needs to be tested while controlling for other columns (not intercept). Its value when calls the function will be used. e.g. testvar = 2 means the second column in the design needs to be tested.
#' @param targetgene the gene names need to be fitted. if NULL, then read in h5 file "feature" of the first sample. 
#' @param maxknotallowed a numeric number. Max number of knots applied in the B-spline fitting.
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param ncores number of cores used when performing the spline fitting for genes.
#' @param model a numeric number (1,2,or 3) indicating which model from the nesting models will be used for fitting.  Model 0 is implemented by the function fitpt.m0().
#' @param knotnum If NULL (default), this function will automatically select the optimal number of knots for each gene. If specified by a numeric vector whose names are gene names, this function will used the specified number of knots for those genes. This argument is used to speed up the fitting if the number of knots has been known.
#' @examples
#' data(mandata)
#' a = fitpt(expr = mandata$expr, pseudotime = mandata$pseudotime, design = mandata$design, maxknotallowed=5, EMmaxiter=10, EMitercutoff=10, ncores=1, model = 1)
fitpt_h5 <- function(expr, 
                     pseudotime, 
                     design, 
                     testvar=testvar, 
                     targetgene=NULL,
                     boot=NULL, 
                     maxknotallowed=10, 
                     EMmaxiter=100, 
                     EMitercutoff=0.05, 
                     ncores=1, 
                     model = 3, 
                     knotnum = NULL) {
  samp <- unname(rownames(design))
  sname <- sapply(samp,function(s) as.vector(h5read(expr,paste0(s,'/barcode')))) ## each sample's cells, h5read is to read the expr in one sample
  if (!is.null(boot)) sname <- sapply(sname,function(i) boot[boot[,1] %in% i,2],simplify = F,USE.NAMES = T) ## use new cell names. boot[,1] old names, boot[,2] new (permuted) cell names
  gn <- h5read(expr,paste0(samp[1],'/feature')) ## read in gene names. since all samples have the same gene names, only need to read in the first sample
  if (is.null(targetgene)) targetgene <- gn
  design = as.matrix(design)
  ## most challenging part
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
  
  philist <- lapply(0:maxknotallowed,function(num.knot) {
    if (num.knot==0) {
      # phi <- cbind(1,bs(pseudotime))
      phi <- bs(pseudotime, intercept = TRUE)
    } else {
      knots = seq(min(pseudotime),max(pseudotime),length.out=num.knot+2)[2:(num.knot+1)]
      # phi <- cbind(1,bs(pseudotime,knots = knots))  
      phi <- bs(pseudotime,knots = knots, intercept = TRUE)
    }
  })
  names(philist) <- as.character(0:maxknotallowed)
  
  maxknot <- 0
  testpos <- 1
  
  while (testpos & maxknot < maxknotallowed) {
    maxknot <- maxknot + 1
    phi <- philist[[as.character(maxknot)]]
    testpos <- mean(sapply(names(sname), function(ss) {
      matrixcalc::is.positive.definite(crossprod(phi[sname[[ss]],,drop=F]))
    })) == 1
  }
  maxknot <- maxknot - 1
  
  ## automatically select knotnum for each gene
  if (is.null(knotnum)){
    bicfunc <- function(num.knot) {
      phi <- philist[[as.character(num.knot)]]
      ll <- sapply(names(sname), function(ss) {
        sexpr <- exprreadfunc(ss,match(targetgene,gn))
        phiss <- phi[sname[[ss]],,drop=F]
        dif2 <- sexpr - sexpr %*% (phiss %*% chol2inv(chol(crossprod(phiss)))) %*% t(phiss)
        dif2 <- rowSums(dif2 * dif2)
        s2 <- dif2/(length(sname[[ss]])-ncol(phi))
        log(2*pi*s2)*nrow(phiss) + dif2/s2
      })
      if (is.vector(ll)){
        sum(ll,na.rm=T) + log(nrow(phi))*((ncol(phi)+1)*sum(!is.na(ll)))
      } else {
        rowSums(ll,na.rm=T) + log(nrow(phi))*((ncol(phi)+1)*rowSums(!is.na(ll)))
      }    
    }
    
    if (ncores!=1) {
      bic <- mclapply(0:maxknot,bicfunc,mc.cores=ncores)
      bic <- do.call(cbind,bic)
    } else {
      bic <- sapply(0:maxknot,bicfunc)
    }
    
    if (is.vector(bic)){
      knotnum <- c(0:maxknot)[which.min(bic)]
    } else {
      knotnum <- c(0:maxknot)[apply(bic,1,which.min)]
    }
    
    names(knotnum) <- targetgene
  } else {
    knotnum <- knotnum[targetgene]
  }
  
  sfit <- function(num.knot) {
    gid <- names(which(knotnum==num.knot))
    phi <- philist[[as.character(num.knot)]]
    phicrossprod <- apply(phi,1,tcrossprod)
    phicrossprod <- sapply(names(sname),function(ss) phicrossprod[,sname[[ss]]],simplify = F)
    phi <- sapply(names(sname),function(ss) phi[sname[[ss]],],simplify = F)
    
    if (model == -1){
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), 1)
      }, simplify = F)
    } else if (model == 1) {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i,-testvar])
      }, simplify = F)
    } else if (model == 2) {
      xs <- sapply(row.names(design), function(i) {  ## change X
        rbind(design[i,testvar],kronecker(diag(num.knot + 4), design[i,-testvar]))
      }, simplify = F)
    } else if (model == 3) {
      xs <- sapply(row.names(design), function(i) {
        kronecker(diag(num.knot + 4), design[i, ])
      }, simplify = F)
    }
    
    as <- names(phi)
    nb <- ncol(phi[[1]])
    cn <- sapply(as,function(s) nrow(phi[[s]]))
    
    phiphi <- sapply(as,function(s) {
      t(phi[[s]]) %*% phi[[s]]
    },simplify = F)
    phiX <- sapply(as, function(s){
      phi[[s]] %*% t(xs[[s]])
    },simplify = F)
    
    phiXTphiX <- sapply(as, function(s){
      t(phiX[[s]]) %*% (phiX[[s]])
    },simplify = F)
    
    ## initialize
    B1 <- Reduce('+',phiXTphiX)
    B2 <- indfit <- s2 <- list()
    for (s in as) {
      sexpr <- exprreadfunc(s,match(gid,gn))
      B2[[s]] <- sexpr %*% phiX[[s]]
      indfit[[s]] <- sexpr %*% (phi[[s]] %*% chol2inv(chol(crossprod(phi[[s]]))))
      tmp <- sexpr-indfit[[s]] %*% t(phi[[s]])
      n <- ncol(tmp)
      s2[[s]] <- (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
    }
    B2 <- Reduce('+',B2)
    s2 <- do.call(cbind,s2)
    
    B <- B2 %*% solve(B1)
    alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
    eta <- (alpha-1)*rowMeans(s2)
    
    diffindfit <- sapply(as,function(s) {
      indfit[[s]] - B %*% xs[[s]]
    },simplify = F)
    
    omega <- Reduce('+',sapply(as,function(s) {
      diffindfit[[s]][,rep(1:nb,nb),drop=FALSE] * diffindfit[[s]][,rep(1:nb,each=nb),drop=FALSE]/max(0.001,s2[,s])
    },simplify = F))
    omega <- omega/length(as)
    for (g in 1:nrow(omega))
      #    if (!is.positive.definite(matrix(omega[g,],nrow=nb))) 
      omega[g,] <- as.vector(diag(pmax(1e-5,diag(matrix(omega[g,],nrow=nb)))))
    
    iter <- 0
    gidr <- gid
    oldpara <- list(beta = B, alpha = alpha, eta = eta, omega = omega)
    while (iter < EMmaxiter && length(gidr) > 0) {
      
      oinv <- sapply(gidr,function(g) {
        chol2inv(chol(matrix(omega[g,],nrow=nb)))
      },simplify = F)
      
      flag <- 0
      for (s in as) {
        sexpr <- exprreadfunc(s,match(gidr,gn))
        sexpr_phibx <- sexpr - B[gidr,] %*% t(phiX[[s]])
        
        L <- rowSums(sexpr_phibx * sexpr_phibx)
        
        Jchol <- sapply(gidr,function(g) {
          chol(phiphi[[s]] + oinv[[g]])
        },simplify = F)
        
        Jsolve <- sapply(gidr,function(g) {
          chol2inv(Jchol[[g]])
        })
        
        K <- tcrossprod(t(phi[[s]]),sexpr_phibx)
        
        L2eKJK <- 2*eta[gidr]+L-colSums(K[rep(1:nb,nb),,drop=FALSE] * K[rep(1:nb,each=nb),,drop=FALSE] * Jsolve)
        A <- log(L2eKJK/2)-digamma(alpha[gidr]+cn[s]/2)
        N <- (2*alpha[gidr]+cn[s])/L2eKJK
        
        JK <- t(rowsum((Jsolve*K[rep(1:nb,nb),,drop=FALSE]),rep(1:nb,each=nb)))
        
        if (flag==0) {
          B1 <- tcrossprod(N,as.vector(phiXTphiX[[s]]))
          rownames(B1) <- gidr
          B2 <- N * ((sexpr - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
          omegalist <- t(Jsolve) + N*JK[,rep(1:nb,nb)] * JK[,rep(1:nb,each=nb)]
          sumA <- A
          sumN <- N
        } else {
          B1 <- B1 + tcrossprod(N,as.vector(phiXTphiX[[s]]))
          B2 <- B2 + N * ((sexpr - t(phi[[s]] %*% t(JK))) %*% phiX[[s]])
          omegalist <- omegalist + t(Jsolve) + N*JK[,rep(1:nb,nb)] * JK[,rep(1:nb,each=nb)]
          sumA <- sumA + A
          sumN <- sumN + N
        }
        flag <- 1
      }
      
      np <- nrow(xs[[1]])
      B[gidr,] <- t(sapply(gidr, function(g){  ## each column is a gene's all betas
        chol2inv(chol(matrix(B1[g,],nrow=np))) %*% B2[g,]
      }))
      
      omega[gidr,] <- omegalist/length(as)
      
      rN <- sumN/length(as)
      rA <- sumA/length(as)
      eta[gidr] <- sapply(gidr, function(g) {
        meanN <- rN[g]
        meanA <- rA[g]
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
      alpha[gidr] <- eta[gidr] * rN
      
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
    }
    oinv <- sapply(gid,function(g) {
      chol2inv(chol(matrix(omega[g,],nrow=nb)))
    },simplify = F)
    omegadet <- sapply(gid,function(g) {
      log(det(matrix(omega[g,],nrow=nb)))/2
    })
    ll <- rowSums(matrix(sapply(as,function(s) {
      sexpr_phibx <- exprreadfunc(s,match(gid,gn))-B[gid,] %*% t(phiX[[s]])
      L <- rowSums(sexpr_phibx * sexpr_phibx)
      K <- tcrossprod(t(phi[[s]]),sexpr_phibx)
      Jchol <- sapply(gid,function(g) {
        chol(phiphi[[s]] + oinv[[g]])
      },simplify = F)
      Jsolve <- sapply(gid,function(g) {
        chol2inv(Jchol[[g]])
      })
      L2eKJK <- 2*eta[gid]+L-colSums(K[rep(1:nb,nb),,drop=FALSE] * K[rep(1:nb,each=nb),,drop=FALSE] * Jsolve)
      logdv <- -2*colSums(log(sapply(Jchol,as.vector)[seq(1,nb*nb,nb+1),,drop=F]))
      alpha[gid]*log(2*eta[gid])+lgamma(cn[s]/2+alpha[gid])-cn[s]*log(pi)/2-lgamma(alpha[gid])-omegadet+logdv/2-(cn[s]/2+alpha[gid])*log(L2eKJK)
    }),nrow=length(gid),dimnames=list(gid,NULL)))
    return(list(beta = B, alpha = alpha, eta = eta, omega = omega,ll=ll))
  }
  
  if (ncores!=1) {
    allres <- mclapply(unique(knotnum),sfit,mc.cores=ncores)
  } else {
    allres <- sapply(unique(knotnum),sfit, simplify = FALSE)
  }
  
  para <- list()
  for (i in 1:length(allres)) {
    for (j in row.names(allres[[i]][[1]])) {
      para[[j]] <- list(beta=allres[[i]][[1]][j,],
                        alpha=allres[[i]][[2]][j],
                        eta=allres[[i]][[3]][j],
                        omega=allres[[i]][[4]][j,],
                        ll=unname(allres[[i]][[5]][j]))
    }
  }
  list(parameter=para[targetgene],knotnum=knotnum[targetgene])
}



