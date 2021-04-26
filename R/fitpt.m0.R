#' Perform the model 0 (beta0 for intercept base only) fitting.
#' 
#' This function is used to perform the model 0 (beta0 for intercept base only) fitting.
#' 
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list containing parameters in the model 0 for genes. 
#' @export
#' @import Matrix splines matrixcalc parallel
#' @param expr gene by cell expression matrix. The expression values should have been library-size-normalized and log-transformed. They can either be imputed or non-imputed.
#' @param cellanno a dataframe where the first column are cell names and second column are sample names. 
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples. 
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param verbose logical. If TRUE, print intermediate information.
#' @examples
#' Data(mandata)
#' a = fitpt.m0(expr = mandata$expr, cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, EMmaxiter=5, EMitercutoff=10, verbose=F)
fitpt.m0 <- function(expr, cellanno, pseudotime, design, EMmaxiter=100, EMitercutoff=0.1, verbose=F) {
  # set.seed(12345)
  pseudotime <- pseudotime[colnames(expr)]
  cellanno <- cellanno[match(colnames(expr),cellanno[,1]),]
  sname <- sapply(unique(cellanno[,2]),function(i) cellanno[cellanno[,2]==i,1],simplify = F)
  as <- names(sname)
  cn <- sapply(sname,length)
  
  ## initialize
  B <- rowMeans(expr)
  
  indfit <- sapply(as,function(s) {
    rowMeans(expr[, cellanno[,2]==s, drop=F])
  },simplify = F)
  
  s2 <- matrix(sapply(as,function(s) {
    tmp <- expr[, cellanno[,2]==s, drop=F]-indfit[[s]]
    n <- ncol(tmp)
    (rowMeans(tmp*tmp)-(rowMeans(tmp))^2)*(n)/(n-1)
  }),nrow=nrow(expr),dimnames=list(rownames(expr),as))
  alpha <- 1/apply(s2,1,var)*rowMeans(s2)^2+2
  eta <- (alpha-1)*rowMeans(s2)
  
  diffindfit <- sapply(as,function(s) {
    indfit[[s]] - B
  })
  
  omega <- apply(diffindfit,1,var)
  
  iter <- 0
  EMitercutoff <- 0
  gidr <- rownames(expr)
  all <- matrix(-Inf, nrow=nrow(expr),ncol=1,dimnames = list(rownames=gidr))
  etalist <- alphalist <- omegalist <- Nlist <- Jslist <- list()
  while (iter < EMmaxiter && length(gidr) > 0) {
    expr_phibx <- sapply(as,function(s) {
      expr[, cellanno[,2]==s, drop=F][gidr,,drop=F]-B[gidr]
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
    
    ll <- rowSums(matrix(sapply(as,function(s) {
      dv <- omega[gidr]/Jsolve[gidr,s,drop=F]
      alpha[gidr]*log(2*eta[gidr])+lgamma(cn[s]/2+alpha[gidr])-cn[s]*log(pi)/2-lgamma(alpha[gidr])-log(dv)/2-(cn[s]/2+alpha[gidr])*log(L2eKJK[[s]])
    }),nrow=length(gidr),dimnames = list(rownames=gidr,colnames=as)))
    
    ## -------------->
    
    
    B1 <- rowSums(matrix(sapply(as, function(s){
      N[,s] * cn[s]
    }), ncol = length(as), dimnames = list(gidr, as)))  ## length(gidr) * length(as) debug here !!!
    B2 <- rowSums(matrix(sapply(as, function(s){   ## debug here !!!
      N[,s,drop=F] * colSums(t(expr[ ,cellanno[,2]==s, drop=F][gidr,,drop=F]) - rep(JK[[s]],each=cn[s]))
    }), ncol = length(as), dimnames = list(gidr, as)))  ## vector of length(gidr) debug here !!!
    B[gidr] <- B2/B1
    
    ## M -step:
    omega[gidr] <- rowMeans(matrix(sapply(as,function(s) {
      Jsolve[,s,drop=F] + N[,s,drop=F]*JK[[s]] * JK[[s]]  ## debug here !!!
    }), ncol=length(as), dimnames = list(gidr, as)))
    
    eta[gidr] <- sapply(gidr,function(g) {
      meanN <- mean(N[g,])
      meanA <- mean(A[g,])
      uniroot(function(eta) {digamma(eta * meanN)-log(eta)+meanA},c(1e-10,1e10))$root
    })
    alpha[gidr] <- eta[gidr] * rowMeans(N)
    
    iter <- iter + 1
    
    llv <- all[,ncol(all)]
    llv[gidr] <- ll
    all <- cbind(all,llv)
    gidr <- names(which(all[,ncol(all)] - all[,ncol(all)-1] > EMitercutoff))
    etalist[[iter]] <- eta
    alphalist[[iter]] <-alpha
    omegalist[[iter]] <-omega
    Nlist[[iter]] <- N
    Jslist[[iter]] <- Jsolve
    rm(list = c('L','Jsolve', 'K'))
  }
  
  allres  <- list(beta = B, alpha = alpha, eta = eta, omega = omega, logL = all)
  para <- list()
  for (j in rownames(expr)) {
    para[[j]] <- list(beta=allres[[1]][j],
                      alpha=allres[[2]][j],
                      eta=allres[[3]][j],
                      omega=allres[[4]][j],
                      ll=allres[[5]][j,ncol(allres[[5]])])
  }
  
  return(list(parameter=para))
}


