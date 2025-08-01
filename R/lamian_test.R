#' Test function for all differential tests in lamian module 2
#'
#' This function is designed for doing different tests.
#'
#' @export
#' @author Wenpin Hou <wh2526@cumc.columbia.edu>
#' @import parallel splines devtools
#' @return a list of results
#' @param expr gene by cell expression matrix. Values are library-size-normalized log-transformed gene expression matrix. They can be either imputed or non-imputed. Zero-expression genes should have been filtered out.
#' @param cellanno a dataframe where the first column are cell names and second column are sample names.
#' @param pseudotime a numeric vector of pseudotime, and the names of this vector are the corresponding cell names.
#' @param design: a matrix. Number of rows should be the same as the number of unique samples. Rownames are sample names. First column is the intercept (all 1), second column is the covariate realization valuels for each of the samples.
#' @param testvar a numeric number indicating the column in the design matrix that needs to be tested while controlling for other columns (not intercept). Default is 2. testvar = 2 means the second column in the design needs to be tested.
#' @param permuiter an integer  indicating the number of permutations performed in the permutation test.
#' @param EMmaxiter an integer indicating the number of iterations in the EM algorithm.
#' @param EMitercutoff a numeric number indicating the log-likelihood cutoff applied to stop the EM algorithm
#' @param verbose.output logical. If TRUE, print intermediate information.
#' @param ncores number of cores for performing the permutation test.
#' @param test.type One of c('Time', 'Variable'). Case insensitive.
#' @param return.all.data logical. If TRUE (default), return all data including inputs.
#' @param overall.only logical. If TRUE (default), only test the overall fdr and skip the trend fdr and mean fdr. Default is FALSE.
#' @param test.method One of c('chisq', 'permutation). If 'permutation' (default), use the permutation test to identify genes. If 'chisq', use the chisq test to identify genes.
#' @param ncores.fit is the ncores for fitpt() or fitfunc()(essentially fitpt()) only. It only works when test.method = 'chisq'.
#' @param fix.all.zero logical. If TRUE (defalt), fix the issue of all zeros in any samples.
#' @param cutoff a numeric number to set the cutoff for the standard deviation of gene expression in any one of the samples. Only useful when fix.all.zero == TRUE.
#' @examples
#' data(mandata)
#' a = lamian_test(expr = mandata$expr[seq(1,3),], cellanno = mandata$cellanno, pseudotime = mandata$pseudotime, design = mandata$design, test.method = 'chisq', permuiter = 10, EMmaxiter = 10, EMitercutoff = 10, verbose.output = FALSE)

lamian_test <- function(expr, cellanno, pseudotime, design=NULL, testvar=2, permuiter=100, maxknotallowed = 10, EMmaxiter=100, EMitercutoff=0.05, verbose.output = FALSE, ncores=detectCores(), test.type='Time', fit.resolution = 1000, return.all.data = TRUE, overall.only = FALSE, test.method = 'permutation', ncores.fit = 1, fix.all.zero = TRUE, cutoff = 1e-5, sd.adjust = 1e-5) { 
  if (test.method == 'permutation') ncores.fit = 1
  set.seed(12345)
  cellanno = data.frame(Cell = as.character(cellanno[,1]), Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  expr <- expr[, names(pseudotime), drop = FALSE]
  cellanno <- cellanno[match(names(pseudotime), cellanno[,1]), ]
  design[,1] <- as.numeric(design[,1])
  design[,2] <- as.numeric(design[,2])
  design = as.matrix(design)
  
  if (fix.all.zero){
    sdm <- sapply(unique(cellanno[,2]),function(us) {
      tmp <- expr[,cellanno[,2]==us, drop=FALSE]
      m <- rowMeans(tmp)
      rowMeans(tmp*tmp)-m*m
    })<cutoff ## first version <20210407 is 0.
    gid <- which(rowSums(sdm) > 0)  ## identify if any genes have sd=0 expression in any one of the samples
    if (length(gid) > 0) { ## if yes, for those genes, add a white-noise with sd=1e-5 on the sample with sd=0.
      mask <- sdm[gid,rep(seq_len(ncol(sdm)),as.vector(table(cellanno[,2])[colnames(sdm)])),drop=F]
      colnames(mask) <- unlist(sapply(colnames(sdm),function(i) cellanno[cellanno[,2]==i,1]))
      expr[gid,] <- expr[gid,] + mask[,colnames(expr),drop=F] * matrix(rnorm(length(mask),sd=sd.adjust), nrow=length(gid)) ## before 20211126, sd.adjust = 1e-5 
      rm('mask')  
    } 
  }
  if (test.method == 'chisq'){
    if (toupper(test.type) == 'TIME'){
      res1 <-
        fitpt(
          expr,
          cellanno,
          pseudotime,
          design = design[, 1, drop = FALSE],
          testvar = testvar,
          maxknotallowed = maxknotallowed,
          maxknotallowed = 10,
          EMmaxiter = EMmaxiter,
          EMitercutoff = EMitercutoff,
          ncores = ncores.fit,
          model = 1
        ) 
      ll1 <- sapply(res1$parameter,function(i) i$ll)
      res0 <-
        fitpt.m0(expr,
                 cellanno,
                 pseudotime,
                 design[, 1, drop = FALSE],
                 EMmaxiter = EMmaxiter,
                 EMitercutoff = EMitercutoff) 
      ll0 <- sapply(res0[[1]], function(i) i$ll)
      paradiff10 <- sapply(res1[[1]], function(i) length(unlist(i[seq_len(4)]))) - sapply(res0[[1]], function(i) length(unlist(i[seq_len(4)])))
      pval.chisq.constantTest <- pchisq(2*(ll1-ll0),df=paradiff10,lower.tail = F)
      fdr.chisq.constantTest <- p.adjust(pval.chisq.constantTest, method='fdr')
      res <- data.frame(fdr.chisq.overall = fdr.chisq.constantTest, 
                        pval.chisq.overall = pval.chisq.constantTest, 
                        llr = ll1-ll0,
                        df.diff= paradiff10,
                        stringsAsFactors = FALSE)
      reslist = list(statistics = res,  parameter = res1$parameter, knotnum = res1$knotnum)  ## function return
    } else if (toupper(test.type) == 'VARIABLE'){
      res1 <-
        fitpt(
          expr,
          cellanno,
          pseudotime,
          design,
          testvar = testvar,
          maxknotallowed = maxknotallowed,
          EMmaxiter = EMmaxiter,
          EMitercutoff = EMitercutoff,
          ncores = ncores.fit,
          model = 1
        )## save 13%
      ll1 <- sapply(res1$parameter, function(i)
        i$ll)
      res2 <-
        fitpt(
          expr,
          cellanno,
          pseudotime,
          design,
          testvar = testvar,
          maxknotallowed = maxknotallowed,
          EMmaxiter = EMmaxiter,
          EMitercutoff = EMitercutoff,
          ncores = ncores.fit,
          model = 2,
          knotnum = res1[[2]]
        )## save 13%
      ll2 <- sapply(res2$parameter, function(i)
        i$ll)
      res3 <-
        fitpt(
          expr,
          cellanno,
          pseudotime,
          design,
          testvar = testvar,
          maxknotallowed = maxknotallowed,
          EMmaxiter = EMmaxiter,
          EMitercutoff = EMitercutoff,
          ncores = ncores.fit,
          model = 3,
          knotnum = res1[[2]]
        )
      ll3 <- sapply(res3$parameter, function(i)
        i$ll)
      
      paradiff31 <- sapply(res3$parameter,function(i) length(i$beta))-sapply(res1$parameter,function(i) length(i$beta))
      paradiff32 <- sapply(res3$parameter,function(i) length(i$beta))-sapply(res2$parameter,function(i) length(i$beta))
      paradiff21 <- sapply(res2$parameter,function(i) length(i$beta))-sapply(res1$parameter,function(i) length(i$beta))
      pval.chisq.overall <- pchisq(2*(ll3-ll1),df=paradiff31,lower.tail = F)
      fdr.chisq.overall <- p.adjust(pval.chisq.overall, method='fdr')
      pval.chisq.trendDiff <- pchisq(2*(ll3-ll2),df=paradiff32,lower.tail = F)
      fdr.chisq.trendDiff <- p.adjust(pval.chisq.trendDiff, method='fdr')
      pval.chisq.meanDiff <- pchisq(2*(ll2-ll1),df=paradiff21,lower.tail = F)
      fdr.chisq.meanDiff <- p.adjust(pval.chisq.meanDiff, method='fdr')
      res <- data.frame(fdr.chisq.overall = fdr.chisq.overall, 
                        pval.chisq.overall = pval.chisq.overall,
                        df.diff.overall = paradiff31,
                        fdr.chisq.trendDiff = fdr.chisq.trendDiff, 
                        pval.chisq.trendDiff = pval.chisq.trendDiff, 
                        df.diff.trendDiff = paradiff32,
                        fdr.chisq.meanDiff = fdr.chisq.meanDiff, 
                        pval.chisq.meanDiff = pval.chisq.meanDiff,
                        df.diff.meanDiff = paradiff21,
                        stringsAsFactors = FALSE)
      reslist = list(statistics = res, ll1 = ll1, ll2 = ll2, ll3 = ll3, parameter = res3$parameter, knotnum = res3$knotnum)  ## function return
    }
    
  } else if (test.method == 'permutation'){
    if (verbose.output) print('fitting model: overall: CovariateTest (Model 3 vs.1) or ConstantTest (Model 1) ...')
    if (ncores == 1){
      fit <-
        lapply(seq_len(permuiter + 1), function(i)
          fitfunc(
            iter = i,
            diffType = 'overall',
            test.type = test.type,
            testvar = testvar,
            maxknotallowed = maxknotallowed,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose.output = verbose.output,
            expr = expr,
            cellanno = cellanno,
            pseudotime = pseudotime,
            design = design
          ))
    } else {
      fit <-
        mclapply(seq_len(permuiter + 1), function(i) {
          set.seed(i)
          fitfunc(
            iter = i,
            diffType = 'overall',
            test.type = test.type,
            testvar = testvar,
            maxknotallowed = maxknotallowed,
            EMmaxiter = EMmaxiter,
            EMitercutoff = EMitercutoff,
            verbose.output = verbose.output,
            expr = expr,
            cellanno = cellanno,
            pseudotime = pseudotime,
            design = design
          )
        }, mc.cores = ncores)
    }
    if (verbose.output) {
      print('The length of fit is ...')  
      print(sapply(fit, length))
      print(summary(sapply(fit,is.null))) 
    }
    fit <- fit[!sapply(fit,is.null)]
    if (verbose.output){
      print('The length of fit after removing null is ...')  
      print(sapply(fit, length))
      print(summary(sapply(fit,is.null))) 
    }
    
    if (length(fit[[1]]) > 1){
      fit <- fit[sapply(fit,length) > 1]
      if (verbose.output){
        print('The length of fit having both null and full model is ...')  
        print(sapply(fit, length))  
      }  
    }
    
    knotnum <- fit[[1]]$fitres.full$knotnum
    parameter <- fit[[1]]$fitres.full$parameter
    ll.full <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    ll.null <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
    llr.overall <- ll.full - ll.null
    pval.overall <- sapply(seq_len(nrow(llr.overall)), function(i) {
      z <- llr.overall[i, seq(2, ncol(llr.overall))]
      z <- z[!is.na(z)]
      den <- density(z)$bw
      mean(pnorm(llr.overall[i, 1], z, sd = den, lower.tail = F))
    })
    log.pval <- sapply(seq_len(nrow(llr.overall)), function(i) {
      z <- llr.overall[i, seq(2, ncol(llr.overall))]
      z <- z[!is.na(z)]
      den <- density(z)$bw
      max(pnorm(llr.overall[i, 1], z, sd = den, lower.tail = F, log.p = T))
    })
    fdr.overall <- p.adjust(pval.overall,method='fdr')
    names(pval.overall) <- names(fdr.overall) <- row.names(llr.overall)
    z.score <- (llr.overall[,1] - rowMeans(llr.overall[, seq(2, (ncol(llr.overall)))]))/apply(llr.overall[, seq(2, (ncol(llr.overall)))],1,sd)
    res.overall <- data.frame(fdr.overall = fdr.overall, pval.overall = pval.overall, z.overall = z.score,
                              log.pval.overall = log.pval, stringsAsFactors = FALSE)
    if (verbose.output) print(paste0('Number of overall XDE genes: ', sum(fdr.overall < 0.05) ))
    
    if (sum(fdr.overall <0.05 ) > 0 & toupper(test.type) == 'VARIABLE' & !overall.only){
      if (verbose.output) print('meanDiff pvalues: Model 2 vs. model 1...')
      if (ncores == 1){
        fit <-
          lapply(seq_len(permuiter + 1), function(i)
            fitfunc(
              iter = i,
              diffType = 'meanDiff',
              gene = names(fdr.overall)[fdr.overall < 0.05],
              test.type = test.type,
              maxknotallowed = maxknotallowed,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose.output = verbose.output,
              expr = expr,
              cellanno = cellanno,
              pseudotime = pseudotime,
              design = design,
              testvar = testvar
            ))
      } else {
        fit <-
          mclapply(seq_len(permuiter + 1), function(i) {
            set.seed(i)
            fitfunc(
              iter = i,
              diffType = 'meanDiff',
              gene = names(fdr.overall)[fdr.overall < 0.05],
              test.type = test.type,
              maxknotallowed = maxknotallowed,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose.output = verbose.output,
              expr = expr,
              cellanno = cellanno,
              pseudotime = pseudotime,
              design = design,
              testvar = testvar
            )
          }, mc.cores = ncores)
      }
      fit <- fit[!sapply(fit,is.null)]
      ll.full <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      ll.null <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      llr <- ll.full - ll.null
      llr <- llr[complete.cases(llr), ]
      if (sum(fdr.overall<0.05) == 1){
        z <- llr[seq(2, length(llr))]
        z <- z[!is.na(z)]
        den <- density(z)$bw
        fdr <- pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F))
        log.pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F,log.p=T))
        names(pval) <- names(fdr) <- names(log.pval) <- names(fdr.overall)[fdr.overall<0.05]
        z.score <- (llr[1] - mean(z))/sd(z)
      } else {
        pval <- sapply(seq_len(nrow(llr)), function(i) {
          z <- llr[i, seq(2, ncol(llr))]
          z <- z[!is.na(z)]
          den <- density(z)$bw
          mean(pnorm(llr[i,1], z, sd=den,lower.tail = F))
        })
        log.pval <- sapply(seq_len(nrow(llr)), function(i) {
          z <- llr[i, seq(2, ncol(llr))]
          z <- z[!is.na(z)]
          den <- density(z)$bw
          max(pnorm(llr[i,1], z, sd=den,lower.tail = F, log.p = T))
        })
        fdr <- p.adjust(pval,method='fdr')
        names(pval) <- names(fdr) <- names(log.pval) <- row.names(llr)
        z.score <- (llr[,1] - rowMeans(llr[, seq(2, (ncol(llr)))]))/apply(llr[, seq(2, (ncol(llr)))],1,sd)
      }
      res.meanDiff <- data.frame(fdr.meanDiff = fdr, pval.meanDiff = pval, z.meanDiff = z.score, log.pval.meanDiff = log.pval, stringsAsFactors = FALSE)
      if (verbose.output) print('trendDiff pvalues: Model 3 vs. model 2...')
      if (ncores == 1){
        fit <-
          lapply(seq_len(permuiter + 1), function(i)
            fitfunc(
              iter = i,
              diffType = 'trendDiff',
              gene = names(fdr.overall)[fdr.overall < 0.05],
              test.type = test.type,
              maxknotallowed = maxknotallowed,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose.output = verbose.output,
              expr = expr,
              cellanno = cellanno,
              pseudotime = pseudotime,
              design = design,
              testvar = testvar
            ))
      } else {
        fit <-
          mclapply(seq_len(permuiter + 1), function(i) {
            set.seed(i)
            fitfunc(
              iter = i,
              diffType = 'trendDiff',
              gene = names(fdr.overall)[fdr.overall < 0.05],
              test.type = test.type,
              maxknotallowed = maxknotallowed,
              EMmaxiter = EMmaxiter,
              EMitercutoff = EMitercutoff,
              verbose.output = verbose.output,
              expr = expr,
              cellanno = cellanno,
              pseudotime = pseudotime,
              design = design,
              testvar = testvar
            )
          }, mc.cores = ncores) ## return a list of (permuiter + 1) where the first is a list of fitres.full and fitres.null
      }
      fit <- fit[!sapply(fit,is.null)]
      ll.full <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.full$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      ll.null <- sapply(seq_len(length(fit)),function(i) sapply(fit[[i]]$fitres.null$parameter,function(j) unname(j$ll),USE.NAMES = F)[row.names(expr)])
      llr <- ll.full - ll.null
      llr <- llr[complete.cases(llr), ]
      if (sum(fdr.overall<0.05) == 1){
        z <- llr[seq(2, length(llr))]
        z <- z[!is.na(z)]
        den <- density(z)$bw
        fdr <- pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F))
        log.pval <- mean(pnorm(llr[1], z, sd=den,lower.tail = F,log.p=T))
        names(pval) <- names(fdr) <- names(log.pval) <- names(fdr.overall)[fdr.overall<0.05]
        z.score <- (llr[1] - mean(z))/sd(z)
      } else {
        pval <- sapply(seq_len(nrow(llr)), function(i) {
          z <- llr[i, seq(2, ncol(llr))]
          z <- z[!is.na(z)]
          den <- density(z)$bw
          mean(pnorm(llr[i,1], z, sd=den,lower.tail = F))
        })
        log.pval <- sapply(seq_len(nrow(llr)), function(i) {
          z <- llr[i, seq(2, ncol(llr))]
          z <- z[!is.na(z)]
          den <- density(z)$bw
          max(pnorm(llr[i,1], z, sd=den,lower.tail = F, log.p = T))
        })
        fdr <- p.adjust(pval,method='fdr')
        names(pval) <- names(fdr) <- names(log.pval) <- row.names(llr)
        z.score <- (llr[,1] - rowMeans(llr[, seq(2, (ncol(llr)))]))/apply(llr[, seq(2, (ncol(llr)))],1,sd)
      }
      
      res.trendDiff <- data.frame(fdr.trendDiff = fdr,  pval.trendDiff = pval, z.trendDiff = z.score, log.pval.trendDiff = log.pval,  stringsAsFactors = FALSE)
      res <- matrix(NA, nrow = nrow(res.overall), ncol = (ncol(res.overall) + ncol(res.meanDiff) + ncol(res.trendDiff)), 
                    dimnames = list(rownames(res.overall),c(colnames(res.overall), colnames(res.meanDiff), colnames(res.trendDiff))))
      res[rownames(res.overall), colnames(res.overall)] <- as.matrix(res.overall)
      res[rownames(res.trendDiff), colnames(res.trendDiff)] <- as.matrix(res.trendDiff)
      res[rownames(res.meanDiff), colnames(res.meanDiff)] <- as.matrix(res.meanDiff)
    } else if (sum(fdr.overall < 0.05) == 0 | (toupper(test.type) == 'VARIABLE' & overall.only) | toupper(test.type) == 'TIME'){
      if (verbose.output) print('Not returning meanDiff and trendDiff: constantTest, user required or no overall XDE genes')
      res <- res.overall
    }
    reslist <- list(statistics = res, 
                    parameter=parameter, 
                    llr.overall = llr.overall,
                    knotnum = knotnum)           ## function return
    
  }
  if (return.all.data){
    return(c(reslist, list(pseudotime = pseudotime[colnames(expr)], design = design, cellanno = cellanno, expr = expr, test.type = test.type, test.method = test.method, testvar = testvar)))
  } else {
    return(c(reslist, list(test.type = test.type, test.method = test.method, testvar = testvar)))
  } 
}
