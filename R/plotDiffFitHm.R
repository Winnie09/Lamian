#' Plot the fitting heatmaps for the sample covariate test (with group difference).
#' 
#' This function is used for plotting the fitting heatmaps for the sample covariate test (with group difference).
#' 
#' @import ggplot2 RColorBrewer gridExtra viridis pheatmap
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param testobj output object from lamian.test(). It is a list.
#' @param showRowName logical. If FALSE (default), row names (i.e. gene names) of the heatmaps will not be shown in the plot.
#' @param cellWidthTotal a numeric number. Total width of each heatmap cell. 
#' @param cellHeightTotal when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
#' @param showCluster (not implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster". 
#' @param colann a data frame. Each column represent the annotation feature for the cells. row names should be the same as the name of the cells in the data.
#' @param rowann a data frame. Each column represent the annotation feature for the genes. row names should be the same as the name of the genes in the data.
#' @param annotation_colors a list. Will be passed onto the annotation_colors input in pheatmap().
#' @param type One of c('Time', 'Variable')
#' @param subsampleCell logical. If TRUE, will subsample cells.
#' @param numSubsampleCell a numeric number indicating the number of cells users want to subsampled. Only useful when subsampleCell == TRUE.
#' @param sep a string of the gene names that need to removed.
#' @param break.0 logical. If TRUE (default), the heatmap color scale will strengthen the difference around 0
#' @examples
#' data(mantestobj)
#' plotDiffFitHm(testobj = mantestobj, type = 'variable')

plotDiffFitHm <- function(testobj, showRowName = FALSE, cellWidthTotal = 250, cellHeightTotal = 400, showCluster = FALSE, colann = NULL, rowann = NULL, annotation_colors = NULL, type = 'time', subsampleCell = TRUE, numSubsampleCell=1e3, sep = NA){
  fit <- testobj$populationFit
  
  if ('DDGType' %in% names(testobj)) {
    DDGType <- testobj$DDGType
  } else {
    DDGType <- getDDGType(testobj)
  }
  DDGType <- DDGType[rownames(testobj$covariateGroupDiff)]
  
  if (subsampleCell){
    if (toupper(type) == 'TIME'){
      id <- round(seq(1, ncol(fit), length.out = numSubsampleCell))
      fit <- fit[, id]
    } else if (toupper(type) == 'VARIABLE'){
      id <- round(seq(1, ncol(fit[[1]]), length.out = numSubsampleCell))
      for (i in 1:length(fit)){
        fit[[i]] <- fit[[i]][, id]
      }
      if (sum(DDGType == 'meanSig') > 0){
        meanid <- which(DDGType == 'meanSig')
        
        FitDiff.scale1 <- scalematrix(testobj$covariateGroupDiff[-meanid,id,drop=F]) ## add FitDiff.scale
        
        
        FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid,id,drop=F])
        FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff),, drop=F]
      } else {
        
        FitDiff.scale <- scalematrix(testobj$covariateGroupDiff[,id,drop=F]) ## add FitDiff.scale  
        
        
      }
      colnames(FitDiff.scale) <- paste0('FitDiff:cell', seq(1, ncol(FitDiff.scale)))
      testobj$pseudotime <- sort(sample(testobj$pseudotime, numSubsampleCell))
    }
    print('subsample done!')
  } else {
    if (toupper(type) == 'VARIABLE'){
      if (sum(DDGType == 'meanSig') > 0){
        meanid <- which(DDGType == 'meanSig')
        FitDiff.scale1 <- scalematrix(testobj$covariateGroupDiff[-meanid,,drop=F]) ## add FitDiff.scale
        
        FitDiff.scale <- rbind(FitDiff.scale1, testobj$covariateGroupDiff[meanid,,drop=F])
        FitDiff.scale <- FitDiff.scale[rownames(testobj$covariateGroupDiff),, drop=F]
      } else {
        FitDiff.scale <- scalematrix(testobj$covariateGroupDiff)
      }
      colnames(FitDiff.scale) <- paste0('FitDiff:cell', seq(1, ncol(FitDiff.scale)))
    }
  }
  
  fit.bak = fit
  clu <- testobj$cluster
  if (type == 'variable'){
    rownames(testobj$cellanno) <- testobj$cellanno[,1]
    testobj$cellanno <- testobj$cellanno[names(testobj$pseudotime), ]
    if ('expr.ori' %in% names(testobj)){
      testobj$expr <- testobj$expr.ori[, names(testobj$pseudotime)]
    } else {
      testobj$expr <- testobj$expr[, names(testobj$pseudotime)]
    }
    
    if ('DDGType' %in% names(testobj)) {
      DDGType <- testobj$DDGType
    } else {
      DDGType <- getDDGType(testobj) 
    }
    
    fit.scale <- do.call(cbind, fit)
    fit.scale <- fit.scale[names(testobj$cluster), ]
    fit.scale <- scalematrix(fit.scale)
    colnames(fit.scale) <- paste0(rep(names(fit), each = ncol(fit.scale)/length(fit)), ';cell', seq(1, ncol(fit.scale)))
  } else {
    fit.scale <- scalematrix(fit)
    dimnames(fit.scale) <- dimnames(fit)
  }
  # res <- data.frame(clu = clu, 
  #                   cor = sapply(names(clu), function(i) cor(fit.scale[i, seq(1, ncol(fit.scale)/2)], seq(1, ncol(fit.scale)/2))),
  #                   changepoint = sapply(names(clu), function(i) which.min(abs(fit.scale[i, seq(1, ncol(fit.scale)/2)]))),
  #                   DDGType = DDGType[names(clu)])
  res <- data.frame(clu = clu, 
                    cor = sapply(names(clu), function(i) cor(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))], seq(1, ncol(FitDiff.scale)))),
                    changepoint = sapply(names(clu), function(i) which.min(abs(FitDiff.scale[i, seq(1, ncol(FitDiff.scale))]))),
                    DDGType = DDGType[names(clu)])
  
  res <- res[order(res$clu, res$DDGType, res$changepoint, res$cor), ]
  fit.scale <- fit.scale[rownames(res), ]
  FitDiff.scale <- FitDiff.scale[rownames(res), ]
  # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
  ## ------------------------
  ## plot original expression 
  ## ------------------------
  cellanno <- testobj$cellanno
  expr = testobj$expr
  expr <- expr[, names(testobj$pseudotime)]
  
  if (type == 'variable'){
    tmp <- lapply(names(fit), function(i){
      expr[rownames(fit.scale), colnames(expr) %in% cellanno[cellanno[, 2] %in% rownames(testobj$design)[testobj$design[, 2] == sub('.*_','', i)], 1]]
    })
    expr.scale <- do.call(cbind, tmp)
  } else if (type == 'time'){
    expr.scale <- expr
  }
  expr.scale <- scalematrix(expr.scale)
  expr.scale <- expr.scale[rownames(fit.scale), ]  

  ## plot ------------------------
  expr.scale[expr.scale > quantile(as.vector(expr.scale), 0.98)] <-
    quantile(as.vector(expr.scale), 0.98)
  expr.scale[expr.scale < quantile(as.vector(expr.scale), 0.02)] <-
    quantile(as.vector(expr.scale), 0.02)
  FitDiff.scale[FitDiff.scale > quantile(as.vector(FitDiff.scale), 0.98)] <-
    quantile(as.vector(FitDiff.scale), 0.98)
  FitDiff.scale[FitDiff.scale < quantile(as.vector(FitDiff.scale), 0.02)] <-
    quantile(as.vector(FitDiff.scale), 0.02)
  fit.scale[fit.scale > quantile(as.vector(fit.scale), 0.98)] <-
    quantile(as.vector(fit.scale), 0.98)
  fit.scale[fit.scale < quantile(as.vector(fit.scale), 0.02)] <-
    quantile(as.vector(fit.scale), 0.02)
  
  ### annotate rows and columns
  
  if (is.null(colann)){
    if (type == 'variable'){
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        group = as.character(testobj$design[cellanno[match(colnames(expr.scale), cellanno[,1]),2],2]),
        expression = 'Original',
        stringsAsFactors = F)
      
      col.group = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$group))+1)
      names(col.group) = c('NA', unique(colann$group))
    } else if (type == 'time'){
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        expression = 'Original',
        stringsAsFactors = F)
    }
  }
  rownames(colann) = colnames(expr.scale)
  col.expression = brewer.pal(n = 8, name = "Pastel2")[1:3]
  names(col.expression) = c('Original', 'ModelFitted', 'ModeledGroupDiff')
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
  names(col.pseudotime) = unique(colann$pseudotime)
  
  if (is.null(rowann)){
    if (type == 'variable'){
      rowann = data.frame(
        cluster = as.character(clu),
        DDGType = as.character(DDGType[names(clu)]),
        stringsAsFactors = F)
    } else if (type == 'time'){
      rowann = data.frame(
        cluster = as.character(clu),
        stringsAsFactors = F)
    }
    rownames(rowann) = names(clu)
  }
  rowann <- rowann[rownames(fit.scale), ,drop=F]
  rowann[,'DDGType'] <- factor(as.character(rowann[,'DDGType']), levels = c('trendSig','meanSig','bothSig','nonDDG','other'))
  
  if (length(unique(clu)) < 8){
    col.clu = brewer.pal(8, 'Set1')[1:length(unique(clu))]
  } else {
    col.clu = colorRampPalette(brewer.pal(8, 'Set1'))[1:length(unique(clu))]
  }
  names(col.clu) = unique(clu)
  
  if (is.null(colann)| is.null(annotation_colors)){
    if (type == 'variable'){
      col.DDGType = brewer.pal(8, 'Set3')[1:5]
      names(col.DDGType) = c('trendSig','meanSig','bothSig','nonDDG','other')
      annotation_colors = list(
        pseudotime = col.pseudotime,
        group = col.group,
        expression = col.expression,
        cluster = col.clu,
        DDGType = col.DDGType)
    } else if (type == 'time'){
      annotation_colors = list(
        pseudotime = col.pseudotime,
        expression = col.expression,
        cluster = col.clu,
        gs = col.gs,
        limmaPb = col.limmaPb)
    }
  }
  col.gs <- c('pink', 'skyblue')
  names(col.gs) <- c('No', 'Yes')
  col.limmaPb <- c('pink', 'skyblue')
  names(col.limmaPb) <- c('nonDiff', 'Diff')
  annotation_colors[['gs']] <- col.gs
  annotation_colors[['limmaPb']] <- col.limmaPb
  
  col.signalType <- brewer.pal(8, 'Set3')[1:3]
  names(col.signalType) <- c('trend only', 'mean only', 'both')
  annotation_colors[['signalType']] <- col.signalType
  
  #### save png
  cpl = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  plist <- list()
  
  if (!is.na(sep)){
    rownames(expr.scale) <-sub(sep, '', rownames(expr.scale))
    rownames(rowann) <- sub(sep, ':.*', rownames(rowann))
  }
  
  p1 <- pheatmap(
    expr.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal / ncol(expr.scale),
    cellheight = cellHeightTotal / nrow(expr.scale),
    border_color = NA, silent = TRUE)
  plist[[1]] <- p1[[4]] 
  
  ## --------------------
  ## plot fitting values
  ## --------------------
  if (type == 'variable'){
    colann.fit1 <-data.frame(pseudotime = rep(1:ncol(fit[[1]]), length(fit)),
                             group = gsub(sub('_.*', '_', names(fit)[1]),'',sub(';.*', '', colnames(fit.scale))), 
                             expression = 'ModelFitted',
                             stringsAsFactors = F)
    colann.fit2 <-data.frame(pseudotime = seq(1, ncol(FitDiff.scale)),
                             group = 'NA', 
                             expression = 'ModeledGroupDiff',
                             stringsAsFactors = F)
    colann.fit <- rbind(colann.fit1, colann.fit2)
    
  } else if (type == 'time'){
    colann.fit <-data.frame(pseudotime = testobj$pseudotime[colnames(fit.scale)],
                            expression = 'ModelFitted',
                            stringsAsFactors = F)
  }
  col.pseudotime = colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
  names(col.pseudotime) = unique(colann.fit$pseudotime)
  annotation_colors$pseudotime <- col.pseudotime
  # col.group <- colorRampPalette(brewer.pal(n=8,'Accent'))(length(unique(colann.fit$group)))
  # names(col.group) = unique(colann.fit$group)
  # annotation_colors$group <-   col.group
  
  fit.scale <- cbind(fit.scale, FitDiff.scale)  ## cbind FitDiff !!!
  rownames(colann.fit) = colnames(fit.scale)
  
  if (!is.na(sep)){
    rownames(fit.scale) <-sub(sep, '', rownames(fit.scale))  
  }
  
  p2 <- pheatmap(
    fit.scale,
    cluster_rows = F,
    cluster_cols = FALSE,
    show_rownames = showRowName,
    show_colnames = FALSE,
    color = cpl,
    annotation_col = colann.fit,
    annotation_row = rowann,
    annotation_colors = annotation_colors,
    cellwidth = cellWidthTotal*1.23 / ncol(fit.scale),
    cellheight = cellHeightTotal / nrow(fit.scale),
    border_color = NA, silent = TRUE)
  plist[[3]] <- p2[[4]] 
  plist[[2]] <- ggplot(data=NULL) + geom_blank() + theme_void()
  grid.arrange(grobs = plist,layout_matrix=matrix(c(1,1,1,1,2,3,3,3,3),nrow=1))
}  




