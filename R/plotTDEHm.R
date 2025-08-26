#' Plot the fitting heatmaps for TDE test.
#'
#' This function is used for plotting the fitting heatmaps.
#'
#' @import ggplot2 RColorBrewer gridExtra viridis pheatmap
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param testobj output object from lamian_test(). It is a list.
#' @param showRowName logical. If FALSE (default), row names (i.e. gene names) of the heatmaps will not be shown in the plot.
#' @param cellWidthTotal a numeric number. Total width of each heatmap cell.
#' @param cellHeightTotal when showRowName = TRUE, cellHeightTotal is suggested to be ten times the number of genes (rows).
#' @param showCluster (no implemented yet). if TRUE, "cluster" should be a slot in testobj, and it will be label in the heatmap. If FALSE, no need to pass in "cluster".
#' @param colann a data frame. Each column represent the annotation feature for the cells. row names should be the same as the name of the cells in the data.
#' @param rowann a data frame. Each column represent the annotation feature for the genes. row names should be the same as the name of the genes in the data.
#' @param annotation_colors a list. Will be passed onto the annotation_colors input in pheatmap().
#' @param subsampleCell logical. If TRUE, will subsample cells.
#' @param numSubsampleCell a numeric number indicating the number of cells users want to subsampled. Only useful when subsampleCell == TRUE.
#' @examples
#' data(mantestobj)
           
plotTDEHm <-
  function(testobj,
           showRowName = FALSE,
           cellWidthTotal = 250,
           cellHeightTotal = 400,
           showCluster = FALSE,
           colann = NULL,
           rowann = NULL,
           annotation_colors = NULL,
           subsampleCell = TRUE,
           numSubsampleCell = 1e3) {
    
    ## extract the objects
    fit <- testobj$populationFit
    clu <- testobj$cluster
    
    # If cell were downsampled in populationFit. Matching the cells in pseudotime, and expr
    if (ncol(fit) < length(testobj$pseudotime)){
      testobj$pseudotime <- testobj$pseudotime[colnames(fit)]
      testobj$expr <- testobj$expr[rownames(fit), colnames(fit)]
      testobj$cellanno <- testobj$cellanno[match(colnames(fit), testobj$cellanno[,1]), ]
    }
    
    ## if users want to further subsampleCells through this plot function
    if (subsampleCell) {
      id <- round(seq(1, ncol(fit), length.out = numSubsampleCell))
      fit <- fit[, id]
      testobj$pseudotime <- testobj$pseudotime[colnames(fit)]
      testobj$expr <- testobj$expr[rownames(fit), colnames(fit)]
      testobj$cellanno <- testobj$cellanno[match(colnames(fit), testobj$cellanno[,1]), ]
      print('Subsample done!')
    }
    
    # ## extract the expression matrix, deprecated in new version (expr.ori was not used anymore)
    # if ('expr.ori' %in% names(testobj)) {
    #   expr <- testobj$expr.ori[, names(testobj$pseudotime)]
    # } else {
    #   expr <- testobj$expr[, names(testobj$pseudotime)]
    # }
    # 
    
    ## standardize the fitting
    fit.bak = fit
    fit.scale <- scalematrix(fit)
    dimnames(fit.scale) <- dimnames(fit)
    
    ## calculate the genes' correlation with pseudotime and changepoint within each cluster  -- will be used for sorting the genes in visualization
    res <- data.frame(
      clu = clu,
      cor = sapply(names(clu), function(i)
        cor(fit.scale[i, seq(1, ncol(fit.scale) / 2)], seq(
          1, ncol(fit.scale) / 2
        ))),
      changepoint = sapply(names(clu), function(i) {
        v <-
          fit.scale[i, seq(round(ncol(fit.scale) * 0.1), round(ncol(fit.scale) *
                                                                 0.9))]
        which(v[-length(v)] * v[-1] < 0)[1]
      })
    )
      # changepoint = sapply(names(clu), function(i) which.min(abs(fit.scale[i, seq(round(ncol(fit.scale)*0.01), round(ncol(fit.scale)*0.99))]))))
    
    ## This is the order of genes in visualization  
    res <- res[order(res$clu, res$changepoint, res$cor),]
    fit.scale <- fit.scale[rownames(res),]
    
    # colnames(fit.scale) <- paste0(colnames(fit.scale), '_', seq(1, ncol(fit.scale)))
    ## ------------------------
    ## plot original expression
    ## ------------------------
    cellanno <- testobj$cellanno
    expr.scale <- testobj$expr[rownames(res), names(testobj$pseudotime)] 
    
    ## <<<<<<<<<<<<<<<<<<<<<<<<<
    ## === sanity check ========
    ## if expr.scale has gene sd = 0, them remove these genes. It is possible that after downsampling the cells, some genes do not have original expression in the remaining cells.
    id = which(rowSds(expr.scale) == 0)
    if (length(id) > 0) {
      expr.scale <- expr.scale[-id, ,drop=FALSE]
      fit.scale <- fit.scale[rownames(expr.scale), ,drop=FALSE]
      clu = clu[names(clu) %in% rownames(expr.scale)]
    }
    ## >>>>>>>>>>>>>>>>>>>>>>>>>>
    
    ## standadize the original expression
    expr.scale <- scalematrix(expr.scale)
    # expr.scale <- expr.scale[rownames(fit.scale),]
    
    ## plot ------------------------
    expr.scale[expr.scale > quantile(as.vector(expr.scale), 0.98, na.rm = TRUE)] <-
      quantile(as.vector(expr.scale), 0.98) ##
    expr.scale[expr.scale < quantile(as.vector(expr.scale), 0.02)] <-
      quantile(as.vector(expr.scale), 0.02)
    fit.scale[fit.scale > quantile(as.vector(fit.scale), 0.98)] <-
      quantile(as.vector(fit.scale), 0.98)
    fit.scale[fit.scale < quantile(as.vector(fit.scale), 0.02)] <-
      quantile(as.vector(fit.scale), 0.02)
    
    ### annotate rows and columns
    if (is.null(colann)) {
      colann <- data.frame(
        # sample = cellanno[match(colnames(expr.scale),cellanno[, 1]), 2],
        pseudotime = testobj$pseudotime[colnames(expr.scale)],
        expression = 'Original',
        stringsAsFactors = FALSE
      )
      
    }
    rownames(colann) = colnames(expr.scale)
    col.expression = brewer.pal(n = 8, name = "Pastel1")[seq_len(2)]
    names(col.expression) = c('Original', 'Model Fitted')
    col.pseudotime = grDevices::colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann$pseudotime)))
    names(col.pseudotime) = unique(colann$pseudotime)
    
    if (is.null(rowann)) {
      rowann = data.frame(cluster = as.character(clu),
                            stringsAsFactors = TRUE)
      
      rownames(rowann) = names(clu)
      # rowann <- rowann[rownames(fit.scale), , drop = FALSE]
      
    }
    
    if (is.null(colann) | is.null(annotation_colors)) {
      if (length(unique(clu)) < 8) {
        col.clu = brewer.pal(8, 'Set1')[seq_len(length(unique(clu)))]
      } else {
        col.clu = grDevices::colorRampPalette(brewer.pal(8, 'Set1'))(length(unique(clu)))
      }
      names(col.clu) = sort(unique(clu))
      
      annotation_colors = list(pseudotime = col.pseudotime,
                               expression = col.expression,
                               cluster = col.clu)
      
    }
    
    #### save png
    cpl = grDevices::colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
    plist <- list()
    
    p1 <- pheatmap::pheatmap(
      expr.scale,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = showRowName,
      show_colnames = FALSE,
      color = cpl,
      annotation_col = colann,
      annotation_row = rowann,
      annotation_colors = annotation_colors,
      cellwidth = cellWidthTotal / ncol(expr.scale),
      cellheight = cellHeightTotal / nrow(expr.scale),
      border_color = NA,
      silent = TRUE
    )
    plist[[1]] <- p1[[4]]
    
    ## --------------------
    ## plot fitting values
    ## --------------------
   
    colann.fit <-
      data.frame(
        pseudotime = testobj$pseudotime[colnames(fit.scale)],
        expression = 'Model Fitted',
        stringsAsFactors = FALSE
      )
    col.pseudotime = grDevices::colorRampPalette(brewer.pal(n = 9, name = "YlGnBu"))(length(unique(colann.fit$pseudotime)))
    names(col.pseudotime) = unique(colann.fit$pseudotime)
    annotation_colors$pseudotime <- col.pseudotime
    rownames(colann.fit) = colnames(fit.scale)
    
    p2 <- pheatmap::pheatmap(
      fit.scale,
      cluster_rows = FALSE,
      cluster_cols = FALSE,
      show_rownames = showRowName,
      show_colnames = FALSE,
      color = cpl,
      annotation_col = colann.fit,
      annotation_row = rowann,
      annotation_colors = annotation_colors,
      cellwidth = cellWidthTotal / ncol(fit.scale),
      cellheight = cellHeightTotal / nrow(fit.scale),
      border_color = NA,
      silent = TRUE
    )
    plist[[3]] <- p2[[4]]
    plist[[2]] <- ggplot(data = NULL) + geom_blank() + theme_void()
    print(grid.arrange(grobs = plist, layout_matrix = matrix(c(
      1, 1, 1, 1, 2, 3, 3, 3, 3
    ), nrow = 1)))
    
  }

