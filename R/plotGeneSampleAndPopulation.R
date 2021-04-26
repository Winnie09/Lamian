#' Plot genes' sample-level and population-level fitting plot for sample covariate test.
#' 
#' This function is used to plot genes' sample-level and population-level fitting plot for sample covariate test.
#' 
#' @param testObj the output object from lamian.test().
#' @param gene a vector of genes that need to do the prediction.
#' @param variable character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
#' @param variable.text a character vector. The text for the legend of the plot, corresponding to each variable values.
#' @param continuous: if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
#' @param expression a character ('demean',or 'original') to define the expression values shown on the plots. if "demean"(default), show demeaned expression. if 'original", show original gene expression.
#' @param ncol only functional when plotting multiple genes. Used to define the number of columns. 
#' @param free.scale logical. If TRUE, the y-axis is on free scale.
#' @param facet.sample logical. If TRUE, each sample is plotted separately.
#' @param plot.point logical. If TRUE, cells will be plotted as ponts.
#' @param line.alpha alpha value for lines.
#' @param line.size size value for lines.
#' @param point.alpha alpha value for points.
#' @param point.size size value for points.
#' @param sep a string in the gene names that needs to replaced with blank.
#' @param palette a RColorBrewer palette name
#' @param ncol number of columns
#' @param axis.text.blank logical. If TRUE, leave axis text as blank.
#' @import splines ggplot2 gridExtra viridis RColorBrewer
#' @export
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @examples
#' data(mantestobj)
#' plotGeneSampleAndPopulation(testobj = mantestobj, gene = rownames(mantestobj$populationFit[[1]])[1], variable = 'gender')

plotGeneSampleAndPopulation <- function(testobj, gene, variable = NULL, variable.text = NULL, free.scale = TRUE, facet.sample = FALSE, plot.point = FALSE, line.alpha = 1, line.size = 1, point.alpha=1, point.size=0.5, continuous = TRUE, sep = NA, palette = 'Dark2', ncol = NULL,  axis.text.blank = F){
  pseudotime <- testobj[['pseudotime']]
  cellanno <- testobj[['cellanno']]
  colnames(cellanno) <- c('Cell', 'Sample')
  if ('expr.ori' %in% names(testobj)) expression <- testobj[['expr.ori']] else expression <- testobj[['expr']]
  if ('populationFit' %in% names(testobj)) fit <- testobj$populationFit else fit = getPopulationFit(testobj, gene = gene, type = type)
  predict.values <- predict_fitting(testobj, gene= gene, test.type = testobj$test.type)
  
  pseudotime = pseudotime[colnames(expression)]
  cellanno <- cellanno[match(colnames(expression), cellanno[,1]), ]
  # predict.values <- predict.values[, colnames(expression),drop=F]
  knotnum <- testobj$knotnum
  knotnum[knotnum==0] <- 1  ## in case the fitting of line would cause bugs
  design <- testobj[['design']]
  
  
  cellanno <- data.frame(Cell = as.character(cellanno[,1]), 
                         Sample = as.character(cellanno[,2]), stringsAsFactors = FALSE)
  variable.d <- if(is.null(variable)) 1 else variable
  if (!is.null(variable.text) & variable.d != 1) {
    design[,variable.d] <- ifelse(design[, variable.d] == 0, variable.text[1], variable.text[2])
  }
  a <- if (free.scale) 'free' else 'fixed'
  if (length(gene) == 1){
    print('plotting one gene ...')
    pd <- data.frame(expr = expression[gene, ], 
                     Sample = cellanno[,2], 
                     Variable = design[match(cellanno[,2], rownames(design)), variable.d],  ##
                     pseudotime = pseudotime[colnames(expression)])
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    linedlist <- lapply(unique(cellanno[,2]), function(p){
      # tmpcell <- cellanno[cellanno[,2]==p,1]
      tmpcellid <- which(cellanno[,2]==p)
      if (toupper(testobj$test.type) == 'TIME'){  ##### add
        tmpdf <- data.frame(expr=predict.values[gene,tmpcellid], 
                            Sample=p, 
                            Variable = design[rownames(design) == p, variable.d], ##
                            pseudotime=pseudotime[tmpcellid])
      } else {
        tmpdf <- data.frame(
          expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][gene, tmpcellid],
          Sample = p,
          Variable = design[rownames(design) == p, variable.d], ##
          pseudotime=pseudotime[tmpcellid])
      }
    })
    ld <- do.call(rbind, linedlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] ## add 20200812
    ### add population Line >>>>>>>>>>>>>>
    pd <- sapply(1:length(fit), function(i){
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      colnames(tmp) <- c('g', 'pseudotime', 'expression')
      tmp <- data.frame(tmp, 
                        type = names(fit)[i],
                        stringsAsFactors = FALSE)
    }, simplify = FALSE)
    pd <- do.call(rbind, pd)
    if (!is.na(sep)) {
      pd$g <- sub(sep, '', pd$g)
      pd$g <- factor(as.character(pd$g), levels = sub(sep, '', gene))
    } else{
      pd$g <- factor(as.character(pd$g), levels = gene)
    }
    pd$type <- sub('.*_', '', pd$type)
    pd2 <- pd
    ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    if (is.null(variable)){
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size) 
        
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample),alpha=line.alpha, size=line.size)   
      }   
    }
    p <- p + geom_line(data = pd2, aes(x = pseudotime, y = expression, group = type, color = type), size = 2, alpha = 1)  ## add population line
    p <- p + 
      theme_classic() +
      # ggtitle(paste0(sub(':.*','',gene),',adj.pvalue=', formatC(testobj$fdr[gene], format = "e", digits = 2))) +
      xlab('Pseudotime') + ylab('Expression') + 
      labs(color = variable) +
      theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm")) +
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
    if (!is.na(sep)){
      p <- p + ggtitle(sub(sep,'',gene)) 
    } else {
      p <- p + ggtitle(gene)
    }
    if (continuous){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1) 
    } else {
      if (length(unique(ld[, 'Variable'])) > 8){
        p <- p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))  
      } else {
        p <- p + scale_color_manual(values = brewer.pal(8, palette)[1:length(unique(ld[, 'Variable']))])
      } 
    }
    if (facet.sample){
      print(p + facet_wrap(~Sample, scales=a))
    } else {
      print(p)
    }
  } else {
    print('plotting multiple genes ...')
    pdlist <- ldlist <- list()
    for (g in gene){
      pd <- data.frame(expr = expression[g, ], 
                       pseudotime = pseudotime[colnames(expression)],
                       Sample = cellanno[match(colnames(expression), cellanno[,1]),2], 
                       Variable = design[match(cellanno[,2], rownames(design)), variable.d],
                       g = g, stringsAsFactors = F)
      pdlist[[g]] <- pd
      linedlist <- lapply(unique(testobj$cellanno[,2]), function(p){
        tmpcellid <- which(cellanno[,2]==p) 
        if (toupper(testobj$test.type) == 'TIME'){  ##### add
          tmpdf <- data.frame(expr=predict.values[g,tmpcellid], 
                              pseudotime=pseudotime[tmpcellid],
                              Sample=p, 
                              Variable = design[rownames(design) == p, variable.d], 
                              g = g, stringsAsFactors = F)
        } else {
          tmpdf <- data.frame(
            expr = predict.values[[which(as.numeric(sub('.*_', '', names(predict.values))) == design[p, variable.d])]][g, tmpcellid],
            pseudotime=pseudotime[tmpcellid],
            Sample = p,
            Variable = design[rownames(design) == p, variable.d], 
            g = g, stringsAsFactors = F)
        }
      })
      ld = do.call(rbind, linedlist)
      ldlist[[g]] <- ld
    }
    pd <- do.call(rbind, pdlist)
    pd[, 'Variable'] <- as.factor(pd[ ,'Variable'])
    ld <- do.call(rbind, ldlist)
    ld[, 'Variable'] <- as.factor(ld[ ,'Variable'])
    ld <- ld[order(ld$pseudotime), ] 
    if (!is.na(sep)) {
      pd$g <- gsub(sep, '', pd$g)
      ld$g <- gsub(sep, '', ld$g)
      pd$g <- factor(pd$g, levels = gsub(sep, '', gene))
      ld$g <- factor(ld$g, levels = gsub(sep, '', gene))
    } else {
      pd$g <- factor(pd$g, levels = gene)
      ld$g <- factor(ld$g, levels = gene)
    }
    pd1 <- pd
    ### add population Line >>>>>>>>>>>>>>
    pd <- sapply(1:length(fit), function(i){
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      colnames(tmp) <- c('g', 'pseudotime', 'expression')
      tmp <- data.frame(tmp, 
                        type = names(fit)[i],
                        stringsAsFactors = FALSE)
    }, simplify = FALSE)
    pd <- do.call(rbind, pd)
    if (!is.na(sep)) {
      pd$g <- sub(sep, '', pd$g)
      pd$g <- factor(as.character(pd$g), levels = sub(sep, '', gene))
    } else{
      pd$g <- factor(as.character(pd$g), levels = gene)
    }
    pd$type <- sub('.*_', '', pd$type)
    pd2 <- pd
    ### <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    if (is.null(variable)){
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Sample), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)
        
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Sample), alpha=line.alpha, size=line.size)   
      }
      
    } else {
      if (plot.point){
        p <- ggplot() + 
          geom_point(data=pd, aes(x=pseudotime, y=expr, color=Variable), alpha=point.alpha, size=point.size) +
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      } else {
        p <- ggplot() + 
          geom_line(data=ld, aes(x=pseudotime, y=expr, color=Variable, group = Sample), alpha=line.alpha, size=line.size) 
      }
    }
    
    p <- p + geom_line(data = pd2, aes(x = pseudotime, y = expression, group = type, color = type), size = 2, alpha = 1)
    
    p <- p + 
      theme_classic() + 
      xlab('Pseudotime') + ylab('Expression') + 
      labs(color = variable) +
      facet_wrap(~g, scales = a, ncol = ncol) +
      theme(legend.spacing.y = unit(0.01, 'cm'), legend.spacing.x = unit(0.01, 'cm'), legend.key.size = unit(0.1, "cm"))+
      guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))
    

    
    if (continuous){
      p <- p + scale_color_viridis(discrete = TRUE, direction = -1) 
    } else {
      if (length(unique(ld[, 'Variable'])) > 8){
        p <- p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))  
      } else {
        p <- p + scale_color_manual(values = brewer.pal(8, palette)[1:length(unique(ld[, 'Variable']))])
      } 
    }
    if (axis.text.blank) {
      p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
    } else {
      p <- p + scale_x_continuous(breaks=c(min(pd$pseudotime),max(pd$pseudotime)))
    }
    print(p)
  }
  
}



