#' Plot the population level fitting curves for genes
#' 
#' This function takes a test output object as input and then plot the population level fitting curves for genes that users specify.
#' 
#' @param testobj the output object from lamian.test().
#' @param gene a vector of genes that need to do the prediction.
#' @param type One of c('Time', 'Variable'). 
#' @param ylim y-axis limits to be passed to ggplot.
#' @param xlim x-axis limits to be passed to ggplot.
#' @param sep a string in the gene names that needs to replaced with blank.
#' @param free.scale logical. If TRUE, the y-axis is on free scale.
#' @param palette a RColorBrewer palette name.
#' @param ncol number of colums for organizing all genes' plot.
#' @param line.size the size of the curves.
#' @param axis.text.blank logical. If TRUE, leave axis text as blank.
#' @export
#' @import ggplot2 RColorBrewer
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @examples
#' data(mantestobj)
#' plotGenePopulation(testobj = mantestobj, type = 'variable', gene = rownames(mantestobj$populationFit[[1]])[seq(1,2)])

plotGenePopulation <- function(testobj,
                               gene = NA,
                               type = 'time',
                               ylim = NA,
                               xlim = NA,
                               sep = NA,
                               free.scale = TRUE,
                               palette = 'Dark2', 
                               ncol = NA,
                               subSampleNumber = NA,
                               line.size = 1, 
                               axis.text.blank = F){
  ## testobj: object returned from lamian.test(). 
  ## gene: a character vector of gene names. It can be of length 1 or > 1.
  ## variable: a character (within the column names in design matrix) to get population pattern. If variable == NA, then return testtime fit. Other wise return testvar fit with the variable.
  
  if (is.na(ncol))  nrow = round(sqrt(length(gene))) else nrow = NA 
  a <- if (free.scale) 'free' else 'fixed'
  if ('populationFit' %in% names(testobj)) fit <- testobj$populationFit else fit = getPopulationFit(testobj, gene = gene, type = type)
  
  if (type == 'time'){
    if (is.na(gene)) gene <- rownames(fit)    
    pd <- sapply(gene, function(g){
      if (!is.na(sep)) g2 <- sub(sep, '', g) else g2 = g  
      tmp <- data.frame(gene = g2, expression = fit[g, ], pseudotime = testobj$pseudotime, stringsAsFactors = FALSE)
    }, simplify = F)
    pd <- do.call(rbind, pd)
    
    if (!is.na(sep)) {
      pd$gene <- factor(as.character(pd$gene), levels = sub(sep, '', gene))
    } else {
      pd$gene <- factor(as.character(pd$gene), levels = gene)
    }
    
    
    p <- ggplot2::ggplot(data= pd, aes(x = pseudotime, y = expression, color = 'red')) + 
      geom_line(size = line.size) +
      theme_classic() +
      xlab('Pseudotime') +
      ylab('Expression') +
      labs(color = '')
    if (is.na(ncol)){
        p <- p + facet_wrap(~gene, nrow = nrow, scales = a)
    } else {
        p <- p + facet_wrap(~gene, ncol = ncol, scales = a)
    }
  } else {
    if (is.na(gene)) gene <- rownames(fit[[1]])
    pd <- sapply(1:length(fit), function(i){
      
      tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      if (!is.na(subSampleNumber)) {
        set.seed(12345)
        id <- sample(1:nrow(fit[[i]]), subSampleNumber)
        tmp <- reshape2::melt(fit[[i]][gene, id ,drop=F])
      } else {
        tmp <- reshape2::melt(fit[[i]][gene, ,drop=F])
      }
      
      colnames(tmp) <- c('gene', 'pseudotime', 'expression')
      if (!is.na(subSampleNumber)) {
        set.seed(12345)
        tmp <- tmp[sample(1:nrow(tmp), subSampleNumber), , drop=F]
      }
      tmp <- data.frame(tmp, 
                        type = names(fit)[i],
                        stringsAsFactors = FALSE)
    }, simplify = FALSE)
    pd <- do.call(rbind, pd)
    if (!is.na(sep)) {
      pd$gene <- sub(sep, '', pd$gene)
      pd$gene <- factor(as.character(pd$gene), levels = sub(sep, '', gene))
    } else{
      
      pd$gene <- factor(as.character(pd$gene), levels = gene)
    
    }
    
    
    p <- ggplot(data= pd, aes(x = pseudotime, y = expression, group = type, color = type)) + 
      geom_line(size = line.size) +
      theme_classic() +
      xlab('Pseudotime') +
      ylab('Expression') +
      labs(color = '') 
    if (is.na(ncol)){
      p <- p + facet_wrap(~gene, nrow = nrow, scales = a)
    } else {
      p <- p + facet_wrap(~gene, ncol = ncol, scales = a)
    }
    if (length(unique(pd$type)) < 8)  {
      p <-  p + scale_color_brewer(palette = palette)
    } else {
      p <- p + scale_color_manual(values = colorRampPalette(brewer.pal(8, palette))(length(unique(pd$type))))
    }
  }
  
  if (!is.na(ylim)[1]) p <- p + ylim(ylim)
  if (!is.na(xlim)[1]) p <- p + xlim(xlim)
  if (axis.text.blank) {
    p <- p + theme(axis.text = element_blank(), axis.ticks = element_blank())
  } else {
    p <- p + scale_x_continuous(breaks=c(min(pd$pseudotime),max(pd$pseudotime)))
  }
  print(p)
}

