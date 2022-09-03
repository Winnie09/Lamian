#' Plot gene(s) by showing both of the original cellcular gene expression and the sample-level, population-level fitting values.
#'
#' This function is used to plot gene(s) by showing both of the original cellcular gene expression and the sample-level, population-level fitting values.
#'
#' @import ggplot2 RColorBrewer splines gridExtra viridis
#' @importFrom grDevices colorRampPalette
#' @return a plot
#' @author Wenpin Hou <whou10@jhu.edu>
#' @export
#' @param testobj object returned from lamian_test().
#' @param gene a character vector of gene names. It can be of length 1 or > 1.
#' @param type One of c('Time', 'Variable').
#' @param variable character, the variable (covariate) to color the samples, should be null or one of the column names of design matrix. Default is NULL, meaning each sample is colored differently. Otherwise, samples are colored by the variable (covariate) values.
#' @param variable.text a character vector. The text for the legend of the plot, corresponding to each variable values.
#' @param facet.sample logical. If TRUE (default), facet_wrap the samples.
#' @param plot.point point size
#' @param line.alpha alpha value of the curves
#' @param continuous if TRUE, samples are colored using viridis continuous colors. If FALSE, RColorBrewer "Dark2" discrete palette.
#' @param cellProp logical. If FALSE (default), plot gene expression. If TRUE, it is cell proportion.
#' @param x.lab a string to indicates x-axis label
#' @param y.lab a string to indicates y-axis label
#' @param point.size the size value of points.
#' @param point.alpha the alpha value of points.
#' @param ylim y-axis limits to be passed to ggplot.
#' @param xlim x-axis limits to be passed to ggplot.
#' @param sep a string in the gene names that needs to replaced with blank.
#' @param free.scale logical. If TRUE, the y-axis is on free scale.
#' @param palette a RColorBrewer palette name.
#' @param ncol number of colums for organizing all genes' plot.
#' @param line.size the size of the curves.
#' @param axis.text.blank logical. If TRUE, leave axis text as blank.
#' @examples
#' data(mantestobj)
#' plotGene(testobj = mantestobj, gene = rownames(mantestobj$populationFit[[1]])[1], variable = 'gender')
plotGene <-
  function(testobj,
           gene,
           variable = NULL,
           variable.text = NULL,
           free.scale = TRUE,
           facet.sample = FALSE,
           plot.point = FALSE,
           line.alpha = 1,
           line.size = 1,
           point.alpha = 1,
           point.size = 0.5,
           continuous = TRUE,
           sep = NA,
           palette = 'Dark2',
           ncol = NULL,
           axis.text.blank = FALSE,
           cellProp = FALSE,
           x.lab = 'Pseudotime',
           y.lab = 'Expression') {
    pseudotime <- testobj[['pseudotime']]
    cellanno <- testobj[['cellanno']]
    colnames(cellanno) <- c('Cell', 'Sample')
    if ('expr.ori' %in% names(testobj)){
      expression <- testobj[['expr.ori']]
    } else{
      expression <- testobj[['expr']]
    }
      
    if (cellProp) {
      ptw <-
        cut(pseudotime,
            seq(min(pseudotime), max(pseudotime), length.out = 100),
            include.lowest = TRUE)
      ptdat <-
        table(ptw, cellanno[match(names(pseudotime), cellanno[, 1]), 2])
      ptdat <-
        t(t(ptdat) / colSums(ptdat)) ## divided by rowsum (rowsum = 1). interval * samples.
      ptdat <- as.data.frame(ptdat)
      colnames(ptdat) <- c('pt', 's', 'prop')
      ptdat[, 1] <- match(ptdat[, 1], levels(ptw))
      
      ptdat$cell <- paste0('cell', seq_len(nrow(ptdat)))
      ptexpr <- t(ptdat[, c('prop', 'prop'), drop = FALSE])
      colnames(ptexpr) <- ptdat$cell
      
      ptpt <- ptdat$pt
      names(ptpt) <- ptdat$cell
      expr = ptexpr
      cellanno = data.frame(cell = ptdat$cell, sample = ptdat$s)
      pseudotime = ptpt
    }
    
    predict.values <-
      predict_fitting(testobj, gene = gene, test.type = testobj$test.type)
    pseudotime = pseudotime[colnames(expression)]
    cellanno <- cellanno[match(colnames(expression), cellanno[, 1]),]
    # predict.values <- predict.values[, colnames(expression),drop=FALSE]
    knotnum <- testobj$knotnum
    knotnum[knotnum == 0] <-
      1  ## in case the fitting of line would cause bugs
    design <- testobj[['design']]
    
    
    cellanno <- data.frame(
      Cell = as.character(cellanno[, 1]),
      Sample = as.character(cellanno[, 2]),
      stringsAsFactors = FALSE
    )
    if (is.null(variable)){
     variable.d <- 1
    } else {
     variable.d <- variable
    }
      
    if (!is.null(variable.text) & variable.d != 1) {
      design[, variable.d] <-
        ifelse(design[, variable.d] == 0, variable.text[1], variable.text[2])
    }
    if (free.scale){
      a <- 'free'
    } else{
      a <- 'fixed'
    }
      
    if (length(gene) == 1) {
      print('plotting one gene ...')
      pd <- data.frame(
        expr = expression[gene,],
        Sample = cellanno[, 2],
        Variable = design[match(cellanno[, 2], rownames(design)), variable.d],
        ##
        pseudotime = pseudotime[colnames(expression)]
      )
      pd[, 'Variable'] <- as.factor(pd[, 'Variable'])
      linedlist <- lapply(unique(cellanno[, 2]), function(p) {
        # tmpcell <- cellanno[cellanno[,2]==p,1]
        tmpcellid <- which(cellanno[, 2] == p)
        if (toupper(testobj$test.type) == 'TIME') {
          ##### add
          tmpdf <- data.frame(
            expr = predict.values[gene, tmpcellid],
            Sample = p,
            Variable = design[rownames(design) == p, variable.d],
            ##
            pseudotime = pseudotime[tmpcellid]
          )
        } else {
          tmpdf <- data.frame(
            expr = predict.values[[which(as.numeric(sub('.*_', '', names(
              predict.values
            ))) == design[p, variable.d])]][gene, tmpcellid],
            Sample = p,
            Variable = design[rownames(design) == p, variable.d],
            ##
            pseudotime = pseudotime[tmpcellid]
          )
        }
      })
      ld <- do.call(rbind, linedlist)
      ld[, 'Variable'] <- as.factor(ld[, 'Variable'])
      ld <- ld[order(ld$pseudotime),] ## add 20200812
      if (is.null(variable)) {
        if (plot.point) {
          p <- ggplot() +
            geom_point(
              data = pd,
              aes(
                x = pd[,4],
                y = pd[,1],
                color = pd[,2]
              ),
              alpha = point.alpha,
              size = point.size
            ) +
            geom_line(
              data = ld,
              aes(
                x = ld[,4],
                y = ld[,1],
                color = ld[,3]
              ),
              alpha = line.alpha,
              size = line.size
            )
          
        } else {
          p <- ggplot() +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        }
      } else {
        if (plot.point) {
          p <- ggplot() +
            geom_point(
              data = pd,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable
              ),
              alpha = point.alpha,
              size = point.size
            ) +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable,
                group = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        } else {
          p <- ggplot() +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable,
                group = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        }
      }
      p <- p +
        theme_classic() +
        # ggtitle(paste0(sub(':.*','',gene),',adj.pvalue=', formatC(testobj$fdr[gene], format = "e", digits = 2))) +
        xlab(x.lab) + ylab(y.lab) +
        labs(color = variable) +
        theme(
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.01, 'cm'),
          legend.key.size = unit(0.1, "cm")
        ) +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
      if (!is.na(sep)) {
        p <- p + ggtitle(sub(sep, '', gene))
      } else {
        p <- p + ggtitle(gene)
      }
      if (continuous) {
        p <- p + scale_color_viridis(discrete = TRUE, direction = -1)
      } else {
        if (length(unique(ld[, 'Variable'])) > 8) {
          p <-
            p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))
        } else {
          p <-
            p + scale_color_manual(values = brewer.pal(8, palette)[seq_len(length(unique(ld[, 'Variable'])))])
        }
      }
      if (axis.text.blank) {
        p <-
          p + theme(axis.text = element_blank(), axis.ticks = element_blank())
      } else {
        p <-
          p + scale_x_continuous(breaks = c(min(pd$pseudotime), max(pd$pseudotime)))
      }
      if (facet.sample) {
        print(p + facet_wrap( ~ Sample, scales = a))
      } else {
        print(p)
      }
    } else {
      print('plotting multiple genes ...')
      pdlist <- ldlist <- list()
      for (g in gene) {
        pd <- data.frame(
          expr = expression[g,],
          pseudotime = pseudotime[colnames(expression)],
          Sample = cellanno[match(colnames(expression), cellanno[, 1]), 2],
          Variable = design[match(cellanno[, 2], rownames(design)), variable.d],
          g = g,
          stringsAsFactors = FALSE
        )
        pdlist[[g]] <- pd
        linedlist <- lapply(unique(testobj$cellanno[, 2]), function(p) {
          tmpcellid <- which(cellanno[, 2] == p)
          if (toupper(testobj$test.type) == 'TIME') {
            ##### add
            tmpdf <- data.frame(
              expr = predict.values[g, tmpcellid],
              pseudotime = pseudotime[tmpcellid],
              Sample = p,
              Variable = design[rownames(design) == p, variable.d],
              g = g,
              stringsAsFactors = FALSE
            )
          } else {
            tmpdf <- data.frame(
              expr = predict.values[[which(as.numeric(sub(
                '.*_', '', names(predict.values)
              )) == design[p, variable.d])]][g, tmpcellid],
              pseudotime = pseudotime[tmpcellid],
              Sample = p,
              Variable = design[rownames(design) == p, variable.d],
              g = g,
              stringsAsFactors = FALSE
            )
          }
        })
        ld = do.call(rbind, linedlist)
        ldlist[[g]] <- ld
      }
      pd <- do.call(rbind, pdlist)
      pd[, 'Variable'] <- as.factor(pd[, 'Variable'])
      ld <- do.call(rbind, ldlist)
      ld[, 'Variable'] <- as.factor(ld[, 'Variable'])
      ld <- ld[order(ld$pseudotime),]
      if (!is.na(sep)) {
        pd$g <- gsub(sep, '', pd$g)
        ld$g <- gsub(sep, '', ld$g)
        pd$g <- factor(pd$g, levels = gsub(sep, '', gene))
        ld$g <- factor(ld$g, levels = gsub(sep, '', gene))
      } else {
        pd$g <- factor(pd$g, levels = gene)
        ld$g <- factor(ld$g, levels = gene)
      }
      if (is.null(variable)) {
        if (plot.point) {
          p <- ggplot() +
            geom_point(
              data = pd,
              aes(
                x = pseudotime,
                y = expr,
                color = Sample
              ),
              alpha = point.alpha,
              size = point.size
            ) +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
          
        } else {
          p <- ggplot() +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        }
        
      } else {
        if (plot.point) {
          p <- ggplot() +
            geom_point(
              data = pd,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable
              ),
              alpha = point.alpha,
              size = point.size
            ) +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable,
                group = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        } else {
          p <- ggplot() +
            geom_line(
              data = ld,
              aes(
                x = pseudotime,
                y = expr,
                color = Variable,
                group = Sample
              ),
              alpha = line.alpha,
              size = line.size
            )
        }
      }
      
      p <- p +
        theme_classic() +
        xlab(x.lab) + ylab(y.lab) +
        labs(color = variable) +
        facet_wrap( ~ g, scales = a, ncol = ncol) +
        theme(
          legend.spacing.y = unit(0.01, 'cm'),
          legend.spacing.x = unit(0.01, 'cm'),
          legend.key.size = unit(0.1, "cm")
        ) +
        guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1)))
      
      if (continuous) {
        p <- p + scale_color_viridis(discrete = TRUE, direction = -1)
      } else {
        if (length(unique(ld[, 'Variable'])) > 8) {
          p <-
            p + scale_color_manual(values = colorRampPalette(rev(brewer.pal(8, palette)))(length(unique(ld[, 'Variable']))))
        } else {
          p <-
            p + scale_color_manual(values = brewer.pal(8, palette)[seq_len(length(unique(ld[, 'Variable'])))])
        }
      }
      if (axis.text.blank) {
        p <-
          p + theme(axis.text = element_blank(), axis.ticks = element_blank())
      } else {
        p <-
          p + scale_x_continuous(breaks = c(min(pd$pseudotime), max(pd$pseudotime)))
      }
      print(p)
    }
    
  }
