#' Infer the pseudotime tree strucutre of the input data
#'
#' This function takes a low dimensional-reduction representation (for example, pca), the cell annotation, etc. as inputs, and then infer the pseudotime tree structures for further assessing the uncertainty of each of the pseudotime tree branches.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a list
#' @export
#' @import TSCAN scattermore RColorBrewer grDevices parallel
#' @importFrom grDevices pdf
#' @importFrom grDevices dev.off
#' @param  pca cell by principal component (pc) matrix. Principal components reduction of the cells.
#' @param cellanno 2-column or 3-column dataframe/matrix, first column is cell name, second column is sample, third column (if exists) is cell type.
#' @param expression only useful when users want to use highly expressed marker genes to determine the starting point of pseudotime. It is a gene by cell expression matrix. The values should be library-size-normalized and log-transformed expression values. They can either be imputed or non-imputed.
#' @param origin.marker: a character or a character vector specifying the highly expressed marker genes used to identify the origin cluster in pseudotime inference.If NA (default), then not used.
#' @param   origin.celltype only useful when users have specify the celltypes for the cells in cellanno[,3], and want to determin the starting point of pseudotime by this celltype.  It is a character of origin celltype. It should be one of the element in cellanno[,3] if it will be used. If NA (default), then not used.
#' @param   number.cluster the number of clusters in cell clustering that will be used in trajectory inference. If NA (default),the number of clusters will be determined automatically by elbow's method.
#' @param  plotdir plot directory for storing the figures of cluster.pdf and pseudotime.pdf. If NA (default), figures will not be generated.
#' @param  xlab the x-axis labels for the figures.
#' @param  ylab the y-axis labels for the figures.
#' @param  max.clunum the maximum number of clusters in the elbew's method.
#' @examples
#' data(hca_bm_pca)
#' data(hca_bm_saver)
#' data(hca_bm_cellanno)
#' res = infer_tree_structure(pca = hca_bm_pca, expression = hca_bm_saver, cellanno = hca_bm_cellanno, origin.marker = c('CD34'), xlab='Principal component 1', ylab = 'Principal component 2')

infer_tree_structure <-
  function(pca,
           cellanno,
           expression,
           origin.marker = NA,
           origin.celltype = NA,
           number.cluster = NA,
           plotdir = NA,
           xlab = 'PC1',
           ylab = 'PC2',
           max.clunum = 50,
           kmeans.seed = 12345,
           ncores = detectCores()) {
    alls <- cellanno[, 2]
    names(alls) <- cellanno[, 1]
    ## set.seed(12345)
    sdev <- apply(pca, 2, sd)
    x <- seq_len(max.clunum)
    optpoint <- which.min(sapply(2:max.clunum, function(i) {
      x2 <- pmax(0, x - i)
      sum(lm(sdev[seq_len(max.clunum)] ~ x + x2)$residuals ^ 2)
    }))
    pcadim = optpoint + 1
    pr <- pca[, seq_len(pcadim)]  # 7
    
    ## clustering
    clu <-
      mykmeans(pr, maxclunum = 50, number.cluster = number.cluster, seed = kmeans.seed, ncores = ncores)$cluster
    table(clu)
    pd = data.frame(x = pr[, 1],
                    y = pr[, 2],
                    cluster = as.factor(clu[rownames(pr)]))
    mypalette = colorRampPalette(brewer.pal(9, 'Set1'))
    if (!is.na(plotdir)) {
      pdf(paste0(plotdir, 'cluster.pdf'),
          width = 3,
          height = 2.1)
        # pdf(paste0(plotdir, 'cluster.pdf'),
        #   width = 6,
        #   height = 4.2)
      print(
        ggplot(data = pd, aes(
          x = pd[,1], y = pd[,2], color = pd[,3]
        )) +
          geom_scattermore() +
          scale_color_manual(values = mypalette(max(clu))) +
          theme_classic() +
          theme(
            legend.spacing.y = unit(0.01, 'cm'),
            legend.spacing.x = unit(0.01, 'cm'),
            legend.key.size = unit(0.1, "cm")
          ) +
          xlab(xlab) + ylab(ylab)
      )
      dev.off()
    }
    ### mclust
    mcl <- exprmclust(t(pr), cluster = clu, reduce = FALSE)
    
    # --------------------
    # construct pseudotime
    # --------------------
    ## find origin
    if (!is.na(origin.celltype)) {
      pd = data.frame(
        x = pr[, 1],
        y = pr[, 2],
        clu = as.factor(mcl$clusterid),
        celltype = cellanno[match(rownames(pr), cellanno[, 1]), 3],
        stringsAsFactors = FALSE
      )
      tab <- table(pd[, 3:4])
      tab <- tab / rowSums(tab)
      pd <- melt(tab)
      pd$clu <-
        factor(as.character(pd$clu), levels = seq(1, max(pd$clu)))
      tmp <- pd[pd$celltype == origin.celltype,]
      origin.cluster <- as.numeric(tmp[which.max(tmp[, 3]), 1])
    } else {
      pd <-
        data.frame(
          x = pr[, 1],
          y = pr[, 2],
          clu = as.factor(mcl$clusterid),
          mark = colMeans(expression[origin.marker, , drop = FALSE]),
          stringsAsFactors = FALSE
        )
      tmp <- tapply(pd[, 4], pd[, 3], mean)
      
      origin.cluster <- as.numeric(which.max(tmp))
    }
    ## construct pseudotime
    ord <-
      TSCANorder(
        mcl,
        startcluster = origin.cluster,
        listbranch = TRUE,
        orderonly = TRUE
      )
    pt <- unlist(sapply(sapply(ord, length), function(i)
      seq(1, i)))
    names(pt) <- unname(unlist(ord))
    
    # ## plot pseudotime
    pd = data.frame(pc1 = pca[, 1],
                    pc2 = pca[, 2],
                    time = as.numeric(pt[rownames(pca)]))
    
    if (!is.na(plotdir)) {
      pdf(paste0(plotdir, 'pseudotime.pdf'),
          width = 3.1,
          height = 2.1)
      print(
        ggplot(data = pd, aes(
          x = pd[,1], y = pd[,2], color = time
        )) +
          geom_scattermore() +
          scale_color_gradient(low = 'yellow', high = 'blue') +
          xlab(xlab) + ylab(ylab) +
          theme_classic()
      )
      dev.off()
    }
    
    # ------------------------------------------------------------
    # get candidate branches to test reproducibility, 20200726 >>
    # ------------------------------------------------------------
    newbranch <-
      findbranch(mst = mcl$MSTtree,
                 order = ord,
                 origin = origin.cluster)
    # -----------------------------------------------------
    # Evaluate robustness of tree branches using resampling
    # -----------------------------------------------------
    # null distribution of Jaccard index, overlap coefficient
    js.null <- lapply(seq(1, length(newbranch)), function(i) {
      b.ori <-
        unlist(sapply(newbranch[[i]], function(c)
          names(mcl$clusterid[mcl$clusterid == c])))
      tmp <- sapply(seq(1, 1e3), function(j) {
        ## set.seed(j)
        b.pm <- sample(rownames(pr), length(b.ori))
        length(intersect(b.pm, b.ori)) / length(union(b.pm, b.ori))
      })
    })
    js.cut <- sapply(js.null, quantile, 0.99)
    
    oc.null <- lapply(seq(1, length(newbranch)), function(i) {
      b.ori <-
        unlist(sapply(newbranch[[i]], function(c)
          names(mcl$clusterid[mcl$clusterid == c])))
      tmp <- sapply(seq(1, 1e3), function(j) {
        ## set.seed(j)
        b.pm <- sample(rownames(pr), length(b.ori))
        length(intersect(b.pm, b.ori)) / min(length(b.pm), length(b.ori))
      })
    })
    oc.cut <- sapply(oc.null, quantile, 0.99)
    mcl$pseudotime <- pt
    mcl$branch <- newbranch
    mcl$js.cut <- js.cut
    mcl$oc.cut <- oc.cut
    mcl$pca <- pr
    mcl$order <- ord
    mcl$allsample <- alls
    return(mcl)
  }
