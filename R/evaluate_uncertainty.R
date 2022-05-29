#' Evaluate the pseudotime tree uncertainty
#'
#' This function is designed to evaluate the pseudotime tree uncertainty, as one of the main functions in lamain module 1.
#'
#' @param inferobj the output object from function infer_tree_structure().
#' @param n.permute: a numeric number of permutation in the permutation test.
#' @param subset.cell a character vector of the names of the selected cells where boostrap will happen on. If NULL, then boostrap from all the cells.
#' @export
#' @return a list of results
#' @author Wenpin Hou <whou10@jhu.edu>
#' @examples
#' data(man_tree_res)
#' a <- evaluate_uncertainty(inferobj = man_tree_res, n.permute = 3)
evaluate_uncertainty <-
  function(inferobj,
           n.permute,
           subset.cell = NULL,
           design = NULL) {
    
    if (is.null(subset.cell)) {
      pr <- inferobj$pca
    } else {
      pr <- inferobj$pca[subset.cell,]
    }
    newbranch <- inferobj$branch
    js.cut <- inferobj$js.cut
    oc.cut <- inferobj$oc.cut
    pt <- inferobj$pseudotime
    ord <- inferobj$order
    alls <- inferobj$allsample
    ctcomplist <- reproduce.js <- reproduce.oc <- corr.score <- list()
    for (pmid in seq(1, n.permute)) {
      ## boostrap cells
      ## set.seed(pmid)
      bstid <- sample(seq_len(nrow(pr)), nrow(pr), replace = TRUE)
      bstid <- unique(bstid)
      pr.pm <- pr[bstid, ]
      
      ## cluster cells
      clu <-
        mykmeans(pr.pm, number.cluster = max(inferobj$clusterid))$cluster ###
      
      ## build pseudotime
      mcl.pm <-
        exprmclust(t(pr.pm), cluster = clu, reduce = FALSE) ###
      
      ## select origin cluster
      pt.pm.mean <-
        tapply(pt[names(mcl.pm[['clusterid']])], list(mcl.pm[['clusterid']]), mean)
      start.cluster <- names(which.min(pt.pm.mean))
      
      ## construct pseudotime
      ord.pm <-
        TSCANorder(
          mcl.pm,
          startcluster = start.cluster,
          listbranch = TRUE,
          orderonly = TRUE
        )
      
      pt.pm <-
        unlist(sapply(sapply(ord.pm, length), function(i)
          seq(1, i)))
      names(pt.pm) <- unname(unlist(ord.pm))
      
      ## plot pseudotime
      pd = data.frame(pc1 = pr[, 1],
                      pc2 = pr[, 2],
                      time = as.numeric(pt.pm[rownames(pr)]))
      
      # get candidate branches
      newbranch.pm <-
        findbranch(mst = mcl.pm$MSTtree,
                   order = ord.pm,
                   origin = start.cluster)
      
      ## compare two MST
      js <- sapply(seq(1, length(newbranch)), function(i) {
        id <-
          which(sapply(paste0(names(ord), ','), function(k)
            grepl(paste0(
              paste0(newbranch[[i]], collapse = ','), ','
            ), k)))[1]
        cells <- ord[[id]]
        b.ori <-
          intersect(unlist(sapply(newbranch[[i]], function(k)
            names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
        sapply(seq(1, length(newbranch.pm)), function(j) {
          id <-
            which(sapply(paste0(names(ord.pm), ','), function(k)
              grepl(paste0(
                paste0(newbranch.pm[[j]], collapse = ','), ','
              ), k)))[1]
          cells <- ord.pm[[id]]
          b.pm <-
            intersect(unlist(sapply(newbranch.pm[[j]], function(k)
              names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
          js <-
            length(intersect(b.pm, b.ori)) / length(union(b.pm, b.ori))
        })
      })
      oc <- sapply(seq(1, length(newbranch)), function(i) {
        id <-
          which(sapply(paste0(names(ord), ','), function(k)
            grepl(paste0(
              paste0(newbranch[[i]], collapse = ','), ','
            ), k)))[1]
        cells <- ord[[id]]
        b.ori <-
          intersect(unlist(sapply(newbranch[[i]], function(k)
            names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
        sapply(seq(1, length(newbranch.pm)), function(j) {
          id <-
            which(sapply(paste0(names(ord.pm), ','), function(k)
              grepl(paste0(
                paste0(newbranch.pm[[j]], collapse = ','), ','
              ), k)))[1]
          cells <- ord.pm[[id]]
          b.pm <-
            intersect(unlist(sapply(newbranch.pm[[j]], function(k)
              names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
          oc <-
            length(intersect(b.pm, b.ori)) / min(length(b.pm), length(b.ori))
        })
      })
      corr <- sapply(seq(1, length(newbranch)), function(i) {
        id <-
          which(sapply(paste0(names(ord), ','), function(k)
            grepl(paste0(
              paste0(newbranch[[i]], collapse = ','), ','
            ), k)))[1]
        cells <- ord[[id]]
        b.ori <-
          intersect(unlist(sapply(newbranch[[i]], function(k)
            names(inferobj$clusterid)[inferobj$clusterid == k])), cells)
        
        sapply(seq(1, length(newbranch.pm)), function(j) {
          id <-
            which(sapply(paste0(names(ord.pm), ','), function(k)
              grepl(paste0(
                paste0(newbranch.pm[[j]], collapse = ','), ','
              ), k)))[1]
          cells <- ord.pm[[id]]
          b.pm <-
            intersect(unlist(sapply(newbranch.pm[[j]], function(k)
              names(mcl.pm$clusterid)[mcl.pm$clusterid == k])), cells)
          ov = intersect(b.ori, b.pm)
          cor(pt[ov], pt.pm[ov])
        })
      })
      corr[is.na(corr)] <- 0
      colnames(corr) <-
        colnames(oc) <-
        colnames(js) <- paste0('original', seq(1, length(newbranch)))
      
      ## get js binary to match branches
      js.binary <- get_binary(js, js.cut)
      corr.score[[pmid]] <- corr * js.binary
      js.melt <- melt(js.binary)
      js.melt <- js.melt[js.melt[, 3] != 0, ]
      colnames(js.melt) <-
        c('permutation.branch', 'original.branch', 'matched')
      reproduce.js[[pmid]] <- as.character(js.melt[, 2])
      
      ## get oc binary to match branches
      oc.binary <- get_binary(oc, oc.cut)
      oc.melt <- melt(oc.binary)
      oc.melt <- oc.melt[oc.melt[, 3] != 0, ]
      reproduce.oc[[pmid]] <- as.character(oc.melt[, 2])
      
      ## samples cell compositions
      ctcomp.new <-
        matrix(0, nrow = length(unique(alls)), ncol = length(newbranch))
      colnames(ctcomp.new) <-
        paste0('origin', seq(1, length(newbranch)))
      rownames(ctcomp.new) <- unique(alls)
      
      if (nrow(js.melt) > 0) {
        ctcomp <-
          sapply(seq_len(nrow(js.melt)), function(i) {
            ## corrected  from 2 to 1. 2020/08/31
            c <- names(clu)[clu %in% newbranch.pm[[js.melt[i, 1]]]]
            ctcomp <- rep(0, length(unique(alls)))
            names(ctcomp) <- unique(alls)
            ctcomp[names(table(alls[c]))] <- table(alls[c])
            ctcomp
          })
        colnames(ctcomp) <- paste0('origin', js.melt[, 2])
        ctcomp <- ctcomp / rowSums(ctcomp)
        
        ctcomp.new[rownames(ctcomp), colnames(ctcomp)] <-
          ctcomp  ## sample by #branch: rowSums = 1  
      }
      ctcomplist[[pmid]] <- t(ctcomp.new)
    }
    
    reproduce.js <- unlist(reproduce.js)
    js.perc <- rep(0, length(newbranch))
    js.perc[as.numeric(names(table(reproduce.js)))] <-
      table(reproduce.js) / n.permute
    names(js.perc) <- newbranch
    
    reproduce.oc <- unlist(reproduce.oc)
    oc.perc <- rep(0, length(newbranch))
    oc.perc[as.numeric(names(table(reproduce.oc)))] <-
      table(reproduce.oc) / n.permute
    names(oc.perc) <- newbranch
    
    corr.score.m <- do.call(rbind, corr.score)
    corr.score.v <- colSums(corr.score.m) / n.permute
    names(corr.score.v) <- newbranch
    
    sort((js.perc + oc.perc) / 2)
    
    detection.rate <-
      data.frame(detection.rate = (js.perc + oc.perc[names(js.perc)]) / 2,
                 stringsAsFactors = FALSE)
    detection.rate[which(detection.rate[, 1] > 1), 1] <- 1
    sample.cellcomp.mean <-
      apply(simplify2array(ctcomplist), seq_len(2), mean)
    sample.cellcomp.sd <-
      apply(simplify2array(ctcomplist), seq_len(2), sd)
    
    rownames(sample.cellcomp.mean) <-
      newbranch[as.numeric(sub('origin', '', rownames(sample.cellcomp.mean)))]
    rownames(sample.cellcomp.sd) <-
      newbranch[as.numeric(sub('origin', '', rownames(sample.cellcomp.sd)))]
    
    if (!is.null(design)) {
      dv <- as.numeric(as.factor(design[colnames(ctcomplist[[1]]), ])) - 1
      sample.cellcomp.pvalue <-
        sapply(ctcomplist, function(i)
          apply(i, 1, function(j) {
            if (length(unique(j)) == 1) {
              1
            } else {
              t.test(j[dv == 1], j[dv == 0])$p.value
            }
          }))
      sample.cellcomp.pvalue <-
        rowMeans(sample.cellcomp.pvalue < 0.05)
      names(sample.cellcomp.pvalue) <-
        newbranch[as.numeric(sub('origin', '', names(sample.cellcomp.pvalue)))]
    }
    
    if (length(newbranch[[length(newbranch)]]) == 2) {
      name <-
        paste0('c(', newbranch[[length(newbranch)]][1], ',', newbranch[[length(newbranch)]][2], ')')
      rownames(detection.rate)[nrow(detection.rate)] <-
        rownames(sample.cellcomp.mean)[nrow(sample.cellcomp.mean)] <-
        rownames(sample.cellcomp.sd)[nrow(sample.cellcomp.sd)] <- name
    }
    
    if (!is.null(design)) {
      result <- list(
        detection.rate = detection.rate,
        sample.cellcomp.mean = sample.cellcomp.mean,
        sample.cellcomp.sd = sample.cellcomp.sd,
        sample.cellcomp.pvalue = sample.cellcomp.pvalue
      )
    } else {
      result <- list(
        detection.rate = detection.rate,
        sample.cellcomp.mean = sample.cellcomp.mean,
        sample.cellcomp.sd = sample.cellcomp.sd
      )
    }
    return(result)
  }
