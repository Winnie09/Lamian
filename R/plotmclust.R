#' Plot the pseudotime tree structure. Overwrite TSCAN::plotmclust().
#' 
#' This function takes the output from infer_tree_structure() to generate a plot. For example, if cells are using pc as dimension reduction, then the plot will be y-axis as PC2, a-axis as PC1 (or otherwise specified by inputs x and y), and dots are cells colored by clusters. A pseudotime tree structure is shown on top of the cells. 
#' 
#' @return A ggplot2 plot
#' @export
#' @import ggplot2 plyr grid igraph
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param x indicates which column in the low-dimrension representation to be used as x-axis coordinates.
#' @param y indicates which column in the low-dimrension representation to be used as y-axis coordinates.
#' @param MSTorder The arbitrary order of cluster to be shown on the plot.
#' @param show_tree Whether to show the links between cells connected in the minimum spanning tree.
#' @param show_full_tree Whether to show the full tree or not. Only useful when show_tree=T. Overrides MSTorder.
#' @param show_cell_names Whether to draw the name of each cell in the plot.
#' @param cell_name_size The size of cell name labels if show_cell_names is TRUE.
#' @param cell_point_size The size of cell point.
#' @param markerexpr The gene expression used to define the size of nodes.
#' @param showcluster logical. indicates whether to show clusters on the plot.
#' @param x.lab x-axis label
#' @param y.lab y-axis label
#' @param subset.cell a character vector of the names of the selected cells to plot, when NULL, plot all cells.
#' @examples
#' data(man_tree_res)
#' plotmclust(man_tree_res)
plotmclust <- function (mclustobj, x = 1, y = 2, MSTorder = NULL, show_tree = T, 
                        show_full_tree = T, show_cell_names = F, cell_name_size = 3, cell_point_size = 3, markerexpr = NULL, showcluster = T, x.lab = 'PC1', y.lab = 'PC2', subset.cell = NULL){   
  
  color_by = "State"
  mypalette = colorRampPalette(brewer.pal(9,'Set1'))
  lib_info_with_pseudo <- data.frame(State = mclustobj$clusterid, 
                                     sample_name = names(mclustobj$clusterid))
  lib_info_with_pseudo$State <- factor(lib_info_with_pseudo$State)
  if (is.null(subset.cell)){
    S_matrix <- mclustobj$pcareduceres
  } else {
    S_matrix <- mclustobj$pcareduceres[subset.cell, , drop = FALSE]  ###
  }
  
  pca_space_df <- data.frame(S_matrix[, c(x, y)])
  colnames(pca_space_df) <- c("pca_dim_1", "pca_dim_2")
  pca_space_df$sample_name <- row.names(pca_space_df)
  edge_df <- merge(pca_space_df, lib_info_with_pseudo, by.x = "sample_name", 
                   by.y = "sample_name")
  edge_df$Marker <- markerexpr[edge_df$sample_name]
  if (!is.null(markerexpr)) {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2, 
                                    size = Marker))
    if (showcluster) {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE) +
        scale_color_manual(values = mypalette(length(unique(edge_df[,color_by]))))
    }
    else {
      g <- g + geom_point(na.rm = TRUE, color = "green")
    }
  }
  else {
    g <- ggplot(data = edge_df, aes(x = pca_dim_1, y = pca_dim_2))
    if (showcluster) {
      g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE, size = cell_point_size) +
        scale_color_manual(values = mypalette(length(unique(edge_df[,color_by]))))
    }
    else {
      g <- g + geom_point(na.rm = TRUE, size = cell_point_size)
    }
  }
  if (show_cell_names) {
    g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
  }
  if (show_tree) {
    clucenter <- mclustobj$clucenter[, c(x, y)]
    clulines <- NULL
    if (show_full_tree) {
      alledges <- as.data.frame(get.edgelist(mclustobj$MSTtree), 
                                stringsAsFactors = F)
      alledges[, 1] <- as.numeric(alledges[, 1])
      alledges[, 2] <- as.numeric(alledges[, 2])
      for (i in 1:nrow(alledges)) {
        clulines <- rbind(clulines, c(clucenter[alledges[i, 
                                                         1], ], clucenter[alledges[i, 2], ]))
      }
    }
    else {
      if (is.null(MSTorder)) {
        clutable <- table(mclustobj$clusterid)
        alldeg <- degree(mclustobj$MSTtree)
        allcomb <- expand.grid(as.numeric(names(alldeg)[alldeg == 
                                                          1]), as.numeric(names(alldeg)[alldeg == 1]))
        allcomb <- allcomb[allcomb[, 1] < allcomb[, 2], 
                           ]
        numres <- t(apply(allcomb, 1, function(i) {
          tmp <- as.vector(get.shortest.paths(mclustobj$MSTtree, 
                                              i[1], i[2])$vpath[[1]])
          c(length(tmp), sum(clutable[tmp]))
        }))
        optcomb <- allcomb[order(numres[, 1], numres[, 
                                                     2], decreasing = T)[1], ]
        MSTorder <- get.shortest.paths(mclustobj$MSTtree, 
                                       optcomb[1], optcomb[2])$vpath[[1]]
      }
      for (i in 1:(length(MSTorder) - 1)) {
        clulines <- rbind(clulines, c(clucenter[MSTorder[i], 
                                                ], clucenter[MSTorder[i + 1], ]))
      }
    }
    clulines <- data.frame(x = clulines[, 1], xend = clulines[, 
                                                              3], y = clulines[, 2], yend = clulines[, 4])
    g <- g + geom_segment(aes_string(x = "x", xend = "xend", 
                                     y = "y", yend = "yend", size = NULL), data = clulines, 
                          size = 1)
    clucenter <- data.frame(x = clucenter[, 1], y = clucenter[, 
                                                              2], id = 1:nrow(clucenter))
    g <- g + geom_text(aes_string(label = "id", x = "x", 
                                  y = "y", size = NULL), data = clucenter, size = 10)
  }
  g <- g + guides(colour = guide_legend(override.aes = list(size = 5))) + 
    # xlab(paste0("PCA_dimension_", x)) + ylab(paste0("PCA_dimension_", y)) + 
    theme(panel.border = element_blank(), axis.line = element_line()) + 
    theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
    theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
    # theme(legend.position = "top", legend.key.size = unit(0.3, 
    #     "in"), legend.text = element_text(size = 20), legend.title = element_text(size = 20), 
    #     legend.box = "vertical") + 
    # theme(legend.key = element_blank()) + 
    theme(legend.position = 'none') +
    theme(panel.background = element_rect(fill = "white")) + 
    theme(axis.text.x = element_text(size = 17, color = "black"), 
          axis.text.y = element_text(size = 17, color = "black"), 
          axis.title.x = element_text(size = 20, vjust = -1), 
          axis.title.y = element_text(size = 20, vjust = 1), 
          plot.margin = unit(c(1, 1, 1, 1), "cm")) +
    xlab(x.lab) + ylab(y.lab)
  g
}


