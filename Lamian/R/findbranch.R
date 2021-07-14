#' Enumerate pseudotime tree branches
#'
#' This function will take the  minimum spanning tree object from the output (a list) of TSCAN::exprmclust() and enumerate all tree branches.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @import stats igraph
#' @importFrom igraph degree
#' @return a list of branch
#' @param mst minimum spanning tree object. It is MSTtree element in the output (a list) of TSCAN::exprmclust().
#' @param order a list of cell orders. Each element is a vector of cells. It is the output (a list) of TSCAN::TSCANorder().
#' @param origin a numeric number indicating the origin cluster.

findbranch <- function(mst, order, origin) {
  deg <- igraph::degree(mst)
  vertex <- names(deg[which(deg > 2 | deg == 1)])
  if (!origin %in% vertex)
    vertex <- c(origin, vertex)
  eg <-
    expand.grid(seq_len(length(vertex)), seq_len(length(vertex)))
  eg <- eg[eg[, 1] < eg[, 2], ]
  eg = data.frame(vertex[eg[, 1]], vertex[eg[, 2]], stringsAsFactors = FALSE)
  
  tmpbranch <- lapply(seq(1, nrow(eg)), function(i) {
    sp <- shortest_paths(mst, from = eg[i, 1], to = eg[i, 2])$vpath[[1]]
    if (sum(vertex %in% sp) == 2)
      as.vector(sp)
  })
  tmpbranch <- tmpbranch[sapply(tmpbranch, length) > 0]
  
  allbranch <-
    gsub('backbone ', '', gsub('branch: ', '', names(order)))
  allbranch <- sapply(allbranch, function(i)
    strsplit(i, ',')[[1]])
  allbranch <- paste0(names(allbranch), collapse = ' ')
  newbranch <- sapply(tmpbranch, function(i) {
    tmp <- paste0(i, collapse = ',')
    if (!grepl(tmp, allbranch)) {
      rev(i)
    } else {
      i
    }
  })
  return(newbranch)
}
