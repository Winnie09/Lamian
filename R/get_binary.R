#' Intermediate function to serve evaluate_uncertainty().
#'
#' This function serves as an intermediate function to serve evaluate_uncertainty(). It matches boostrap and origin branches.
#'
#' @return a binary matrix
#' @param  matrix #boostrap.branch * #origin.branch, values are js or oc
#' @param matrix.cut js or oc null distribution cutoff
#' @author Wenpin Hou <whou10@jhu.edu>

get_binary <- function(matrix, matrix.cut) {
  matrix.binary <- sapply(seq(1, ncol(matrix)), function(c) {
    (matrix[, c] > matrix.cut[c]) + 0
  })
  while (length(which(rowSums(matrix.binary) > 1)) > 0 |
         length(which(colSums(matrix.binary) > 1)) > 0) {
    dup.id <- which(rowSums(matrix.binary) > 1)
    if (length(dup.id) == 1) {
      addid <- which.max(matrix[dup.id,])
      matrix.binary[dup.id,] <- 0
      matrix.binary[dup.id, addid] <- 1
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id) {
        addid <- which.max(matrix[dup.i,])
        matrix.binary[dup.i,] <- 0
        matrix.binary[dup.i, addid] <- 1
      }
    }
    
    dup.id <- which(colSums(matrix.binary) > 1)
    if (length(dup.id) == 1) {
      addid <- which.max(matrix[, dup.id])
      matrix.binary[, dup.id] <- 0
      matrix.binary[addid, dup.id] <- 1
    } else if (length(dup.id) > 1) {
      for (dup.i in dup.id) {
        addid <- which.max(matrix[, dup.i])
        matrix.binary[, dup.i] <- 0
        matrix.binary[addid, dup.i] <- 1
      }
    }
  }
  return(matrix.binary)
}
