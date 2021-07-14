#' Standardize the values by rows
#'
#' This function standardize the values for a matrix by its rows
#'
#' @return a gene by cell (or pseudotime) expression matrix
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param data a matrix.

scalematrix <- function(data) {
  cm <- rowMeans(data)
  csd <- apply(data, 1, sd)
  (data - cm) / csd
}
