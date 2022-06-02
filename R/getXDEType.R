#' Obtain the XDE types for the XDE significant genes in sample covariate test.
#'
#' This function is used to obtain the XDE types for the XDE genes in sample covariate test.
#'
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a vector of XDE Types, including meanSig, trendSig, bothSig, other, and nonXDE. The names of this vector are the XDE genes.
#' @export
#' @param testobj output object of the function tespt().
#' @param cutoff the FDR cutoff to select "statistically significant" genes.
#' @examples
#' data(mantestobj)
#' a = getXDEType(mantestobj)

getXDEType <- function(testobj, cutoff = 0.05) {
  res <- testobj$statistics
  if (toupper(testobj$test.type) == 'VARIABLE') {
    diffType <- sapply(rownames(res), function(g) {
      if (res[g, grep('^fdr.*overall$', colnames(res))] < cutoff) {
        if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff &
            res[g, grep('^fdr.*trendDiff$', colnames(res))] > cutoff) {
          'meanSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] > cutoff &
                   res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff) {
          'trendSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff &
                   res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff) {
          'bothSig'
        } else {
          'other'
        }
      } else {
        'nonXDE'
      }
    })
  } else if (toupper(testobj$test.type) == 'TIME') {
    print('ConstantTimeTest does not lead to XDEType!')
  }
  return(diffType)
}
