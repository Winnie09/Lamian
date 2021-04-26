#' Obtain the DDG types for the DDG in sample covariate test.
#' 
#' This function is used to obtain the DDG types for the DDG in sample covariate test.
#' 
#' @author Wenpin Hou <whou10@jhu.edu>
#' @return a vector of DDGType. The names of this vector are the DDG. 
#' @export
#' @param testobj output object of the function tespt().
#' @param cutoff the FDR cutoff to select "statistically significant" genes. 
#' @examples
#' data(mantestobj)
#' a = getDDGType(mantestobj)

getDDGType <- function(testobj, cutoff = 0.05){
  res <- testobj$statistics
  if (testobj$test.type == 'Variable' | testobj$test.type == 'variable'){
    diffType <- sapply(rownames(res), function(g){
      if (res[g, grep('^fdr.*overall$', colnames(res))] < cutoff){
        if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] > cutoff){
          'meanSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] > cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff){
          'trendSig'
        } else if (res[g, grep('^fdr.*meanDiff$', colnames(res))] < cutoff & res[g, grep('^fdr.*trendDiff$', colnames(res))] < cutoff){
          'bothSig'
        } else {
          'other'
        }
      } else {
        'nonDDG'
      }
    })
  } else if (testobj$test.type == 'Time' | testobj$test.type == 'time'){
    print('ConstantTimeTest does not lead to DDGType!')
  }
  return(diffType)
}

