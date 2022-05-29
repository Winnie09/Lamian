#' Calculate FDR difference and AUC.
#'
#' This function is used to calculate the FDR difference (Fdr.Diff) and
#' Area Under Sensitifity-Real_FDR Curve (AUC).
#'
#' @export
#' @import stats
#' @return a data frame. Columns are c('Fdr.Diff','Area').
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param sensfdr the output object or SensFdr() which contains the information
#' of sensitivity, real FDR, and reported FDR.
#' @param cutoff the reported FDR cutoff for calculating the Fdr.Diff and AUC.
#' Ranges of reported FDR greater than this cutoff will not be considered.
#' @examples
#' a = AreaUnderSensFdr(sensfdr = data.frame(Sensitivity = abs(rnorm(10)), Real_FDR = abs(rnorm(10)), Reported_FDR = abs(rnorm(10))))
AreaUnderSensFdr <- function(sensfdr, cutoff = 0.25){
  ## [1:3] "Sensitivity" "Real_FDR" "Reported_FDR"
  if (cutoff != 1){
    tmp <- sensfdr
    bound <- approx(x=tmp[,3],y=tmp[,2],xout=cutoff)$y
    tmp <- rbind(tmp[tmp[,3] < cutoff,2:3],c(bound,cutoff))
    tmp <- unique(tmp)
    diff <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)-cutoff*cutoff/2   ## (area under Real_FDR ~ Reported_FDR)-cutoff*cutoff/2
    tmp <- sensfdr
    bound <- approx(x=tmp[,2],y=tmp[,1],xout=cutoff)$y
    tmp <- rbind(tmp[tmp[,2] < cutoff,1:2],c(bound,cutoff))
    area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm=T)/cutoff
  } else if (cutoff == 1){
    tmp <- sensfdr
    bound <- tmp[nrow(tmp),2]
    tmp <- rbind(tmp[,2:3],c(bound,1))
    tmp <- unique(tmp)
    diff <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm = T)-1/2   ## (area under Real_FDR ~ Reported_FDR)-cutoff*cutoff/2
    tmp <- sensfdr
    bound <- tmp[nrow(tmp),1]
    tmp <- rbind(tmp[,1:2],c(bound,1))
    area <- sum(sapply(2:nrow(tmp),function(i) (tmp[i-1,1]+tmp[i,1])*(tmp[i,2]-tmp[i-1,2])/2),na.rm=T)
  }
  
  res <- c(diff, area)
  names(res) = c('Fdr.Diff','Area')
  res
}
