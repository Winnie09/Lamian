## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# load in Lamian
options(warn=-1)
suppressMessages(library(Lamian))

## -----------------------------------------------------------------------------
data(man_tree_data)

## -----------------------------------------------------------------------------
str(man_tree_data)

## -----------------------------------------------------------------------------
set.seed(12345)
res = infer_tree_structure(pca = man_tree_data[['pca']], 
                           expression = man_tree_data[['expression']], 
                           cellanno = man_tree_data[['cellanno']], 
                           origin.marker = c('CD34'), 
                           number.cluster = 5,
                           xlab='Principal component 1', 
                           ylab = 'Principal component 2')

## -----------------------------------------------------------------------------
names(res)

## ----fig_plotmclust, fig.height = 5, fig.width = 5.5, fig.align = "center", out.width = '40%'----
Lamian::plotmclust(res, cell_point_size = 0.5, 
                   x.lab = 'Principal component 1', 
                   y.lab = 'Principal component 2')

## -----------------------------------------------------------------------------
result <- evaluate_uncertainty(res, n.permute=5)
names(result)

## -----------------------------------------------------------------------------
result[[1]]

## -----------------------------------------------------------------------------
result[[2]]

## -----------------------------------------------------------------------------
result[[3]]

## -----------------------------------------------------------------------------
data = result[[2]]
rownames(data) <-
  c('HSC->lymphocyte', 'HSC->myeloid', 'HSC->erythroid')
design = data.frame(sample = paste0('BM', 1, 8), sex = c(0, rep(1, 4), rep(0, 3)))
branchPropTest(data = data, design = design)

## -----------------------------------------------------------------------------
data(expdata)

## -----------------------------------------------------------------------------
str(expdata$expr)

## -----------------------------------------------------------------------------
head(expdata$cellanno)

## -----------------------------------------------------------------------------
str(expdata$pseudotime)

## -----------------------------------------------------------------------------
print(expdata$design)

## ----eval=FALSE---------------------------------------------------------------
#  saveh5(expr = expdata$expr, pseudotime = expdata$pseudotime, cellanno = expdata$cellanno, path = 'data/multi.h5')

## -----------------------------------------------------------------------------
Res <- lamian_test(
  expr = expdata$expr,
  cellanno = expdata$cellanno,
  pseudotime = expdata$pseudotime,
  design = expdata$design,
  test.type = 'variable',
  testvar = 2,
  permuiter = 5,
  ## this is for permutation test only. Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
  ncores = 1
)

## ----eval=FALSE---------------------------------------------------------------
#  Res <- lamian_test_h5(
#    expr = 'data/multi.h5',
#    cellanno = expdata$cellanno,
#    pseudotime = expdata$pseudotime,
#    design = expdata$design,
#    test.type = 'variable',
#    testvar = 2,
#    permuiter = 5,
#    ## this is for permutation test only. Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.
#    ncores = 1
#  )

## -----------------------------------------------------------------------------
## get differential dynamic genes statistics
stat <- Res$statistics
stat <- stat[order(stat[, 1],-stat[, 3]),]
## identify XDE genes with FDR.overall < 0.05 cutoff
diffgene <-
  rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05, ])

## -----------------------------------------------------------------------------
## population fit
Res$populationFit <-
  getPopulationFit(Res, gene = diffgene, type = 'variable')

## -----------------------------------------------------------------------------
## clustering
Res$covariateGroupDiff <-
  getCovariateGroupDiff(testobj = Res, gene = diffgene)
Res$cluster <-
  clusterGene(Res, gene = diffgene, type = 'variable', k = 5)

## ----fig_sct_plotDiffFitHm3, fig.height = 5, fig.width = 8, fig.align = "center", eval=TRUE----
colnames(Res$populationFit[[1]]) <-
  colnames(Res$populationFit[[2]]) <- colnames(Res$expr)
plotXDEHm(
  Res,
  cellWidthTotal = 180,
  cellHeightTotal = 350,
  subsampleCell = FALSE,
  sep = ':.*'
)

## ----fig_plotClusterMeanAndDiff, fig.height = 8, fig.width = 3.2, fig.align = "center", eval = TRUE----
## plotClusterMeanAndDiff
plotClusterMeanAndDiff(Res)

## -----------------------------------------------------------------------------
Res <- lamian_test(
  expr = expdata$expr,
  cellanno = expdata$cellanno,
  pseudotime = expdata$pseudotime,
  design = expdata$design,
  test.type = 'time',
  permuiter = 5
)

## ----eval=FALSE---------------------------------------------------------------
#  Res <- lamian_test_h5(
#    expr = 'data/multi.h5',
#    cellanno = expdata$cellanno,
#    pseudotime = expdata$pseudotime,
#    design = expdata$design,
#    test.type = 'time',
#    permuiter = 5
#  ) ## this is for permutation test only. Alternatively, we can use test.method = 'chisq' to swich to the chi-square test.

## -----------------------------------------------------------------------------
head(Res$statistics)

## -----------------------------------------------------------------------------
diffgene <- rownames(Res$statistics)[Res$statistics[, 1] < 0.05]

## -----------------------------------------------------------------------------
Res$populationFit <-
  getPopulationFit(Res, gene = diffgene, type = 'time')

## -----------------------------------------------------------------------------
Res$cluster <-
  clusterGene(Res, gene = diffgene, type = 'time', k = 3)

## ----fig_sct_plotTDEHm, fig.height = 4.5, fig.width = 10, fig.align = "center", eval = TRUE----
plotTDEHm(
  Res,
  subsampleCell  = FALSE,
  showCluster = TRUE,
  type = 'time',
  cellWidthTotal = 200,
  cellHeightTotal = 200
)

## ----fig_sct_plotclustermean, fig.height = 5, fig.width = 5, fig.align = "center", out.width = '40%', eval = TRUE----
plotClusterMean(testobj = Res,
                cluster = Res$cluster,
                type = 'time')

## ----eval = TRUE--------------------------------------------------------------
Res <-
  cellPropTest(
    cellanno = expdata$cellanno,
    pseudotime = expdata$pseudotime,
    design = expdata$design[, 1, drop = F],
    ncores = 4,
    test.type = 'Time'
  )

## -----------------------------------------------------------------------------
head(Res$statistics)

## ----eval = TRUE--------------------------------------------------------------
Res <-
  cellPropTest(
    cellanno = expdata$cellanno,
    pseudotime = expdata$pseudotime,
    design = expdata$design,
    ncores = 4,
    test.type = 'Variable',
    testvar = 2
  )

## -----------------------------------------------------------------------------
sessionInfo()

