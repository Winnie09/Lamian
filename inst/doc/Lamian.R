## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## -----------------------------------------------------------------------------
# load in Lamian
options(warn=-1)
suppressMessages(library(Lamian))

## -----------------------------------------------------------------------------
data(hca_bm_saver)
data(hca_bm_pca)
data(hca_bm_cellanno)

## -----------------------------------------------------------------------------
str(hca_bm_saver)

## -----------------------------------------------------------------------------
str(hca_bm_pca)

## -----------------------------------------------------------------------------
str(hca_bm_cellanno)

## -----------------------------------------------------------------------------
res = infer_tree_structure(pca = hca_bm_pca, expression = hca_bm_saver, cellanno = hca_bm_cellanno, origin.marker = c('CD34'), xlab='Principal component 1', ylab = 'Principal component 2')

## -----------------------------------------------------------------------------
names(res)

## ----fig_plotmclust, fig.height = 2.6, fig.width = 3, fig.align = "center"----
plotmclust(res, cell_point_size = 0.1, x.lab = 'Principal component 1', y.lab = 'Principal component 2')

## -----------------------------------------------------------------------------
result <- evaluate_uncertainty(res, n.permute=3)
names(result)

## -----------------------------------------------------------------------------
result[[1]]

## -----------------------------------------------------------------------------
result[[2]]

## -----------------------------------------------------------------------------
result[[3]]

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

## -----------------------------------------------------------------------------
Res <- lamian.test(expr = expdata$expr, cellanno = expdata$cellanno, pseudotime = expdata$pseudotime, design = expdata$design, test.type = 'variable')

## -----------------------------------------------------------------------------
## get differential dynamic genes statistics 
stat <- Res$statistics
head(stat)
stat <- stat[order(stat[,1], -stat[,3]), ]
diffgene <- rownames(stat[stat[, grep('^fdr.*overall$', colnames(stat))] < 0.05,])
str(diffgene)

## -----------------------------------------------------------------------------
## population fit
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'variable')

## -----------------------------------------------------------------------------
## clustering
Res$covariateGroupDiff <- getCovariateGroupDiff(testobj = Res, gene = diffgene)
Res$cluster <- clusterGene(Res, gene = diffgene, type = 'variable', k=5)
table(Res$cluster)

## ----fig_plotClusterMeanAndDiff, fig.height = 8, fig.width = 3.5, fig.align = "center"----
## plotClusterMeanAndDiff
plotClusterMeanAndDiff(Res)

## ---- fig_plotGOEnrich, fig.height = 2.6, fig.width = 3, fig.align = "center", eval=FALSE----
#  ## GO analysis
#  goRes <- GOEnrich(testobj = Res, type = 'variable')
#  plotGOEnrich(goRes = goRes, fdr.cutoff = 0.05, fc.cutoff = 2)

## ----fig_plotFitHm, fig.height = 8, fig.width = 12, fig.align = "center"------
plotFitHm(Res, type = 'variable', cellWidthTotal = 200, cellHeightTotal = 350)

## ----fig_sct_plotDiffFitHm3, fig.height = 9, fig.width = 16, fig.align = "center"----
plotDiffFitHm3(Res, cellWidthTotal = 180, cellHeightTotal = 350)

## -----------------------------------------------------------------------------
Res <- lamian.test(expr = expdata$expr, cellanno = expdata$cellanno, pseudotime = expdata$pseudotime, design = expdata$design, test.type = 'time', ncores = 1)
names(Res)

## -----------------------------------------------------------------------------
head(Res$statistics)

## -----------------------------------------------------------------------------
diffgene <- rownames(Res$statistics)[Res$statistics[,1] < 0.05]
str(diffgene)

## -----------------------------------------------------------------------------
Res$populationFit <- getPopulationFit(Res, gene = diffgene, type = 'time')

## -----------------------------------------------------------------------------
Res$cluster <- clusterGene(Res, gene = diffgene, type = 'time', k=3)

## ----fig_sct_plotclustermean, fig.height = 2.3, fig.width = 3.3, fig.align = "center"----
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')

## ----eval = FALSE-------------------------------------------------------------
#  goRes <- GOEnrich(testobj = Res, type = 'time', sep = ':.*')
#  plotGOEnrich(goRes = goRes)

## ----fig_sct_plotFitHm, fig.height = 4.5, fig.width = 10, fig.align = "center"----
plotFitHm(Res, subsampleCell  = F, showCluster = T, type = 'time', cellWidthTotal = 200, cellHeightTotal = 200)

## -----------------------------------------------------------------------------
Res <- cell_prop_test(cellanno = expdata$cellanno, pseudotime = expdata$pseudotime, design = expdata$design, ncores = 4)

## -----------------------------------------------------------------------------
Res$statistics

## -----------------------------------------------------------------------------
sessionInfo()

