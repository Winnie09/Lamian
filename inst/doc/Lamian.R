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
res = infer_tree_structure(pca = man_tree_data[['pca']], 
                           expression = man_tree_data[['expression']], 
                           cellanno = man_tree_data[['cellanno']], 
                           origin.marker = c('CD34'), 
                           xlab='Principal component 1', 
                           ylab = 'Principal component 2')

## -----------------------------------------------------------------------------
names(res)

## ----fig_plotmclust, fig.height = 5, fig.width = 6, fig.align = "center"------
plotmclust(res, cell_point_size = 0.5, 
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
Res <- lamian.test(expr = expdata$expr, 
                   cellanno = expdata$cellanno, 
                   pseudotime = expdata$pseudotime, 
                   design = expdata$design, 
                   test.type = 'variable', 
                   permuiter = 5)

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

## ----fig_plotClusterMeanAndDiff, fig.height = 14, fig.width = 7, fig.align = "center"----
## plotClusterMeanAndDiff
plotClusterMeanAndDiff(Res)

## ---- eval=FALSE--------------------------------------------------------------
#  ## GO analysis
#  goRes <- GOEnrich(testobj = Res, type = 'variable')

## ----fig_sct_plotDiffFitHm3, fig.height = 6, fig.width = 9, fig.align = "center"----
colnames(Res$populationFit[[1]]) <- colnames(Res$populationFit[[1]]) <- colnames(Res$expr) 
plotXDEHm(Res, cellWidthTotal = 180, cellHeightTotal = 350, subsampleCell = FALSE, sep = ':.*')

## -----------------------------------------------------------------------------
Res <- lamian.test(expr = expdata$expr, 
                   cellanno = expdata$cellanno, 
                   pseudotime = expdata$pseudotime, 
                   design = expdata$design, 
                   test.type = 'time', 
                   permuiter = 5)
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

## ----fig_sct_plotclustermean, fig.height = 7, fig.width = 8, fig.align = "center"----
plotClusterMean(testobj=Res, cluster = Res$cluster, type = 'time')

## ----eval = FALSE-------------------------------------------------------------
#  goRes <- GOEnrich(testobj = Res, type = 'time', sep = ':.*')
#  plotGOEnrich(goRes = goRes)

## ----fig_sct_plotTDEHm, fig.height = 4.5, fig.width = 10, fig.align = "center"----
plotTDEHm(
  Res,
  subsampleCell  = FALSE,
  showCluster = TRUE,
  type = 'time',
  cellWidthTotal = 200,
  cellHeightTotal = 200
)

## -----------------------------------------------------------------------------
Res <-
  cellPropTest(
    cellanno = expdata$cellanno,
    pseudotime = expdata$pseudotime,
    design = expdata$design,
    ncores = 4
  )

## -----------------------------------------------------------------------------
Res$statistics

## -----------------------------------------------------------------------------
sessionInfo()

