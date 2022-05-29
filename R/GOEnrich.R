#' Obtain enrich Gene Ontology (GO) terms
#'
#' This function is used to obtain enrich Gene Ontology (GO) terms (currently only for human genes)
#'
#' @import  topGO org.Hs.eg.db igraph grid
#' @importFrom topGO graph
#' @importFrom topGO algorithm
#' @importFrom topGO depth
#' @export
#' @return a list. Each element is a data frame about the GO terms and their statistics.
#' @author Wenpin Hou <whou10@jhu.edu>
#' @param testobj output object of lamian_test()
#' @param fdr.cutoff FDR cutoff to select "statistically significant" GO terms.
#' @param k number of clusters to be used for clustering, if the gene clustering results have not yet been included in testobj.
#' @param use.clusters If TRUE (default), use the clusters in testobj.
#' @param type the type of the differential test. One of c('Time', 'Variable'). Use for selecting the columns in statistics of FDR. Only useful when "cluster" in not in testobj.
#' @param species currently only work for "human". will include "mouse".
#' @examples
#' data(mantestobj)
GOEnrich <-
  function(testobj,
           fdr.cutoff = 0.05,
           k = 5,
           use.clusters = TRUE,
           type = 'Variable',
           species = 'human',
           sep = ':.*') {
    fdr.col.id <- grep('^fdr.*overall$', colnames(testobj$statistics))
    if (toupper(type) == 'VARIABLE') {
      fdr <- testobj$statistics[, fdr.col.id]
    } else if (toupper(type) == 'TIME') {
      fdr <- testobj$statistics[, 'fdr.overall']
    }
    if (sum(fdr < fdr.cutoff) == 0) {
      print('There is no differential genes! GoEnrich stopped.')
      break
    } else {
      if (use.clusters) {
        if ('cluster' %in% names(testobj)) {
          clu <- testobj$cluster
        } else {
          clu <-
            clusterGene(
              testobj,
              gene = rownames(testobj$statistics)[testobj$statistics[, fdr.col.id] < fdr.cutoff],
              type = type,
              k = k
            )
        }
        diffgeneList <-
          sapply(sort(unique(clu)), function(i) {
            #########
            names(clu)[clu == i]
          })
      } else {
        diffgeneList <-
          list(diffgene = rownames(testobj$statistics)[fdr < fdr.cutoff])
      }
      
      
      resList <- lapply(diffgeneList, function(diffgene) {
        allgene <- rownames(testobj$expr)
        gl <- sub(sep, '', diffgene)
        back <- sub(sep, '', allgene)
        geneList <- factor(as.integer(back %in% gl))
        names(geneList) <- back
        suppressMessages({
          GOdata <-
             new(
              "topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = function(a) {
                a
              },
              annot = annFUN.org,
              mapping = "org.Hs.eg.db",
              ID = "Symbol"
            )
          resultFisher <-
            topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
          sigres <-
            topGO::GenTable(
              GOdata,
              classicFisher = resultFisher,
              topNodes = length(resultFisher@score),
              orderBy = "classicFisher",
              numChar = 1000
            )
        })
        sigres$classicFisher[sigres$classicFisher == "< 1e-30"] <- 0
        sigres <- sigres[sigres$Annotated >= 10,]
        sigres$FDR <-
          stats::p.adjust(sigres$classicFisher, method = "fdr")
        ptcount <- 0
        fc <-
          ((sigres[, "Significant"] + ptcount) / (sum(GOdata@allScores[GOdata@feasible] ==
                                                        1) + ptcount)) / ((sigres[, "Annotated"] + ptcount) / (sum(GOdata@feasible) +
                                                                                                                 ptcount))
        sigres <- data.frame(sigres, FC = fc)
        sigres <- sigres[order(sigres$FDR,-sigres$FC),]
      })
      names(resList) <- sort(unique(clu))
      return(resList)
    }
    
  }
