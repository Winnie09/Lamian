\name{man_tree_res}
\docType{data}
\alias{man_tree_res.rda}
\title{
'man_tree_res.rda.Rd' data
}
\description{
a list of data as the example input for the function infer_tree_structure().
}
\details{
a list containing gene by cell expression matrix, low-dimension representation (PCA) of the cells, and the cell annotation where the first column are the cell names, the second column are the sample names, and the row names are the cell names. 
}
\examples{
man_tree_res <- system.file("data", "man_tree_res.rda", package = "Lamian")
str(man_tree_res)
}
\keyword{datasets}
