\name{man_tree_data.rda}
\docType{data}
\alias{man_tree_data.rda}
\title{
'man_tree_data.rda.Rd' data
}
\description{
A list of input files: expression matrix, cellanno dataframe, pseudotime, and design matrix.
}
\details{
This data was subsetted from the bone marrow samples in Human Cell Atlas that were used in the manuscript.
}
\examples{
man_tree_data <- system.file("data", "man_tree_data.rda", package = "Lamian")
str(man_tree_data)
}
\keyword{datasets}
