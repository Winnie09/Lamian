\name{expdata.rda}
\docType{data}
\alias{expdata.rda}
\title{
'expdata.rda.Rd'data
}
\description{
A list of input files: expression matrix, cellanno dataframe, pseudotime, and design matrix.
}
\details{
This data was subsetted from the bone marrow samples in Human Cell Atlas that were used in the manuscript.
}
\examples{
expdata <- system.file("data", "expdata.rda", package = "Lamian")
str(expdata)
}
\keyword{datasets}
