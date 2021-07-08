\name{mantestobj}
\docType{data}
\alias{mantestobj.rda}
\title{
'mantestobj.rda.Rd' data
}
\description{
a list of data as the example output from the functions lamian.test(). This data can be used as the input of functions such as plotting (plotGene(), plotXDEHm(), plotTDEHm, etc..) or downstream analysis (getDDGType(), etc..)
}
\details{
a list containing gene by cell expression matrix, the cell annotation dataframe, pseudotime, parameter estimates from Lamian model, etc.
}
\examples{
mantestobj <- system.file("data", "mantestobj.rda", package = "Lamian")
str(mantestobj)
}
\keyword{datasets}
