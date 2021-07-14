\name{mandata}
\docType{data}
\alias{mandata.rda}
\title{
'mandata.rda.Rd' data
}
\description{
a list of data as the example input for the functions such as lamian.test() or cellPropTest().
}
\details{
a list containing gene by cell expression matrix, low-dimension representation (PCA) of the cells, and the cell annotation where the first column are the cell names, the second column are the sample names, and the row names are the cell names. 
}
\examples{
mandata <- system.file("data", "mandata.rda", package = "Lamian")
str(mandata)
}
\keyword{datasets}
