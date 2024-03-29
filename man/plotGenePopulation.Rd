% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotGenePopulation.R
\name{plotGenePopulation}
\alias{plotGenePopulation}
\title{Plot the population level fitting curves for genes}
\usage{
plotGenePopulation(
  testobj,
  gene = NA,
  type = "time",
  ylim = NA,
  xlim = NA,
  sep = NA,
  free.scale = TRUE,
  palette = "Dark2",
  ncol = NA,
  subSampleNumber = NA,
  line.size = 1,
  axis.text.blank = FALSE
)
}
\arguments{
\item{testobj}{the output object from lamian_test().}

\item{gene}{a vector of genes that need to do the prediction.}

\item{type}{One of c('Time', 'Variable').}

\item{ylim}{y-axis limits to be passed to ggplot.}

\item{xlim}{x-axis limits to be passed to ggplot.}

\item{sep}{a string in the gene names that needs to replaced with blank.}

\item{free.scale}{logical. If TRUE, the y-axis is on free scale.}

\item{palette}{a RColorBrewer palette name.}

\item{ncol}{number of colums for organizing all genes' plot.}

\item{line.size}{the size of the curves.}

\item{axis.text.blank}{logical. If TRUE, leave axis text as blank.}
}
\value{
a plot
}
\description{
This function takes a test output object as input and then plot the population level fitting curves for genes that users specify.
}
\examples{
data(mantestobj)
plotGenePopulation(testobj = mantestobj, type = 'variable', gene = rownames(mantestobj$populationFit[[1]])[seq(1,2)])
}
\author{
Wenpin Hou <whou10@jhu.edu>
}
