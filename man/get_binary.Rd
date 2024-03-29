% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_binary.R
\name{get_binary}
\alias{get_binary}
\title{Intermediate function to serve evaluate_uncertainty().}
\usage{
get_binary(matrix, matrix.cut)
}
\arguments{
\item{matrix}{#boostrap.branch * #origin.branch, values are js or oc}

\item{matrix.cut}{js or oc null distribution cutoff}
}
\value{
a binary matrix
}
\description{
This function serves as an intermediate function to serve evaluate_uncertainty(). It matches boostrap and origin branches.
}
\author{
Wenpin Hou <whou10@jhu.edu>
}
