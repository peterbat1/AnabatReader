% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AnabatTools.R
\name{print.Anabat}
\alias{print.Anabat}
\title{Generic function to print summary data for an object of class Anabat}
\usage{
\method{print}{Anabat}(AnabatData)
}
\arguments{
\item{AnabatData}{An object of class Anabat as produced by the function \code{\link{readAnabat}}}
}
\description{
This is convenience function which calls \code{\link{summary.Anabat}} and produces the same output to the console.
}
\examples{
d <- readAnabat("blah")
print(d)
}

