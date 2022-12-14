% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLEce.R
\name{summary.MLEce}
\alias{summary.MLEce}
\alias{print.summary.MLEce}
\title{Summarizing MLEce function}
\usage{
\method{summary}{MLEce}(object, ...)

\method{print}{summary.MLEce}(x, digits = max(3, getOption("digits") - 3), ...)
}
\arguments{
\item{object}{an object of class "MLEce" made by the function \code{MLEce}.}

\item{...}{not used, but exists because of the compatibility.}

\item{x}{an object of class "MLEce".}

\item{digits}{a numeric number of significant digits.}
}
\value{
\code{summary} describes information about MLEce. (quantile statistics, correlation, estimates)
}
\description{
\code{summary} method for a class "MLEce".
}
\examples{
datt = rBiGam(100, c(1,4,5))
res = MLEce(datt, "BiGam", boots = 50)
summary(res)

}
