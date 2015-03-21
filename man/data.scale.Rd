\name{data.scale}
\alias{data.scale}
\docType{package}
\title{
Scale a matrix }
\description{
Scale the data so each column has mean 0 and variance 1. This function is used as a pre-processing step to prep the data for analysis in all functions of the \code{mht} package.}
\usage{data.scale(data,warning)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{warning}{Logical value. A warning message is printed if the intercept is added. Default is TRUE.}
}


\details{Scale the data so each column has mean 0 and variance 1. If we note \code{x} a column of the output scaled matrix -except the first one which is the intercept, we have \code{sum(x)=0} and \code{sum(x^2)/n=1}.
}

\value{
\item{data}{Scaled data.}
\item{intercept}{Logical value. TRUE if the intercept was already included in the input data.}
\item{means.data}{Vector of means of the input data matrix.}
\item{sigma.data}{Vector of variances of the input data matrix.}
}

\references{
Multiple hypotheses testing for variable selection; F. Rohart 2011\cr}


\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
res=data.scale(x)
x.scaled=res$data
means.x=res$means.data
sigma.x=res$sigma.data
}
}