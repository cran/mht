\name{decompbaseortho}
\alias{decompbaseortho}
\docType{package}
\title{
Decomposition of Gram-schmidt
}
\description{
Orthonormalization of an input matrix with the Gram-Schmidt algorithm.}
\usage{
decompbaseortho(data)
}

\details{Performed an orthonormalization of the input matrix, recording the useless columns.
}
\arguments{
\item{data}{Input matrix of dimension n * p; each column is a parameter.}
}


\value{
\item{U}{The orthonormal basis obtained from \code{data}.}
\item{nonind}{Set of variables with no contribution. }
}


\examples{
\dontrun{
x1=rnorm(100)
x2=rnorm(100)
x3=cbind(x1,x2,x1+x2)
dec=decompbaseortho(x3)
}
}