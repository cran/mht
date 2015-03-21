\name{decompbaseortho}
\alias{decompbaseortho}
\docType{package}
\title{
Gram-Schmidt algorithm
}
\description{
Orthonormalization of an input matrix with the Gram-Schmidt algorithm.}
\usage{
decompbaseortho(data)
}

\details{Performs an orthonormalization of the input matrix, recording the columns that are linear combination of the previous ones.
}
\arguments{
\item{data}{Input matrix of dimension n * p; each column is a variable.}
}


\value{
\item{U}{The orthonormal basis obtained from \code{data}.}
\item{nonind}{Set of variables with no contribution. }
\item{trueind}{Set of variables with contribution. }
\item{rank}{Rank of the input matrix, calculated with the function \code{\link[Matrix]{rankMatrix}} .}

}


\examples{
\dontrun{
x1=rnorm(100)
x2=rnorm(100)
x3=cbind(x1,x2,x1+x2)

dec=decompbaseortho(x3)
dec$nonind
dec$trueind
dec$rank
}
}