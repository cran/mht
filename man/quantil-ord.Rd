\name{quantil_ord}
\alias{quantil_ord}
\docType{package}
\title{Calculation of the type I error to use in the proc_ord procedure
}
\description{
Calculation of the type I error to use in the proc_ord procedure}
\usage{
quantil_ord(n,p,k,alpha,IT)
}

\arguments{
\item{n}{Positive integer.}
\item{p}{Positive integer, can be lower or larger than n.}
\item{k}{Positive integer, has to be lower than n and p.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 20000.}

}


\details{Calculation of the estimated type I error to be used in the second step of the procedure proc_ord in order to have a final test of level alpha.
}

\value{
\item{quantile}{The estimated type I error to be used in the second step of the procedure. It is of size the length of alpha.}
}

