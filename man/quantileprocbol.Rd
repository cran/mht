\name{quantileprocbol}
\alias{quantileprocbol}
\docType{package}
\title{Calculation of the quantiles of the procbol procedure
}
\description{
Calculation of the quantiles of the procbol procedure}
\usage{
quantileprocbol(data,k,alpha,IT,maxq,sigma)
}

\arguments{
\item{data}{Input matrix of dimension n*p; each row is an observation vector.}
\item{k}{Positive integer.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 1000.}
\item{maxq}{Number of maximum multiple hypotheses testing to do. Default is log(min(n,p)-1,2).}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}

}


\details{Calculation of the quantile of the statistic that tests the null hypothesis Hk: there is no more signal remaining in the variables (k+1,...,p).
}

\value{
\item{quantile}{The estimated quantiles.}
\item{nbrprob}{Number of times there was not enough simulation to estimate the alpha-quantile.}
}

