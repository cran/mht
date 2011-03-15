\name{quantiletest}
\alias{quantiletest}
\docType{package}
\title{Call to the C++ function
}
\section{Warning}{You must not call this function unless you are sure. Otherwise you can have problem with the C++ code.}

\description{
Simulation for the calculation of the quantiles of the procbol procedure}
\usage{
quantiletest(k,a,b,d,n,p,F,I,maxq)}

\arguments{
\item{k}{Positive integer.}
\item{a}{Input matrix of size n*p.}
\item{b}{Matrix of the k first columns of data.}
\item{d}{Orthonormal matrix obtained with the Gram-Schmidt algorithm from the k first columns of data.}
\item{n}{Positive integer.}
\item{p}{Positive integer.}
\item{F}{Matrix of size IT*maxq.}
\item{I}{Number of simulations in the calculation of the quantile.}
\item{maxq}{Number of maximum multiple hypotheses testing to do.}
}

\details{IT simulations stored in F.}
