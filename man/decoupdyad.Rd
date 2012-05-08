\name{dyadiqueordre}
\alias{dyadiqueordre}
\docType{package}
\title{Dyadic algorithm to order the variables
}
\description{
Dyadic algorithm using the Bolasso technique to order the variables }
\usage{
dyadiqueordre(data,Y,m,maxordre,var_nonselect,maxit,showtest,showordre)
}

\arguments{
\item{data}{Input matrix of dimension n*p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{m}{Number of bootstrapped iteration of the Lasso. Default is \code{m}=100.}
\item{maxordre}{Number of variables one wishes to order. Default is min(n-1,p-1).}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{maxit}{Number of times the algorithm is restarted if unabled to order \code{maxordre} variables. Default is \code{maxit}=5}
\item{showtest}{Logical value. If TRUE, show the number of regularization parameters tested to show the steps of the algorithm. Default is FALSE.}
\item{showordre}{Logical value. If TRUE, show the ordered variables at each step of the algorithm. Default is TRUE.}
}


\details{The algorithm starts from a large regularization parameter given by one run of Lasso. It proceeds by dyadic splitting until one variable is isolated; e.g one variable alone achieve a frequency of 1; it is the first ordered variable. And so on until \code{maxordre} variables are ordered.
}
\seealso{\code{\link{bolasso}}}

\value{
\item{data}{The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{ordre}{The order obtained on the variables.}
\item{prob}{Number of times the algorithm was not able to ordered the desired number of variables.}
}

