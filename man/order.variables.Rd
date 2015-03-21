\name{order.variables}
\alias{order.variables}
\docType{package}
\title{
Rank the variables}
\description{
Gives an order to the variables and rearrange the input matrix following that order.}
\usage{order.variables(data,Y,maxordre,ordre=c("bolasso","pval","pval_hd","FR"),
var_nonselect,m,showordre)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{maxordre}{Number of variables to be ordered. Default is min(n/2-1,p/2-1).}
\item{ordre}{Several possible algorithms to order the variables, ordre=c("bolasso","pval","pval_hd","FR"). "bolasso" uses the dyadic algorithm with the Bolasso technique \code{\link{dyadiqueordre}}, "pval" uses the p-values obtained with a regression on the full set of variables (only when p<n), "pval_hd" uses marginal regression, "FR" uses Forward Regression. Default is "bolasso".}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of \code{data}. Default is 1, the selection is not performed on the intercept.}
\item{m}{Number of bootstrapped iteration of the Lasso. Only use if the algorithm is set to "bolasso". Default is m=100.}
\item{showordre}{If showordre=TRUE, show the variables being ordered at each step of the algorithm.}
}


\details{Rank the variables of \code{data} taking into account the vector of observations \code{Y} and rearrange the input matrix following that order.}

\value{
\item{data}{A list containing:
\itemize{\item{X - The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{Y - the input response vector}
\item{means.X - Vector of means of the input data matrix.}
\item{sigma.X - Vector of variances of the input data matrix.}
}}

\item{data_ord}{Input data matrix rearranged by ORDREBETA}
\item{ORDRE}{Gives the \code{maxordre} most important variables of the data matrix.}
\item{ORDREBETA}{Gives the order on all the variables of the data matrix (either arbitrary completion of ORDRE -`Bolasso' and `FR', or the true order -`pval' and `pval_hd').}
}

\references{
Multiple hypotheses testing for variable selection; F. Rohart 2011\cr}


\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

res.bolasso=order.variables(x,y,maxordre=15,ordre="bolasso")
res.pval=order.variables(x,y,ordre="pval")
res.FR=order.variables(x,y,maxordre=15,ordre="FR")
res.pval.hd=order.variables(x,y,maxordre=15,ordre="pval_hd")

}
}