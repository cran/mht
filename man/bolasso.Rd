\name{bolasso}
\alias{bolasso}
\docType{package}
\title{
Bolasso: Bootstrapped Lasso}
\description{
Perform a bootstrapped Lasso on some random subsamplings of the input data}
\usage{
bolasso(data,Y,mu,m,probaseuil,penalty.factor,random)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{mu}{Positive regularization sequence to be used for the Lasso.}
\item{m}{Number of bootstrap iteration of the Lasso. Default is m=100.}
\item{probaseuil}{A frequency threshold for selecting the most stable variables over the \code{m} boostrap iteration of the Lasso. Default is 1.}
\item{penalty.factor}{Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables except the intercept.}
\item{random}{optionnal parameter, matrix of size n*m. If \code{random} is provided, the \code{m} bootstrap samples are constructed from its m columns.}
}


\details{The Lasso from the \code{\link[glmnet]{glmnet}} package is performed with the regularization parameter mu over m bootstrap samples. An appearance frequency is obtained which shows the predictive power of each variable. It is calculated as the number of times a variables has been selected by the Lasso over the \code{m} bootstrap iteration.
}

\value{A 'bolasso' object is returned for which the method \code{plot} is available.

\item{data}{A list containing:
\itemize{\item{X - The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{Y - the input response vector}
\item{means.X - Vector of means of the input data matrix.}
\item{sigma.X - Vector of variances of the input data matrix.}
}}

\item{ind}{Set of selected variables for the regularization \code{mu} and the threshold \code{probaseuil}.}
\item{frequency}{Appearance frequency of each variable; number of times each variables is selected over the m bootstrap iterations.}
}

\seealso{\code{\link{plot.bolasso}}, \code{\link{dyadiqueordre}}}


\references{
Model-consistent sparse estimation through the bootstrap; F. Bach 2009 \cr}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(1,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

mod=bolasso(x,y,mu=seq(1.5,0.1,-0.1))
mod
}
}
