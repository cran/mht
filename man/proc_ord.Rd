\name{proc_ord}
\alias{proc_ord}
\docType{package}
\title{
Pre-ordered procedure}
\description{
Performed a multiple hypotheses testing for ordered variable selection.}
\usage{proc_ord(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{ordre}{Vector that gave the order on the variables, could be a partial order. If absent, data is considers to be already ordered; Default is (1,2,..,p).}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{alpha}{A user supplied type I error sequence. Default is \code{alpha}=(0.1,0.05)}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 20000.}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
}

\details{The details of the procedure are in 'Multiple hypotheses testing for variable selection; F. Rohart 2011'. If \code{showresult}=TRUE, the statistics and quantile are given. If the statistics is greater than the quantile, the test is rejected and takes the value 1. The procedure stops when all tests are accepted (with values 0).}

\value{A 'proctest_ord' object is returned for which the methods \code{predict} and \code{refit} are available.
\item{data}{List of the user-data: the scaled matrix used in the algorithm -the first column being (1,...,1)- and Y.}
\item{coefficients}{Matrix of the estimated coefficients. Each row concerns a specific user level alpha.}
\item{relevant_var}{Set of the relevant variables.}
\item{fitted.values}{Matrix of the fitted values, each column concerns a specific user level alpha.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for the matrix containing the intercept, for each values of alpha.}
\item{quantile}{The estimated type I error to be used in the second step of the procedure in order to have a test of level alpha, each column stands for one test.}
\item{call}{The call that has been used.}
}

\references{
Adaptive tests of linear hypotheses by model selection; Baraud & al 2002\cr
Multiple hypotheses testing for variable selection; F. Rohart 2011\cr}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod=proc_ord(x,y,ordre=5:1,alpha=c(0.1,0.05))
}
}