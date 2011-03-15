\name{proc_ord}
\alias{proc_ord}
\docType{package}
\title{
Pre-ordered procedure}
\description{
Performed a multiple hypotheses testing when the variables are ordered regardless of the response variable.}
\usage{
proc_ord(data,Y,ordre,alpha,IT,showresult)}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{ordre}{Vector that gave the order on the variables, could be a partial order. If absent, data is considers to be already ordered; Default is (1,2,..,p).}
\item{alpha}{A user supplied type I error sequence.}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 20000.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
}


\value{
\item{data}{The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{relevant_var}{Set of the relevant variables of the input matrix.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for the matrix containing the intercept, for each values of alpha.}
\item{quantile}{The estimated type I error to be used in the second step of the procedure in order to have a test of level alpha, each column stands for one test.}
}

\references{
Adaptive tests of linear hypotheses by model selection; Baraud & al \cr}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod=proc_ord(x,y,ordre=5:1,alpha=c(0.1,0.05))
}
}