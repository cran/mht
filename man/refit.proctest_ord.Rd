\name{refit.proctest_ord}
\alias{refit.proctest_ord}
\title{
Refit a \code{proctest_ord} object}
\description{
Refit a \code{proctest_ord} object for a new observation Ynew}
\usage{
\method{refit}{proctest_ord}(object,Ynew,ordrenew,IT,var_nonselect,sigma,showresult,...)
}

\arguments{
\item{object}{Object of class "proctest_ord".}
\item{Ynew}{Response variable of length n.}
\item{ordrenew}{Vector that gave the order on the variables, could be a partial order. If absent, data is considers to be already ordered; Default is (1,2,..,p).}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 20000.}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
\item{\dots}{not used.}
}


\details{see \code{\link{proc_ord}} for details.
}

\value{A 'proctest_ord' object is returned.
\item{data}{List of the user-data: the scaled matrix used in the algorithm -the first column being (1,...,1)- and Y.}
\item{coefficients}{Matrix of the estimated coefficients. Each row concerns a specific user level alpha.}
\item{relevant_var}{Set of the relevant variables.}
\item{fitted.values}{Matrix of the fitted values, each column concerns a specific user level alpha.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for the matrix containing the intercept, for each values of alpha.}
\item{quantile}{The estimated type I error to be used in the second step of the procedure in order to have a test of level alpha, each column stands for one test.}
\item{call}{The call that has been used.}
\item{call.old}{The call that produced the initial 'object'.}
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
ynew=x\%*\%beta+rnorm(100)
mod=procbol(x,y,maxordre=15)
mod2=refit(mod,ynew,maxordre=15)
}
}