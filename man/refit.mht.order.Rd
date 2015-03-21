\name{refit.mht.order}
\alias{refit.mht.order}
\title{
Refit a \code{mht.order} object}
\description{
Refit a \code{mht.order} object for a new observation \code{Ynew} and/or a new order \code{ordrenew}}
\usage{
\method{refit}{mht.order}(object,Ynew,ordrenew,IT,var_nonselect,sigma,showresult,...)
}

\arguments{
\item{object}{Object of class "mht.order"  as obtained from \code{\link{mht}}.}
\item{Ynew}{Response variable of length n.}
\item{ordrenew}{Vector from which the varibles are to be ordered, it can be a partial order. If absent, data is considers to be already ordered; Default is (1,2,..,p).}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 10000.}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
\item{\dots}{not used.}
}


\details{see \code{\link{mht.order}} for details.
}

\value{A 'mht.order' object is returned.
\item{data}{A list containing:
\itemize{\item{X - The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{Y - the input response vector}
\item{means.X - Vector of means of the input data matrix.}
\item{sigma.X - Vector of variances of the input data matrix.}
}}

\item{coefficients}{Matrix of the estimated coefficients. Each row concerns a specific user level alpha.}
\item{residuals}{Matrix of the residuals. Each row concerns a specific user level \code{alpha}.}
\item{relevant_var}{Set of the relevant variables. Each row concerns a specific user level \code{alpha}}
\item{fitted.values}{Matrix of the fitted values, each column concerns a specific user level \code{alpha}.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for each values of \code{alpha}.}
\item{quantile}{The estimated type I error to be used in the second step of the procedure in order to have a test of level alpha, each column stands for one test. See F.Rohart (2011) for details.}
\item{call}{The call that has been used.}
\item{call.old}{The call that produced the initial 'object'.}
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
ynew=x\%*\%beta+rnorm(100)

# mht.order
mod.order=mht.order(x,y)

# refit mht.order on a new vector of observation
mod.refit=refit(mod.order,ynew)
}
}