\name{mht.order}
\alias{mht.order}
\docType{package}
\title{
Multiple testing procedure for ordered variable selection}
\description{
Performs multiple hypotheses testing for ordered variable selection.}
\usage{mht.order(data,Y,ordre,var_nonselect,alpha,IT,sigma,showresult)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{ordre}{Vector from which the varibles are to be ordered, it can be a partial order. If absent, data is considers to be already ordered; Default is (1,2,..,p).}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of \code{data}. Default is 1, the selection is not performed on the intercept.}
\item{alpha}{A user supplied type I error sequence. Default is \code{alpha}=(0.1,0.05)}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 10000.}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}
\item{showresult}{Logical value. if TRUE, shows the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
}

\details{The details of the procedure are in 'Multiple hypotheses testing for variable selection; F. Rohart 2011'. If \code{showresult}=TRUE, the statistics and quantile are printed through the algorithm. If the statistic is greater than the quantile, the test is rejected (takes the value 1). The procedure stops when the null huypothesis is accepted (all alternative hypotheses are 0).\cr
    The statistics to test the null hypotheses are different whether the variance \code{sigma} is known.

}

\value{A 'mht.order' object is returned for which the methods \code{predict}, \code{refit} and \code{plot} are available.
\item{data}{A list containing:
\itemize{\item{X - The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{Y - the input response vector}
\item{means.X - Vector of means of the input data matrix.}
\item{sigma.X - Vector of variances of the input data matrix.}
}}

\item{coefficients}{Matrix of the estimated coefficients. Each row concerns a specific user level \code{alpha}.}
\item{residuals}{Matrix of the residuals. Each row concerns a specific user level \code{alpha}.}
\item{relevant_var}{Set of the relevant variables. Each row concerns a specific user level \code{alpha}}
\item{fitted.values}{Matrix of the fitted values, each column concerns a specific user level \code{alpha}.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for each values of \code{alpha}.}
\item{quantile}{The estimated type I error to be used in the second step of the procedure in order to have a test of level alpha, each column stands for one test. See F.Rohart (2011) for details.}
\item{call}{The call that has been used.}
}

\seealso{\code{\link{predict.mht.order}},  \code{\link{refit.mht.order}}, \code{\link{plot.mht.order}}}


\references{
Adaptive tests of linear hypotheses by model selection; Baraud & al 2002\cr
Multiple hypotheses testing for variable selection; F. Rohart 2011\cr}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

mod.order=mht.order(x,y,ordre=5:1,alpha=c(0.1,0.05))
mod.order
}
}