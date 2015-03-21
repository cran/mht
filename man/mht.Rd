\name{mht}
\alias{mht}
\docType{package}
\title{
Multiple testing procedure for non-ordered variable selection}
\description{
Performs multiple hypotheses testing in a linear model}
\usage{mht(data,Y,var_nonselect,alpha,sigma,maxordre,ordre,m,show,IT,maxq)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of \code{data}. Default is 1, the selection is not performed on the intercept.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{sigma}{Value of the variance if it is known; 0 otherwise. Default is 0.}
\item{maxordre}{Number of variables to be ordered. Default is min(n/2-1,p/2-1).}
\item{ordre}{Several possible algorithms to order the variables, ordre=c("bolasso","pval","pval_hd","FR"). "bolasso" uses the dyadic algorithm with the Bolasso technique \code{\link{dyadiqueordre}}, "pval" uses the p-values obtained with a regression on the full set of variables (only when p<n), "pval_hd" uses marginal regression, "FR" uses Forward Regression. Default is "bolasso".}
\item{m}{Number of bootstrap iteration of the Lasso. Only used if the algorithm is set to "bolasso". Default is m=100.}
\item{show}{Vector of logical values, show=(showordre,showtest,showresult). Default is (1,0,1). If showordre==TRUE, show the ordered variables at each step of the algorithm.
If showtest==TRUE, show the number of regularization parameters tested to show the advancement of the dyadic algorithm. Only use if the algorithm is set to "bolasso".
if showresult==TRUE, show the value of the statistics and the estimated quantile at each step of the procedure.}
\item{IT}{Number of simulations for the calculation of the quantile. Default is 1000.}
\item{maxq}{Number of maximum multiple hypotheses testing to perform. Default is log(min(n,p)-1,2).}
}


\details{ \code{mht} is a two-step procedure that performs variable selection in high dimensional linear model. The first step orders the variables taking into account the vector of observations \code{Y}. The second step finds a cut-off between the relevant variables (high rank) and the irrelevant ones (low rank) through multiple hypotheses testing.\cr
The input maxordre is not to be forgotten: the more variables to order, the more difficult for the algorithm to distinguish which noisy variable is more important that another noisy variable. It is advised to limit maxordre to \code{p/2} or \code{n/2} if they are large. The parameter maxq can be useful for large value of \code{n}, it is advised to limit it to 5-6 in order to minimize computational time (for the calculation of the quantile).
}

\value{A 'mht' object is returned for which the methods \code{refit}, \code{predict} and  \code{plot} are available.
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
\item{ordre}{Order obtained on the \code{maxordre} variables.}
\item{ordrebeta}{The full order on all the variables.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for each values of \code{alpha}.}
\item{quantile}{The estimated quantiles used in the second step of the procedure.}
\item{call}{The call that produced this object.}
}

\seealso{\code{\link{predict.mht}}, \code{\link{refit.mht}}, \code{\link{plot.mht}}}

\references{
Multiple hypotheses testing for variable selection; F. Rohart 2011\cr}


\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

mod=mht(x,y,alpha=c(0.1,0.05),maxordre=15)
mod
}
}