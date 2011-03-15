\name{procbol}
\alias{procbol}
\docType{package}
\title{
Multiple testing procedure for variable selection}
\description{
Performed a multiple hypotheses testing}
\usage{
procbol(data,Y,var_nonselect,alpha,maxordre,choix_ordre=c("bolasso","pval","pval_hd"),m,maxit,showordre,IT,maxq,showtest,showresult)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{var_nonselect}{Number of variables to not select. They have to be in the first columns. Default is 1, the selection is not performed on the intercept.}
\item{alpha}{A user supplied type I error sequence. Default is (0.1,0.05).}
\item{maxordre}{Number of variables one wishes to order. Default is min(n-1,p-1).}
\item{choix_ordre}{Three possible algorithms to order the variables. "bolasso" uses the dyadic algorithm with the Bolasso, "pval" use the p-values obtained with a regression on the full set of variables (only when p<n), the last "pval_hd" makes regression on one variable at a time. Default is "bolasso".}
\item{m}{Number of bootstrapped iteration of the Lasso. Only use if the algorithm is set to "bolasso". Default is m=100.}
\item{maxit}{Number of times the algorithm is restarted if unabled to order maxordre variables. Default is 5}
\item{showordre}{Logical value. If TRUE, show the ordered variables at each step of the algorithm. Default is TRUE.}
\item{IT}{Number of simulations in the calculation of the quantile. Default is 20000.}
\item{maxq}{Number of maximum multiple hypotheses testing to do. Default is log(min(n,p)-1,2).}
\item{showtest}{Logical value. Only use if the algorithm is set to "bolasso". if TRUE, show the number of regularization parameters tested to show the advancement of the dyadic algorithm. Default is FALSE.}
\item{showresult}{Logical value. if TRUE, show the value of the statistics and the estimated quantile at each step of the procedure. Default is TRUE.}
}


\details{Procbol is a two-step procedure that order the variables taking into account the response variable and performs a multiple hypotheses testing to cut-off the variables in two set: the relevant and the irrelevant variables. The input maxordre is not to be forgotten, indeed the more they are variables to order, the more difficult it will be for the algorithm to distinguish which noise variable is more important that another noise variable. It is advised to limit maxordre to p/2 or n/2 if they are large. The parameter maxq can be usefull for large value of n, it is advised to limit it to 5-6 in order to minimizes computational time in the calculation of the quantile.
}

\value{
\item{data}{The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{relevant_var}{Set of the relevant variables of the input matrix.}
\item{ordre}{Order obtained on the maxordre variable.}
\item{ordrebeta}{The full order on all the variables.}
\item{kchap}{Vector containing the length of the estimated set of relevant variables, for the matrix containing the intercept, for each values of alpha.}
\item{quantile}{The estimated quantiles used in the second step of the procedure.}
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod=procbol(x,y,alpha=c(0.1,0.05),maxordre=15)
}
}