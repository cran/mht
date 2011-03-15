\name{bolasso}
\alias{bolasso}
\docType{package}
\title{
Bolasso: Bootstrapped Lasso}
\description{
Performed a bootstrapped Lasso}
\usage{
bolasso(data,Y,mu,m,probaseuil,var_nonselect)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{mu}{Positive regularization sequence to be used for the Lasso.}
\item{m}{Number of bootstrapped iteration of the Lasso. Default is m=100.}
\item{probaseuil}{A frequency threshold for selecting the set of relevant variables. Default is 0.5}
\item{var_nonselect}{Number of variables to not select. They have to be in the first columns. Default is 1, the selection is not performed on the intercept.}
}


\details{The Lasso is performed with the regularization parameter mu over the m bootstrap samples. An appearance frequency is obtained which shows the predictive power of each variable. 
}

\value{
\item{data}{The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{ind}{Set of selected variables for the regularization mu and the threshold probaseuil}
\item{frequency}{Appearance frequency of each variable; number of times each variables is selected over the m bootstrap iterations.}
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
y=rnorm(100)
mod=bolasso(x,y,mu=0.1,probaseuil=0.8)
}
}
