\name{bolasso}
\alias{bolasso}
\docType{package}
\title{
Bolasso: Bootstrapped Lasso}
\description{
Performed a bootstrapped Lasso}
\usage{
bolasso(data,Y,mu,m,probaseuil,penalty.factor,random)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each row is an observation vector. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{mu}{Positive regularization sequence to be used for the Lasso.}
\item{m}{Number of bootstrapped iteration of the Lasso. Default is m=100.}
\item{probaseuil}{A frequency threshold for selecting the set of relevant variables. Default is 1}
\item{penalty.factor}{Separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables except the intercept.}
\item{random}{optionnal parameter. Matrix of size n*m, the m bootstrap samples are constructed from the m columns.}
}


\details{The Lasso is performed with the regularization parameter mu over m bootstrapped samples. An appearance frequency is obtained which shows the predictive power of each variable. 
}

\value{
\item{data}{The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{ind}{Set of selected variables for the regularization \code{mu} and the threshold \code{probaseuil}.}
\item{frequency}{Appearance frequency of each variable; number of times each variables is selected over the m bootstrapped iterations.}
}

\references{
Model-consistent sparse estimation through the bootstrap; F. Bach 2009 \cr}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(1,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod=bolasso(x,y,mu=seq(1.5,0.1,-0.1))}
}
