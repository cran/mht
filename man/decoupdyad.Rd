\name{dyadiqueordre}
\alias{dyadiqueordre}
\docType{package}
\title{Dyadic algorithm to order variables
}
\description{
Dyadic algorithm using the Bolasso technique to order the variables }
\usage{
dyadiqueordre(data,Y,m,maxordre,var_nonselect,showtest,showordre,random)
}

\arguments{
\item{data}{Input matrix of dimension n * p; each of the n rows is an observation vector of p variables. The intercept should be included in the first column as (1,...,1). If not, it is added.}
\item{Y}{Response variable of length n.}
\item{m}{Number of bootstrap iteration of the Lasso. Default is \code{m}=100.}
\item{maxordre}{Number of variables to order. Default is min(n/2-1,p/2-1).}
\item{var_nonselect}{Number of variables that don't undergo feature selection. They have to be in the first columns of data. Default is 1, the selection is not performed on the intercept.}
\item{showtest}{Logical value. If TRUE, show the number of regularization parameters tested to show the steps of the algorithm. Default is FALSE.}
\item{showordre}{Logical value. If TRUE, shows the ordered variables at each step of the algorithm. Default is TRUE.}
\item{random}{optionnal parameter. Matrix of size n*m, the m bootstrap samples are constructed from the m columns.}
}


\details{The algorithm starts from a large regularization parameter given by one run of Lasso. It proceeds by dyadic splitting until one variable is isolated; e.g one variable alone achieve a frequency of 1; it is the first ordered variable. And so on until \code{maxordre} variables are ordered.
}

\value{A 'bolasso' object is returned for which the method \code{plot} is available.

\item{data}{A list containing:
\itemize{\item{X - The scaled matrix used in the algorithm, the first column being (1,...,1).}
\item{Y - the input response vector}
\item{means.X - Vector of means of the input data matrix.}
\item{sigma.X - Vector of variances of the input data matrix.}
}}

\item{ordre}{The order obtained on the variables.}
\item{mu}{Vector of the positive regularization sequence that was used in the algorithm.}
\item{frequency}{Matrix of p rows. Appearance frequency of each variable for the regularization parameter in \code{mu}.}
}

\seealso{\code{\link{bolasso}}, \code{\link{plot.bolasso}}}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod=dyadiqueordre(x,y,maxordre=15)
}
}
