\name{mht-package}
\alias{mht-package}
\docType{package}
\title{
Multiple hypotheses testing for variable selection
}
\description{
This package performed variable selection with multiple hypotheses testing, either for ordered variables or non ordered variables.
For this doing, calculation of the quantiles of the unknown statistics is done. It is a sequential procedure: it starts to test the null hypothesis: none variable is relevant, if this hypothesis is rejected, it then test "only the first variable is relevant", and so on until the null hypothesis is accepted.
}
\details{
\tabular{ll}{
Package: \tab mht\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-03-15\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
Two major functions: proc_ord and procbol. The first gives the set of relevant variables when the variables are pre-ordered, the last does the same when the variables are not pre-ordered.}
\author{
Florian Rohart \cr
Maintainer: florian.rohart@gmail.com
}
\references{
Multiple hypotheses testing for variable selection; f.Rohart \cr
Model-consistent sparse estimation through the bootstrap; F.Bach \cr
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod1=proc_ord(x,y,ordre=5:1)
mod2=procbol(x,y,alpha=c(0.1,0.05),maxordre=15)
}
}
