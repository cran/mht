\name{mht-package}
\alias{mht-package}
\docType{package}
\title{
Multiple hypotheses testing for variable selection
}
\description{
This package performs variable selection with multiple hypotheses testing, either for ordered variable selection or non-ordered variable selection.
Each is a sequential procedure: it starts to test the null hypothesis "none variable is relevant"; if this hypothesis is rejected, it then test "only the first variable is relevant", and so on until the null hypothesis is accepted.
}
\details{
\tabular{ll}{
Package: \tab mht\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-05-08\cr
License: \tab  GPL-3\cr
LazyLoad: \tab yes\cr
}
Two major functions: \code{\link{mht.order}} and \code{\link{mht}} (formerly proc_ord and procbol, it was changed to give more clarity and flexibility). The first gives the set of relevant variables for ordered variable selection, the last does the same for non-ordered variable selection.}
\author{
Florian Rohart \cr
Maintainer: florian.rohart@gmail.com
}
\references{
Multiple hypotheses testing for variable selection; F. Rohart 2011 \cr
Model-consistent sparse estimation through the bootstrap; F. Bach 2009 \cr
Adaptive tests of linear hypotheses by model selection; Baraud & al 2002\cr
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)
mod1=mht.order(x,y,ordre=5:1)
mod2=mht(x,y,alpha=c(0.1,0.05),maxordre=15)
}
}
