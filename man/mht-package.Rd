\name{mht-package}
\alias{mht-package}
\docType{package}
\title{
Multiple hypothesis testing for variable selection
}
\description{
Multiple hypothesis testing for variable selection in high dimensional linear models.\cr
This package performs variable selection with multiple hypothesis testing, either for ordered variable selection or non-ordered variable selection. In both cases, a sequential procedure is performed. It starts to test the null hypothesis "no variable is relevant"; if this hypothesis is rejected, it then tests "only the first variable is relevant", and so on until the null hypothesis is accepted. \cr
    More details are available in the paper `Multiple hypothesis testing for variable selection', Rohart F. (2011).
}
\details{
\tabular{ll}{
Package: \tab mht\cr
Type: \tab Package\cr
Version: \tab 3.1.2\cr
License: \tab  GPL-3\cr
date: 20-03-2015\cr
}
Two major functions: \code{\link{mht.order}} and \code{\link{mht}} (proc_ord and procbol in version <3.00, it was changed to give more clarity and flexibility). The first estimates the set of relevant variables for ordered variable selection, e.g. if an apriori of the importance of the variables is known; the last does the same for non-ordered variable selection.}
\author{
Florian Rohart \cr
Maintainer: florian.rohart@gmail.com
}
\references{
Multiple hypothesis testing for variable selection; F. Rohart 2011 \cr
Model-consistent sparse estimation through the bootstrap; F. Bach 2009 \cr
Adaptive tests of linear hypotheses by model selection; Baraud & al 2002\cr
}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

# mht.order
mod.order=mht.order(x,y,ordre=5:1)
mod.order

# mht
mod=mht(x,y,alpha=c(0.1,0.05),maxordre=15)
mod
}
}
