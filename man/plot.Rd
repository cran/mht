\name{plot}
\alias{plot.mht.order}
\alias{plot.mht}
\docType{package}
\title{
graphical output for a \code{mht} or \code{mht.order} object}
\description{
Graphical output for a \code{mht} or \code{mht.order} object. Four plots (selectable by which.plot) are currently available: a plot of the fitted values against the true values, a plot of the residuals against the fitted values, a Normal Q-Q plot, and a barplot showing the coefficients.
}
\usage{
\method{plot}{mht}(x,which.plot=1:4,\dots)
\method{plot}{mht.order}(x,which.plot=1:4,\dots)
}

\arguments{
\item{x}{Object of class "mht" as obtained from \code{\link{mht}} or "mht.order" as obtained from \code{\link{mht.order}}.}
\item{which.plot}{if a subset of the plots is required, specify a subset of the numbers 1:4.}
\item{\dots}{not used.}
}


\details{ The four plots are shown for each level \code{alpha} of the \code{mht} or \code{mht.order} object.\cr

For the Normal Q-Q plot, the standardized residuals are used. They are calculated as R[i]/sqrt(var(R)).
}


\seealso{\code{\link{mht}}, \code{\link{mht.order}}}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

# mht
mod=mht(x,y,alpha=c(0.1,0.05),maxordre=15)
plot(mod)

# mht.order
mod.order=mht.order(x,y,ordre=5:1,alpha=c(0.1,0.05))
plot(mod.order)
}
}
