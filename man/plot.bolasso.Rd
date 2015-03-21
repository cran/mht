\name{plot.bolasso}
\alias{plot.bolasso}
\docType{package}
\title{
graphical output for a \code{bolasso} object}
\description{
graphical output for a \code{bolasso} object. Plot the frequency of selection of each variable depending on the regularization parameter \code{mu} from the "bolasso" object.}
\usage{
\method{plot}{bolasso}(x,\dots)
}

\arguments{
\item{x}{Object of class "bolasso". As obtained from \code{\link{bolasso}}.}
\item{\dots}{not used.}
}

\details{Plot the frequency of selection of all variables depending on the regularization parameter.
}


\seealso{\code{\link{bolasso}}, \code{\link{dyadiqueordre}}}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

# BoLasso
mod=bolasso(x,y,mu=seq(1.5,0.1,-0.1))
plot(mod)
}
}

