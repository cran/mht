\name{predict}
\alias{predict.mht}
\alias{predict.mht.order}
\title{
Predict a \code{mht} or \code{mht.order} object}
\description{
Predict a \code{mht} or \code{mht.order} object for new data newx}
\usage{
\method{predict}{mht}(object,newx,level,\dots)
\method{predict}{mht.order}(object,newx,level,\dots)
}

\arguments{
\item{object}{Object of class "mht.order" as obtained from \code{\link{mht.order}} or "mht" as obtained from \code{\link{mht}}.}
\item{newx}{Data matrix of size n*p.}
\item{level}{Level of the prediction interval. Default is 0.95.}
\item{\dots}{not used.}
}


\details{The prediction is available for each level \code{alpha} of the object. The prediction values and the prediction interval are derived from the \code{\link{predict.lm}} function. If \code{newx} is missing, the fitted values of the object are returned.
}

\value{
{Array of predicted values and prediction interval. The third dimension is relative to the type I error alpha -from the initial object-. For each alpha, a matrix with column names fit, lwr, and upr (from \code{\link{predict.lm}})}
}

\seealso{\code{\link{mht}}, \code{\link{mht.order}}, \code{\link{predict.lm}}}

\examples{
\dontrun{
x=matrix(rnorm(100*20),100,20)
beta=c(rep(2,5),rep(0,15))
y=x\%*\%beta+rnorm(100)

# mht
mod=mht(x,y,alpha=c(0.1,0.05),maxordre=15)

# predict without new data: gives the fitted values
pred=predict(mod)

# predict with new data
pred=predict(mod,newx=matrix(rnorm(40*20),40,20))
pred
}
}

