\name{predict.mht.order}
\alias{predict.mht.order}
\title{
Predict a \code{mht.order} object}
\description{
Predict a \code{mht.order} object for new data newx}
\usage{
\method{predict}{mht.order}(object,newx,\dots)
}

\arguments{
\item{object}{Object of class "mht.order".}
\item{newx}{Data matrix of size n*p.}
\item{\dots}{not used.}
}


\details{see \code{\link{mht.order}} for details.
}

\value{
{Matrix of predicted values, each columns concerns a particular type I error alpha -from the initial object-.}
}
