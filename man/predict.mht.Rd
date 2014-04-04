\name{predict.mht}
\alias{predict.mht}
\title{
Predict a \code{mht} object}
\description{
Predict a \code{mht} object for new data newx}
\usage{
\method{predict}{mht}(object,newx,\dots)
}

\arguments{
\item{object}{Object of class "mht".}
\item{newx}{Data matrix of size n*p.}
\item{\dots}{not used.}
}


\details{see \code{\link{mht}} for details.
}

\value{
{Matrix of predicted values, each columns concerns a particular type I error alpha -from the initial object-.}
}
