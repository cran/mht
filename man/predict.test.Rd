\name{predict.proctest}
\alias{predict.proctest}
\title{
Predict a \code{proctest} object}
\description{
Predict a \code{proctest} object for new data newx}
\usage{
\method{predict}{proctest}(object,newx,\dots)
}

\arguments{
\item{object}{Object of class "proctest".}
\item{newx}{Data matrix of size n*p.}
\item{\dots}{not used.}
}


\details{see \code{\link{procbol}} for details.
}

\value{
{Matrix of predicted values, each columns concerns a particular type I error alpha -from the initial object-.}
}
