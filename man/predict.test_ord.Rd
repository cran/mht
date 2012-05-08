\name{predict.proctest_ord}
\alias{predict.proctest_ord}
\title{
Predict a \code{proctest_ord} object}
\description{
Predict a \code{proctest_ord} object for new data newx}
\usage{
\method{predict}{proctest_ord}(object,newx,\dots)
}

\arguments{
\item{object}{Object of class "proctest_ord".}
\item{newx}{Data matrix of size n*p.}
\item{\dots}{not used.}
}


\details{see \code{\link{proc_ord}} for details.
}

\value{
{Matrix of predicted values, each columns concerns a particular type I error alpha -from the initial object-.}
}
