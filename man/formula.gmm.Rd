\name{formula.gmm}
\alias{formula.gmm}
\title{Formula method for gmm objects}
\description{
 Method to extract the formula from \code{gmm} objects produced by \code{\link{gmm}}.
}
\usage{
\method{formula}{gmm}(x, ...)
}
\arguments{
 \item{x}{An object of class \code{gmm} returned by the function \code{\link{gmm}}}
\item{...}{Other arguments to pass to other methods}
}

\examples{

# GMM is like GLS for linear models without endogeneity problems

set.seed(345)
n = 200
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
y <- 10+5*rnorm(n) + x

res <- gmm(y~x,x)
formula(res)
}

