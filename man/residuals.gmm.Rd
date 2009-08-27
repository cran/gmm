\name{residuals.gmm}
\alias{residuals.gmm}
\title{Residuals of GMM}
\description{
 Method to extract the residuals of the model estimated by \code{gmm}.
}
\usage{
\method{residuals}{gmm}(object, ...)
}

\arguments{
 \item{object}{An object of class \code{gmm} returned by the function \code{\link{gmm}}}
\item{...}{Other arguments when \code{residuals} is applied to an other classe object}
}

\value{
It returns the  matrix of residuals \eqn{(y-\hat{y})} in \code{g=y~x} as it is done by  \code{residuals.lm}.
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
plot(x,residuals(res), main="Residuals of an estimated model with GMM")
}

