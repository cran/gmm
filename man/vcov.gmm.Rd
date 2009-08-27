\name{vcov.gmm}
\alias{vcov.gmm}
\title{Variance-covariance matrix of GMM}
\description{
 It extracts the matrix of variances and covariances from \code{gmm} objects.
}
\usage{
\method{vcov}{gmm}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{gmm} returned by the function \code{\link{gmm}}}
\item{...}{Other arguments when \code{vcov} is applied to another classe object}
}

\value{
A matrix of variances and covariances
}


\examples{
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
g <- y~ym1+ym2
x <- H

res <- gmm(g,x)
vcov(res)
}

