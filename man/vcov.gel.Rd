\name{vcov.gel}
\alias{vcov.gel}
\title{Variance-covariance matrix of GEL}
\description{
 It extracts the matrix of variances and covariances from \code{gel} objects.
}
\usage{
\method{vcov}{gel}(object, lambda=FALSE, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} returned by the function \code{\link{gel}}}
\item{lambda}{If set to TRUE, the covariance matrix of the Lagrange multipliers is produced.}
\item{...}{Other arguments when \code{vcov} is applied to an other classe object}
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
t0 <- c(0,.5,.5)

res <- gel(g,x,t0)
vcov(res)
vcov(res,lambda=TRUE)
}

