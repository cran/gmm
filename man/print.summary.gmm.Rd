\name{print.summary.gmm}
\alias{print.summary.gmm}
\title{Printing the summary of gmm}
\description{
 It is a printing method for \code{summary.gmm} objects.
}
\usage{
\method{print}{summary.gmm}(x, digits=5, ...)
}
\arguments{
 \item{x}{An object of class \code{summary.gmm} returned by the function \code{\link{summary.gmm}}}
\item{digits}{The number of digits to be printed}
\item{...}{Other arguments when \code{print} is applied to an other classe object}
}

\value{
It prints some results from \code{\link{summary.gmm}}
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

print(summary(res))

}


