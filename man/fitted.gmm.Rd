\name{fitted.gmm}
\alias{fitted.gmm}
\title{Fitted values of GMM}
\description{
 Method to extract the fitted values of the model estimated by \code{\link{gmm}}.
}
\usage{
\method{fitted}{gmm}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{gmm} returned by the function \code{\link{gmm}}}
\item{...}{Other arguments when \code{fitted} is applied to an other classe object}
}

\value{
It returns a matrix of the estimated mean \eqn{\hat{y}} in \code{g=y~x} as it is done by \code{fitted.lm}.
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
plot(x,y, main="Fitted model with GMM")
lines(x,fitted(res),col=2,)
legend("topright",c("Y","Yhat"),col=1:2,lty=c(1,1))
}

