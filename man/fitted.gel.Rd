\name{fitted.gel}
\alias{fitted.gel}
\title{Fitted values of GEL}
\description{
 Method to extract the fitted values of the model estimated by \code{\link{gel}}.
}
\usage{
\method{fitted}{gel}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} returned by the function \code{\link{gel}}}
\item{...}{Other arguments when \code{fitted} is applied to an other classe object}
}

\value{
It returns a matrix of the estimated mean \eqn{\hat{y}} in \code{g=y~x} as it is done by \code{fitted.lm}.
}

\examples{

# GEL can deal with endogeneity problems

n = 200
phi<-c(.2,.7)
thet <- 0.2
sd <- .2
set.seed(123)
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)

y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]
H <- cbind(x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
g <- y~ym1+ym2
x <- H

res <- gel(g,x,c(0,.3,.6))
plot(y, main="Fitted ARMA with GEL")
lines(fitted(res),col=2)
}

