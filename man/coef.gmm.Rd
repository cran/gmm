\name{coef.gmm}
\alias{coef.gmm}
\title{Coefficients of GMM}
\description{
 It extracts the coefficients from \code{gmm} objects.
}
\usage{
\method{coef}{gmm}(object, ...)
}
\arguments{
 \item{object}{An object of class gmm returned by the function \code{\link{gmm}}}
\item{...}{Other arguments when \code{coef} is applied to another classe object}
}

\value{
Vector of coefficients
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
coef(res)
}

