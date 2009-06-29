\name{summary.gmm}
\alias{summary.gmm}
\title{Method for object of class gmm}
\description{
 It presents the results from the \code{gmm} estimation in the same fashion as summary does for the lm class objects for example. It also compute the J-test for overidentifying restriction.
}
\usage{
\method{summary}{gmm}(object, interval=FALSE, ...)
}
\arguments{
 \item{object}{An object of class gmm returned by the function \code{gmm}}
\item{interval}{Should the results include the confidence intervals of \eqn{\hat{\theta}}. If so, "interval" should be eqal to the confidence level.}
\item{...}{Other arguments when summary is applied to an other classe object}
}

\value{
It returns a list with the parameter estimates and theirs standard deviations, t-stat and p-values. It also returns the J-test and p-value for the null hypothesis that \eqn{E(g(\theta,X)=0} 
}

\references{
  Hansen, L.P. (1982),
  Large Sample Properties of Generalized Method of Moments Estimators.
  \emph{Econometrica}, \bold{50},
  1029-1054,

  Hansen, L.P. and Heaton, J. and Yaron, A.(1996),
  Finit-Sample Properties of Some Alternative GMM Estimators.
  \emph{Journal of Business and Economic Statistics}, \bold{14}
  262-280.
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

summary(res)
summary(res,interval=0.95)

}

