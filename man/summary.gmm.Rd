\name{summary.gmm}
\alias{summary.gmm}
\title{Method for object of class gmm}
\description{
 It presents the results from the \code{gmm} estimation in the same fashion as summary does for the lm class objects for example. It also compute the J-test for overidentifying restriction.
}
\usage{
\method{summary}{gmm}(object, ...)
}
\arguments{
 \item{object}{An object of class gmm returned by the function \code{gmm}}
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

g <- function(tet,x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u,u*x[5:(n-2)],u*x[4:(n-3)],u*x[3:(n-4)])
	#f <- cbind(u,u*x[5:(n-2)],u*x[4:(n-3)],u*x[3:(n-4)],u*x[2:(n-5)],u*x[1:(n-6)])
	#f <- cbind(u,u*x[4:(n-3)],u*x[3:(n-4)],u*x[2:(n-5)],u*x[1:(n-6)])
	return(f)
	}

Dg <- function(tet,x)
	{
	n <- nrow(x)
	xx <- cbind(rep(1,(n-6)),x[6:(n-1)],x[5:(n-2)])
	 H  <- cbind(rep(1,(n-6)),x[5:(n-2)],x[4:(n-3)],x[3:(n-4)])
	# H  <- cbind(rep(1,(n-6)),x[5:(n-2)],x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
	#H  <- cbind(rep(1,(n-5)),x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
	f <- -crossprod(H,xx)/(n-6)
	return(f)
	}

resgmm <- gmm(g,c(0,.3,.6),x,grad=Dg)
summary(resgmm)

}

