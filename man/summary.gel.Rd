\name{summary.gel}
\alias{summary.gel}
\title{Method for object of class gel}
\description{
 It presents the results from the \code{gel} estimation in the same fashion as \code{summary} does for the \code{lm} class objects. It also compute the J-test, LM and LR tests for overidentifying restriction.
}
\usage{
\method{summary}{gel}(object, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} returned by the function \code{\link{gel}}}
\item{...}{Other arguments when summary is applied to an other classe object}
}

\value{
It returns a list with the parameter estimates and their standard deviations, t-stats and p-values. It also returns the three tests (J, LM and LR) and p-values for the null hypothesis that \eqn{E(g(\theta,X)=0} and several convergence codes. 
}

\references{

Anatolyev, S. (2005), GMM, GEL, Serial Correlation, and Asymptotic Bias. \emph{Econometrica}, \bold{73}, 983-1002.

Kitamura, Yuichi (1997), Empirical Likelihood Methods With Weakly Dependent Processes.
\emph{The Annals of Statistics}, \bold{25}, 2084-2102.

Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.


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
summary(res)

}

