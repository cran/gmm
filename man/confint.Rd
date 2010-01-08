\name{confint}
\alias{confint.gel}
\alias{confint.gmm}
\title{Confidence intervals for GMM or GEL}
\description{
It produces confidence intervals for the coefficients from \code{gel} or \code{gmm} estimation.
}
\usage{
\method{confint}{gel}(object, parm, level = 0.95, lambda = FALSE, ...)
\method{confint}{gmm}(object, parm, level = 0.95, ...)
}
\arguments{
 \item{object}{An object of class \code{gel} or \code{gmm} returned by the function \code{\link{gel}} or \code{\link{gmm}}}
\item{parm}{A specification of which parameters are to be given confidence intervals, either a vector of numbers or a vector
          of names.  If missing, all parameters are considered.}
\item{level}{The confidence level}
\item{lambda}{If set to TRUE, the confidence intervals for the Lagrange multipliers are produced.}
\item{...}{Other arguments when \code{confint} is applied to another classe object}
}

\value{
It returns a matrix with the first column being the lower bound and the second the upper bound.} 


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
#################
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2,0,1), ar = phi, ma = thet, sd = sd)), ncol = 1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)], x[3:(n-4)], x[2:(n-5)], x[1:(n-6)])
g <- y ~ ym1 + ym2
x <- H
t0 <- c(0,.5,.5)

res <- gel(g, x, t0)

confint(res)
confint(res, level = 0.90)
confint(res, lambda = TRUE)

########################

res <- gmm(g, x)

confint(res)
confint(res, level = 0.90)

}

