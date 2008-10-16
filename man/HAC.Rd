\name{HAC}
\alias{HAC}
\title{Covariance matrix of weakly dependent time series}
\description{
 Function to compute a consistent covariance matrix of the sample
 mean of a random vector of time series}
\usage{
HAC(x, weights = weightsAndrews2, prewhite = FALSE, ar.method = "ols",kernel=c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning")) 
}
\arguments{
 \item{x}{a \eqn{n\times q} matrix of time series, where n is the sample size.}
 \item{weights}{The method to compute the kernel weights. For now, weightsAndrews2 is the only one possible. I leave the option there because I am planning to give more choices in futur versions of the package.}
\item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}.}
\item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening.}
\item{kernel}{The choice of kernel}
}

\details{
\code{HAC} is simply a modified version of \code{meatHAC} from the package sandwich. The modifications have been made so that the argument x can be a matrix instead of an object of class lm or glm. The details on how is works can be found on the sandwich manual.
}

\value{A \eqn{q \times q} matrix containing an estimator of the asymptotic variance of 
\eqn{\sqrt{n} \bar{x}}, where \eqn{\bar{x}} is \eqn{q\times 1} vector with 
typical element \eqn{\bar{x}_i = \frac{1}{n}\sum_{j=1}^nx_{ji}}. This function is called by gmm() but can also be used by itself.}


\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.
}

\examples{
x <- arima.sim(n=200,list(order=c(1,0,1),ar=.5,ma=.3,sd=.5))
y <- .4*x+rnorm(200)
x <- cbind(x,y)
vcov <- HAC(x)
}


