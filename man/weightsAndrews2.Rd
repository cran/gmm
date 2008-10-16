\name{weightsAndrews2}
\alias{weightsAndrews2}
\alias{bwAndrews2}
\alias{kweights}
\title{Kernel weights}
\description{
 Function to compute the kernel weights used to compute the HAC covariance matrix}
\usage{
weightsAndrews2(x, bw = bwAndrews2, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol = 1e-7, verbose = FALSE)

bwAndrews2(x, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols")  
}
\arguments{
  \item{x}{A \eqn{n\times q} matrix of time series from which we want to compute the covariance matrix.}
 \item{bw}{The method to compute the bandwith parameter. For now, bwAndrews2 is the only one possible. I leave the option there because I am planning to give more choices in futur versions of the package.}
 \item{prewhite}{logical or integer. Should the estimating functions
    be prewhitened? If \code{TRUE} or greater than 0 a VAR model of
    order \code{as.integer(prewhite)} is fitted via \code{ar} with
    method \code{"ols"} and \code{demean = FALSE}.}
\item{ar.method}{character. The \code{method} argument passed to
   \code{\link{ar}} for prewhitening.}
 \item{verbose}{logical. Should the bandwidth parameter used be
    printed?}
\item{approx}{a character specifying the approximation method if the
    bandwidth has to be chosen by \code{bwAndrews2}.}
\item{tol}{numeric. Weights that exceed \code{tol} are used for computing
   the covariance matrix, all other weights are treated as 0.}
\item{kernel}{The choice of kernel}
}

\details{\code{weightsAndrews2} and \code{bwAndrews2} are simply modified version of \code{weightsAndrews} and \code{bwAndrews} from the package sandwich. The modifications have been made so that the argument x can be a matrix instead of an object of class lm or glm. The details on how is works can be found on the sandwich manual. \code{kweights} is the same as the one included in the package sandwich.
}

\value{
\code{weightsAndrews} returns a vector of weights.

\code{bwAndrews} returns the selected bandwidth parameter.
}

\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.
}

\examples{
x <- arima.sim(n=200,list(ordre=c(1,0,1),ar=.9,ma=.4))
y <- 2*x + rnorm(200)
x = cbind(x,y)
w <- weightsAndrews2(x, bw = bwAndrews2, kernel = "Quadratic")
plot(w,type='l')
}


