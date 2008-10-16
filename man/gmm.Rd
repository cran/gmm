\name{gmm}

\alias{gmm}

\title{Generalized method of moment estimation}

\description{
 Function to estimate a vector of parameters based on moment conditions using the GMM method of Hansen(82). It is a preliminary version. It has been tested with monte carlo experiment but it may have bugs. Comments are welcome.
}
\usage{
gmm(g,t0,x,grad=NULL,type=c("twoStep","cue","iterative"),kernel=c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),iid=FALSE,crit=10e-7,
     itermax=100,algo=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), vcov=HAC, ...)
}

\arguments{
 \item{g}{A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix of time series. This series is than used to build the q sample moment conditions. }
 \item{t0}{A \eqn{k \times 1} vector of starting values}
\item{x}{The matrix or vector of data from which the function \eqn{g(\theta,x)} is computed}
\item{grad}{A function of the form \eqn{G(\theta,x)} which returns a \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)}. By default, the numerical algorithm \code{numericDeriv} is used.}
\item{type}{The GMM method: "twostep" is the two step GMM proposed by Hansen(1982) and the cue and iterative are respectively the continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron (1996)}
\item{kernel}{type of kernel used to compute the covariance matrix of the vector of moment conditions. See \code{\link{HAC}} for more details}
\item{iid}{Hypothesis on the properties of x. By default, x is a weakly dependant time series}
\item{crit}{The level of precision required for the iterative GMM}
\item{itermax}{The maximum number of iterations for the iterative GMM}
\item{algo}{The numerical algorithm for the optimization problem. See \code{\link{optim}} for more details.}
\item{vcov}{The method used to compute de covariance matrix. By default it is the \code{\link{HAC}}. For now it is the only option. But in a futur version, more choices will be available.}
\item{...}{More option to give to the \code{\link{HAC}} procedure.}
}

\details{
\code{weightsAndrews2} and \code{bwAndrews2} are simply modified version of \code{weightsAndrews} and \code{bwAndrews} from the package sandwich. The modifications have been made so that the argument x can be a matrix instead of an object of class lm or glm. The details on how is works can be found on the sandwich manual.
}

\value{
'gmm' returns an object of 'class' '"gmm"' 

     The functions 'summary' is used to obtain and print a
     summary of the results.  It also compute the J-test of overidentying
     restriction
The object of class "gmm" is a list containing:

par: \eqn{k\times 1} vector of parameters

vcov: the covariance matrix

objective: the value of the objective function \eqn{\bar{g}'var(\bar{g})^{-1}\bar{g}}
 }


\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

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
g <- function(tet,x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u,u*x[5:(n-2)],u*x[4:(n-3)],u*x[3:(n-4)])
	return(f)
	}

Dg <- function(tet,x)
	{
	n <- nrow(x)
	xx <- cbind(rep(1,(n-6)),x[6:(n-1)],x[5:(n-2)])
        H  <- cbind(rep(1,(n-6)),x[5:(n-2)],x[4:(n-3)],x[3:(n-4)])
	f <- -crossprod(H,xx)/(n-6)
	return(f)
	}
n = 500
phi<-c(.2,.7)
thet <- 0
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
resgmm <- gmm(g,c(0,.3,.6),x,grad=Dg)
}

