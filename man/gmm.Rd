\name{gmm}

\alias{gmm}

\title{Generalized method of moment estimation}

\description{
 Function to estimate a vector of parameters based on moment conditions using the GMM method of Hansen(82). 
}
\usage{
gmm(g,x,t0=NULL,gradv=NULL, type=c("twoStep","cue","iterative"), 
wmatrix = c("optimal","ident"),  vcov=c("HAC","iid"), 
kernel=c("Quadratic Spectral","Truncated", "Bartlett", 
"Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews2, 
prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7, 
itermax=100,intercept=TRUE,optfct=c("optim","optimize"), ...)
}
\arguments{
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build the q sample moment conditions. It can also be a formula if the model is linear (see details below).}

\item{x}{The matrix or vector of data from which the function \eqn{g(\theta,x)} is computed. If "g" is a formula, it is an \eqn{n \times Nh} matrix of instruments (see details below).}

\item{t0}{A \eqn{k \times 1} vector of starting values. It is required only when "g" is a function because only then a numerical algorithm is used to minimize the objective function. If the dimension of \eqn{\theta} is one, see the argument "optfct".}

\item{gradv}{A function of the form \eqn{G(\theta,x)} which returns a \eqn{q\times k} matrix of derivatives of \eqn{\bar{g}(\theta)} with respect to \eqn{\theta}. By default, the numerical algorithm \code{numericDeriv} is used. It is of course strongly suggested to provide this function when it is possible. This gradiant is used compute the asymptotic covariance matrix of \eqn{\hat{\theta}}. If "g" is a formula, the gradiant is not required (see the details below).}

\item{type}{The GMM method: "twostep" is the two step GMM proposed by Hansen(1982) and the "cue" and "iterative" are respectively the continuous updated and the iterative GMM proposed by Hansen, Eaton et Yaron (1996)}

\item{wmatrix}{Which weighting matrix should be used in the objective function. By default, it is the inverse of the covariance matrix of \eqn{g(\theta,x)}. The other choice is the identity matrix which is usually used to obtain a first step estimate of \eqn{\theta} }

\item{vcov}{Assumption on the properties of the random vector x. By default, x is a weakly dependant process. The "iid" option will avoid using the HAC matrix which will accelerate the estimation if one is ready to make that assumption.}

\item{kernel}{type of kernel used to compute the covariance matrix of the vector of sample moment conditions (see \code{\link{HAC}} for more details)}

\item{crit}{The stoping rule for the iterative GMM. It can be reduce to increase the precision.}

\item{bw}{The method to compute the bandwidth parameter. By default it is \code{\link{bwAndrews2}} which is proposed by Andrews (1991). The alternative is \code{\link{bwNeweyWest2}} of Newey-West(1994).}

\item{prewhite}{logical or integer. Should the estimating functions be prewhitened? If \code{TRUE} or greater than 0 a VAR model of order \code{as.integer(prewhite)} is fitted via \code{ar} with method \code{"ols"} and \code{demean = FALSE}.}

\item{ar.method}{character. The \code{method} argument passed to \code{\link{ar}} for prewhitening.}

\item{approx}{A character specifying the approximation method if the bandwidth has to be chosen by \code{bwAndrews2}.}

\item{tol}{Weights that exceed \code{tol} are used for computing the covariance matrix, all other weights are treated as 0.}

\item{itermax}{The maximum number of iterations for the iterative GMM. It is unlikely that the algorithm does not converge but we keep it as a safety.}

\item{intercept}{If "g" is a formula, should the model include a constant? It should always be the case but the choice is yours.}

\item{optfct}{Only when the dimension of \eqn{\theta} is 1, you can choose between the algorithm \code{\link{optim}} or \code{\link{optimize}}. In that case, the former is unreliable. If \code{\link{optimize}} is chosen, "t0" must be \eqn{1\times 2} which represents the interval in which the algorithm seeks the solution.}


\item{...}{More options to give to \code{\link{optim}}.}
}

\details{
\code{weightsAndrews2} and \code{bwAndrews2} are simply modified version of \code{weightsAndrews} and \code{bwAndrews} from the package sandwich. The modifications have been made so that the argument x can be a matrix instead of an object of class lm or glm. The details on how is works can be found on the sandwich manual.

If we want to estimate a model like \eqn{Y_t = \theta_1 + X_{2t}\theta_2 + \cdots + X_{k}\theta_k + \epsilon_t} using the moment conditions \eqn{Cov(\epsilon_tH_t)=0}, where \eqn{H_t} is a vector of \eqn{Nh} instruments, than we can define "g" like we do for \code{\link{lm}}. We would have \code{g = y~x2+x3+\cdots+xk} and the argument "x" above would become the matrix H of instruments. As for \code{\link{lm}}, \eqn{Y_t} can be a \eqn{Ny \times 1} vector which would imply that \eqn{k=Nh \times Ny}. The intercept is included by default so you do not have to add a column of ones to the matrix \eqn{H}. You do not need to provide the gradiant in that case since in that case it is embedded in \code{\link{gmm}}.
}

\value{
'gmm' returns an object of 'class' '"gmm"' 

The functions 'summary' is used to obtain and print a summary of the results.  It also compute the J-test of overidentying restriction

The object of class "gmm" is a list containing:

par: \eqn{k\times 1} vector of parameters

vcov: the covariance matrix of the parameter

objective: the value of the objective function \eqn{\| var(\bar{g})^{-1/2}\bar{g}\|^2}
 }


\references{
  Zeileis A (2006), Object-oriented Computation of Sandwich Estimators.
  \emph{Journal of Statistical Software}, \bold{16}(9), 1--16.
  URL \url{http://www.jstatsoft.org/v16/i09/}.

 Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.

 Newey WK & West KD (1987), A Simple, Positive Semi-Definite,
 Heteroskedasticity and Autocorrelation Consistent Covariance
 Matrix. \emph{Econometrica}, \bold{55}, 703--708.

 Newey WK & West KD (1994), Automatic Lag Selection in Covariance
 Matrix Estimation. \emph{Review of Economic Studies}, \bold{61}, 631-653.

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
# First, an exemple with the fonction g()

g <- function(tet,x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u,u*x[4:(n-3)],u*x[3:(n-4)],u*x[2:(n-5)],u*x[1:(n-6)])
	return(f)
	}

Dg <- function(tet,x)
	{
	n <- nrow(x)
	xx <- cbind(rep(1,(n-6)),x[6:(n-1)],x[5:(n-2)])
        H  <- cbind(rep(1,(n-6)),x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
	f <- -crossprod(H,xx)/(n-6)
	return(f)
	}
n = 500
set.seed(123)
phi<-c(.2,.7)
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)

res_2s <- gmm(g,x,c(0,.3,.6),gradv=Dg)
summary(res_2s)

res_iter <- gmm(g,x,c(0,.3,.6),gradv=Dg,type="iterative")
summary(res_iter)

# The same model but with g as a formula....  much simpler in that case

y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
g <- y~ym1+ym2
x <- H

res_2s <- gmm(g,x)
summary(res_2s)

res_iter <- gmm(g,x,type="iterative")
summary(res_iter)

}

