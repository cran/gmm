\name{getLamb}

\alias{getLamb}

\title{Solving for the Lagrange multipliers of Generalized Empirical Likelihood (GEL) }

\description{
 It computes the vector of Lagrange multipliers, which maximizes the GEL objective function, using an iterative Newton method.
}
\usage{
getLamb(g, tet, x, type = c('EL','ET','CUE'), tol_lam = 1e-12, maxiterlam = 1000, tol_obj = 1e-7, k = 1)
}
\arguments{
\item{g}{A function of the form \eqn{g(\theta,x)} and which returns a \eqn{n \times q} matrix with typical element \eqn{g_i(\theta,x_t)} for \eqn{i=1,...q} and \eqn{t=1,...,n}. This matrix is then used to build the q sample moment conditions.}

\item{tet}{A \eqn{k \times 1} vector of parameters at which the function \eqn{g(\theta,x)} has to be evaluated}

\item{x}{The matrix or vector of data from which the function \eqn{g(\theta,x)} is computed.}

\item{type}{"EL" for empirical likelihood, "ET" for exponential tilting and "CUE" for continuous updated estimator.}

\item{tol_lam}{Tolerance for \eqn{\lambda} between two iterations. The algorithm stops when \eqn{\|\lambda_i -\lambda_{i-1}\|} reaches \code{tol_lam} }

\item{maxiterlam}{The algorithm stops if there is no convergence after "maxiterlam" iterations.}

\item{tol_obj}{Tolerance for the gradiant of the objective function. The algorithm returns a non-convergence message if \eqn{\max(|gradiant|)} does not reach \code{tol_obj}. It helps the \code{gel} algorithm to select the right space to look for \eqn{\theta}}

\item{k}{It represents the ratio k1/k2, where \eqn{k1=\int_{-\infty}^{\infty} k(s)ds} and \eqn{k2=\int_{-\infty}^{\infty} k(s)^2 ds}. See Smith(2004).}
}
\details{
It solves the problem \eqn{\frac{1}{n}\sum_{t=1}^n D\rho(<g(\theta,x_t),\lambda>)g(\theta,x_t)=0}.
}

\value{
lambda: A \eqn{q\times 1} vector of Lagrange multipliers which solve the system of equations given above.
singular: 0 for a normal solution, 1 if the algorithm does not converge and 2 if the algorithm produces a singular system, NaN or Inf values. 
\code{conv_mes}: A message with details about the convergence.
 }

\references{
Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

Smith, R.J. (2004), GEL Criteria for Moment Condition Models. \emph{Working paper, CEMMAP}.

}

\examples{
g <- function(tet,x)
	{
	n <- nrow(x)
	u <- (x[7:n] - tet[1] - tet[2]*x[6:(n-1)] - tet[3]*x[5:(n-2)])
	f <- cbind(u, u*x[4:(n-3)], u*x[3:(n-4)], u*x[2:(n-5)], u*x[1:(n-6)])
	return(f)
	}
n = 500
phi<-c(.2, .7)
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n = n, list(order = c(2, 0, 1), ar = phi, ma = thet, sd = sd)), ncol = 1)
getLamb(g, c(0, phi), x, type = "EL")
}


