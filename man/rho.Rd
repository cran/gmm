\name{rho}

\alias{rho}

\title{Objective function of Generalized Empirical Likelihood (GEL)}

\description{
 It computes the objective function of GEL, its first and second analytical derivatives. It is called by \code{\ref{gel}.}
}
\usage{
rho(x,lamb,derive=0,type=c("EL","ET","CUE"),drop=TRUE)
}
\arguments{
\item{x}{A \eqn{n\times q} matrix with typical element \eqn{(t,i)}, \eqn{g_i(\theta,x_t)}}

\item{lamb}{A \eqn{q \times 1} vector of lagrange multipliers}

\item{derive}{0 for the objective function, 1 for the first derivative with respect to \eqn{\lambda} and 2 for the second derivative with respect to \eqn{\lambda}.}

\item{type}{"EL" for empirical likelihood, "ET" for exponential tilting and "CUE" for continuous updated estimator.}

\item{drop}{Because the solution may not be in the domain of \eqn{\rho(v)} \eqn{\forall t} in small sample, we can drop those observations to avoid the return of NaN}
}
\details{
The objective function is \eqn{\frac{1}{n}\sum_{t=1}^n \rho(<g(\theta,x_t),\lambda>)}, where \eqn{\rho(v)=\log{(1-v)}} for empirical likelihood, \eqn{-e^v} for exponential tilting and \eqn{(-v-0.5v^2)} for continuous updated estimator. 
}

\value{
'rho' returns a scalar if "derive=0", a \eqn{q\time 1} vector if "derive" = 1 and a \eqn{q\times q} matrix if derive = 2.
 }

\references{
Newey, W.K. and Smith, R.J. (2004), Higher Order Properties of GMM and 
Generalized Empirical Likelihood Estimators. \emph{Econometrica}, \bold{72}, 219-255.

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
	f <- cbind(u,u*x[4:(n-3)],u*x[3:(n-4)],u*x[2:(n-5)],u*x[1:(n-6)])
	return(f)
	}
n = 500
phi<-c(.2,.7)
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
gt <- g(thet,x)
lamb <- rep(0,5)
rho(gt,lamb,type='EL',derive=0)
}


