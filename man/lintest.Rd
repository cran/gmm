\name{lintest}
\alias{lintest}
\title{Test of linear restrictions for "gmm" class objects}
\description{
 It tests the null hypothesis \eqn{R\theta =c}
}
\usage{
lintest(object,R,c)
}
\arguments{
 \item{object}{An object of class "gmm" or "gel" returned by the function \code{\link{gmm}} or \code{\link{gel}}.}
\item{R}{A \eqn{r\times k} matrix (see details).}
\item{c}{A \eqn{r\times 1} matrix (see details).}
}
\details{ The class object \code{\ref{gmm}} or \code{\link{gel}} returns a \eqn{k\time 1} vector of estimates. The null hypothesis is \eqn{H0:R\theta=c} which tests r linearly independent restrictions. }
\value{
It returns the wald test and the p-value of the test
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
thet <- 0.2
sd <- .2
x <- matrix(arima.sim(n=n,list(order=c(2,0,1),ar=phi,ma=thet,sd=sd)),ncol=1)
y <- x[7:n]
ym1 <- x[6:(n-1)]
ym2 <- x[5:(n-2)]

H <- cbind(x[4:(n-3)],x[3:(n-4)],x[2:(n-5)],x[1:(n-6)])
g <- y~ym1+ym2
x <- H

res <- gmm(g,x)
# The test is Thets(ym1)=.3 and Theta(ym2)=.6
R <- cbind(c(0,0),diag(2))
c <- c(.3,.6)
lintest(res,R,c)


}

