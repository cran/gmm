\name{plot.gel}
\alias{plot.gel}
\title{Plot Diagnostics for a gel Object}
\description{
 It is a plot method for \code{gel} objects.
}
\usage{
\method{plot}{gel}(x, which = c(1L:4),
	    main = list("Residuals vs Fitted values", "Normal Q-Q",
	    "Response variable and fitted values","Implied probabilities"),
	    panel = if(add.smooth) panel.smooth else points,
	    ask = prod(par("mfcol")) < length(which) && dev.interactive(), ...,
	    add.smooth = getOption("add.smooth"))
}
\arguments{
  \item{x}{\code{gel} object, typically result of \code{\link{gel}}.}
  \item{which}{if a subset of the plots is required, specify a subset of
    the numbers \code{1:4}.}
  \item{main}{Vector of titles for each plot.
  }
  \item{panel}{panel function.  The useful alternative to
    \code{\link{points}}, \code{\link{panel.smooth}} can be chosen
    by \code{add.smooth = TRUE}.}
  \item{ask}{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, see \code{\link{par}(ask=.)}.}
  \item{\dots}{other parameters to be passed through to plotting
    functions.}
  \item{add.smooth}{logical indicating if a smoother should be added to
    most plots; see also \code{panel} above.}
}
\details{
It is a beta version of a plot method for \code{gel} objects. It is a modified version of \code{plot.lm}. For now, it is available only for linear models expressed as a formula. Any suggestions are welcome regarding plots or options to include. 
The first two plots are the same as the ones provided by \code{plot.lm}, the third is the dependant variable \eqn{y} with its mean \eqn{\hat{y}} (the fitted values) and the last plots the implied probabilities with the empirical density \eqn{1/T}.
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

plot(res,which=3)
plot(res,which=4)
}
