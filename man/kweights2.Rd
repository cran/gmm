\name{kweights2}
\alias{kweights2}
\title{Kernel Weights}
\description{
 Kernel weights for kernel-based heteroskedasticity
 and autocorrelation consistent (HAC) covariance matrix estimators
 as introduced by Andrews (1991). This function comes from the \code{sandwich} package.
}

\usage{
kweights2(x, kernel = c("Truncated", "Bartlett", "Parzen",
  "Tukey-Hanning", "Quadratic Spectral"), normalize = FALSE)
}

\arguments{
 \item{x}{numeric.}
 \item{kernel}{a character specifying the kernel used. All kernels used
    are described in Andrews (1991).}
 \item{normalize}{logical. If set to \code{TRUE} the kernels are
   normalized as described in Andrews (1991).}    
}

\value{
Value of the kernel function at \code{x}.
}

\references{
  Andrews DWK (1991),
  Heteroskedasticity and Autocorrelation Consistent Covariance Matrix Estimation.
  \emph{Econometrica}, \bold{59},
  817--858.
}

\keyword{regression}
\keyword{ts}
