#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

kweights2 <- function(x, kernel = c("Truncated", "Bartlett", "Parzen",
                     "Tukey-Hanning", "Quadratic Spectral"), normalize = FALSE)
{
  kernel <- match.arg(kernel)
  if(normalize) {
    ca <- switch(kernel,  
      "Truncated" = 2,
      "Bartlett" = 2/3,
      "Parzen" = .539285,
      "Tukey-Hanning" = 3/4,
      "Quadratic Spectral" = 1)
  } else ca <- 1

  switch(kernel,  
  "Truncated" = { ifelse(ca * x > 1, 0, 1) },
  "Bartlett" = { ifelse(ca * x > 1, 0, 1 - abs(ca * x)) },
  "Parzen" = { 
    ifelse(ca * x > 1, 0, ifelse(ca * x < 0.5,
      1 - 6 * (ca * x)^2 + 6 * abs(ca * x)^3, 2 * (1 - abs(ca * x))^3))
  },
  "Tukey-Hanning" = {
    ifelse(ca * x > 1, 0, (1 + cos(pi * ca * x))/2)
  },
  "Quadratic Spectral" = {
    y <- 6 * pi * x/5
    ifelse(x < 1e-4, 1, 3 * (1/y)^2 * (sin(y)/y - cos(y)))
  })
}



bwNeweyWest2 <- function (x, kernel = c("Bartlett", "Parzen", 
    "Quadratic Spectral", "Truncated", "Tukey-Hanning"), 
    prewhite = 1, ar.method = "ols",...) 
{
    kernel <- match.arg(kernel)
    if (kernel %in% c("Truncated", "Tukey-Hanning")) 
        stop(paste("Automatic bandwidth selection only available for ", 
            dQuote("Bartlett"), ", ", dQuote("Parzen"), " and ", 
            dQuote("Quadratic Spectral"), " kernel. Use ", sQuote("bwAndrews2"), 
            " instead.", sep = ""))
    prewhite <- as.integer(prewhite)
    n <- nrow(x)
    k <- ncol(x)
    weights <- rep(1, k)
    if (length(weights) < 2) 
        weights <- 1
    mrate <- switch(kernel, Bartlett = 2/9, Parzen = 4/25, "Quadratic Spectral" = 2/25)
    m <- floor(ifelse(prewhite > 0, 3, 4) * (n/100)^mrate)
    if (prewhite > 0) {
        x <- as.matrix(na.omit(ar(x, order.max = prewhite, 
            demean = FALSE, aic = FALSE, method = ar.method)$resid))
        n <- n - prewhite
    }
    hw <- x %*% weights
    sigmaj <- function(j) sum(hw[1:(n - j)] * hw[(j + 1):n])/n
    sigma <- sapply(0:m, sigmaj)
    s0 <- sigma[1] + 2 * sum(sigma[-1])
    s1 <- 2 * sum(1:m * sigma[-1])
    s2 <- 2 * sum((1:m)^2 * sigma[-1])
    qrate <- 1/(2 * ifelse(kernel == "Bartlett", 1, 2) + 1)
    rval <- switch(kernel, Bartlett = {
        1.1447 * ((s1/s0)^2)^qrate
    }, Parzen = {
        2.6614 * ((s2/s0)^2)^qrate
    }, "Quadratic Spectral" = {
        1.3221 * ((s2/s0)^2)^qrate
    })
    rval <- rval * (n + prewhite)^qrate
    return(rval)
}


HAC <- function (x, weights = weightsAndrews2, bw = bwAndrews2, prewhite = FALSE, ar.method = "ols", kernel=c("Quadratic Spectral", "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx="AR(1)",tol = 1e-7) 
{
    n.orig <- n <- nrow(x)
    k <- ncol(x)
    kernel=match.arg(kernel)	
    if(prewhite > 0) 
	{
    	var.fit <- ar(x, order.max = prewhite, demean = FALSE, aic = FALSE, method = ar.method)
    	if(k > 1) D <- solve(diag(ncol(x)) - apply(var.fit$ar, 2:3, sum))
     	else D <- as.matrix(1/(1 - sum(var.fit$ar)))
	x <- as.matrix(na.omit(var.fit$resid))
	n <- n - prewhite
	}
    weights <- weights(x, ar.method = ar.method,kernel=kernel,bw=bw, approx = approx, prewhite = 1, tol = tol)
    if (length(weights) > n) 
	{
        warning("more weights than observations, only first n used")
        weights <- weights[1:n]
	}
    utu <- 0.5 * crossprod(x) * weights[1]
    wsum <- n * weights[1]/2
    w2sum <- n * weights[1]^2/2
    if (length(weights) > 1) {
        for (ii in 2:length(weights)) {
            utu <- utu + weights[ii] * crossprod(x[1:(n - 
                ii + 1), , drop = FALSE], x[ii:n, , drop = FALSE])
            wsum <- wsum + (n - ii + 1) * weights[ii]
            w2sum <- w2sum + (n - ii + 1) * weights[ii]^2
        }
    }
    utu <- utu + t(utu)
    
    if(prewhite > 0) {
    utu <- crossprod(t(D), utu) %*% t(D)
     }
    wsum <- 2 * wsum
    w2sum <- 2 * w2sum
    bc <- n^2/(n^2 - wsum)
    df <- n^2/w2sum
    rval <- utu/n.orig
    
    return(rval)
}

weightsAndrews2 <- function (x, bw = bwAndrews2, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol = 1e-7, verbose = FALSE)
{
    kernel <- match.arg(kernel)
    approx=match.arg(approx)

    if (is.function(bw)) 
        bw <- bw(x, kernel = kernel, prewhite = prewhite, ar.method = ar.method, approx=approx)
    n <- NROW(x) 
    weights <- kweights2(0:(n - 1)/bw, kernel = kernel)
    weights <- weights[1:max(which(abs(weights) > tol))]
    return(weights)
}


bwAndrews2 <- function (x, kernel = c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"), approx = c("AR(1)", 
    "ARMA(1,1)"), prewhite = 1, ar.method = "ols") 
{
    kernel <- match.arg(kernel)
    approx <- match.arg(approx)
    n <- nrow(x)
    k <- ncol(x)

    if (approx == "AR(1)") {
        fitAR1 <- function(x) {
            rval <- ar(x, order.max = 1, aic = FALSE, method = "ols")
            rval <- c(rval$ar, sqrt(rval$var.pred))
            names(rval) <- c("rho", "sigma")
            return(rval)
        }
        ar.coef <- apply(x, 2, fitAR1)
        denum <- sum((ar.coef["sigma", ]/(1 - ar.coef["rho", 
            ]))^4)
        alpha2 <- sum(4 * ar.coef["rho", ]^2 * ar.coef["sigma", 
            ]^4/(1 - ar.coef["rho", ])^8)/denum
        alpha1 <- sum(4 * ar.coef["rho", ]^2 * ar.coef["sigma", 
            ]^4/((1 - ar.coef["rho", ])^6 * (1 + ar.coef["rho", 
            ])^2))/denum
    }
    else {
        fitARMA11 <- function(x) {
            rval <- arima(x, order = c(1, 0, 1), include.mean = FALSE)
            rval <- c(rval$coef, sqrt(rval$sigma2))
            names(rval) <- c("rho", "psi", "sigma")
            return(rval)
        }
        arma.coef <- apply(x, 2, fitARMA11)
        denum <- sum(((1 + arma.coef["psi", ]) * arma.coef["sigma", 
            ]/(1 - arma.coef["rho", ]))^4)
        alpha2 <- sum(4 * ((1 + arma.coef["rho", ] * 
            arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", 
            ]))^2 * arma.coef["sigma", ]^4/(1 - arma.coef["rho", 
            ])^8)/denum
        alpha1 <- sum(4 * ((1 + arma.coef["rho", ] * 
            arma.coef["psi", ]) * (arma.coef["rho", ] + arma.coef["psi", 
            ]))^2 * arma.coef["sigma", ]^4/((1 - arma.coef["rho", 
            ])^6 * (1 + arma.coef["rho", ])^2))/denum
    }
    rval <- switch(kernel, Truncated = {
        0.6611 * (n * alpha2)^(1/5)
    }, Bartlett = {
        1.1447 * (n * alpha1)^(1/3)
    }, Parzen = {
        2.6614 * (n * alpha2)^(1/5)
    }, "Tukey-Hanning" = {
        1.7462 * (n * alpha2)^(1/5)
    }, "Quadratic Spectral" = {
        1.3221 * (n * alpha2)^(1/5)
    })
   return(rval)
}

