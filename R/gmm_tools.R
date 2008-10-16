HAC <- function (x, weights = weightsAndrews2, prewhite = FALSE, ar.method = "ols",kernel=c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning")) 
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
    if (is.function(weights)) 
        weights <- weights(x, ar.method = ar.method,kernel=kernel)
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
    if (is.function(bw)) 
        bw <- bw(x, kernel = kernel, prewhite = prewhite, ar.method = ar.method, approx=approx)
    n <- NROW(x) 
    weights <- kweights(0:(n - 1)/bw, kernel = kernel)
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

summary.gmm <- function(object, ...)
	{
	z <- object
	se <- sqrt(diag(z$vcov))
	par <- z$par
	tval <- par/se
	j <- z$objective*z$n
	ans <- list(met=z$met,kernel=z$kernel,algo=z$algo)
	names(ans$met) <- "GMM method"
	names(ans$kernel) <- "kernel for cov matrix"
	names(ans$algo) <- "Numerical algorithm"
	
	
	ans$par <- round(cbind(par,se, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)

    	dimnames(ans$par) <- list(names(z$par), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$j <- cbind(j,ifelse(z$df>0,pchisq(j,z$df,lower.tail = FALSE),"*******"))
	dimnames(ans$j) <- list("Test E(g)=0   ",c("J-test","Pz(>j)"))
	
	class(ans) <- "summary.gmm"
	ans
	}




