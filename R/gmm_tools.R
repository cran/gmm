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

summary.gmm <- function(object, ...)
	{
	z <- object
	se <- sqrt(diag(z$vcov))
	par <- z$coefficients
	tval <- par/se
	j <- z$objective*z$n
	ans <- list(met=z$met,kernel=z$kernel,algo=z$algo,call=z$call)
	names(ans$met) <- "GMM method"
	names(ans$kernel) <- "kernel for cov matrix"
		
	ans$coefficients <- round(cbind(par,se, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)

    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$J_test <- noquote(paste("Test-J degrees of freedom is ",z$df,sep=""))
	ans$j <- noquote(cbind(j,ifelse(z$df>0,pchisq(j,z$df,lower.tail = FALSE),"*******")))
	dimnames(ans$j) <- list("Test E(g)=0:  ",c("J-test","Pz(>j)"))
	class(ans) <- "summary.gmm"
	ans
	}

confint.gmm <- function(object, parm, level=0.95, ...)
		{
		z <- object
		se <- sqrt(diag(z$vcov))
		par <- z$coefficients
			
		zs <- qnorm((1-level)/2,lower.tail=FALSE)
		ch <- zs*se
		ans <- cbind(par-ch,par+ch)
		dimnames(ans) <- list(names(par),c((1-level)/2,0.5+level/2))
		if(!missing(parm))
			ans <- ans[parm,]
		ans
		}


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

		
residuals.gmm <- function(object,...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gmm <- function(object,...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

print.gmm <- function(x, digits=5, ...)
	{
	cat("Method\n", x$met,"\n\n")
	cat("Objective function value: ",x$objective,"\n\n")
	print.default(format(coef(x), digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	invisible(x)
	}

coef.gmm <- function(object,...) object$coefficients

vcov.gmm <- function(object,...) object$vcov


print.summary.gmm <- function(x, digits = 5, ...)
	{
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
	cat("\nMethod: ", x$met,"\n\n")
	cat("Kernel: ", x$kernel,"\n\n")
	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)

	cat("\nJ-test:\n")
	print.default(format(x$j, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	invisible(x)
	}

charStable <- function(theta,tau,pm=0)
	{
	# pm is the type parametrization as described by Nolan(2009)
	# It takes the value 0 or 1 

	# const can fixe parameters. It is NULL for no constraint or
	# a matrix in which case the constraint is theta[const[,1]]=const[,2]

	a <- theta[1]
	b <- theta[2]
	g <- theta[3]
	d <- theta[4]
	if(pm == 0)
		{
		if(a == 1)
			{
			if(g == 0)
				{
				the_car <- exp(complex(ima=d*tau)) 
				}
			else
				{
				re_p <- -g*abs(tau)
				im_p <- d*tau
				im_p[tau!=0] <- im_p[tau!=0] + re_p[tau!=0]*2/pi*b*sign(tau[tau!=0])*log(g*abs(tau[tau!=0]))
				the_car <- exp(complex(re=re_p,ima=im_p))
				}
			}
		else
			{
			if(g == 0)
				{
				the_car <- exp(complex(ima=d*tau)) 
				}
			else
				{
				phi <- tan(pi*a/2)
				re_p <- -g^a*abs(tau)^a
				im_p <- d*tau*1i
				im_p[tau!=0] <- im_p[tau!=0] + re_p*( b*phi*sign(tau[tau!=0])*(abs(g*tau[tau!=0])^(1-a)-1) )
				the_car <- exp(complex(re=re_p,ima=im_p))
				}
			}
		}

	if(pm == 1)
		{
		if(a == 1)
			{
			re_p <- -g*abs(tau)
			im_p <- d*tau
			im_p[tau!=0] <- im_p[tau!=0]+re_p*(b*2/pi*sign(tau[tau!=0])*log(abs(tau[tau!=0])))			
			the_car <- exp(complex(re=re_p,ima=im_p))
			}
		else
			{
			phi <- tan(pi*a/2)
			re_p <- -g^a*abs(tau)^a
			im_p <- re_p*(-phi*b*sign(tau))+d*tau
			the_car <- exp(complex(re=re_p,ima=im_p))
			}
		}
	return(the_car)
	}










		


