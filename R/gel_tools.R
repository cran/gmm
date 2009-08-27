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


smooth_g <- function (x, bw = bwAndrews2, prewhite = 1, ar.method = "ols",weights=weightsAndrews2,
			kernel=c("Bartlett","Parzen","Truncated","Tukey-Hanning"), approx = c("AR(1)","ARMA(1,1)"),
			tol = 1e-7) 
	{
	kernel <- match.arg(kernel)
	approx <- match.arg(approx)
		
	n <- nrow(x)
	if (is.function(weights))
		{
			w <- weights(x, bw = bw, kernel = kernel,  
			prewhite = prewhite, ar.method = ar.method, tol = tol, 
			verbose = FALSE, approx = approx)
		}
		else
			w <- weights


	rt <- length(w)
	if (rt >= 2)
		{
		w <- c(w[rt:2],w)
		w <- w / sum(w)
		rt <- rt-1
		sgt <- function(t) crossprod(x[(t-rt):(t+rt),],w)
		x[(rt+1):(n-rt),] <- t(sapply((rt+1):(n-rt),sgt))
		sx <- list("smoothx"=x,"kern_weights"=w)
		return(sx)		
		}
	else
		sx <- list("smoothx"=x,"kern_weights"=1)
		return(sx)		
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

confint.gel <- function(object, parm, level=0.95, lambda=FALSE, ...)
		{
		z <- object	
		n <- nrow(z$gt)
		
		se_par <- sqrt(diag(z$vcov_par))
		par <- z$coefficients
		tval <- par/se_par

		se_parl <- sqrt(diag(z$vcov_lambda))
		lamb <- z$lambda

		zs <- qnorm((1-level)/2,lower.tail=FALSE)
		ch <- zs*se_par

		if(!lambda)
			{
			ans <- cbind(par-ch,par+ch)
			dimnames(ans) <- list(names(par),c((1-level)/2,0.5+level/2))
			}
		if(lambda)
			{
			chl <- zs*se_parl
			ans <- cbind(lamb-chl,lamb+chl)
			dimnames(ans) <- list(names(lamb),c((1-level)/2,0.5+level/2))
			}		
		if(!missing(parm))
			ans <- ans[parm,]
		ans
		}

coef.gel <- function(object, lambda=FALSE, ...) 
	{
	if(!lambda)
		object$coefficients
	else
		object$lambda
	}

vcov.gel <- function(object, lambda=FALSE, ...) 
	{
	if(!lambda)
		object$vcov_par
	else
		object$vcov_lambda
	}

print.gel <- function(x, digits=5, ...)
	{
	cat("Type de GEL: ", x$type,"\n\n")
	cat("Coefficients:\n")
	print.default(format(coef(x), digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	cat("Lambdas:\n")
	print.default(format(coef(x,lambda=TRUE), digits=digits),
                      print.gap = 2, quote = FALSE)
	invisible(x)
	}

print.summary.gel <- function(x, digits = 5, ...)
	{
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
	cat("\nType of GEL: ", x$type,"\n\n")
	cat("Kernel: ", x$kernel,"\n\n")
	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)

	cat("\nLambdas:\n")
	print.default(format(x$lambda, digits=digits),
                      print.gap = 2, quote = FALSE)

	cat("\nTests of overidentifying restrictions:\n")
	print.default(format(x$test, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\nConvergence code for the coefficients: ",x$conv_par,"\n")
	cat("\nConvergence code for the lambdas: ",x$conv_lambda,"\n")
	
	invisible(x)
	}


summary.gel <- function(object, ...)
	{
	z <- object
	n <- nrow(z$gt)
	khat <- crossprod(z$gt)/n
	gbar <- colMeans(z$gt)
	
	se_par <- sqrt(diag(z$vcov_par))
	par <- z$coefficients
	tval <- par/se_par

	se_parl <- sqrt(diag(z$vcov_lambda))
	lamb <- z$lambda
	tvall <- lamb/se_parl

	LR_test <- 2*z$objective*n
	LM_test <- n*crossprod(z$lambda,crossprod(khat,z$lambda))
	J_test <- n*crossprod(gbar,solve(khat,gbar))
	test <- c(LR_test,LM_test,J_test)
	vptest <- pchisq(test,(ncol(z$gt)-length(z$par)),lower.tail=FALSE)
	ans <- list(type=z$type,call=z$call)
	names(ans$type) <-"Type of GEL"
	
	ans$coefficients <- round(cbind(par,se_par, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)
	ans$lambda <- round(cbind(lamb,se_parl, tvall, 2 * pnorm(abs(tvall), lower.tail = FALSE)),5)

    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    	dimnames(ans$lambda) <- list(names(z$lambda), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$test <- cbind(test,vptest)
	dimnames(ans$test) <- list(c("LR test","LM test","J test"),c("statistics","p-value"))	

	if (z$type == "EL")
		ans$badrho <- z$badrho
	if (!is.null(z$weights))
		{
		ans$weights <- z$weights
		}
	ans$conv_par <- z$conv_par
	ans$conv_pt <- z$conv_pt
	ans$conv_moment <- cbind(z$conv_moment)
	ans$conv_lambda <- z$conv_lambda
	names(ans$conv_par) <- "Convergence_code_theta"
	names(ans$conv_pt) <- "Sum_of_pt"
	names(ans$conv_lambda) <- "Convergence_code_for_lambda"
	dimnames(ans$conv_moment) <- list(names(z$gt),"Sample_moment_with_pt")
	class(ans) <- "summary.gel"
	ans	
}

get_dat <- function (formula,h) 
{
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	if (!is.matrix(h))
		h <- cbind(rep(1,length(h)),h)
	else	
		h <- cbind(rep(1,nrow(h)),h)
	colnames(h) <- c("(Intercept)",paste("h",1:(ncol(h)-1),sep=""))
	y <- as.matrix(model.response(mf, "numeric"))
	x <- as.matrix(model.matrix(mt, mf, NULL))
	if (attr(mt,"intercept")==0)
		{
		h <- as.matrix(h[,2:ncol(h)])
		}
	ny <- ncol(y)
	k <- ncol(x)
	nh <- ncol(h)
	if (nrow(y) != nrow(x) | nrow(x) != nrow(h) | nrow(y)!=nrow(h))
		stop("The number of observations of X, Y and H must be the same")
	if (nh<k)
		stop("The number of moment conditions must be at least equal to the number of coefficients to estimate")
	if (is.null(colnames(y)))
		{
		if (ny>1) 
			colnames(y) <- paste("y",1:ncol(y),sep="")
		if (ny == 1) 
			colnames(y) <- "y"
		}
	x <- cbind(y,x,h)
	return(list(x=x,nh=nh,ny=ny,k=k,mf=mf,mt=mt,cl=cl))
}


residuals.gel <- function(object,...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gel <- function(object,...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

formula.gel <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gel object was not created by a formula")
    else
	formula(x$terms)
}

