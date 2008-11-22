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

summary.gel <- function(object,interval=FALSE, ...)
	{
	z <- object
	n <- nrow(z$gt)
	khat <- crossprod(z$gt)/n
	gbar <- colMeans(z$gt)
	
	se_par <- sqrt(diag(z$vcov_par))
	par <- z$par
	tval <- par/se_par

	se_parl <- sqrt(diag(z$vcov_lambda))
	lamb <- z$lambda
	tvall <- lamb/se_parl

	LR_test <- 2*z$objective*n
	LM_test <- n*crossprod(z$lambda,crossprod(khat,z$lambda))
	J_test <- n*crossprod(gbar,solve(khat,gbar))
	test <- c(LR_test,LM_test,J_test)
	vptest <- pchisq(test,(ncol(z$gt)-length(z$par)),lower.tail=FALSE)
	ans <- list(type=z$type)
	names(ans$type) <-"Type of GEL"
	
	ans$par <- round(cbind(par,se_par, tval, 2 * pnorm(abs(tval), lower.tail = FALSE)),5)
	ans$lambda <- round(cbind(lamb,se_parl, tvall, 2 * pnorm(abs(tvall), lower.tail = FALSE)),5)

    	dimnames(ans$par) <- list(names(z$par), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    	dimnames(ans$lambda) <- list(names(z$lambda), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))

	ans$test <- cbind(test,vptest)
	dimnames(ans$test) <- list(c("LR test","LM test","J test"),c("statistics","p-value"))	

	if (interval != FALSE)
		{
		zs <- qnorm((1-interval)/2,lower.tail=FALSE)
		ch <- zs*se_par
		ans$interval_par <- cbind(par-ch,par+ch)
		dimnames(ans$interval_par) <- list(names(par),c("Theta_lower","Theta_upper"))
		chl <- zs*se_parl
		ans$interval_lam <- cbind(lamb-chl,lamb+chl)
		dimnames(ans$interval_lam) <- list(names(lamb),c("Lambda_lower","Lambda_upper"))
		}


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

get_dat <- function (formula,h,intercept=TRUE) 
{
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")
	
	h <- cbind(rep(1,nrow(h)),h)
	colnames(h) <- c("Intercept",paste("h",1:(ncol(h)-1),sep=""))
	y <- as.matrix(model.response(mf, "numeric"))
	x <- as.matrix(model.matrix(mt, mf, NULL))
	if (!intercept)
		{
		x <- as.matrix(x[,2:ncol(x)])
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
	return(list(x=x,nh=nh,ny=ny,k=k))
}

