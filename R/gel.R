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

rho <- function(x, lamb, derive = 0, type = c("EL", "ET", "CUE"), drop = TRUE, k = 1)
	{

	type <- match.arg(type)
	lamb <- matrix(lamb, ncol = 1)
	gml <- x%*%lamb*k
	ch <- 0
	if (derive == 0)
		{
		if (type == "EL")
			{
			ch <- sum(gml >= 1)
			if (drop)
				{				
				gml <- (gml < 1)*gml
				rhomat <- log(1 - gml) 
				}
			else
				{
				if (ch > 0)
					rhomat <- NaN
				else
					rhomat <- log(1 - gml) 
				}
			}
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -gml -0.5*gml^2
		}
	if (derive==1)
		{
		if (type == "EL")
			rhomat <- -1/(1 - gml) 
			
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -1 - gml
		}
	if (derive==2)
		{
		if (type == "EL")
			rhomat <- -1/(1 - gml)^2 
			
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -1
		}
	rhom <-list(ch = ch, rhomat = rhomat) 
	return(rhom)
	}

getLamb <- function(g, tet, x, type = c('EL', 'ET', 'CUE'), tol_lam = 1e-12, maxiterlam = 1000, tol_obj = 1e-7, k = 1)
	{
	type <- match.arg(type)	
	gt <- g(tet, x)

	n <- nrow(gt)
	tol_cond=1e-12
	gb <- colMeans(gt)
	khat <- crossprod(gt)/n
	lamb0 <- -solve(khat,gb)

	conv_mes <- "Normal convergence" 
	singular <-0
	crit <-1e30
	crit0 <- crit
	dcrit <- 10
	dgblam <- -10
	gblam0 <- NULL

	j <- 1
	while ((crit > tol_lam*( 1 + sqrt( crossprod(lamb0) ) ) ) & (j <= maxiterlam))
		{ 
		rho2 <- as.numeric(rho(gt, lamb0, derive = 2, type = type, k = k)$rhomat)
		rho1 <- as.numeric(rho(gt, lamb0, derive = 1, type = type, k = k)$rhomat)
		gblam <- colMeans(rho1*gt)
		klam <- crossprod(rho2*gt, gt)/n
		chklam <- sum(abs(klam))
		if (!is.null(gblam0))
			dgblam <- crossprod(gblam) - crossprod(gblam0)
		
		#if (is.na(chklam) | chklam == 0 | chklam == Inf |  dgblam > 0 | dgblam == Inf | is.na(dgblam) | dcrit < 0)
		if (is.na(chklam) | chklam == 0 | chklam == Inf | dgblam == Inf | is.na(dgblam))
			{
			lamb1 <- rep(0, length(lamb0))
			crit <- 0
			singular=2
			conv_mes <- "The algorithm produced singular system,  NaN or Inf" 
			}
		else
			{
			if (rcond(klam) > tol_cond)
				{
				lamb1 <- lamb0 - solve(klam, gblam)
                                crit <- sqrt(crossprod(lamb0 - lamb1))
				lamb0 <- lamb1
				}
			else
				{
				lamb1 <- rep(0 , length(lamb0))
				crit <- 0
				singular <- 2
				conv_mes <- "The algorithm produced singular system" 
				}
			}
		gblam0 <- gblam
		j <- j + 1
		dcrit<- crit0 - crit
		crit0 <- crit
		}
	z <- list("lambda" = lamb1, singular = singular, conv_mes = conv_mes)
	if (j > maxiterlam)
		{
		singular <- 1
		conv_mes <- "No convergence after 'maxiterlam' iterations"
		z$singular <- singular
                stop("Maxiterlam reached.\n Increase it, try other starting values \n or use the option optlam=\"numeric\".")		
		}
		z$obj <- crossprod(gblam)
	return(z)
	}

smoothG <- function (x, bw = bwAndrews, prewhite = 1, ar.method = "ols", weights = weightsAndrews,
			kernel = c("Bartlett", "Parzen", "Truncated", "Tukey-Hanning"), approx = c("AR(1)", "ARMA(1,1)"),
			tol = 1e-7) 
	{
	kernel <- match.arg(kernel)
	approx <- match.arg(approx)
		
	n <- nrow(x)
	if (is.function(weights))
		{
                        class(x) <- "gmmFct"
			w <- weights(x, bw = bw, kernel = kernel,  
			prewhite = prewhite, ar.method = ar.method, tol = tol, 
			verbose = FALSE, approx = approx)
		}
		else
			w <- weights


	rt <- length(w)
	if (rt >= 2)
		{
		w <- c(w[rt:2], w)
		w <- w / sum(w)
		rt <- rt - 1
		sgt <- function(t) crossprod(x[(t-rt):(t+rt),], w)
		x[(rt+1):(n-rt),] <- t(sapply((rt + 1):(n - rt), sgt))
		sx <- list("smoothx" = x, "kern_weights" = w)
		return(sx)		
		}
	else
		sx <- list("smoothx" = x,"kern_weights" = 1)
		return(sx)		
	}


gel <- function(g, x, tet0, gradv = NULL, smooth = FALSE, type = c("EL", "ET", "CUE", "ETEL"), 
                kernel = c("Truncated", "Bartlett"), bw = bwAndrews, approx = c("AR(1)", 
    		"ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol_weights = 1e-7, tol_lam = 1e-9, tol_obj = 1e-9, 
		tol_mom = 1e-9, maxiterlam = 100, constraint = FALSE, optfct = c("optim", "optimize", "nlminb"), 
                optlam = c("iter", "numeric"), model = TRUE, X = FALSE, Y = FALSE, TypeGel = "baseGel", ...)
	{

	type <- match.arg(type)
	optfct <- match.arg(optfct)
	optlam <- match.arg(optlam)
	weights <- weightsAndrews
	approx <- match.arg(approx)
	kernel <- match.arg(kernel)

	all_args <- list(g = g, x = x, tet0 = tet0, gradv = gradv, smooth = smooth, type = type,
                kernel = kernel, bw = bw, approx = approx, prewhite = prewhite, ar.method = ar.method, 
		tol_weights = tol_weights, tol_lam = tol_lam, tol_obj = tol_obj, tol_mom = tol_mom, 
		maxiterlam = maxiterlam, constraint = constraint, optfct = optfct, weights = weights,
                optlam = optlam, model = model, X = X, Y = Y, TypeGel = TypeGel, call = match.call())

	class(all_args)<-TypeGel
	Model_info<-getModel(all_args)
	z <- momentEstim(Model_info, ...)

	class(z) <- "gel"
	return(z)
	}


  .thetf <- function(tet, P)
    {
    if(!is.null(P$gform))
      {
      dat <- P$dat
      x <- dat$x
      }
    else
      x <- P$x

    if (P$optlam == "iter")
      {
      lamblist <- getLamb(P$g, tet, x, type = P$typel, tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, tol_obj = P$tol_obj, k = P$k1/P$k2)
      lamb <- lamblist$lambda
      gt <- P$g(tet, x)
      pt <- -rho(gt, lamb, type = P$typet, derive = 1, k = P$k1/P$k2)$rhomat/nrow(gt)
      checkmom <- sum(as.numeric(pt)*gt)
      if (lamblist$singular == 0)		
        p <- sum(rho(gt, lamb, type = P$typet, k = P$k1/P$k2)$rhomat) + abs(checkmom)/P$tol_mom
      if (lamblist$singular == 1)		
        p <- sum(rho(gt, lamb, type = P$typet, k = P$k1/P$k2)$rhomat) + abs(checkmom)/P$tol_mom + lamblist$obj/P$tol_mom
      if (lamblist$singular == 2)		
        p <- 1e50*proc.time()[3]
      }
    else
      {
      gt <- P$g(tet, x)
      rhofct <- function(lamb)
        {
        rhof <- -sum(rho(gt, lamb, type = P$typel, k = P$k1/P$k2)$rhomat)
        return(rhof)
        }
      if (ncol(gt) > 1)
        rlamb <- optim(rep(0, ncol(gt)), rhofct, control = list(maxit = 1000))
      else
        {
        rlamb <- optimize(rhofct, c(-1,1))
        rlamb$par <- rlamb$minimum
        rlamb$value <- rlamb$objective
        }
      lamb <- rlamb$par
      pt <- -rho(gt, lamb, type = P$typet, derive = 1, k = P$k1/P$k2)$rhomat/nrow(gt)
      checkmom <- sum(as.numeric(pt)*gt)
      p <- -rlamb$value + (checkmom)^2/P$tol_mom + (sum(as.numeric(pt)) - 1)^2/P$tol_mom
      }
    return(p)
    }


