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

rho <- function(x,lamb,derive=0,type=c("EL","ET","CUE"),drop=TRUE)
	{

	type <- match.arg(type)
	lamb <- matrix(lamb,ncol=1)
	gml <- x%*%lamb
	ch <- 0
	if (derive==0)
		{
		if (type == "EL")
			{
			ch <- sum(gml>=1)
			if (drop)
				{				
				gml <- (gml<1)*gml
				rhomat <- log(1-gml) 
				}
			else
				{
				if (ch>0)
					rhomat <- NaN
				else
					rhomat <- log(1-gml) 
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
			rhomat <- -1/(1-gml) 
			
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -1 -gml
		}
	if (derive==2)
		{
		if (type == "EL")
			rhomat <- -1/(1-gml)^2 
			
		if (type == "ET")
			rhomat <- -exp(gml)
		
		if (type == "CUE")
			rhomat <- -1
		}
	rhom <-list(ch = ch, rhomat=rhomat) 
	return(rhom)
	}

get_lamb <- function(g,tet,x,type=c('EL','ET','CUE'),tol_lam=1e-12,maxiterlam=1000,tol_obj=1e-7)
	{
	type <- match.arg(type)	
	gt <- g(tet,x)
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
	while ((crit > tol_lam*( 1+sqrt( crossprod(lamb0) ) ) ) & (j<=maxiterlam))
		{ 
		rho2 <- as.numeric(rho(gt,lamb0,derive=2,type=type)$rhomat)
		rho1 <- as.numeric(rho(gt,lamb0,derive=1,type=type)$rhomat)
		gblam <- colMeans(rho1*gt)
		klam <- crossprod(rho2*gt,gt)/n
		chklam <- sum(abs(klam))
		if (!is.null(gblam0))
			dgblam <- crossprod(gblam)-crossprod(gblam0)
		
		if (is.na(chklam) | chklam == 0 | chklam == Inf |  dgblam>0 | dgblam == Inf | is.na(dgblam) | dcrit < 0)
			{
			lamb1 <- rep(sqrt(1/n),length(lamb0))
			crit <- 0
			singular=2
			conv_mes <- "The algorithm produced singular system,  NaN or Inf" 
			}
		else
			{
			if (rcond(klam)>tol_cond)
				{
				lamb1 <- lamb0-solve(klam,gblam)
				crit <- sqrt(crossprod(lamb0-lamb1))
				lamb0 <- lamb1
				}
			else
				{
				lamb1 <- rep(sqrt(1/n),length(lamb0))
				crit <- 0
				singular=2
				conv_mes <- "The algorithm produced singular system" 
				}
			}
		gblam0 <- gblam
		j <- j+1
		dcrit<- crit0-crit
		crit0 <- crit
		}
	z <- list("lambda"=lamb1,singular=singular,conv_mes=conv_mes)
	if (j>maxiterlam | max(abs(gblam))>tol_obj)
		{
		singular <- 1
		conv_mes <- "No convergence after 'maxiterlam' iterations"
		z$singular <- singular		
		}
		z$obj <- crossprod(gblam)
	return(z)
	}

gel <- function(g,x,tet0,gradv=NULL,smooth=FALSE,type=c("EL","ET","CUE","ETEL"), vcov=c("HAC","iid"), kernel = c("Bartlett", "Parzen", 
  		"Truncated", "Tukey-Hanning"), bw=bwAndrews2, approx = c("AR(1)", 
    		"ARMA(1,1)"), prewhite = 1, ar.method = "ols", tol_weights = 1e-7, tol_lam=1e-9, tol_obj = 1e-9, 
		tol_mom = 1e-9,maxiterlam=1000, constraint=FALSE,optfct=c("optim","optimize","nlminb"),optlam=c("iter","numeric"),
		model=TRUE, X=FALSE, Y=FALSE,...)
	{
	vcov=match.arg(vcov)
	type <- match.arg(type)
	optfct <- match.arg(optfct)
	optlam <- match.arg(optlam)
	weights=weightsAndrews2
	if (type == "ETEL")
		{
		typel <- "ET"
		typet <- "EL"	
		}
	else
		{
		typel <- type
		typet <- type
		}
	approx <- match.arg(approx)
	kernel <- match.arg(kernel)
	if(optfct=="optim")
		k <- length(tet0)
	else
		k <- 1

	typeg=0
	if (is(g,"formula"))
		{
		typeg=1
		dat <- get_dat(g,x)
		x <- dat$x
		
		g <- function(tet,x,ny=dat$ny,nh=dat$nh,k=dat$k)
			{
			tet <- matrix(tet,ncol=k)
			e <- x[,1:ny] -  x[,(ny+1):(ny+k)]%*%t(tet)
			gt <- e*x[,ny+k+1]
			if (nh > 1)
				{	
				for (i in 2:nh)
					{
					gt <- cbind(gt,e*x[,(ny+k+i)])
					}
				}
			return(gt)
			}
		gradv <- function(tet,x,ny=dat$ny,nh=dat$nh,k=dat$k)
			{
			tet <- matrix(tet,ncol=k)
			dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1,ny))/nrow(x)
			return(dgb)
			}
		}	
		if (typeg)
			n <- nrow(x)
		else
			n = nrow(g(tet0,x))

	if (smooth)
		{
		g1 <- g
		rgmm <- gmm(g,x,tet0,wmatrix="ident")

		if (is.function(weights))
			w <- weights(g(rgmm$coefficients,x),kernel=kernel,bw=bw,prewhite = prewhite,ar.method=ar.method,approx=approx,tol=tol_weights)
		else
			w <- weights
		sg <- function(thet,x)
			{
			gf <- g1(thet,x)
			gt <- smooth_g(gf, weights=w)$smoothx 
			return(gt)
			}
		g <- sg
		}	
	lll <- 1
	thetf <- function(tet)
		{
		if (optlam=="iter")
			{
			lamblist <- get_lamb(g,tet,x,type=typel,tol_lam=tol_lam,maxiterlam=maxiterlam,tol_obj=tol_obj)
			lamb <- lamblist$lambda
			gt <- g(tet,x)
			pt <- -rho(gt,lamb,type=typet,derive=1)$rhomat/nrow(gt)
			checkmom <- sum(as.numeric(pt)*gt)
			if (lamblist$singular==0)		
				p <- sum(rho(gt,lamb,type=typet)$rhomat) + abs(checkmom)/tol_mom
			if (lamblist$singular==1)		
				p <- sum(rho(gt,lamb,type=typet)$rhomat) + abs(checkmom)/tol_mom + lamblist$obj/tol_mom
			if (lamblist$singular==2)		
				p <- 1e50*lll
			lll <- lll+1
			}
		else
			{
			gt <- g(tet,x)
			rhofct <- function(lamb)
				{
				rhof <- -sum(rho(gt,lamb,type=typel)$rhomat)
				return(rhof)
				}
			if (ncol(gt)>1)
				rlamb <- optim(rep(0,ncol(gt)),rhofct,control=list(maxit=1000))
			else
				{
				rlamb <- optimize(rhofct,c(-1,1))
				rlamb$par <- rlamb$minimum
				rlamb$value <- rlamb$objective
				}
			lamb <- rlamb$par
			pt <- -rho(gt,lamb,type=typet,derive=1)$rhomat/nrow(gt)
			checkmom <- sum(as.numeric(pt)*gt)
			p <- -rlamb$value + (checkmom)^2/tol_mom + (sum(as.numeric(pt))-1)^2/tol_mom
			}
		return(p)
		}

	if (!constraint)
		{
		if (optfct == "optim")
			res <- optim(tet0,thetf,...)
		if (optfct == "nlminb")
			res <- nlminb(tet0,thetf,...)
		
		if (optfct == "optimize")
			{
			res <- optimize(thetf,tet0,...)
			res$par <- res$minimum
			res$convergence <- "There is no convergence code for optimize"
			}
		}
	if(constraint)
		res <- constrOptim(tet0,thetf,grad=NULL,...)


	if (optlam=="iter")
		{
		rlamb <- get_lamb(g,res$par,x,type=typel,tol_lam=tol_lam,maxiterlam=maxiterlam,tol_obj=tol_obj)
		z <- list(coefficients=res$par,lambda=rlamb$lam,conv_lambda=rlamb$conv_mes,conv_par=res$convergence)
		z$foc_lambda <- rlamb$obj
		}
	if (optlam=="numeric")
		{
		gt<-g(res$par,x)
		rhofct <- function(lamb)
			{
			rhof <- -sum(rho(gt,lamb,type=typel)$rhomat)
			return(rhof)
			}
		rlamb <- optim(rep(0,ncol(gt)),rhofct,control=list(maxit=1000))
		z <- list(coefficients=res$par,conv_par=res$convergence,lambda=rlamb$par)
		z$conv_lambda=paste("Lambda by optim. Conv. code = ",rlamb$convergence,sep="")
		rho1 <- as.numeric(rho(gt,z$lambda,derive=1,type=typel)$rhomat)
		z$foc_lambda <- crossprod(colMeans(rho1*gt))
		}
	
	z$type <- type
	z$gt <- g(z$coefficients,x)
	rhom <- rho(z$gt,z$lambda,type=typet)
	z$pt <- -rho(z$gt,z$lambda,type=typet,derive=1)$rhomat/n
	z$conv_moment <- colSums(as.numeric(z$pt)*z$gt)
	z$conv_pt <- sum(as.numeric(z$pt))
	z$objective <- sum(as.numeric(rhom$rhomat)-rho(1,0,type=typet)$rhomat)/n

	if (type == "EL")	
		{
		z$badrho <- rhom$ch
		names(z$badrho) <- "Number_of_bad_rho"
		}

	Gf <- function(thet)
		{
		myenv <- new.env()
		assign('x',x,envir=myenv)
		assign('thet',thet,envir=myenv)
		barg <- function(x,thet)
			{
			gt <- g(thet,x)
			gbar <- as.vector(colMeans(gt))
			return(gbar)
			}
		G <- attr(numericDeriv(quote(barg(x,thet)),"thet",myenv),"gradient")
		return(G)
		}
	if (!is.function(gradv)) 
		G <- Gf(z$coefficients)
	else
		G <- gradv(z$coefficients,x)
	if (vcov == "iid")
		khat <- crossprod(z$gt)
	else
		khat <- HAC(g(z$coefficients,x), kernel=kernel, bw=bw,prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol_weights)

	kg <- solve(khat,G)
	z$vcov_par <- solve(crossprod(G,kg))/n
	z$vcov_lambda <- ((solve(khat)-kg%*%z$vcov_par%*%t(kg)))/n
	
	if (smooth) z$weights<-w

	if (typeg ==0)
		{
		names(z$coefficients) <- paste("Theta[",1:k,"]",sep="")
		colnames(z$gt) <- paste("gt[",1:ncol(z$gt),"]",sep="")
		names(z$lambda) <- paste("Lambda[",1:ncol(z$gt),"]",sep="")
		}
	if (typeg == 1)
		{
		namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
		nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
		if (dat$ny > 1)
			{
			namey <- colnames(dat$x[,1:dat$ny])
			names(z$coefficients) <- paste(rep(namey,dat$k),"_",rep(namex,rep(dat$ny,dat$k)),sep="")
			colnames(z$gt) <- paste(rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),sep="")
			names(z$lambda) <- paste("Lam(",rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),")",sep="")
			}
		if (dat$ny == 1)
			{
			names(z$coefficients) <- namex
			colnames(z$gt) <- nameh
			names(z$lambda) <- nameh
			}
		}
	dimnames(z$vcov_par) <- list(names(z$coefficients),names(z$coefficients))
	dimnames(z$vcov_lambda) <- list(names(z$lambda),names(z$lambda))
	if (typeg==1)
		{
		b <- z$coefficients
		y <- as.matrix(model.response(dat$mf, "numeric"))
		ny <- dat$ny
		b <- t(matrix(b,nrow=ny))
		x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
		yhat <- x%*%b
		z$fitted.values<-yhat	
		z$residuals<-y-yhat	
		z$terms<- dat$mt
		if(model) z$model<-dat$mf
		if(X) z$x<-x
		if(Y) z$y<-y
		}
	else
		if(X) z$x<-x
	
	z$call <- match.call()
	class(z) <- "gel"
	return(z)
	}
