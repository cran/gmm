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

gmm <- function(g,x,t0=NULL,gradv=NULL, type=c("twoStep","cue","iterative"), wmatrix = c("optimal","ident"),  vcov=c("HAC","iid"), 
	      kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews2, 
	      prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7, itermax=100,optfct=c("optim","optimize","nlminb"), model=TRUE, X=FALSE, 		      Y=FALSE, ...)
{
type <- match.arg(type)
kernel <- match.arg(kernel)
vcov <- match.arg(vcov)
wmatrix <- match.arg(wmatrix)
optfct <- match.arg(optfct)
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
		gradv <- function(x,ny=dat$ny,nh=dat$nh,k=dat$k)
			{
			dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1,ny))/nrow(x)
			return(dgb)
			}
		tetlin <- function(x,w,ny=dat$ny,nh=dat$nh,k=dat$k)
			{
			n <- nrow(x)
			ym <- as.matrix(x[,1:ny])
			xm <- as.matrix(x[,(ny+1):(ny+k)])
			hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
			whx <- solve(w,(crossprod(hm,xm)%x%diag(ny)))	
			wvecyh <- solve(w,matrix(crossprod(ym,hm),ncol=1))
			dg <- gradv(x)
			xx <- crossprod(dg,whx)
			par <- solve(xx,crossprod(dg,wvecyh))
			gb <- matrix(colSums(g(par,x))/n,ncol=1)
			value <- crossprod(gb,solve(w,gb)) 
			res <- list(par=par,value=value)
			return(res)
			}
		}
if (optfct == "optimize")
	{
	n = nrow(g(t0[1],x))
	q = ncol(g(t0[1],x))
	k = 1
	k2 <- k
	df <- q-k
	}
else
	{
	if (typeg)
		{
		k <- dat$k
		k2 <- k*dat$ny
		n <- nrow(x)
		q <- dat$ny*dat$nh
		df <- q-k*dat$ny
		}
	else
		{
		n = nrow(g(t0,x))
		q = ncol(g(t0,x))
		k = length(t0)
		k2 <- k
		df <- q-k
		}
	}
obj1 <- function(thet)
	{
	gt <- g(thet,x)
	gbar <- as.vector(colMeans(gt))
	obj <- crossprod(gbar,solve(w,gbar))
	return(obj)
	}
iid <- function(thet)
	{
	gt <- g(thet,x)
	v <- crossprod(gt,gt)/n
	return(v)
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

if (q == k2 | wmatrix == "ident")
	{
	if (typeg)
		{
		w <- diag(q)
		res <- tetlin(x,w)
		z = list(coefficients=res$par,objective=res$value)
		}
	else
		{
		w=diag(rep(1,q))
		if (optfct == "optim")
			res <- optim(t0,obj1, ...)
		if (optfct == "nlminb")
			{
			res <- nlminb(t0,obj1, ...)
			res$value <- res$objective
			}
		if (optfct == "optimize")
					{
			res <- optimize(obj1,t0, ...)
			res$par <- res$minimum
			res$value <- res$objective
			}	
		z = list(coefficients=res$par,objective=res$value)	
		}
	}
else
	{
	if (type=="twoStep")
		{
		w=diag(rep(1,q))
		if (typeg)
			{
			res1 <- tetlin(x,w)
			}
		else
			{
			if (optfct == "optim")
				res1 <- optim(t0,obj1, ...)
			if (optfct == "nlminb")
				{
				res1 <- nlminb(t0,obj1, ...)
				res1$value <- res1$objective
				}
			if (optfct == "optimize")
				{
				res1 <- optimize(obj1,t0, ...)
				res1$par <- res1$minimum
				res1$value <- res1$objective
				}	
			}
		if (vcov == "iid")
			w <- iid(res1$par)
		if (vcov == "HAC")
			w <- HAC(g(res1$par,x), kernel=kernel, bw=bw,prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
		if (typeg)
			{
			res2 <- tetlin(x,w)
			}
		else
			{
			if (optfct == "optim")
				res2 <- optim(res1$par,obj1, ...)
			if (optfct == "nlminb")
				{
				res2 <- nlminb(res1$par,obj1, ...)
				res2$value <- res2$objective
				}
			if (optfct == "optimize")
				{
				res2 <- optimize(obj1,t0, ...)
				res2$par <- res2$minimum
				res2$value <- res2$objective
				}	
			}
		z = list(coefficients=res2$par,objective=res2$value)	
		}
	if (type=="cue")
		{
		obj_cue <- function(thet)
			{
			gt <- g(thet,x)
			gbar <- as.vector(colMeans(gt))
			if (vcov == "iid")
				w2 <- iid(thet)
			if (vcov == "HAC")
				w2 <- HAC(g(thet,x), kernel=kernel, bw=bw,prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
			obj <- crossprod(gbar,solve(w2,gbar))
			return(obj)
			}	
		if (typeg)
			{
			if (is.null(t0))
				t0 <- tetlin(x,diag(rep(1,q)))$par
			if (optfct == "optim")
				res2 <- optim(t0,obj_cue, ...)
			if (optfct == "nlminb")
				{
				res2 <- nlminb(t0,obj_cue, ...)
				res2$value <- res2$objective
				}
			if (optfct == "optimize")
				{
				res2 <- optimize(obj_cue,t0, ...)
				res2$par <- res2$minimum
				res2$value <- res2$objective
				}
			}
		else
			{
			if (optfct == "optim")
				res2 <- optim(t0,obj_cue, ...)
			if (optfct == "nlminb")
				{
				res2 <- nlminb(t0,obj_cue, ...)
				res2$value <- res2$objective
				}
			if (optfct == "optimize")
				{
				res2 <- optimize(obj_cue,t0, ...)
				res2$par <- res2$minimum
				res2$value <- res2$objective
				}	
			}
		z = list(coefficients=res2$par,objective=res2$value)	
		}
	if (type=="iterative")
		{
		w=diag(rep(1,q))
		if (typeg)
			res <- tetlin(x,w)
		else
			{
			if (optfct == "optim")
				res <- optim(t0,obj1, ...)
			if (optfct == "nlminb")
				{
				res <- nlminb(t0,obj1, ...)
				res$value <- res$objective
				}
			if (optfct == "optimize")
				{
				res <- optimize(obj1,t0, ...)
				res$par <- res$minimum
				res$value <- res$objective
				}	
			}
		ch <- 100000
		j <- 1
		while(ch>crit)
		{
			tet <- res$par
			if (vcov == "iid")
				w <- iid(tet)
			if (vcov == "HAC")
				w <- HAC(g(tet,x), kernel=kernel, bw=bw,prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)
			if (typeg)
				res <- tetlin(x,w)
			else
				{	
				if (optfct == "optim")
					res <- optim(tet,obj1, ...)
				if (optfct == "nlminb")
					{
					res <- nlminb(tet,obj1, ...)
					res$value <- res$objective
					}
				if (optfct == "optimize")
					{
					res <- optimize(obj1,t0, ...)
					res$par <- res$minimum
					res$value <- res$objective
					}	
				}
			ch <- crossprod(abs(tet-res$par)/tet,abs(tet-res$par)/tet)
			if (j>itermax)
				{
				cat("No convergence after ",itermax," iterations")
				ch <- crit
				}
			j <- j+1	
		}

		z = list(coefficients=res$par,objective=res$value)	
		}
	}

	if (!is.function(gradv)) 
		G <- Gf(z$coefficients)
	else
		if (typeg)
			G <- gradv(x)
		else	
			G <- gradv(z$coefficients,x)

	if (vcov == "iid")
		v <- iid(z$coefficients)/n
	else
		v <- HAC(g(z$coefficients,x), kernel=kernel, bw=bw,prewhite=prewhite,ar.method=ar.method,approx=approx,tol=tol)/n
	
	if (wmatrix == "optimal")
		{
		z$vcov <- solve(crossprod(G,solve(v,G)))
		}
	else
		{
		GGG <- solve(crossprod(G),t(G))
		z$vcov <- GGG%*%v%*%t(GGG)
		}

z$gt <- g(z$coefficients,x)
if (typeg==0)
	names(z$coefficients) <- paste("Theta[",1:k,"]",sep="")
if (typeg == 1)
	{
	namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
	nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
	if (dat$ny > 1)
		{
		namey <- colnames(dat$x[,1:dat$ny])
		names(z$coefficients) <- paste(rep(namey,dat$k),"_",rep(namex,rep(dat$ny,dat$k)),sep="")
		colnames(z$gt) <- paste(rep(namey,dat$nh),"_",rep(nameh,rep(dat$ny,dat$nh)),sep="")
		}
	if (dat$ny == 1)
		{
		names(z$coefficients) <- namex
		colnames(z$gt) <- nameh
		}
	}

dimnames(z$vcov) <- list(names(z$coefficients),names(z$coefficients))
z$df <- df
z$k <- k
z$n <- n
if(typeg==1)
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
if(typeg==0)
	if(X) z$x <- x

z$call<-match.call()
z$met <- type
z$kernel <- kernel
z$coefficients <- c(z$coefficients)
class(z) <- "gmm"
z
}

formula.gmm <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gmm object was not created by a formula")
    else
	formula(x$terms)
}
