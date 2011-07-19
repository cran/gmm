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

gmm <- function(g,x,t0=NULL,gradv=NULL, type=c("twoStep","cue","iterative"), wmatrix = c("optimal","ident"),  vcov=c("HAC","iid","TrueFixed"), 
	      kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews, 
	      prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7, itermax=100,optfct=c("optim","optimize","nlminb"),
	      model=TRUE, X=FALSE, Y=FALSE, TypeGmm = "baseGmm", centeredVcov = TRUE, weightsMatrix = NULL, data, ...)
{

type <- match.arg(type)
kernel <- match.arg(kernel)
vcov <- match.arg(vcov)
wmatrix <- match.arg(wmatrix)
optfct <- match.arg(optfct)

if(vcov=="TrueFixed" & is.null(weightsMatrix))
	stop("TrueFixed vcov only for fixed weighting matrix")
if(!is.null(weightsMatrix))
	wmatrix <- "optimal"

if(missing(data))
	data<-NULL
all_args<-list(data = data, g = g, x = x, t0 = t0, gradv = gradv, type = type, wmatrix = wmatrix, vcov = vcov, kernel = kernel,
                   crit = crit, bw = bw, prewhite = prewhite, ar.method = ar.method, approx = approx, 
                   weightsMatrix = weightsMatrix, centeredVcov = centeredVcov,
                   tol = tol, itermax = itermax, optfct = optfct, model = model, X = X, Y = Y, call = match.call())
class(all_args)<-TypeGmm
Model_info<-getModel(all_args)
z <- momentEstim(Model_info, ...)

z <- FinRes(z, Model_info)
z
}

getDat <- function (formula, h, data) 
{
	cl <- match.call()
	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula", "data"), names(mf), 0L)
	mf <- mf[c(1L, m)]
	mf$drop.unused.levels <- TRUE
	mf$na.action <- "na.pass"	
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	mt <- attr(mf, "terms")

	if (inherits(h,'formula'))
		{
		mfh <- match.call(expand.dots = FALSE)
		mh <- match(c("h", "data"), names(mfh), 0L)
		mfh <- mfh[c(1L, mh)]
		mfh$formula <- mfh$h
		mfh$h <- NULL
		mfh$drop.unused.levels <- TRUE
		mfh$na.action <- "na.pass"
		mfh[[1L]] <- as.name("model.frame")
		mfh <- eval(mfh, parent.frame())
		mth <- attr(mfh, "terms")
		h <- as.matrix(model.matrix(mth, mfh, NULL))
		}
	else
		{		
		if (!is.matrix(h))
			h <- cbind(rep(1,length(h)),h)
		else	
			h <- cbind(rep(1,nrow(h)),h)
		h <- as.matrix(h)	
		if(is.null(colnames(h)))
			colnames(h) <- c("h.(Intercept)",paste("h",1:(ncol(h)-1),sep=""))
		else
			attr(h,'dimnames')[[2]][1] <- "h.(Intercept)"
		if (attr(mt,"intercept")==0)
			{
			h <- as.matrix(h[,2:ncol(h)])
			}
		}
	y <- as.matrix(model.response(mf, "numeric"))
	xt <- as.matrix(model.matrix(mt, mf, NULL))
	ny <- ncol(y)
	k <- ncol(xt)
	nh <- ncol(h)

	if (nrow(y) != nrow(xt) | nrow(xt) != nrow(h) | nrow(y)!=nrow(h))
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
	rownames(xt) <- rownames(y)
	rownames(h) <- rownames(y)
	x <- cbind(y,xt,h)
	if(any(is.na(x)))
		{
		warning("There are missing values. Associated observations have been removed")
		x <- na.omit(x)
		if (nrow(x)<=k)
			stop("The number of observations must be greater than the number of coefficients")		
		}
	colnames(x)<-c(colnames(y),colnames(xt),colnames(h))
	return(list(x=x,nh=nh,ny=ny,k=k,mf=mf,mt=mt,cl=cl))
}


.tetlin <- function(x, w, ny, nh, k, gradv, g, type=NULL, inv=TRUE)
  {
  n <- nrow(x)
  ym <- as.matrix(x[,1:ny])
  xm <- as.matrix(x[,(ny+1):(ny+k)])
  hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
  if (!is.null(type))
  	{
  	if(type=="2sls")
	  	{
  		fsls <- lm(xm~hm-1)$fitted
  	     par <- lm(ym~fsls-1)$coef
		if (ny == 1)
		{
  	     	e2sls <- ym-xm%*%par
 	     	v2sls <- lm(e2sls~hm-1)$fitted
  	     	value <- sum(v2sls^2)/sum(e2sls^2)
  	     }
  	     else
  	     {
  	     	par <- c(t(par))	
  	     	g2sls <- g(par, x, ny, nh, k)
  	     	w <- crossprod(g2sls)/n
  	     	gb <- matrix(colMeans(g2sls), ncol = 1)
   			value <- crossprod(gb, solve(w, gb)) 
  	     }
	  	}
  	}
  else
  	{
  if (ny>1)
  	{
     if (inv) 
	{
	whx <- solve(w, (crossprod(hm, xm) %x% diag(ny)))
	wvecyh <- solve(w, matrix(crossprod(ym, hm), ncol = 1))	
	}
     else
        {
	whx <- w%*% (crossprod(hm, xm) %x% diag(ny))
	wvecyh <- w%*%matrix(crossprod(ym, hm), ncol = 1)
        }
     dg <- gradv(NULL,x, ny, nh, k)
     xx <- crossprod(dg, whx)
     par <- solve(xx, crossprod(dg, wvecyh))
     }
  else
  	{   
     if (nh>k)
     	{
	if(inv)
           xzwz <- crossprod(xm,hm)%*%solve(w,t(hm))	
	else
     	   xzwz <- crossprod(xm,hm)%*%w%*%t(hm)
     	par <- solve(xzwz%*%xm,xzwz%*%ym)	
	     }
	else
		par <- solve(crossprod(hm,xm),crossprod(hm,ym))  	}
	}
  gb <- matrix(colSums(g(par, x, ny, nh, k))/n, ncol = 1)
  if(inv)
	  value <- crossprod(gb, solve(w, gb)) 
  else
	  value <- crossprod(gb, w%*%gb) 

  res <- list(par = par, value = value)
  return(res)
  }


.obj1 <- function(thet, x, w, gf, INV = TRUE)
  {
  gt <- gf(thet, x)
  gbar <- as.vector(colMeans(gt))
  if (INV)		
  	obj <- crossprod(gbar, solve(w, gbar))
  else
	obj <- crossprod(gbar,w)%*%gbar
  return(obj)
  }

.Gf <- function(thet, x, g)
  {
  myenv <- new.env()
  assign('x', x, envir = myenv)
  assign('thet', thet, envir = myenv)
  barg <- function(x, thet)
    {
    gt <- g(thet, x)
    gbar <- as.vector(colMeans(gt))
    return(gbar)
    }
  G <- attr(numericDeriv(quote(barg(x, thet)), "thet", myenv), "gradient")
  return(G)
  }

.objCue <- function(thet, x, P)
  {
  gt <- P$g(thet,x)
  gbar <- as.vector(colMeans(gt))
  if (P$vcov == "iid")
    w2 <- P$iid(thet, x, P$g, P$centeredVcov)
  if (P$vcov == "HAC")
    {
    if(P$centeredVcov)
        gmat <- lm(P$g(thet,x)~1)
    else
      {
        gmat <- P$g(thet,x)
        class(gmat) <- "gmmFct"
      }
    w2 <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
            ar.method = P$ar.method, approx = P$approx, tol = P$tol, sandwich = FALSE)
   }
  obj <- crossprod(gbar,solve(w2,gbar))
  return(obj)
}	


