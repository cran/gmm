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
	      kernel=c("Quadratic Spectral","Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),crit=10e-7,bw = bwAndrews, 
	      prewhite = FALSE, ar.method = "ols", approx="AR(1)",tol = 1e-7, itermax=100,optfct=c("optim","optimize","nlminb"),
	      model=TRUE, X=FALSE, Y=FALSE, TypeGmm = "baseGmm", centeredVcov = TRUE, weightsMatrix = NULL, ...)
{

type <- match.arg(type)
kernel <- match.arg(kernel)
vcov <- match.arg(vcov)
wmatrix <- match.arg(wmatrix)
optfct <- match.arg(optfct)
all_args<-list(g = g, x = x, t0 = t0, gradv = gradv, type = type, wmatrix = wmatrix, vcov = vcov, kernel = kernel,
                   crit = crit, bw = bw, prewhite = prewhite, ar.method = ar.method, approx = approx, 
                   weightsMatrix = weightsMatrix, centeredVcov = centeredVcov,
                   tol = tol, itermax = itermax, optfct = optfct, model = model, X = X, Y = Y, call = match.call())
class(all_args)<-TypeGmm
Model_info<-getModel(all_args)
z <- momentEstim(Model_info, ...)
z <- FinRes(z, Model_info)
z
}

getDat <- function (formula,h) 
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
	colnames(h) <- c("h.(Intercept)",paste("h",1:(ncol(h)-1),sep=""))
	y <- as.matrix(model.response(mf, "numeric"))
	xt <- as.matrix(model.matrix(mt, mf, NULL))
	if (attr(mt,"intercept")==0)
		{
		h <- as.matrix(h[,2:ncol(h)])
		}
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
	x <- cbind(y,xt,h)
	colnames(x)<-c(colnames(y),colnames(xt),colnames(h))
	return(list(x=x,nh=nh,ny=ny,k=k,mf=mf,mt=mt,cl=cl))
}


.tetlin <- function(x, w, ny, nh, k, gradv, g)
  {
  n <- nrow(x)
  ym <- as.matrix(x[,1:ny])
  xm <- as.matrix(x[,(ny+1):(ny+k)])
  hm <- as.matrix(x[,(ny+k+1):(ny+k+nh)])
  whx <- solve(w, (crossprod(hm, xm) %x% diag(ny)))
  wvecyh <- solve(w, matrix(crossprod(ym, hm), ncol = 1))
  dg <- gradv(NULL,x, ny, nh, k)
  xx <- crossprod(dg, whx)
  par <- solve(xx, crossprod(dg, wvecyh))
  gb <- matrix(colSums(g(par, x, ny, nh, k))/n, ncol = 1)
  value <- crossprod(gb, solve(w, gb)) 
  res <- list(par = par, value = value)
  return(res)
  }


.obj1 <- function(thet, x, w, gf)
  {
  gt <- gf(thet, x)
  gbar <- as.vector(colMeans(gt))
  obj <- crossprod(gbar, solve(w, gbar))
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

