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

getModel <- function(object, ...)
  {
  UseMethod("getModel")
  }

getModel.baseGmm <- function(object, ...)
  {
  if(is(object$g, "formula"))
    {
    object$gradvf <- FALSE
    dat <- getDat(object$g, object$x)
    if(is.null(object$weightsMatrix))
      clname <- paste(class(object), ".", object$type, ".formula", sep = "")
    else
      {
      clname <- "fixedW.formula"
      object$type <- "One step GMM with fixed W"
      }
    object$gform<-object$g
    g <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] - x[,(ny+1):(ny+k)] %*% t(tet)
      gt <- e * x[, ny+k+1]
      if(nh > 1)
	for (i in 2:nh)	  gt <- cbind(gt, e*x[, (ny+k+i)])
      return(gt)
      }
    gradv <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k, g = NULL)
      {
      a <- g
      tet <- NULL
      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)]) %*% x[,(ny+1):(ny+k)]) %x% diag(rep(1,ny))/nrow(x)
      return(dgb)
      }
    object$g <- g
    }
  else
    {
    if(is.null(object$weightsMatrix))
      clname <- paste(class(object), "." ,object$type, sep = "")
    else
      clname <- "fixedW"
    if (!is.function(object$gradv))
      { 
      gradv <- .Gf
      object$gradvf <- FALSE
      }
    else
      {
      gradv <- object$gradv
      object$gradvf <- TRUE
      }
    }
	
  iid <- function(thet, x, g, centeredVcov)
    {
    gt <- g(thet,x)
    if(centeredVcov) gt <- residuals(lm(gt~1))
    n <- ifelse(is.null(nrow(x)), length(x), nrow(x))
    v <- crossprod(gt,gt)/n
    return(v)
    }
 
  object$iid<-iid
  object$TypeGmm <- class(object)
  object$gradv <- gradv	
 
  class(object)  <- clname
  return(object)
  }

getModel.baseGel <- function(object, ...)
  {
  P <- object
  if (P$type == "ETEL")
    {
    P$typel <- "ET"
    P$typet <- "EL"	
    }
  else
    {
    P$typel <- P$type
    P$typet <- P$type
    }
  if(P$optfct == "optim")
    P$k <- length(P$tet0)
  else
    P$k <- 1
  
  if (is(P$g, "formula"))
    {
    clname <- paste(class(P), ".modFormula", sep = "")
    dat <- getDat(P$g, P$x)
    x <- dat$x

    g <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      e <- x[,1:ny] -  x[, (ny+1):(ny+k)]%*%t(tet)
      gt <- e*x[, ny+k+1]
      if (nh > 1)
        {	
        for (i in 2:nh)
          {
          gt <- cbind(gt, e*x[,(ny+k+i)])
          }
        }
      return(gt)
      }

    gradv <- function(tet, x, ny = dat$ny, nh = dat$nh, k = dat$k)
      {
      tet <- matrix(tet, ncol = k)
      dgb <- -(t(x[,(ny+k+1):(ny+k+nh)])%*%x[,(ny+1):(ny+k)])%x%diag(rep(1, ny))/nrow(x)
      return(dgb)
      }
    P$dat <- dat
    P$gform <- P$g
    P$g <- g
    P$gradv <- gradv
    }	
  else
    {
    clname <- paste(class(P), ".mod", sep = "")
    P$gform <- NULL
    if (!is.function(object$gradv))
      { 
      P$gradv <- .Gf
      P$gradvf <- FALSE
      }
    else
      {
      P$gradvf <- TRUE
      }

    }

  if (P$smooth)
    {
    if(P$kernel == "Truncated")
        {
        P$wkernel <- "Bartlett"
        P$k1 <- 2
        P$k2 <- 2
        }
    if(P$kernel == "Bartlett")
        {
        P$wkernel <- "Parzen"
        P$k1 <- 1
        P$k2 <- 2/3
        }
    P$g1 <- P$g
    rgmm <- gmm(P$g, x, P$tet0, wmatrix = "ident")


    P$bwVal <- P$bw(lm(P$g(rgmm$coefficients, x)~1), kernel = P$wkernel, prewhite = P$prewhite, 
               ar.method = P$ar.method, approx = P$approx)
    P$w <- P$weights(lm(P$g(rgmm$coefficients, x)~1), kernel = P$kernel, bw = P$bwVal, prewhite = P$prewhite, 
               ar.method = P$ar.method, approx = P$approx, tol = P$tol_weights)

    P$g <- function(thet, x, g1 = P$g1, w = P$w)
      {
      gf <- g1(thet, x)
      gt <- smoothG(gf, weights = w)$smoothx
      return(gt)
      }
    }
  else
   {
   P$k1 <- 1
   P$k2 <- 1
   P$w <- 1
   P$bwVal <- 1
   }	
  class(P) <- clname
  return(P)
  }

