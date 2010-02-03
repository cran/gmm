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


FinRes <- function(z, object, ...)
  {
# object is computed by the getModel method #
  UseMethod("FinRes")
  }

FinRes.baseGmm.res <- function(z, object, ...)
  {
  P <- object
  if(!is.null(object$gform))
    {
    dat <- z$dat
    x <- dat$x
    }
  else
    x <- z$x

  n <- z$n
  gradv <- z$gradv
  iid <- z$iid 	

  if(P$gradvf)
    G <- gradv(z$coefficients, x)
  else
    G <- gradv(z$coefficients, x, g = object$g)

  if (P$vcov == "iid")
    v <- iid(z$coefficients, x, z$g, P$centeredVcov)/n
  else
    {
    if(P$centeredVcov) 
	gmat <- lm(z$gt~1)
    else
       {
       gmat <- z$gt
       class(gmat) <- "gmmFct"
       }
    v <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol, sandwich = FALSE)/n
    }

  z$vcov <- try(solve(crossprod(G, solve(v, G))), silent = TRUE)
  if(class(z$vcov) == "try-error")
    {
    z$vcov <- matrix(Inf,length(z$coef),length(z$coef))
    warning("The covariance matrix of the coefficients is singular")
    }

  dimnames(z$vcov) <- list(names(z$coefficients), names(z$coefficients))
  z$call <- P$call
  
  
  if(is.null(P$weightsMatrix))
    {
    if(P$wmatrix == "ident")
      z$w <- diag(ncol(z$gt))
    else
      {
      z$w <- try(solve(v), silent = TRUE)
      if(class(z$w) == "try-error")
         warning("The covariance matrix of the moment function is singular")
      }
    }
  else
    z$w <- P$weightsMatrix

  z$G <- G
  z$met <- P$type
  z$kernel <- P$kernel
  z$coefficients <- c(z$coefficients)
  class(z) <- "gmm"
  return(z)
  }


