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
    v <- iid(z$coefficients, x, z$g)/n
  else
    {
    v <- HAC(z$gt, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol)/n
    }
  if (P$wmatrix == "optimal")
    {
    z$vcov <- solve(crossprod(G, solve(v, G)))
    }
  else
    {
    GGG <- solve(crossprod(G), t(G))
    z$vcov <- GGG %*% v %*% t(GGG)
    }
  dimnames(z$vcov) <- list(names(z$coefficients), names(z$coefficients))
  z$call <- P$call
  z$met <- P$type
  z$kernel <- P$kernel
  z$coefficients <- c(z$coefficients)
  class(z) <- "gmm"
  return(z)
  }


