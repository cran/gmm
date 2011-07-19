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


summary.gmm <- function(object, ...)
	{
	z <- object
	se <- sqrt(diag(z$vcov))
	par <- z$coefficients
	tval <- par/se
	ans <- list(met=z$met,kernel=z$kernel,algo=z$algo,call=z$call)
	names(ans$met) <- "GMM method"
	names(ans$kernel) <- "kernel for cov matrix"
		
	ans$coefficients <- cbind(par,se, tval, 2 * pnorm(abs(tval), lower.tail = FALSE))
    	dimnames(ans$coefficients) <- list(names(z$coefficients), 
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
	ans$stest <- specTest(z)
        ans$algoInfo <- z$algoInfo
	class(ans) <- "summary.gmm"
	ans
	}

print.summary.gmm <- function(x, digits = 5, ...)
	{
	cat("\nCall:\n")
	cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
	cat("\nMethod: ", x$met,"\n\n")
	cat("Kernel: ", x$kernel,"\n\n")
	cat("Coefficients:\n")
	print.default(format(x$coefficients, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	cat(x$stest$ntest,"\n")
	print.default(format(x$stest$test, digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	if(!is.null(x$algoInfo))
		{	
		cat("#############\n")
 		cat("Information related to the numerical optimization\n")
		}
	if(!is.null(x$algoInfo$convergence))
		cat("Convergence code = ", x$algoInfo$convergence,"\n")
	if(!is.null(x$algoInfo$counts))
		{	
		cat("Function eval. = ",x$algoInfo$counts[1],"\n")
		cat("Gradian eval. = ",x$algoInfo$counts[2],"\n")
		}	
	if(!is.null(x$algoInfo$message))
		cat("Message: ",x$algoInfo$message,"\n")
	invisible(x)
	}


formula.gmm <- function(x, ...)
{
    if(is.null(x$terms))
	stop("The gmm object was not created by a formula")
    else
	formula(x$terms)
}

confint.gmm <- function(object, parm, level=0.95, ...)
		{
		z <- object
		se <- sqrt(diag(z$vcov))
		par <- z$coefficients
			
		zs <- qnorm((1-level)/2,lower.tail=FALSE)
		ch <- zs*se
		ans <- cbind(par-ch,par+ch)
		dimnames(ans) <- list(names(par),c((1-level)/2,0.5+level/2))
		if(!missing(parm))
			ans <- ans[parm,]
		ans
		}
		
residuals.gmm <- function(object,...) 
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$residuals
	}

fitted.gmm <- function(object,...)
	{
	if(is.null(object$model))
		stop("The residuals method is valid only for g=formula")
	object$fitted.value
	}

print.gmm <- function(x, digits=5, ...)
	{
	cat("Method\n", x$met,"\n\n")
	cat("Objective function value: ",x$objective,"\n\n")
	print.default(format(coef(x), digits=digits),
                      print.gap = 2, quote = FALSE)
	cat("\n")
	if(!is.null(x$algoInfo$convergence))
		cat("Convergence code = ", x$algoInfo$convergence,"\n")
	invisible(x)
	}

coef.gmm <- function(object,...) object$coefficients

vcov.gmm <- function(object,...) object$vcov

estfun.gmmFct <- function(x, y = NULL, theta = NULL, ...)
	{
	if (is(x, "function"))
		{
		gmat <- x(y, theta)
		return(gmat)
		}
	else
		return(x)
	}

estfun.gmm <- function(x, ...)
  {
  foc <- x$gt %*% x$w %*% x$G
  return(foc)
  }

bread.gmm <- function(x, ...)
  {
  GWG <- crossprod(x$G, x$w %*% x$G)
  b <- try(solve(GWG), silent = TRUE)
  if (class(b) == "try-error")
    stop("The bread matrix is singular")
  return(b)
  }





		


