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

momentEstim <- function(object, ...)
  {
  UseMethod("momentEstim")
  }

momentEstim.baseGmm.twoStep <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    {
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    }
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	
  if (q == k2 | P$wmatrix == "ident")
    z = list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)	
  else
    {
    if (P$vcov == "iid")
      w <- P$iid(res$par, P$x, P$g, P$centeredVcov)

    if (P$vcov == "HAC")
     {
     if(P$centeredVcov) 
      	gmat <- lm(P$g(res$par, P$x)~1)
     else
       {
       gmat <- P$g(res$par, P$x)
       class(gmat) <- "gmmFct"
       }
      w <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol, sandwich = FALSE)

     }

    if (P$optfct == "optim")
      res2 <- optim(res$par, .obj1, x = P$x, w = w, gf = P$g, ...)

    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(res$par, .obj1, x = P$x, w = w, gf = P$g, ...)
      res2$value <- res2$objective
      }

    if (P$optfct == "optimize")
      {
      res2 <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }	

     z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df)	
    }

  if(is.null(names(P$t0)))
    names(z$coefficients) <- paste("Theta[" ,1:k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$t0)

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

momentEstim.baseGmm.twoStep.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  if (is.null(P$data))
    	dat <- getDat(P$gform, P$x)
    else
    	dat <- getDat(P$gform, P$x, P$data)
  
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
    
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, P$g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    if (P$vcov == "iid")
    	{
      res2 <- .tetlin(x, diag(q), dat$ny, dat$nh, dat$k, P$gradv, P$g, type="2sls")
      }
    if (P$vcov == "HAC")
      {
      w=diag(q)
      res1 <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, P$g)
      if(P$centeredVcov) 
       	 gmat <- lm(g(res1$par, x)~1)
      else
        {
        gmat <- g(res1$par, x)
        class(gmat) <- "gmmFct"
        }
      w <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, 
		ar.method = P$ar.method, approx = P$approx, tol = P$tol, sandwich = FALSE)
	 res2 <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
      }
    
    z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df)	
    }
  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.iterative.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  if (is.null(P$data))
    	dat <- getDat(P$gform, P$x)
    else
    	dat <- getDat(P$gform, P$x, P$data)
  
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    w=diag(rep(1, q))
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "iid")
        w <- P$iid(tet, x, g, P$centeredVcov)
      if (P$vcov == "HAC")
	{
        if (P$centeredVcov)
          gmat <- lm(g(tet, x)~1)
        else
          {
          gmat <- g(tet, x)
          class(gmat) <- "gmmFct"
          }
        w <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, ar.method = P$ar.method, 
                 approx = P$approx, tol = P$tol, sandwich = FALSE)
        }
      res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
      ch <- crossprod(abs(tet- res$par)/tet)^.5
      if (j>P$itermax)
        {
        cat("No convergence after ", P$itermax, " iterations")
        ch <- P$crit
        }
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df)	
   }
  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.iterative <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
    }	
  else
    {
    ch <- 100000
    j <- 1
    while(ch > P$crit)
      {
      tet <- res$par
      if (P$vcov == "iid")
        w <- P$iid(tet, P$x, P$g, P$centeredVcov)
      if (P$vcov == "HAC")
	{
        if (P$centeredVcov)
          gmat <- lm(P$g(tet, P$x)~1)
        else
          {
          gmat <- P$g(tet, P$x)
          class(gmat) <- "gmmFct"
          }
        w <- kernHAC(gmat, kernel = P$kernel, bw = P$bw, prewhite = P$prewhite, ar.method = P$ar.method, 
                 approx = P$approx, tol = P$tol, sandwich = FALSE)
        }

      if (P$optfct == "optim")
        res <- optim(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
      if (P$optfct == "nlminb")
        {
        res <- nlminb(tet, .obj1, x = P$x, w = w, gf = P$g, ...)
        res$value <- res$objective
        }
      if (P$optfct == "optimize")
        {
        res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
        res$par <- res$minimum
        res$value <- res$objective
        }	
        ch <- crossprod(abs(tet-res$par)/tet)^.5	
        if (j>P$itermax)
          {
          cat("No convergence after ", P$itermax, " iterations")
          ch <- P$crit
          }
        j <- j+1	
      }
    z = list(coefficients = res$par, objective = res$value,k=k, k2=k2, n=n, q=q, df=df)	
    }

  if(is.null(names(P$t0)))
    names(z$coefficients) <- paste("Theta[" ,1:k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$t0)

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

momentEstim.baseGmm.cue.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  if (is.null(P$data))
    	dat <- getDat(P$gform, P$x)
    else
    	dat <- getDat(P$gform, P$x, P$data)
  
  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  
  if (q == k2 | P$wmatrix == "ident")
    {
    w <- diag(q)
    res <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
    z = list(coefficients = res$par, objective = res$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }
  else
    {
    if (is.null(P$t0))
      P$t0 <- .tetlin(x,diag(q), dat$ny, dat$nh, dat$k, P$gradv, g)$par
    if (P$optfct == "optim")
      res2 <- optim(P$t0,.objCue, x = x, P = P, ...)
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0,.objCue, x = x, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = x, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients = res2$par, objective = res2$value, dat = dat, k = k, k2 = k2, n = n, q = q, df = df)
    }

  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.baseGmm.cue <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w=diag(q)
  if (P$optfct == "optim")
    res <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
  if (P$optfct == "nlminb")
    {
    res <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res$value <- res$objective
    }
  if (P$optfct == "optimize")
    {
    res <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res$par <- res$minimum
    res$value <- res$objective
    }	

  if (q == k2 | P$wmatrix == "ident")
    {
    z <- list(coefficients = res$par, objective = res$value, k=k, k2=k2, n=n, q=q, df=df)
    }	
  else
    {
    if (P$optfct == "optim")
      res2 <- optim(P$t0, .objCue, x = x, P = P, ...)
    if (P$optfct == "nlminb")
      {
      res2 <- nlminb(P$t0, .objCue, x = x, P = P, ...)
      res2$value <- res2$objective
      }
    if (P$optfct == "optimize")
      {
      res2 <- optimize(.objCue,P$t0, x = x, P = P, ...)
      res2$par <- res2$minimum
      res2$value <- res2$objective
      }
    z = list(coefficients=res2$par,objective=res2$value, k=k, k2=k2, n=n, q=q, df=df)	
    }

  if(is.null(names(P$t0)))
    names(z$coefficients) <- paste("Theta[" ,1:k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$t0)

  z$x <- P$x
  z$gradv <- P$gradv
  z$gt <- P$g(z$coefficients, P$x)
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm, ".res", sep = "")	
  return(z)
  }

momentEstim.baseGel.modFormula <- function(object, ...)
  {
  P <- object
  g <- P$g
  dat <- getDat(P$gform, P$x)
  x <- dat$x
  n <- nrow(x)

  if (!P$constraint)
    {
    if (P$optfct == "optim")
      res <- optim(P$tet0, .thetf, P = P, ...)
    if (P$optfct == "nlminb")
      res <- nlminb(P$tet0, .thetf, P = P, ...)
	
    if (P$optfct == "optimize")
      { 
      res <- optimize(.thetf, P$tet0, P = P, ...)
      res$par <- res$minimum
      res$convergence <- "There is no convergence code for optimize"
      }
    }

  if(P$constraint)
    res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P, ...)


  if (P$optlam == "iter")
    {
    rlamb <- getLamb(P$g, res$par, x, type = P$typel, tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, tol_obj = P$tol_obj, k = P$k1/P$k2)
    z <- list(coefficients = res$par, lambda = rlamb$lam, conv_lambda = rlamb$conv_mes, conv_par = res$convergence)
    z$foc_lambda <- rlamb$obj
    }

  if (P$optlam == "numeric")
    {
    gt <- P$g(res$par, x)
    rhofct <- function(lamb)
      {
      rhof <- -sum(rho(gt, lamb, type = P$typel, k = P$k1/P$k2)$rhomat)
      return(rhof)
      }
    rlamb <- optim(rep(0, ncol(gt)), rhofct, control = list(maxit = 1000))
    z <- list(coefficients = res$par, conv_par = res$convergence, lambda = rlamb$par)
    z$conv_lambda = paste("Lambda by optim. Conv. code = ", rlamb$convergence, sep = "")
    rho1 <- as.numeric(rho(gt, z$lambda, derive = 1, type = P$typel, k = P$k1/P$k2)$rhomat)
    z$foc_lambda <- crossprod(colMeans(rho1*gt))
    }
  z$type <- P$type
  z$gt <- P$g(z$coefficients, x)
  rhom <- rho(z$gt, z$lambda, type = P$typet, k = P$k1/P$k2)
  z$pt <- -rho(z$gt, z$lambda, type = P$typet, derive = 1, k = P$k1/P$k2)$rhomat/n
  z$conv_moment <- colSums(as.numeric(z$pt)*z$gt)
  z$conv_pt <- sum(as.numeric(z$pt))
  z$objective <- sum(as.numeric(rhom$rhomat) - rho(1, 0, type = P$typet)$rhomat)/n

  namex <- colnames(dat$x[, (dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[, (dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])

  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[, 1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny,dat$nh)), sep = "")
    names(z$lambda) <- paste("Lam(",rep(namey,dat$nh), "_", rep(nameh, rep(dat$ny,dat$nh)), ")", sep = "")
    }
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    names(z$lambda) <- nameh
    }

  if (P$type == "EL")	
    {
    z$badrho <- rhom$ch
    names(z$badrho) <- "Number_of_bad_rho"
    }

  G <- P$gradv(z$coefficients, x)

  khat <- crossprod(z$gt)/(n*P$k2)*P$bwVal
  G <- G/P$k1 

  kg <- solve(khat, G)
  z$vcov_par <- solve(crossprod(G, kg))/n
  p_temp <- solve(khat,G)
  z$vcov_lambda <- solve(khat, ( diag(ncol(khat)) - G %*% (z$vcov_par*n) %*% t(p_temp) ))/n*P$bwVal^2

  z$weights <- P$w
  z$bwVal <- P$bwVal
  names(z$bwVal) <- "Bandwidth"

  dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
  dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x%*%b
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$call <- P$call
  z$k1 <- P$k1
  z$k2 <- P$k2
  z$khat <- khat

  class(z) <- paste(P$TypeGel, ".res", sep = "")
  return(z)
  }

momentEstim.baseGel.mod <- function(object, ...)
  {
  P <- object
  x <- P$x
  n <- ifelse(is.null(dim(x)),length(x),nrow(x))
  if (!P$constraint)
    {
    if (P$optfct == "optim")
      res <- optim(P$tet0, .thetf, P = P, ...)
    if (P$optfct == "nlminb")
      res <- nlminb(P$tet0, .thetf, P = P, ...)
	
    if (P$optfct == "optimize")
      { 
      res <- optimize(.thetf, P$tet0, P = P, ...)
      res$par <- res$minimum
      res$convergence <- "There is no convergence code for optimize"
      }
    }

  if(P$constraint)
    res <- constrOptim(P$tet0, .thetf, grad = NULL, P = P, ...)


  if (P$optlam == "iter")
    {
    rlamb <- getLamb(P$g, res$par, x, type = P$typel, tol_lam = P$tol_lam, maxiterlam = P$maxiterlam, tol_obj = P$tol_obj, k = P$k1/P$k2)
    z <- list(coefficients = res$par, lambda = rlamb$lam, conv_lambda = rlamb$conv_mes, conv_par = res$convergence)
    z$foc_lambda <- rlamb$obj
    }

  if (P$optlam == "numeric")
    {
    gt <- P$g(res$par, x)
    rhofct <- function(lamb)
      {
      rhof <- -sum(rho(gt, lamb, type = P$typel, k = P$k1/P$k2)$rhomat)
      return(rhof)
      }
    rlamb <- optim(rep(0, ncol(gt)), rhofct, control = list(maxit = 1000))
    z <- list(coefficients = res$par, conv_par = res$convergence, lambda = rlamb$par)
    z$conv_lambda = paste("Lambda by optim. Conv. code = ", rlamb$convergence, sep = "")
    rho1 <- as.numeric(rho(gt, z$lambda, derive = 1, type = P$typel, k = P$k1/P$k2)$rhomat)
    z$foc_lambda <- crossprod(colMeans(rho1*gt))
    }
  z$type <- P$type
  z$gt <- P$g(z$coefficients, x)
  rhom <- rho(z$gt, z$lambda, type = P$typet, k = P$k1/P$k2)
  z$pt <- -rho(z$gt, z$lambda, type = P$typet, derive = 1, k = P$k1/P$k2)$rhomat/n
  z$conv_moment <- colSums(as.numeric(z$pt)*z$gt)
  z$conv_pt <- sum(as.numeric(z$pt))
  z$objective <- sum(as.numeric(rhom$rhomat) - rho(1, 0, type = P$typet, k = P$k1/P$k2)$rhomat)/n

  if (P$type == "EL")	 
    {
    z$badrho <- rhom$ch
    names(z$badrho) <- "Number_of_bad_rho"
    }

  if(P$gradvf)
    G <- P$gradv(z$coefficients, x)
  else
    G <- P$gradv(z$coefficients, x, g = P$g)
  
  khat <- crossprod(z$gt)/(n*P$k2)*P$bwVal
  G <- G/P$k1 

  kg <- solve(khat, G)
  z$vcov_par <- solve(crossprod(G, kg))/n
  p_temp <- solve(khat,G)
  z$vcov_lambda <- solve(khat, ( diag(ncol(khat)) - G %*% (z$vcov_par*n) %*% t(p_temp) ))/n*P$bwVal^2
	
  z$weights <- P$w
  z$bwVal <- P$bwVal
  names(z$bwVal) <- "Bandwidth"
 
  if(is.null(names(P$tet0)))
    names(z$coefficients) <- paste("Theta[" ,1:P$k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$tet0)

  colnames(z$gt) <- paste("gt[",1:ncol(z$gt),"]", sep = "")
  names(z$lambda) <- paste("Lambda[",1:ncol(z$gt),"]", sep = "")
  dimnames(z$vcov_par) <- list(names(z$coefficients), names(z$coefficients))
  dimnames(z$vcov_lambda) <- list(names(z$lambda), names(z$lambda))
  if(P$X) z$x <- x
  z$call <- P$call
  z$k1 <- P$k1
  z$k2 <- P$k2
  z$khat <- khat

  class(z) <- paste(P$TypeGel, ".res", sep = "")
  return(z)
  }

momentEstim.fixedW.formula <- function(object, ...)
  {
  P <- object
  g <- P$g
  if (is.null(P$data))
    	dat <- getDat(P$gform, P$x)
    else
    	dat <- getDat(P$gform, P$x, P$data)

  x <- dat$x
  k <- dat$k
  k2 <- k*dat$ny
  n <- nrow(x)
  q <- dat$ny*dat$nh
  df <- q-k*dat$ny
  w <- P$weightsMatrix
  if(!all(dim(w) == c(q,q)))
    stop("The matrix of weights must be qxq")
  if(!is.real(eigen(w)$values))
    stop("The matrix of weights must be strictly positive definite")
  if(is.real(eigen(w)$values))
    {
    if(sum(eigen(w)$values<=0)!=0)
      stop("The matrix of weights must be strictly positive definite")
    }
  
  res2 <- .tetlin(x, w, dat$ny, dat$nh, dat$k, P$gradv, g)
  z = list(coefficients = res2$par, objective = res2$value, dat=dat, k=k, k2=k2, n=n, q=q, df=df)	

  z$gt <- g(z$coefficients, x) 
  b <- z$coefficients
  y <- as.matrix(model.response(dat$mf, "numeric"))
  ny <- dat$ny
  b <- t(matrix(b, nrow = dat$ny))
  x <- as.matrix(model.matrix(dat$mt, dat$mf, NULL))
  yhat <- x %*% b
  z$dat <- dat 
  z$fitted.values <- yhat	
  z$residuals <- y - yhat	
  z$terms <- dat$mt
  if(P$model) z$model <- dat$mf
  if(P$X) z$x <- x
  if(P$Y) z$y <- y
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
  
  namex <- colnames(dat$x[,(dat$ny+1):(dat$ny+dat$k)])
  nameh <- colnames(dat$x[,(dat$ny+dat$k+1):(dat$ny+dat$k+dat$nh)])
 
  if (dat$ny > 1)
    {
    namey <- colnames(dat$x[,1:dat$ny])
    names(z$coefficients) <- paste(rep(namey, dat$k), "_", rep(namex, rep(dat$ny, dat$k)), sep = "")
    colnames(z$gt) <- paste(rep(namey, dat$nh), "_", rep(nameh, rep(dat$ny, dat$nh)), sep = "")
    }
 
  if (dat$ny == 1)
    {
    names(z$coefficients) <- namex
    colnames(z$gt) <- nameh
    }
  class(z) <- paste(P$TypeGmm,".res",sep="")
  return(z)	
  }

momentEstim.fixedW <- function(object, ...)
  {
  P <- object
  x <- P$x
  if (P$optfct == "optimize")
    {
    n = nrow(P$g(P$t0[1], x))
    q = ncol(P$g(P$t0[1], x))
    k = 1
    }
  else
    {
    n = nrow(P$g(P$t0, x))
    q = ncol(P$g(P$t0, x))
    k = length(P$t0)
    }
  k2 <- k
  df <- q - k
  w <- P$weightsMatrix
  if(!all(dim(w) == c(q,q)))
    stop("The matrix of weights must be qxq")
  if(!is.real(eigen(w)$values))
    stop("The matrix of weights must be strictly positive definite")
  if(is.real(eigen(w)$values))
    {
    if(sum(eigen(w)$values<=0)!=0)
      stop("The matrix of weights must be strictly positive definite")
    }

  if (P$optfct == "optim")
    res2 <- optim(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)

  if (P$optfct == "nlminb")
    {
    res2 <- nlminb(P$t0, .obj1, x = P$x, w = w, gf = P$g, ...)
    res2$value <- res2$objective
    }
  if (P$optfct == "optimize")
    {
    res2 <- optimize(.obj1, P$t0, x = P$x, w = w, gf = P$g, ...)
    res2$par <- res2$minimum
    res2$value <- res2$objective
    }	
  z = list(coefficients = res2$par, objective = res2$value, k=k, k2=k2, n=n, q=q, df=df)	
  
  if(is.null(names(P$t0)))
    names(z$coefficients) <- paste("Theta[" ,1:k, "]", sep = "")
  else
    names(z$coefficients) <- names(P$t0)

  z$x <- P$x
  z$gt <- P$g(z$coefficients, P$x)
  z$gradv <- P$gradv
  z$iid <- P$iid
  z$g <- P$g
 
  class(z) <- paste(P$TypeGmm,".res",sep="")	
  return(z)
  }

