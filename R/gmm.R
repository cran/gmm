gmm <- function(g,t0,x,grad=NULL,type=c("twoStep","cue","iterative"),kernel=c("Quadratic Spectral", 
    "Truncated", "Bartlett", "Parzen", "Tukey-Hanning"),iid=FALSE,crit=10e-7,itermax=100,algo=c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"), vcov=HAC, ...)
{
type <- match.arg(type)
kernel <- match.arg(kernel)
algo <- match.arg(algo)

n = nrow(g(t0,x))
q = ncol(g(t0,x))
k = length(t0)
df <- q-k

obj1 <- function(thet)
	{
	gt <- g(thet,x)
	gbar <- as.vector(colMeans(gt))
	obj <- crossprod(gbar,solve(w,gbar))
	return(obj)
	}
iidf <- function(thet)
	{
	gt <- g(thet,x)
	v <- crossprod(gt,gt)/n
	return(v)
	}
Gf <- function(thet)
	{
	warning("It is strongly suggested to provide the gradient function of g_bar") 
	myenv <- new.env()
	assign('x',x,envir=myenv)
	assign('thet',res$par,envir=myenv)
	barg <- function(x,thet)
		{
		gt <- g(thet,x)
		gbar <- as.vector(colMeans(gt))
		return(gbar)
		}
	G <- attr(numericDeriv(quote(barg(x,thet)),"thet",myenv),"gradient")
	return(G)
	}

if (q == k)
	{
	w=diag(rep(1,q))
	res <- optim(t0,obj1,method=algo)
	z = list(par=res$par,objective=res$value)	
	}
else
	{
	if (type=="twoStep")
		{
		w=diag(rep(1,q))
		res1<-optim(t0,obj1,method=algo)
		if (iid)
			w <- iidf(res1$par)
		else
			w <- vcov(g(res1$par,x), kernel=kernel, ...)

		res2<-optim(res1$par,obj1,method=algo)
		z = list(par=res2$par,objective=res2$value)	
		}
	if (type=="cue")
		{
		obj_cue <- function(thet)
			{
			gt <- g(thet,x)
			gbar <- as.vector(colMeans(gt))
			if (iid)
				w2 <- iidf(thet)
			else
				w2 <- vcov(g(thet,x), kernel=kernel, ...)
			obj <- crossprod(gbar,solve(w2,gbar))
			return(obj)
			}		
		w=diag(rep(1,q))
		res1<-optim(t0,obj1,method=algo)
		res2<-optim(res1$par,obj_cue,method=algo)
		z = list(par=res2$par,objective=res2$value)	
		}
	if (type=="iterative")
		{
		w=diag(rep(1,q))
		res<-optim(t0,obj1,method=algo)
		ch <- 100000
		j <- 1
		while(ch>crit)
		{
			tet <- res$par
			if (iid)
				w <- iidf(tet)
			else
				w <- vcov(g(tet,x), kernel=kernel, ...)
	
			res<-optim(tet,obj1,method=algo)
			ch <- crossprod(abs(tet-res$par)/tet,abs(tet-res$par)/tet)
			if (j>itermax)
				{
				cat("No convergence after ",itermax," iterations")
				ch <- crit
				}
			j <- j+1	
		}
		z = list(par=res$par,objective=res$value)	
		}
	if (!is.function(grad)) 
		G <- Gf(z$par)
	else
		G <- grad(z$par,x)
	
	if (iid)
		v <- iidf(z$par)/n
	else
		v <- vcov(g(z$par,x), kernel=kernel, ...)/n

	z$vcov <- solve(crossprod(G,solve(v,G)))
	}
names(z$par) <- paste("Theta[",1:k,"]",sep="")
z$df <- df
z$k <- k
z$n <- n
z$met <- type
z$kernel <- kernel
z$algo=algo

class(z) <- "gmm"
z
}
