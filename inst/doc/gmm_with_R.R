## ----echo=FALSE---------------------------------------------------------------
library(knitr)
opts_chunk$set(size='footnotesize', fig.height=5, out.width='70%')

## -----------------------------------------------------------------------------
library(gmm)

## -----------------------------------------------------------------------------
g1 <- function(tet,x)
        {
        m1 <- (tet[1]-x)
        m2 <- (tet[2]^2 - (x - tet[1])^2)
        m3 <- x^3-tet[1]*(tet[1]^2+3*tet[2]^2)
        f <- cbind(m1,m2,m3)
        return(f)
        }

## -----------------------------------------------------------------------------
Dg <- function(tet,x)
        {
        G <- matrix(c( 1,
                              2*(-tet[1]+mean(x)),
                              -3*tet[1]^2-3*tet[2]^2,0,
                              2*tet[2],-6*tet[1]*tet[2]),
                           nrow=3,ncol=2)
        return(G)
        }

## -----------------------------------------------------------------------------
set.seed(123)
n <- 200
x1 <- rnorm(n, mean = 4, sd = 2)

## -----------------------------------------------------------------------------
print(res <- gmm(g1,x1,c(mu = 0, sig = 0), grad = Dg))

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
specTest(res)

## -----------------------------------------------------------------------------
sim_ex <- function(n,iter)
    {
    tet1 <- matrix(0,iter,2)
    tet2 <- tet1 
    for(i in 1:iter)
      {
      x1 <- rnorm(n, mean = 4, sd = 2)
      tet1[i,1] <- mean(x1)
      tet1[i,2] <- sqrt(var(x1)*(n-1)/n)
      tet2[i,] <- gmm(g1,x1,c(0,0),grad=Dg)$coefficients
      }
    bias <- cbind(rowMeans(t(tet1)-c(4,2)),rowMeans(t(tet2)-c(4,2)))
    dimnames(bias)<-list(c("mu","sigma"),c("ML","GMM"))
    Var <- cbind(diag(var(tet1)),diag(var(tet2)))
    dimnames(Var)<-list(c("mu","sigma"),c("ML","GMM"))
    MSE <- cbind(rowMeans((t(tet1)-c(4,2))^2),rowMeans((t(tet2)-c(4,2))^2))
    dimnames(MSE)<-list(c("mu","sigma"),c("ML","GMM"))
    return(list(bias=bias,Variance=Var,MSE=MSE))
    }

## -----------------------------------------------------------------------------
g2 <- function(theta,x)
	{
	tau <- seq(1,5,length.out=10)
	pm <- 1
	x <- matrix(c(x),ncol=1)
	x_comp <- x%*%matrix(tau,nrow=1)
	x_comp <- matrix(complex(ima=x_comp),ncol=length(tau))
	emp_car <- exp(x_comp)
	the_car <- charStable(theta,tau,pm)
	gt <- t(t(emp_car) - the_car)
	gt <- cbind(Im(gt),Re(gt))
	return(gt)
	}

## -----------------------------------------------------------------------------
library(stabledist)
set.seed(345)
x2 <- rstable(500,1.5,.5,pm=1)
t0 <- c(alpha = 2, beta = 0, gamma = sd(x2)/sqrt(2), delta = 0)
print(res <- gmm(g2,x2,t0))

## -----------------------------------------------------------------------------
summary(res)

## -----------------------------------------------------------------------------
res2 <- gmm(g2,x2,t0,optfct="nlminb",lower=c(0,-1,0,-Inf),upper=c(2,1,Inf,Inf))
summary(res2)

## ----warning=FALSE------------------------------------------------------------
data(Finance)
x3 <- Finance[1:1500,"WMK"]
t0<-c(alpha = 1.8, beta = 0.1, gamma = sd(x3)/sqrt(2),delta = 0)
res3 <- gmm(g2,x3,t0,optfct="nlminb")
summary(res3)

## -----------------------------------------------------------------------------
library(car)
linearHypothesis(res3,cbind(diag(2),c(0,0),c(0,0)),c(2,0))

## -----------------------------------------------------------------------------
library(mvtnorm)
set.seed(112233)
sig <- matrix(c(1,.5,.5,1),2,2)
n <- 400
e <- rmvnorm(n,sigma=sig)
x4 <- rnorm(n)
w <- exp(-x4^2) + e[,1]
y <- 0.1*w + e[,2]

## -----------------------------------------------------------------------------
h <- cbind(x4, x4^2, x4^3)
g3 <- y~w

## -----------------------------------------------------------------------------
summary(res <- gmm(g3,x=h))

## -----------------------------------------------------------------------------
res2 <- gmm(g3,x=h,type='iterative',crit=1e-8,itermax=200)
coef(res2)

## -----------------------------------------------------------------------------
res3 <- gmm(g3,x=h,res2$coef,type='cue')
coef(res3)

## -----------------------------------------------------------------------------
confint(res3,level=.90)

## -----------------------------------------------------------------------------
plot(w,y,main="LS vs GMM estimation")
lines(w,fitted(res),col=2)
lines(w,fitted(lm(y~w)),col=3,lty=2)
lines(w,.1*w,col=4,lty=3)
legend("topleft",c("Data","Fitted GMM","Fitted LS","True line"),pch=c(1,NA,NA,NA),col=1:3,lty=c(NA,1,2,3))

## -----------------------------------------------------------------------------
t <- 400
set.seed(345)
x5 <- arima.sim(n=t,list(ar=c(1.4,-0.6),ma=c(0.6,-0.3)))
x5t<-cbind(x5)
for (i in 1:6) x5t<-cbind(x5t,lag(x5,-i))
x5t<-na.omit(x5t)
g4<-x5t[,1]~x5t[,2]+x5t[,3]
res<-gmm(g4,x5t[,4:7])
summary(res)

## -----------------------------------------------------------------------------
res2 <- gmm(g4,x=x5t[,4:7],kernel="Truncated")
coef(res2)
res3 <- gmm(g4,x=x5t[,4:7],kernel="Bartlett")
coef(res3)
res4 <- gmm(g4,x=x5t[,4:7],kernel="Parzen")
coef(res4)
res5<- gmm(g4,x=x5t[,4:7],kernel="Tukey-Hanning")
coef(res5)

## -----------------------------------------------------------------------------
diag(vcov(res2))^.5
diag(vcov(res3))^.5
diag(vcov(res4))^.5
diag(vcov(res5))^.5

## -----------------------------------------------------------------------------
plot(res,which=2)

## -----------------------------------------------------------------------------
plot(res,which=3)

## -----------------------------------------------------------------------------
data(Finance)
r <- Finance[1:500,1:5]
rm <- Finance[1:500,"rm"]
rf <- Finance[1:500,"rf"]
z <- as.matrix(r-rf)
zm <- as.matrix(rm-rf)
res <- gmm(z~zm,x=zm)
coef(res)
R <- cbind(diag(5),matrix(0,5,5))
c <- rep(0,5)
linearHypothesis(res,R,c,test = "Chisq")

## ----eval=FALSE---------------------------------------------------------------
# test <- paste(names(coef(res)[1:5])," = 0",sep="")
# linearHypothesis(res,test)

## -----------------------------------------------------------------------------
res2<-gmm(z~zm-1,cbind(1,zm))
specTest(res2)

## -----------------------------------------------------------------------------
g5 <- function(tet, x) {
     gmat <- (tet[1] + tet[2] * (1 + c(x[, 1]))) * (1 + x[, 2:6]) -  1
     return(gmat)
 }
res_sdf <- gmm(g5, x = as.matrix(cbind(rm, r)), c(0, 0))
specTest(res_sdf)

## -----------------------------------------------------------------------------
g6 <- function(theta, x) {
     t <- length(x)
     et1 <- diff(x) - theta[1] - theta[2] * x[-t]
     ht <- et1^2 - theta[3] * x[-t]^(2 * theta[4])
     g <- cbind(et1, et1 * x[-t], ht, ht * x[-t])
     return(g)
 }

## -----------------------------------------------------------------------------
rf <- Finance[,"rf"]
rf <- ((1 + rf/100)^(365) - 1) * 100
dr <- diff(rf)
res_0 <- lm(dr ~ rf[-length(rf)])
tet0 <- c(res_0$coef, var(residuals(res_0)), 0)
names(tet0) <- c("alpha", "beta", "sigma^2", "gamma")
res_rf <- gmm(g6, rf, tet0, control = list(maxit = 1000, reltol = 1e-10))
coef(res_rf)

## ----eval=FALSE---------------------------------------------------------------
# y <- rbind(y1-mean(y1),y2-mean(y2),y3-mean(y3))
# x <- rbind(x1-mean(x1),x2-mean(x2),x3-mean(x3))
# res <- gmm(y~x,h)

## ----eval=FALSE---------------------------------------------------------------
# y <- rbind(y1,y2,y3)
# x <- rbind(x1,x2,x3)
# res <- gmm(y~x,h)

## ----eval=FALSE---------------------------------------------------------------
# gt <- g(t0, x)
# V <- kernHAC(lm(gt~1),sandwich = FALSE)
# W <- solve(V)

## -----------------------------------------------------------------------------
print(res<-gmm(g4,x5t[,4:7],wmatrix="ident"))
diag(vcovHAC(res))^.5

## -----------------------------------------------------------------------------
diag(vcov(res))^.5

## -----------------------------------------------------------------------------
print(res<-gmm(g4,x5t[,4:7], weightsMatrix = diag(5)))

## -----------------------------------------------------------------------------
tet0 <- c(mu = mean(x1), sig = sd(x1))
res_el <- gel(g1,x1,tet0)
summary(res_el)

## -----------------------------------------------------------------------------
res_et <- gel(g1,x1,tet0,type="ET")
coef(res_et)

## -----------------------------------------------------------------------------
res_cue <- gel(g1,x1,tet0,type="CUE")
coef(res_cue)

## -----------------------------------------------------------------------------
res_etel <- gel(g1,x1,c(mu=1,sig=1),type="ETEL")
coef(res_etel)

## ----warning=FALSE------------------------------------------------------------
tet0 <- gmm(g4,x=x5t[,4:7],wmatrix="ident")$coef
res <- gel(g4,x=x5t[,4:7],tet0,smooth=TRUE,kernel="Truncated")
summary(res)

## -----------------------------------------------------------------------------
specTest(res)

## -----------------------------------------------------------------------------
plot(res,which=4)

## ----eval=FALSE---------------------------------------------------------------
# ui=cbind(0,-1,-1)
# ci <- -1

## ----eval=FALSE---------------------------------------------------------------
# res <- gel(g4,x=dat5[,4:7],tet0,smooth=TRUE,kernel="Truncated",
# constraint=TRUE, ui=ui,ci=ci)

## ----eval=FALSE---------------------------------------------------------------
# res <- gel(g4,x=dat5[,4:7],tet0,smooth=TRUE, optlam="optim")
# res <- gel(g4,x=dat5[,4:7],tet0,smooth=TRUE, optlam="optim",
# LambdaControl=list(trace=TRUE, parscale=rep(.1,5)))

