\name{Finance}
\docType{data}
\alias{Finance}
\title{Returns on selected stocks}
\description{
Returns on selected stocks, the Market portfolio and factors of Fama and French for CAPM and APT analysis 
}
\usage{data(Finance)}
\format{
A list of six daily time series from 1993-01-05 to 2009-01-30 in matrix format. Dates are reported as rownames()
 \describe{
    \item{r}{Daily returns on these following 10 stocks: (WMK, UIS, ORB, MAT, ABAX, T, EMR, JCS, VOXX, ZOOM),}
    \item{r2}{Daily returns on these following 10 stocks: (TDW, ROG, GGG, PC, GCO, EBF, F, FNM, NHP, AA),}
    \item{rf}{Risk-free rate,}
    \item{rm}{Return of the market portfolio of Fama-French,}
    \item{hml}{Factor High-Minus-Low of Fama-French,}
    \item{smb}{Factor Small-Minus-Big of Fama-French.}
 }
}

\source{\url{http://ca.finance.yahoo.com/} and \url{http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/}}

\examples{

data(Finance)

## CAPM test with GMM
r <- Finance$r
rm <- Finance$rm
rf <- Finance$rf

z <- r-rf
t <- nrow(z)
zm <- rm-rf
h <- matrix(zm,t,1)
res <- gmm(z~zm,x=h)
R <- cbind(diag(10),matrix(0,10,10))
c <- matrix(0,10,1)
lintest(res,R,c)

## APT test with Fama-French factors and GMM
r <- Finance$r
rm <- Finance$rm
rf <- Finance$rf

z <- r-rf
f1 <- rm-rf
f2 <- Finance$hml-rf
f3 <- Finance$smb-rf
h <- cbind(f1,f2,f3)
res2 <- gmm(z~f1+f2+f3,x=h)
R <- cbind(diag(10),matrix(0,10,30))
c <- matrix(0,10,1)
lintest(res2,R,c)
}

\keyword{datasets}
