## ECON675: ASSIGNMENT 2
## Erin Markiewitz
## 10/10/2018
######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
library(Matrix)             #fast matrix calcs
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Q3 (a): data generation, ploynomial basis
######################################################################
d   = 5
N   = 500
M   = 1000

DGP    = function(n=N){
  x      = t(as.matrix(replicate(n,runif(d,-1,1))))
  v      = rnorm(n)
  x.norm = sapply(1:n,function(i) t(x[i,])%*%x[i,])
  e      = 0.3637899*(1+x.norm)*v
  g0.x   =exp(x.norm)
  u      = rnorm(n)
  tt     = as.numeric((sqrt(x.norm)+u)>1)
  y      = tt + g0.x + e
  
  return(list(y=y, x=x, tt=tt))
}

# generate the polynomial basis
gen.P = function(Z,K) {
  if (K==0)   out = NULL;
  if (K==1)   out = poly(Z,degree=1,raw=TRUE);
  if (K==2)  {out = poly(Z,degree=1,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^2);}
  if (K==2.5) out = poly(Z,degree=2,raw=TRUE);
  if (K==3)  {out = poly(Z,degree=2,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^3);}
  if (K==3.5) out = poly(Z,degree=3,raw=TRUE);
  if (K==4)  {out = poly(Z,degree=3,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^4);}
  if (K==4.5) out = poly(Z,degree=4,raw=TRUE);
  if (K==5)  {out = poly(Z,degree=4,raw=TRUE); for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^5);}
  if (K==5.5) out = poly(Z,degree=5,raw=TRUE);
  if (K>=6)  {out = poly(Z,degree=5,raw=TRUE); for (k in 6:K) for (j in 1:ncol(Z)) out = cbind(out,Z[,j]^k);}
  ## RETURN POLYNOMIAL BASIS
  return(out)
}

######################################################################
# Q3 (b): monte carlo simulation
######################################################################
K   <- c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 8, 9, 10)
K.r <- c(6, 11, 21, 26, 56, 61, 126, 131, 252, 257, 262, 267, 272, 277)
nK  <- length(K)
theta.hat <- matrix(NaN, ncol=nK, nrow=M)
se.hat    <- theta.hat


for (m in 1:M) {
  data <- DGP(N)
  X    <- data$x
  Y    <- data$y 
  TT   <- data$tt
  
  for (k in 1:nK) {
    X.pol <- cbind(1, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    
    # Compute annihalator matrix
    MP     <- diag(rep(1,N)) - X.Q %*% t(X.Q)
    
    # Pre-multiplly by MP
    Y.M <- MP %*% Y
    TT.M <- MP %*% TT
    
    # Get theta.hat using partition regression
    theta.hat[m, k] <- (t(TT.M) %*% Y.M) / (t(TT.M) %*% TT.M)
    
    # Get standard errors
    Sigma <- diag((as.numeric((Y.M - TT.M*theta.hat[m, k])))^2)
    se.hat[m, k] <- sqrt(t(TT.M) %*% Sigma %*% TT.M) / (t(TT.M) %*% TT.M)
  }
}

# Tabulate results
table <- matrix(NaN, ncol=6, nrow=nK)
for (k in 1:nK) {
  table[k, 1] <- K.r[k]
  table[k, 2] <- mean(theta.hat[, k])                               # point estimate
  table[k, 3] <- mean(theta.hat[, k]) - 1                           # bias
  table[k, 4] <- sd(theta.hat[, k])                                 # standard deviation
  table[k, 5] <- mean(se.hat[, k])                                  # mean standard error
  table[k, 6] <- mean((theta.hat[, k] - 1.96 * se.hat[, k] > 1) | 
                        (theta.hat[, k] + 1.96 * se.hat[, k] < 1))  # rejection rate
}
write.table(round(table,3), "PhD_Coursework/econ675/hw2/partial_linear.txt", sep=" & ", eol="\\\\ \n", col.names = FALSE, row.names = FALSE)

######################################################################
# Q3 (c): cross-validation
######################################################################

# cross validation function
K.CV <- function(tt, X, Y) {
  temp <- rep(NaN, nK)
  for (k in 1:nK) {
    X.pol <- cbind(1, tt, gen.P(X, K[k]))
    X.Q   <- qr.Q(qr(X.pol))
    XX <- X.Q %*% t(X.Q)
    Y.hat <- XX %*% Y
    W <- diag(XX)
    temp[k] <- mean(((Y-Y.hat) / (1-W))^2)
  }
  return(which.min(temp))
}

theta.hat2 <- rep(NaN, M)
se.hat2    <- theta.hat2
K.hat2     <- theta.hat2



for (m in 1:M) {
  data <- DGP(N)
  X    <- data$x
  Y    <- data$y 
  tt   <- data$tt
  
  k.opt <- K.CV(tt, X, Y)
  K.hat2[m]     <- K.r[k.opt]
}

table(K.hat2)