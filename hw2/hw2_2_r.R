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
#Q2 (a): generate random data
######################################################################
N   = 1000
M   = 1000

# Generate x's
x.mat   = replicate(M,runif(N,-1,1))

# Generate chi-squared r.v
chi = replicate(M,rchisq(N,5))

# Generate errors
u.mat   = (chi-5)*x.mat

# Generate y's
y.mat   = exp(-0.1*(4*x.mat-1)^2)*sin(5*x.mat) + u.mat

######################################################################
#Q2 (b): cross-validation series estimator
######################################################################

# Write function to compute CV errors for K in 1:20
cross.val <- function(i){
  
  data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
  temp = rep(NaN, 20)
  
  for (k in 1:20) {
    
    data[, temp := x^k]
    setnames(data, "temp", paste0("x_", k))
    
    X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
    Y <- as.matrix(data[,y])
    
    # Compute projection matrix using QR decomp
    X.Q   <- qr.Q(qr(X))
    XX <- X.Q %*% t(X.Q)
    Y.hat <- XX %*% Y
    W <- diag(XX)
    
    temp[k] <- mean(((Y-Y.hat) / (1-W))^2)
    
  }
  return(temp)
}

# RUN SIMULATION -- RUNTIME 5 MINS
results     <- sapply(1:M,function(i) cross.val(i))

# Get average CV errors across simulations
results.avg <- rowMeans(results)

# Get the optimal K
K.hat       <- which.min(results.avg)


# #Plot CV
g <- as.data.frame(cbind(1:20,results.avg))
colnames(g) <- c("K", "CV")


ggplot(g,aes(x=K, y=CV)) +
  geom_line(linetype = "dashed")+
  geom_point()+
  labs(title="Simulated Cross-Validation Errors for M=1000 Simulations") +theme(plot.title = element_text(hjust = 0.5))

######################################################################
#Q2 (c): diagnostics
######################################################################

# Generate grid of x-values for plot
x.grid = seq(-1,1,0.1)

# Write function to compute optimal beta's (i.e. for K=7)
cv.beta <- function(i){
  
  data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
  
  for (k in 1:7) {
    data[, temp := x^k]
    setnames(data, "temp", paste0("x_", k))
  }  
  
  X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
  Y <- as.matrix(data[,y])
  
  beta <- solve(crossprod(X))%*%crossprod(X,Y) 
  
  return(t(beta)) 
  
}

# Write function to compute optimal standard errors (i.e for K=7)
cv.se <- function(i){
  
  data = data.table(y=y.mat[,i],x=x.mat[,i],const=1)
  
  for (k in 1:7) {
    data[, temp := x^k]
    setnames(data, "temp", paste0("x_", k))
  }  
  
  X <- as.matrix(data[, c("const",grep("x_", colnames(data), value = TRUE)), with = FALSE])
  Y <- as.matrix(data[,y])
  
  X.Q   <- qr.Q(qr(X))
  XX <- X.Q %*% t(X.Q)
  Y.hat <- XX %*% Y
  W <- diag(XX)
  
  s <- (1/(N-1))*sum((Y-Y.hat)^2)
  
  V <- s*(t(W)%*%W)
  
  se <- sqrt(V)
  
  return(se) 
  
}

# Get optimal betas for each m
results.beta <- sapply(1:M,function(i) cv.beta(i))

# Get associated standard errors
results.se <- sapply(1:M,function(i) cv.se(i))

# Compute averages over M
opt.beta  <- rowMeans(results.beta)
opt.se    <- rowMeans(results.se)

# Compute regressors for each number in the x.grid
X.new     <- t(sapply(x.grid, function(x) return(cbind(1,x,x^2,x^3,x^4,x^5,x^6,x^7))))

# Compute y.hats
y.hats    <- as.numeric(X.new%*%as.vector(opt.beta))

# Write the true regression function and compute for x.grid values
f.true <- function(x) exp(-0.1*(4*x-1)^2)*sin(5*x)
y.true <- f.true(x.grid)


# MAKE PLOT

# Get data in right format for ggplot
plot.data = melt(as.data.frame(cbind(x.grid,y.hats,y.true)),id="x.grid")

ggplot(plot.data,aes(x=x.grid,y=value,color=variable))+
  geom_line(linetype = "dashed")+geom_point()+
  labs(title="True and Series Estimate of the Regression function")+
  labs(y=expression(paste(mu(x))),x=expression(paste(x))) +theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "blue"),labels = c(expression(paste(hat(mu))),expression(paste(mu))))

######################################################################
#Q2 (d): derivative of the regression function
######################################################################


# Compute regressors for each number in the x.grid
X.der     <- t(sapply(x.grid, function(x) return(cbind(0,1,2*x,3*x^2,4*x^3,5*x^4,6*x^5,7*x^6))))

# Compute y.hats
dy.hats    <- as.numeric(X.der%*%as.vector(opt.beta))

# Write the true derivative function and compute for x.grid values
df.true <- function(x) exp(-0.1*(4*x-1)^2)*(5*cos(5*x)-0.8*sin(5*x)*(4*x-1))
dy.true <- df.true(x.grid)


# MAKE PLOT
dplot.data = melt(as.data.frame(cbind(x.grid,dy.hats,dy.true)),id="x.grid")

ggplot(dplot.data,aes(x=x.grid,y=value,color=variable))+
  geom_line(linetype = "dashed")+geom_point()+
  labs(title="True and Series Estimate of the Derivative Regression Function")+
  labs(y=expression(paste(mu(x))),x=expression(paste(x))) +theme(plot.title = element_text(hjust = 0.5))+
  scale_color_manual(values=c("black", "blue"),labels = c(expression(paste(d*hat(mu)/dx)),expression(paste(d*mu/dx))))

