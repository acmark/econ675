## ECON675: PS 2
## Erin Markiewitz 
## 10/12/2018
######################################################################
# Load packages, clear workspace
######################################################################
rm(list = ls())             #clear workspace
library(foreach)
library(dplyr)              #for data manipulation
library(data.table)         #for data manipulation
library(ggplot2)            #for pretty plots
library(boot)               #for bootstrapping
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
#Q1 (a): compute theoretically optimal bandwidth
######################################################################
# NB. This code only makes sense with the associated tex file...

# Write function to compute second derivative of normal density
d2norm  <- function(x, mu=0, v=1) {
  dnorm(x,mu,sqrt(v))*(((x-mu)/v)^2-1/v)
}

# Second derivative, squared of given Gaussian mixture
myf     <- function(x){
  (0.5*d2norm(x,-1.5,1.5)+0.5*d2norm(x,1,1))^2
}

# Compute required integral
k1     <-integrate(myf, -Inf, Inf)$val

# Compute optimal bandwidth
n      <- 1000
k2     <- 1/5
k3     <- 3/5
P      <- 2

h_aimse <- ((1/(2*P*n))*(factorial(P)/k2)^2*(k3/k1))^(1/(1+2*P))


######################################################################
#Q1 (b): monte carlo
######################################################################

# Function for EP kernel
K.ep    <- function(x){
  y <- .75 * (1-x^2) * (abs(x) <= 1)
  return(y)
}

# Function to compute true density value
f.true  <- function(x){
  y<-0.5*dnorm(x,-1.5,sqrt(1.5))+0.5*dnorm(x,1,1)
  return(y)
}

# Create vector of bandwidths
h.list = h_aimse*seq(0.5,1.5,0.1)

# Generate big matrix of random draws from the given Gaussian DGP
N          <- 1000
M          <- 1000

components <- sample(1:2,prob=c(0.5,0.5),size=n,replace=TRUE)
mu.vec     <- c(-1.5,1)
sd.vec     <- sqrt(c(1.5,1))

set.seed(5290)
X.mat      <- replicate(M,rnorm(n=N,mean=mu.vec[components],sd=sd.vec[components]))

# Function for computing LOO imse for a given bandwidth and random sample
imse.lo         <- function(x.rand=randx, h=h_aimse){
  
  # Compute leave-one-out fhats for each x_i
  y   = sapply(1:N,function(i) 1/(1000*h)*sum(K.ep((as.matrix(x.rand)[-i,]-x.rand[i])/h)))
  
  # Convert y to data.table for easy manipulation
  y   = as.data.table(y)
  
  # Add true density values
  y[, y.true := f.true(x.rand)]
  
  # Compute squared errors
  y[, sq_er.lo := (y - y.true)^2]
  
  # Compute imse.lo
  imse.lo <- y[, mean(sq_er.lo)]
  
  output <- imse.lo
  
  return(output)
}  

# Function for computing full-sample imse for a given bandwidth and random sample
imse.li         <- function(x.rand, h=h_aimse){
  
  # First compute vector of density estimates at each x_i
  y   = sapply(x.rand,function(x) 1/(1000*h)*sum(K.ep((x.rand-x)/h)))
  
  # Convert y to data.table for easy manipulation
  y   = as.data.table(y)
  
  # Add true density values
  y[, y.true := f.true(x.rand)]
  
  # Compute squared errors
  y[, sq_er.li := (y - y.true)^2]
  
  # Compute imse.li
  imse.li <- y[, mean(sq_er.li)]
  
  output <- imse.li
  
  return(output)
}  

#RUN SIMULATIONS - TOTAL RUNTIME APPROX 13-15 MINS
IMSE_LI <- foreach(h=h.list, .combine='cbind') %:%
   foreach(i=1:1000, .combine='c') %do% {
     imse.li(X.mat[,i],h)
   }

 IMSE_LO <- foreach(h=h.list, .combine='cbind') %:%
   foreach(i=1:1000, .combine='c') %do% {
     imse.lo(X.mat[,i],h)
   }

#Plot IMSEs
 IMSE_comb <- as.data.frame(cbind(h.list,colMeans(IMSE_LI),colMeans(IMSE_LO)))
 colnames(IMSE_comb) <- c("h", "IMSE_li","IMSE_lo")
 g <- melt(IMSE_comb, id="h")

 ggplot(g) +
   geom_line(aes(x=h, y=value, colour=variable)) +
   scale_colour_manual(values=c("blue","green")) +
   labs(title="Simulated IMSEs for Full and LOO Samples",y="IMSE") +theme(plot.title = element_text(hjust = 0.5))


######################################################################
#Q1 (d): rule-of-thumb bandwidth
######################################################################

# Write function to compute squared second derivative of normal density
d2normsq  <- function(x, mu=0, v=1) {
  (dnorm(x,mu,sqrt(v))*(((x-mu)/v)^2-1/v))^2
}


# Write function to compute ROT bandwidth for random sample
h.rot <- function(x.rand){
  
  # Compute sample mean and variance
  mu = mean(x.rand)
  v  = var(x.rand)
  
  # Compute second derivative of normal density
  k1     <-integrate(d2normsq,mu=mu,v=v, -Inf, Inf)$val
  
  # Compute ROT bandwidth
  h <- ((1/N)*(1/k2)^2*(k3/k1))^(1/5)
  
}

# Run simulation using foreach
h.rot.vec <- foreach(i=1:1000, .combine='c') %do% h.rot(X.mat[,i])

# Run simulation using sapply - FASTER!
h.rot.vec2<- sapply(1:M,function(i) h.rot(X.mat[,i]))

# Compute mean h.rot.vec
mean(h.rot.vec2)