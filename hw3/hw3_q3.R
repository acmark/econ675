###################################################################
# ECON 675, Assignment 3
# Fall 2018
# University of Michigan
# Latest update: Oct 22, 2018
###################################################################


rm(list=ls(all=TRUE))
library(foreign); library(MASS);
library(boot)
library(data.table)
library(foreach)
library(data.table)
library(Matrix)
library(ggplot2)
library(sandwich)
library(xtable)
library(gmm)

######################################################################
# Q3 1 
######################################################################
set.seed(123)
# Set up Enviorment 
N = 1000
X = runif(N,0,1)
x.max = max(X)

# Write function for bootrap statistic
boot.stat = function(data, i){
  N*(x.max -max(data[i]))
}

# Run bootsrap with 599 replications
boot.results = boot(data = X, R = 599, statistic = boot.stat)

# Make frequency plot
h         = hist(boot.results$t,plot=FALSE)
h$density = h$counts/sum(h$counts)
plot(h,freq=FALSE,main="Distribution of Bootstrap Statistic",xlab="Bootstrap statistic")
dev.copy(png,'hw3_q3_2_r.png')
dev.off()

######################################################################
# Q3: 2
######################################################################

# Generate parametric bootstrap samples
X.boot = replicate(599,runif(N,0,x.max))

# Compute maximums for each replications
x.max.boot = sapply(1:599,function(i) max(X.boot[,i]))

# Compute bootstrap statistic
t.boot     = N*(x.max -x.max.boot)

x.quant = range(c(0, 1, 100))
x.exp = dexp(x.quant, rate = 1, log = FALSE)

# Make frequency plot
h2         = hist(t.boot,plot=FALSE)
h2$density = h2$counts/sum(h2$counts)
plot(h2,freq=FALSE,main="Distribution of Parametric Bootstrap Statistic",xlab="Parametric bootstrap statistic",ylim=c(0,0.4),xlim=c(0,8))
dev.copy(png,'hw3_q3_3_r.png')
dev.off()


