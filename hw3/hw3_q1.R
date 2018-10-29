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


setwd("/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw3")

# load the data
pisofirme <- read.csv("pisofirme.csv", header = TRUE)
complete  <- complete.cases(pisofirme[, 5:27])
pisofirme <- pisofirme[complete, ]
# s_i: non-missing indicator
pisofirme$log_inc <- log(pisofirme$S_incomepc+1)
pisofirme$nmissing <- 1 - pisofirme$dmissing
pisoframe = as.data.frame(pisofirme)


# Get Piso Firme data
pisoframe <- as.data.table(read.csv('pisofirme.csv'))

# Create dependent variable for logistic regression
pisoframe[,nmissing:= 1-dmissing]

# Create income regressor
pisoframe[,log_inc:= log(S_incomepc+1)]

# Create income regressor
pisoframe[,log_inc:= log(S_incomepc+1)]

###################################################################
# Q 1_9 a
###################################################################

#estimate logit model
logit_q1 <- glm(nmissing ~ S_age + S_HHpeople + log_inc,
    family = "binomial", data = pisoframe)

#extract point estimates and calculate standard errors
b.hat   <- as.data.table(logit_q1["coefficients"])
V.hat   <- vcovHC(logit_q1, type = "HC1")
se.hat  <- as.data.table(sqrt(diag(V.hat)))
V.out <- diag(V.hat)
#compute t-stats and p values
t.stat <- b.hat/se.hat
n = nrow(pisofirme)
d = 4
p = round(2*pt(abs(t.stat[[1]]),df=n-d,lower.tail=FALSE),3)

#compute CI
CI.lower = b.hat - qnorm(0.975)*se.hat
CI.upper = b.hat + qnorm(0.975)*se.hat

results.a  = as.data.frame(cbind(b.hat,V.out,t.stat,p,CI.lower,CI.upper))
colnames(results.a) = c("Coef.","V","t-stat","p-val","CI.lower","CI.upper")
rownames(results.a) = c("Const.", "S_age","S_HHpeople","log_inc")

# Get latex table output
xtable(results.a,digits=5)
print(xtable(results.a, type = "latex"), file = "hw3_q1_9a_r.tex")

###################################################################
# Q 1_9 b
###################################################################


# set up logistic bootstrap
boot.logit <- function(data, i){
  logit  <- glm(nmissing ~ S_age + S_HHpeople + +I(log(S_incomepc+1)),
                data = data[i, ], family = "binomial")
  V      <- vcovHC(logit, type = "HC1")
  se     <- sqrt(diag(V.hat))
  t.boot <- (coef(logit)-coef(logit_q1))/se

  return(t.boot)
}

# run logistic bootstrap
set.seed(123)
boot.out <- boot(data=pisofirme, R=499, statistic = boot.logit, stype = "i")

# back out quantiles of boot t-dist. for CIs
boot.quant <- sapply(1:4, function (i) quantile(boot.out$t[,i], c(0.025, 0.975)))

#CIs
boot.ci.lower = b.hat + t(boot.quant)[,1]*se.hat
boot.ci.upper = b.hat + t(boot.quant)[,2]*se.hat

boot.p = sapply(1:4,function(i) 1/499*sum(boot.out$t[,i]>=t.stat[i]))

# Tabulate bootstrap results
results.b  = as.data.frame(cbind(b.hat,boot.ci.lower,boot.ci.upper,boot.p))
colnames(results.b) = c("Coef.","CI.lower","CI.upper","p-val")
rownames(results.b) = c("Const.", "S_age","S_HHpeople","log_inc")

# Get latex table output
xtable(results.b,digits=4)
print(xtable(results.b, type = "latex"), file = "hw3_q1_9b_r.tex")

###################################################################
# Q 1_9 C
###################################################################

# subset data
X  = pisoframe[, c("S_age","S_HHpeople","log_inc")]
X$const = 1
setcolorder(X,c("const","S_age","S_HHpeople","log_inc"))
b.hat = coef(logit_q1)

# Construct link function
mu = function(u){(1+exp(-u))^(-1)}

# Construct vector of x_i'*beta.hats
XB = as.matrix(X)%*%b.hat
# Compute predicted probabilities
mu.hat = mu(XB)
X[,mu.hat:=mu.hat]

#Make plot
plot(density(mu.hat,kernel="e", adjust = 5, bw="ucv",na.rm=TRUE),main="Kernel Density Estimate")
dev.copy(png,'hw3_q1_9c_r.png')
dev.off()


