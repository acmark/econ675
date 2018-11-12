###################################################################
# ECON 675, Assignment 4
# Erin Markiewitz
# Fall 2018
# University of Michigan
# Latest update: Nov 9, 2018
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
library(mvtnorm)
set.seed(12345)
setwd("/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw4")

####
# Generate data
####

obs = 50
reps = 1000 
covmatrix = matrix(c(1,0.85,.85,1),2,2)


#Generate X and Z, e, and y
W = replicate(reps,rmvnorm(obs, mean = c(0,0), sigma = covmatrix, method = "chol"))
e = replicate(reps,rnorm(50))
Y = sapply(1:reps,function(i) rep(1,obs) + W[,,i]%*%c(.5,1)+e[,i])

#estimate beta hat and gamma_tstats
beta_hat = sapply(1:reps, function(i) lm(Y[,i]~W[,,i])$coefficients[2])
gamma_tstat = sapply(1:reps,function(i) summary(lm(Y[,i]~W[,,i]))[["coefficients"]][, "t value"][3])
beta_hat_se = sapply(1:reps,function(i) summary(lm(Y[,i]~W[,1,i]))[["coefficients"]][, "Std. Error"][2])
lol = mean(beta_hat_se)
  
  
#estimate beta tilde
beta_tilde = sapply(1:reps, function(i) lm(Y[,i]~W[,2,i])$coefficients[2])

#estimate beta check
beta_check= ifelse(gamma_tstat>=1.96,beta_hat,beta_tilde)

# Summary statistics 
sum_beta = rbind(summary(beta_hat),summary(beta_tilde),summary(beta_check))
print(xtable(sum_beta, type = "latex"), file = "hw4_q3_1_r.tex")

plot_dat = data.frame(beta = c(beta_hat,beta_tilde,beta_check),Estimator=rep(c("hat", "tilde","check"), each = reps))

densplot = ggplot(plot_dat,aes(x=beta,fill=Estimator))+ 
  geom_density(alpha=0.5, kernel="e",bw="ucv")+
  ggtitle("Kernel Density Plots")+
  xlab("Point Estimator")+
  ylab("Density")+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_discrete( 
    name="Estimator",
    breaks=c("hat", "tilde", "check"),
    labels=c("hat", "tilde", "check" ))+
  theme(legend.justification = c(0.05, 0.98), legend.position = c(0.05, 0.98)) +  stat_function(fun = dnorm, n = 5000, args = list(mean = 0.5, sd = lol)) 
ggsave("hw4_q3_1_r.png")

#Coverage rates

#Compute coverage rate for beta_hat
beta_hat_cis     = cbind(beta_hat-1.96*beta_hat_se,beta_hat+1.96*beta_hat_se)
beta_hat_covdum  = ifelse(0.5>=beta_hat_cis[,1]&0.5<=beta_hat_cis[,2],1,0)
beta_hat_cr        = mean(beta_hat_covdum)

# Compute coverage rate for beta_tilde

beta_tilde_se       = sapply(1:reps,function(i) summary(lm(Y[,i]~W[,1,i]))[["coefficients"]][, "Std. Error"][2])
beta_tilde_cis     = cbind(beta_tilde-1.96*beta_tilde_se,beta_tilde+1.96*beta_tilde_se)
beta_tilde_covdum  = ifelse(0.5>=beta_tilde_cis[,1]&0.5<=beta_tilde_cis[,2],1,0)
beta_tilde_cr        = mean(beta_tilde_covdum)

# Compute coverage rate for beta_check
beta_check_cip    = ifelse(beta_hat==beta_check,beta_hat-1.96*beta_hat_se,beta_tilde-1.96*beta_tilde_se)
beta_check_cil    = ifelse(beta_hat==beta_check,beta_hat+1.96*beta_hat_se,beta_tilde+1.96*beta_tilde_se)
beta_check_cis    = cbind(beta_check_cil,beta_check_cip)
beta_check_covdum = ifelse(0.5>=beta_check_cis[,1]&0.5<=beta_check_cis[,2],1,0)
beta_check_cr     = mean(beta_check_covdum)

# Put results together
cr_table           = rbind(beta_hat_cr,beta_tilde_cr,beta_check_cr)
rownames(cr_table) = c("beta_hat","beta_tilde","beta_check")
colnames(cr_table) = c("Coverage Rate")

print(xtable(cr_table, type = "latex"), file = "hw4_q3_2_r.tex")


