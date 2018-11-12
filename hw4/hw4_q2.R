###################################################################
# ECON 675, Assignment 4
# Erin Markiewitz
# Fall 2018
# University of Michigan
# Latest update: Nov 9, 2018
###################################################################
rm(list = ls())             #clear workspace
library(foreach)            #for looping
library(data.table)         #for data manipulation
library(Matrix)             #fast matrix calcs
library(ggplot2)            #for pretty plots
library(sandwich)           #for variance-covariance estimation 
library(xtable)             #for latex tables
library(boot)               #for bootstrapping
library(CausalGAM)          
options(scipen = 999)       #forces R to use normal numbers instead of scientific notation

######################################################################
# Input data, subset/organize data
######################################################################
setwd("/Users/erinmarkiewitz/Dropbox/Phd_Coursework/Econ675/hw4")
data <- as.data.table(read.csv('LaLonde_all.csv'))

data = data[,log.re74:=log(re74+1)]
data = data[,log.re75:=log(re75+1)]
data = data[,age.sq:=age^2]
data = data[,educ.sq:=educ^2]
data = data[,age.cu:=age^3]
data = data[,black.u74:=black*u74]
data = data[,educ.logre74:=educ*log.re74]

#subset the data into subsets for LaLonde and PSID controls
X.l = data[treat ==1 | treat ==0]
Y.l = data[treat ==1 | treat ==0,.(re78)]

X.p = data[treat ==1 | treat ==2]
Y.p = data[treat ==1 | treat ==2,.(re78)]

#recode psid treatment indiator
X.p = X.p[,treat:=as.numeric(treat==1)]

#specify covariate sets
X.l.0  = X.l[,.(treat)]
X.p.0  = X.p[,.(treat)]

X.l.a  = X.l[,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75")]
X.p.a  = X.p[,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75")]

X.l.b  =  X.l[,-c("age.cu","black.u74","educ.logre74","re78","re74","re75")]
X.p.b  =  X.p[,-c("age.cu","black.u74","educ.logre74","re78","re74","re75")]

X.l.c  =  X.l[,-c("re78","re74","re75")]
X.p.c  =  X.p[,-c("re78","re74","re75")]

######################################################################
# 1) difference in means 
######################################################################
# lalonde
dm.l    = lm(as.matrix(Y.l)~as.matrix(X.l.0))
dm.l.se = sqrt(diag(vcovHC(dm.l, type = "HC1")))
dm.l.ciu = dm.l$coefficients + 1.96*dm.l.se
dm.l.cil = dm.l$coefficients - 1.96*dm.l.se
dm.l.out = cbind(dm.l$coefficients,dm.l.se,dm.l.cil,dm.l.ciu)

#psid
dm.p = lm(as.matrix(Y.p)~as.matrix(X.p.0))
dm.p.se = sqrt(diag(vcovHC(dm.p, type = "HC1")))
dm.p.ciu = dm.p$coefficients + 1.96*dm.p.se
dm.p.cil = dm.p$coefficients - 1.96*dm.p.se
dm.p.out = cbind(dm.p$coefficients,dm.p.se,dm.p.cil,dm.p.ciu)

######################################################################
# 2) OLS
######################################################################
# lalonde
ols.l.a    = lm(as.matrix(Y.l)~as.matrix(X.l.a))
ols.l.se.a = sqrt(diag(vcovHC(ols.l.a, type = "HC1")))
ols.l.ciu.a = ols.l.a$coefficients + 1.96*ols.l.se.a
ols.l.cil.a = ols.l.a$coefficients - 1.96*ols.l.se.a

ols.l.b    = lm(as.matrix(Y.l)~as.matrix(X.l.b))
ols.l.se.b = sqrt(diag(vcovHC(ols.l.b, type = "HC1")))
ols.l.ciu.b = ols.l.b$coefficients + 1.96*ols.l.se.b
ols.l.cil.b = ols.l.b$coefficients - 1.96*ols.l.se.b

ols.l.c    = lm(as.matrix(Y.l)~as.matrix(X.l.c))
ols.l.se.c = sqrt(diag(vcovHC(ols.l.c, type = "HC1")))
ols.l.ciu.c = ols.l.c$coefficients + 1.96*ols.l.se.c
ols.l.cil.c = ols.l.c$coefficients - 1.96*ols.l.se.c

# PSID
ols.p.a    = lm(as.matrix(Y.p)~as.matrix(X.p.a))
ols.p.se.a = sqrt(diag(vcovHC(ols.p.a, type = "HC1")))
ols.p.ciu.a = ols.p.a$coefficients + 1.96*ols.p.se.a
ols.p.cil.a = ols.p.a$coefficients - 1.96*ols.p.se.a

ols.p.b    = lm(as.matrix(Y.p)~as.matrix(X.p.b))
ols.p.se.b = sqrt(diag(vcovHC(ols.p.b, type = "HC1")))
ols.p.ciu.b = ols.p.b$coefficients + 1.96*ols.p.se.b
ols.p.cil.b = ols.p.b$coefficients - 1.96*ols.p.se.b

ols.p.c    = lm(as.matrix(Y.p)~as.matrix(X.p.c))
ols.p.se.c = sqrt(diag(vcovHC(ols.p.c, type = "HC1")))
ols.p.ciu.c = ols.p.c$coefficients + 1.96*ols.p.se.c
ols.p.cil.c = ols.p.c$coefficients - 1.96*ols.p.se.c

#out
ols.l.out = cbind(c(ols.l.a$coefficients[2],ols.l.b$coefficients[2],ols.l.c$coefficients[2]),c(ols.l.se.a[2],ols.l.se.b[2],ols.l.se.c[2]),c(ols.l.cil.a[2],ols.l.cil.b[2],ols.l.cil.c[2]),c(ols.l.ciu.a[2],ols.l.ciu.b[2],ols.l.ciu.c[2]))
ols.p.out = cbind(c(ols.p.a$coefficients[2],ols.p.b$coefficients[2],ols.p.c$coefficients[2]),c(ols.p.se.a[2],ols.p.se.b[2],ols.p.se.c[2]),c(ols.p.cil.a[2],ols.l.cil.b[2],ols.l.cil.c[2]),c(ols.l.ciu.a[2],ols.l.ciu.b[2],ols.l.ciu.c[2]))

######################################################################
# 3) RI
######################################################################

################
# covariates a #
################
# subset the outcome and coveriate series 
Y.treat       = data[treat==1,.(re78)]
Y.control.l   = data[treat==0,.(re78)]
Y.control.p   = data[treat==2,.(re78)]

#subset covariates for imputation
X.treat.a      = data[treat==1,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.l.a  = data[treat==0,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.p.a  = data[treat==2,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]

#estimate ols coefficients for imputation
ols.treat.a         = lm(as.matrix(Y.treat)~as.matrix(X.treat.a))
ols.control.l.a     = lm(as.matrix(Y.control.l)~as.matrix(X.control.l.a))
ols.control.p.a     = lm(as.matrix(Y.control.p)~as.matrix(X.control.p.a))

#insert constants
X.treat.a[,const:=1]
setcolorder(X.treat.a,c("const"))
X.control.l.a[,const:=1]
setcolorder(X.control.l.a,c("const"))
X.control.p.a[,const:=1]
setcolorder(X.control.p.a,c("const"))

#impute individual treatment effects 
tvec.ri.treat.l.a      = as.matrix(X.treat.a)%*%(as.vector(ols.treat.a$coefficients)-as.vector(ols.control.l.a$coefficients))
tvec.ri.treat.p.a      = as.matrix(X.treat.a)%*%(as.vector(ols.treat.a$coefficients)-as.vector(ols.control.p.a$coefficients))

tvec.ri.control.l.a    = as.matrix(X.control.l.a)%*%(as.vector(ols.treat.a$coefficients)-as.vector(ols.control.l.a$coefficients))  
tvec.ri.control.p.a    = as.matrix(X.control.p.a)%*%(as.vector(ols.treat.a$coefficients)-as.vector(ols.control.p.a$coefficients))  


#ATE
ate.ri.l.a       = mean(c(tvec.ri.treat.l.a,tvec.ri.control.l.a))
ate.ri.p.a       = mean(c(tvec.ri.treat.p.a,tvec.ri.control.p.a))


#ATT
att.ri.a          = mean(tvec.ri.treat.l.a)


################
# covariates b #
################
Y.treat       = data[treat==1,.(re78)]
Y.control.l   = data[treat==0,.(re78)]
Y.control.p   = data[treat==2,.(re78)]

#subset covariates for imputation
X.treat.b      = data[treat==1,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.l.b  = data[treat==0,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.p.b  = data[treat==2,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]

#estimate ols coefficients for imputation
ols.treat.b         = lm(as.matrix(Y.treat)~as.matrix(X.treat.b))
ols.control.l.b     = lm(as.matrix(Y.control.l)~as.matrix(X.control.l.b))
ols.control.p.b     = lm(as.matrix(Y.control.p)~as.matrix(X.control.p.b))

#insert constants
X.treat.b[,const:=1]
setcolorder(X.treat.b,c("const"))
X.control.l.b[,const:=1]
setcolorder(X.control.l.b,c("const"))
X.control.p.b[,const:=1]
setcolorder(X.control.p.b,c("const"))

#impute individual treatment effects 
tvec.ri.treat.l.b      = as.matrix(X.treat.b)%*%(as.vector(ols.treat.b$coefficients)-as.vector(ols.control.l.b$coefficients))
tvec.ri.treat.p.b      = as.matrix(X.treat.b)%*%(as.vector(ols.treat.b$coefficients)-as.vector(ols.control.p.b$coefficients))

tvec.ri.control.l.b    = as.matrix(X.control.l.b)%*%(as.vector(ols.treat.b$coefficients)-as.vector(ols.control.l.b$coefficients))  
tvec.ri.control.p.b    = as.matrix(X.control.p.b)%*%(as.vector(ols.treat.b$coefficients)-as.vector(ols.control.p.b$coefficients))  


#ATE
ate.ri.l.b       = mean(c(tvec.ri.treat.l.b,tvec.ri.control.l.b))
ate.ri.p.b       = mean(c(tvec.ri.treat.p.b,tvec.ri.control.p.b))


#ATT
att.ri.b          = mean(tvec.ri.treat.l.b)

################
# covariates c #
################

Y.treat       = data[treat==1,.(re78)]
Y.control.l   = data[treat==0,.(re78)]
Y.control.p   = data[treat==2,.(re78)]

#subset covariates for imputation
X.treat.c      = data[treat==1,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.l.c  = data[treat==0,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]
X.control.p.c  = data[treat==2,-c("age.sq","educ.sq","age.cu","black.u74","educ.logre74","u74","u75","re78","re74","re75","treat")]

#estimate ols coefficients for imputation
ols.treat.c         = lm(as.matrix(Y.treat)~as.matrix(X.treat.c))
ols.control.l.c     = lm(as.matrix(Y.control.l)~as.matrix(X.control.l.c))
ols.control.p.c     = lm(as.matrix(Y.control.p)~as.matrix(X.control.p.c))

#insert constants
X.treat.c[,const:=1]
setcolorder(X.treat.c,c("const"))
X.control.l.c[,const:=1]
setcolorder(X.control.l.c,c("const"))
X.control.p.c[,const:=1]
setcolorder(X.control.p.c,c("const"))

#impute individual treatment effects 
tvec.ri.treat.l.c      = as.matrix(X.treat.c)%*%(as.vector(ols.treat.c$coefficients)-as.vector(ols.control.l.c$coefficients))
tvec.ri.treat.p.c      = as.matrix(X.treat.c)%*%(as.vector(ols.treat.c$coefficients)-as.vector(ols.control.p.c$coefficients))

tvec.ri.control.l.c    = as.matrix(X.control.l.c)%*%(as.vector(ols.treat.c$coefficients)-as.vector(ols.control.l.c$coefficients))  
tvec.ri.control.p.c    = as.matrix(X.control.p.c)%*%(as.vector(ols.treat.c$coefficients)-as.vector(ols.control.p.c$coefficients))  


#ATE
ate.ri.l.c       = mean(c(tvec.ri.treat.l.c,tvec.ri.control.l.c))
ate.ri.p.c       = mean(c(tvec.ri.treat.p.c,tvec.ri.control.p.c))

#ATT
att.ri.c          = mean(tvec.ri.treat.l.c)

######################################################################
# IPW and Doubly Robust  (ATT)
######################################################################

# Generate treatment outcome variables
T.l = data[treat==1|treat==0,.(treat)]
T.p = data[treat==1|treat==2,.(treat)]

#Recode 2's to 0's for PSID sample
T.p = T.p[,treat:=as.numeric(treat==1)]

# Get propensity scores using logit regression
prop.l.a = glm(as.matrix(T.l) ~ as.matrix(X.l.a[,-c("treat")]),family = "binomial")
prop.l.b = glm(as.matrix(T.l) ~ as.matrix(X.l.b[,-c("treat")]),family = "binomial")
prop.l.c = glm(as.matrix(T.l) ~ as.matrix(X.l.c[,-c("treat")]),family = "binomial")

prop.p.a = glm(as.matrix(T.p) ~ as.matrix(X.p.a[,-c("treat")]),family = "binomial")
prop.p.b = glm(as.matrix(T.p) ~ as.matrix(X.p.b[,-c("treat")]),family = "binomial")
prop.p.c = glm(as.matrix(T.p) ~ as.matrix(X.p.c[,-c("treat")]),family = "binomial")

# Add prop scores to the data matrices for easy computing of treatment effects
X.l.ipw = X.l
X.l.ipw[,ps.a:=prop.l.a$fitted.values]
X.l.ipw[,ps.b:=prop.l.b$fitted.values]
X.l.ipw[,ps.c:=prop.l.c$fitted.values]

X.p.ipw = X.p
X.p.ipw[,ps.a:=prop.p.a$fitted.values]
X.p.ipw[,ps.b:=prop.p.b$fitted.values]
X.p.ipw[,ps.c:=prop.p.c$fitted.values]

# Create variables for computing ATEs 
X.l.ipw[,t1.a:=treat*re78/ps.a]
X.l.ipw[,t0.a:=(1-treat)*re78/(1-ps.a)]
X.l.ipw[,t1.b:=treat*re78/ps.b]
X.l.ipw[,t0.b:=(1-treat)*re78/(1-ps.b)]
X.l.ipw[,t1.c:=treat*re78/ps.c]
X.l.ipw[,t0.c:=(1-treat)*re78/(1-ps.c)]

# Compute proportion of treated respondents
p.l           = mean(X.l[,treat])

# Create additional variables for computing ATTs
X.l.ipw[,t1.att:=treat*re78/p.l]
X.l.ipw[,t0.a2:=(1-treat)*re78/(1-ps.a)*(ps.a/p.l)]
X.l.ipw[,t0.b2:=(1-treat)*re78/(1-ps.b)*(ps.b/p.l)]
X.l.ipw[,t0.c2:=(1-treat)*re78/(1-ps.c)*(ps.c/p.l)]

# Compute ATTs
att.ipw.l.a  = mean(X.l.ipw[,t1.att])-mean(X.l.ipw[,t0.a2])
att.ipw.l.b  = mean(X.l.ipw[,t1.att])-mean(X.l.ipw[,t0.b2])
att.ipw.l.c  = mean(X.l.ipw[,t1.att])-mean(X.l.ipw[,t0.c2])




######################################################################
# [4.b] Inverse Probability Weighting, PSID control
######################################################################

# Create variables for computing ATEs 
X.p.ipw[,t1.a:=treat*re78/ps.a]
X.p.ipw[,t0.a:=(1-treat)*re78/(1-ps.a)]
X.p.ipw[,t1.b:=treat*re78/ps.b]
X.p.ipw[,t0.b:=(1-treat)*re78/(1-ps.b)]
X.p.ipw[,t1.c:=treat*re78/ps.c]
X.p.ipw[,t0.c:=(1-treat)*re78/(1-ps.c)]

# Compute proportion of treated respondents
p.p           = mean(X.p[,treat])

# Create additional variables for computing ATTs
X.p.ipw[,t1.att:=treat*re78/p.p]
X.p.ipw[,t0.a2:=(1-treat)*re78/(1-ps.a)*(ps.a/p.p)]
X.p.ipw[,t0.b2:=(1-treat)*re78/(1-ps.b)*(ps.b/p.p)]
X.p.ipw[,t0.c2:=(1-treat)*re78/(1-ps.c)*(ps.c/p.p)]

# Compute ATTs
att.ipw.p.a  = mean(X.p.ipw[,t1.att])-mean(X.p.ipw[,t0.a2])
att.ipw.p.b  = mean(X.p.ipw[,t1.att])-mean(X.p.ipw[,t0.b2])
att.ipw.p.c  = mean(X.p.ipw[,t1.att])-mean(X.p.ipw[,t0.c2])



######################################################################
# IPW and Doubly Robust  (ATE)
######################################################################

# Covariates A
ATE.l.a <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.l),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
) 

# Covariates B
ATE.l.b <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.l),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
) 

# Covariates C 
ATE.l.c <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.l),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
)

 #Covariates A, PSID control #can't calculate standard error
 ATE.p.a <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                          pscore.family = binomial,
                          outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                          outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75,
                          outcome.family = gaussian,
                          treatment.var = "treat",
                          data=as.data.frame(X.p),
                          divby0.action="t",
                          divby0.tol=0.001,
                          var.gam.plot=FALSE,
                          nboot=0,
                          variance.smooth.span = 0
 )

# Covariates B, PSID control
ATE.p.b <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.p),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
) 

# Covariates C, PSID control
ATE.p.c <- estimate.ATE(pscore.formula = treat ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         pscore.family = binomial,
                         outcome.formula.t = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.formula.c = re78 ~ age + educ + black + hisp + married + nodegr + log.re74 + log.re75 + age.sq + educ.sq + u74 + u75 + age.cu + black.u74 + educ.logre74,
                         outcome.family = gaussian,
                         treatment.var = "treat",
                         data=as.data.frame(X.p),
                         divby0.action="t",
                         divby0.tol=0.001,
                         var.gam.plot=FALSE,
                         nboot=0
)

######################################################################
#  construct table 
######################################################################

## lalonde controls
# Mean Diff + OLS results
a  =    rbind(dm.l.out[2,],ols.l.out)

# Reg imputation results
b1  =   c(ATE.l.a$ATE.reg.hat,ATE.l.a$ATE.reg.asymp.SE,ATE.l.a$ATE.reg.hat-1.96*ATE.l.a$ATE.reg.asymp.SE,ATE.l.a$ATE.reg.hat+1.96*ATE.l.a$ATE.reg.asymp.SE)
b2  =   c(ATE.l.b$ATE.reg.hat,ATE.l.b$ATE.reg.asymp.SE,ATE.l.b$ATE.reg.hat-1.96*ATE.l.b$ATE.reg.asymp.SE,ATE.l.b$ATE.reg.hat+1.96*ATE.l.b$ATE.reg.asymp.SE)
b3  =   c(ATE.l.c$ATE.reg.hat,ATE.l.c$ATE.reg.asymp.SE,ATE.l.c$ATE.reg.hat-1.96*ATE.l.c$ATE.reg.asymp.SE,ATE.l.c$ATE.reg.hat+1.96*ATE.l.c$ATE.reg.asymp.SE)

# IPW results
c1  =   c(ATE.l.a$ATE.IPW.hat,ATE.l.a$ATE.IPW.asymp.SE,ATE.l.a$ATE.IPW.hat-1.96*ATE.l.a$ATE.IPW.asymp.SE,ATE.l.a$ATE.IPW.hat+1.96*ATE.l.a$ATE.IPW.asymp.SE)
c2  =   c(ATE.l.b$ATE.IPW.hat,ATE.l.b$ATE.IPW.asymp.SE,ATE.l.b$ATE.IPW.hat-1.96*ATE.l.b$ATE.IPW.asymp.SE,ATE.l.b$ATE.IPW.hat+1.96*ATE.l.b$ATE.IPW.asymp.SE)
c3  =   c(ATE.l.c$ATE.IPW.hat,ATE.l.c$ATE.IPW.asymp.SE,ATE.l.c$ATE.IPW.hat-1.96*ATE.l.c$ATE.IPW.asymp.SE,ATE.l.c$ATE.IPW.hat+1.96*ATE.l.c$ATE.IPW.asymp.SE)

# Doubly robust results
d1  =   c(ATE.l.a$ATE.AIPW.hat,ATE.l.a$ATE.AIPW.asymp.SE,ATE.l.a$ATE.AIPW.hat-1.96*ATE.l.a$ATE.AIPW.asymp.SE,ATE.l.a$ATE.AIPW.hat+1.96*ATE.l.a$ATE.AIPW.asymp.SE)
d2  =   c(ATE.l.b$ATE.AIPW.hat,ATE.l.b$ATE.AIPW.asymp.SE,ATE.l.b$ATE.AIPW.hat-1.96*ATE.l.b$ATE.AIPW.asymp.SE,ATE.l.b$ATE.AIPW.hat+1.96*ATE.l.b$ATE.AIPW.asymp.SE)
d3  =   c(ATE.l.c$ATE.AIPW.hat,ATE.l.c$ATE.AIPW.asymp.SE,ATE.l.c$ATE.AIPW.hat-1.96*ATE.l.c$ATE.AIPW.asymp.SE,ATE.l.c$ATE.AIPW.hat+1.96*ATE.l.c$ATE.AIPW.asymp.SE)

## PSID control

# Mean Diff + OLS results
e  =    rbind(dm.p.out[2,],ols.p.out)

# Reg imputation results
f1  =   c(ATE.p.a$ATE.reg.hat,0,0,0)
f2  =   c(ATE.p.b$ATE.reg.hat,ATE.p.b$ATE.reg.asymp.SE,ATE.p.b$ATE.reg.hat-1.96*ATE.p.b$ATE.reg.asymp.SE,ATE.p.b$ATE.reg.hat+1.96*ATE.p.b$ATE.reg.asymp.SE)
f3  =   c(ATE.p.c$ATE.reg.hat,ATE.p.c$ATE.reg.asymp.SE,ATE.p.c$ATE.reg.hat-1.96*ATE.p.c$ATE.reg.asymp.SE,ATE.p.c$ATE.reg.hat+1.96*ATE.p.c$ATE.reg.asymp.SE)

# IPW results
g1  =   c(ATE.p.a$ATE.IPW.hat,0,0,0)
g2  =   c(ATE.p.b$ATE.IPW.hat,ATE.p.b$ATE.IPW.asymp.SE,ATE.p.b$ATE.IPW.hat-1.96*ATE.p.b$ATE.IPW.asymp.SE,ATE.p.b$ATE.IPW.hat+1.96*ATE.p.b$ATE.IPW.asymp.SE)
g3  =   c(ATE.p.c$ATE.IPW.hat,ATE.p.c$ATE.IPW.asymp.SE,ATE.p.c$ATE.IPW.hat-1.96*ATE.p.c$ATE.IPW.asymp.SE,ATE.p.c$ATE.IPW.hat+1.96*ATE.p.c$ATE.IPW.asymp.SE)

# Doubly robust results
h1  =   c(ATE.p.a$ATE.AIPW.hat,0,0,0)
h2  =   c(ATE.p.b$ATE.AIPW.hat,ATE.p.b$ATE.AIPW.asymp.SE,ATE.p.b$ATE.AIPW.hat-1.96*ATE.p.b$ATE.AIPW.asymp.SE,ATE.p.b$ATE.AIPW.hat+1.96*ATE.p.b$ATE.AIPW.asymp.SE)
h3  =   c(ATE.p.c$ATE.AIPW.hat,ATE.p.c$ATE.AIPW.asymp.SE,ATE.p.c$ATE.AIPW.hat-1.96*ATE.p.c$ATE.AIPW.asymp.SE,ATE.p.c$ATE.AIPW.hat+1.96*ATE.p.c$ATE.AIPW.asymp.SE)


## PUT RESULTS TOGETHER
l.out  = rbind(a,b1,b2,b3,c1,c2,c3,d1,d2,d3)
p.out  = rbind(e,f1,f2,f3,g1,g2,g3,h1,h2,h3)
ate.out = round(cbind(l.out,p.out),2)

# EXPORT RESULTS AS CSV
write.table(ate.out, file = "Table1_ATE_resultq.csv",row.names=FALSE,col.names=FALSE,sep=",")



