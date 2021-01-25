rm(list = ls())

### use cf.R file as source
source("C:/Users/owner/Dropbox/Taehyeon/Violence-Data/cf.R")


library(mvtnorm)
library(AER)
library(foreign)

### simulation study
n <- 10000
set.seed(0125)
#setting 1

z1 <- rnorm(n); z2 <- rnorm(n)
mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
err <- mvrnorm(n,mu=mu,Sigma = V)
u1 <- err[,1]; v2 <- err[,2]
y2 <- 1+z1/8+z2/3+z2^2/8+v2
y1 <- 1+z1+10*y2+10*y2^2+u1
true.coef <- c(1,1,10,10)

iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
iv.fit$coefficients

cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
cf.fit$coefficients
pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
pret.fit

#setting 2 (drastically violated)
z2 <- rnorm(n)
beta2 <- 1; beta3 <- 0.2; beta4 <- 1
gamma1 <- 1; gamma2 <- 0.2
delta <- 0.5
err2 <- mvrnorm(n,mu=mu,Sigma = diag(2))
u1 <- rnorm(n); v2 <- rnorm(n)
temp <- rnorm(n)

y2 <- -gamma2+gamma1*z2+gamma2*z2^2+v2
w <- delta*v2^2+temp
y1 <- beta2*y2+beta3*y2^2+beta4*w+u1

iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
summary(iv.fit)
cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
cf.fit$coefficients
pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
pret.fit$coefficients



### real data analysis

rm(list = ls())
source("C:/Users/owner/Dropbox/Taehyeon/Violence-Data/cf.R")
library(mvtnorm)
library(AER)
library(foreign)



setwd("C:/Users/owner/Dropbox/Taehyeon/Violence-Data/2010_0015_data") # set directory

tempdata<-read.dta("All-Tables---Respondent-Data.dta") # our data
attach(tempdata)

### baseline covariates and preprocessing
X=cbind(lit,age,sex,totland_ha,dethnic,landgini_sc,dmarket,confland,ethnic_homo,soc_homo,ldensity08,ltotexp_ae)
mydata=data.frame(d,dead9303_bypop100,ldistance,laltitude,X)
cleandata=na.omit(mydata) # delete na
baseline=data.matrix(cleandata[,-c(1,2,3,4)]) # measured covariate
y=data.matrix(cleandata[,1]) # outcome
treatment=data.matrix(cleandata[,2]) # endogenous
z1=data.matrix(cleandata[,3]) # instrument variable
z2=data.matrix(cleandata[,4]) # instrument variable


outcome.formula <- y~treatment+I(treatment^2)+baseline # object formula
treatment.formula <- treatment~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline # object formula
cf.fit <- cf(outcome.formula,treatment.formula)
# abs(cf.fit$coefficients/sqrt(diag(cf.fit$vcov)))>pnorm(0.975)


### comparing to calculation with error.iv
fstg<-lm(treatment~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline)
fstg.quad<-lm(I(treatment^2)~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline+resid(fstg))
error.iv<-resid(fstg.quad)
ivmodel.cf=ivreg(y~treatment+I(treatment^2)+baseline|error.iv+z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)++baseline)
summary(ivmodel.cf)

### using pretest estimator
pretest.fit <- pretest(outcome.formula,treatment.formula)
pretest.fit
