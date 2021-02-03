rm(list = ls())

### use cf.R file as source
source("C:/Users/owner/Dropbox/Taehyeon/Violence-Data/cf.R")

library(MASS)
library(mvtnorm)
library(AER)
library(foreign)

### simulation study
n <- 10000
iter <- 1000
set.seed(2021)
#setting 1

mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
mse <- matrix(0,nrow = iter, ncol = 2)

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  u1 <- err[,1]; v2 <- err[,2]
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1

  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  iv.mse <- mean((iv.coef[i,]-true.coef)^2)

  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  cf.mse <- mean((cf.coef[i,]-true.coef)^2)

  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients

  ifelse(iv.mse<cf.mse,mse[i,1] <-1,mse[i,2] <- 1)
}

apply(mse,2,sum) # comparing mse of 2sls vs cf


#setting 3 (drastically violated)

beta2 <- 1; beta3 <- 0.2; beta4 <- 1
gamma1 <- 1; gamma2 <- 0.2
delta <- 0.5
true.coef <- c(beta2,beta3)

iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
mse <- matrix(0,nrow = iter, ncol = 2)

for (i in 1:iter) {
  # generating data
  z2 <- rnorm(n)
  err2 <- mvrnorm(n,mu=mu,Sigma = diag(2))
  u1 <- rnorm(n); v2 <- rnorm(n)
  temp <- rnorm(n)

  y2 <- -gamma2+gamma1*z2+gamma2*z2^2+v2
  w <- delta*v2^2+temp
  y1 <- beta2*y2+beta3*y2^2+beta4*w+u1
  
  iv.fit <- ivreg(y1~0+y2+I(y2^2)+w|z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients[-3]
  iv.mse <- mean((iv.coef[i,]-true.coef)^2)

  cf.fit <- cf(y1~0+y2+I(y2^2)+w,y2~z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients[-3]
  cf.mse <- mean((cf.coef[i,]-true.coef)^2)

  pret.fit <- pretest(y1~0+y2+I(y2^2)+w,y2~z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients[-3]

  ifelse(iv.mse<cf.mse,mse[i,1] <-1,mse[i,2] <- 1)

}

apply(mse,2,sum) # comparing mse of 2sls vs cf

## non-normal distribution

### 1. double exponential

library(nimble) # to get function which can generate double exponential

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
mse <- matrix(0,nrow = iter, ncol = 2)

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  eps1 <- rdexp(n); eps2 <- rdexp(n)
  u1 <- exp(eps1); v2 <- eps1/2+sqrt(3)/2*eps2
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1


  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  iv.mse <- mean((iv.coef[i,]-true.coef)^2)

  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  cf.mse <- mean((cf.coef[i,]-true.coef)^2)

  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients

  ifelse(iv.mse<cf.mse,mse[i,1] <-1,mse[i,2] <- 1)
}

apply(mse,2,sum) # comparing mse of 2sls vs cf



### 2. bivariate log-normal

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
mse <- matrix(0,nrow = iter, ncol = 2)

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  eps1 <- err[,1]; eps2 <- err[,2]
  u1 <- exp(eps1); v2 <- exp(eps2)
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1


  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  iv.mse <- mean((iv.coef[i,]-true.coef)^2)

  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  cf.mse <- mean((cf.coef[i,]-true.coef)^2)

  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients

  ifelse(iv.mse<cf.mse,mse[i,1] <-1,mse[i,2] <- 1)
}

apply(mse,2,sum) # comparing mse of 2sls vs cf

### 3. bivariate absolutely normal

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
mse <- matrix(0,nrow = iter, ncol = 2)

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  eps1 <- err[,1]; eps2 <- err[,2]
  u1 <- abs(eps1)-sqrt(2/pi); v2 <- abs(eps2)-sqrt(2/pi)
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1


  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  iv.mse <- mean((iv.coef[i,]-true.coef)^2)

  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  cf.mse <- mean((cf.coef[i,]-true.coef)^2)

  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients

  ifelse(iv.mse<cf.mse,mse[i,1] <-1,mse[i,2] <- 1)
}

apply(mse,2,sum) # comparing mse of 2sls vs cf


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
z1=data.matrix(cleandata[,3]) # instrument variable 1
z2=data.matrix(cleandata[,4]) # instrument variable 2

iv.fit <- ivreg(y~treatment+I(treatment^2)+baseline|z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline)
summary(iv.fit)
outcome.formula <- y~treatment+I(treatment^2)+baseline # object formula
treatment.formula <- treatment~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline # object formula

cf.fit <- cf(outcome.formula,treatment.formula)
# abs(cf.fit$coefficients/sqrt(diag(cf.fit$vcov)))>pnorm(0.975)
cf.fit$coefficients
sqrt(cf.fit$vcov[2,2])
sqrt(cf.fit$vcov[3,3]) # results are same

### comparing to calculation with error.iv reg
fstg<-lm(treatment~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline)
fstg.quad<-lm(I(treatment^2)~z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)+baseline+resid(fstg))
error.iv<-resid(fstg.quad)
ivmodel.cf=ivreg(y~treatment+I(treatment^2)+baseline|error.iv+z1+z2+I(z1^2)+I(z2^2)+I(z1*z2)++baseline)
summary(ivmodel.cf)

### using pretest estimator
pretest.fit <- pretest(outcome.formula,treatment.formula)
pretest.fit
