### Obtain control-function code from github 
source("https://raw.githubusercontent.com/zijguo/Control-function/main/cf.R")

# Define winsorize and rmse function for simulation result
winsorize <- function(x,frac=.05){
  if(length(frac) != 1 || frac < 0 ||
     frac > 0.5) {
    stop("bad value for 'fraction'")
  }
  lim <- quantile(x, probs=c(frac, 1-frac))
  x[ x < lim[1] ] <- lim[1]
  x[ x > lim[2] ] <- lim[2]
  x
}
rmse <- function(coef,true.coeff){
  coef.diff <- t(apply(coef,1,function(x) x-true.coeff))
  temp <- sqrt(apply(coef.diff^2,2,mean))
  return(temp)
}

### R Packages to Load ###
### The AER package is only needed to run the working example.
library(MASS)
library(mvtnorm)
library(AER)
library(foreign)
                       
########################
### simulation study ###
########################
                       
n <- 10000
iter <- 10000
set.seed(2021)

### Setting 1 : Satisfied model

mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0


for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  u1 <- err[,1]; v2 <- err[,2]
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2),verbose=FALSE) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}


### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR

### Setting 2.1: Cubic model

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0


for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  c <- sd(1+1/8*z1+1/3*z2+1/8*z2^2+z2^3)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  u1 <- err[,1]; v2 <- err[,2]
  y2 <- 1/(2*c)*(1+1/8*z1+1/3*z2+1/8*z2^2+z2^3)+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2),verbose=FALSE) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}


### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR


### Setting 2.2: Exp model

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0


for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  c <- sd(1+1/8*z1+1/3*z2+1/8*z2^2+exp(z2))
  err <- mvrnorm(n,mu=mu,Sigma = V)
  u1 <- err[,1]; v2 <- err[,2]
  y2 <- 1/(2*c)*(1+1/8*z1+1/3*z2+1/8*z2^2+exp(z2))+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2),verbose=FALSE) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}

### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR


### setting 3 (drastically violated)

beta2 <- 1; beta3 <- 0.2; beta4 <- 1
gamma1 <- 1; gamma2 <- 0.2
delta <- 0.5
true.coef <- c(1,1,beta2,beta3)
ERR <- 0

iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0

for (i in 1:iter) {
  # generating data
  z2 <- rnorm(n)
  u1 <- rnorm(n); v2 <- rnorm(n)
  temp <- rnorm(n)
  
  y2 <- -gamma2+gamma1*z2+gamma2*z2^2+v2
  w <- delta*v2^2+temp
  y1 <- beta2*y2+beta3*y2^2+beta4*w+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2),verbose=FALSE) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}

### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 1st and 2nd components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR



###################################################
### Sensitivity of joint distribution of errors ###
###################################################

### 1. double exponential

library(nimble) # to get function which can generate double exponential

true.coef <- c(1,1,10,10)
iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  # variance of double exponential distribution is 2*scale^2
  eps1 <- rdexp(n,scale = 1/sqrt(2)); eps2 <- rdexp(n,scale = 1/sqrt(2)) 
  u1 <- eps1; v2 <- eps1/2+sqrt(3)/2*eps2
  y2 <- 1+1/8*z1+1/3*z2+1/8*z2^2+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}

### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3

iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR


### 2. bivariate log-normal


iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0



for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  eps1 <- err[,1]; eps2 <- err[,2]
  u1 <- exp(eps1); v2 <- exp(eps2)
  y2 <- 1+1/8*z1+1/3*z2+1/8*z2^2+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}

### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR


### 3. bivariate absolutely normal

iv.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
cf.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
pret.coef <- matrix(NA,nrow = iter,ncol = length(true.coef))
ERR <- 0

for (i in 1:iter) {
  # generating data
  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  eps1 <- err[,1]; eps2 <- err[,2]
  u1 <- abs(eps1)-sqrt(2/pi); v2 <- abs(eps2)-sqrt(2/pi)
  y2 <- 1+1/8*z1+1/3*z2+1/8*z2^2+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1
  
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  iv.coef[i,] <- iv.fit$coefficients
  
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.coef[i,] <- cf.fit$coefficients
  
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.coef[i,] <- pret.fit$coefficients
  if (pret.fit$p_value<0.05) {
    ERR <- ERR+1
  }
  if (i==iter) {
    ERR <- ERR/iter
  }
}

### calculation
iv.wins.coef <- apply(iv.coef, 2, winsorize)
iv.wins.bias <- (apply(iv.wins.coef, 2, mean)-true.coef)/true.coef
iv.wins.rmse <- rmse(iv.wins.coef,true.coef)
iv.bias <- (apply(iv.coef, 2, mean)-true.coef)/true.coef
iv.rmse <- rmse(iv.coef,true.coef)

cf.wins.coef <- apply(cf.coef, 2, winsorize)
cf.wins.bias <- (apply(cf.wins.coef, 2, mean)-true.coef)/true.coef
cf.wins.rmse <- rmse(cf.wins.coef,true.coef)
cf.bias <- (apply(cf.coef, 2, mean)-true.coef)/true.coef
cf.rmse <- rmse(cf.coef,true.coef)


pret.wins.coef <- apply(pret.coef, 2, winsorize)
pret.wins.bias <- (apply(pret.wins.coef, 2, mean)-true.coef)/true.coef
pret.wins.rmse <- rmse(pret.wins.coef,true.coef)
pret.bias <- (apply(pret.coef, 2, mean)-true.coef)/true.coef
pret.rmse <- rmse(pret.coef,true.coef)


### 3rd and 4th components of terms are related to the value of beta2 and beta3
iv.wins.bias[3:4]
cf.wins.bias[3:4]
pret.wins.bias[3:4]
(cf.wins.rmse/iv.wins.rmse)[3:4]
(pret.wins.rmse/iv.wins.rmse)[3:4]
ERR



##########################
### real data analysis ###
##########################                       

setwd("~/2010_0015_data") # set directory

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
