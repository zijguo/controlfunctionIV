# Control-function
--------------------
The methods in these R functions are based on work by Guo and Small (2016).

## Installation
To load these functions into R, run the following commands. (Indeed, we need to change this repo to public or change a setting about token)
```{r}
source("https://raw.githubusercontent.com/zijguo/Control-function/main/cf.R")
```
## Examples
We wil introduce a simulation example to show how we can use. The code example.R has additional working examples.

```{r}
rm(list = ls())

### Use cf.R file as source
source("C:/Users/owner/Dropbox/Taehyeon/Violence-Data/cf.R")


library(mvtnorm)
library(AER) # to use ivreg function which have two stage least square
library(foreign)

## Simulation study

### Generate data
n <- 10000

set.seed(2021)

### setting 1

mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
true.coef <- c(1,1,10,10) # our true beta coefficient

### generating data

  z1 <- rnorm(n); z2 <- rnorm(n)
  err <- mvrnorm(n,mu=mu,Sigma = V)
  u1 <- err[,1]; v2 <- err[,2]
  y2 <- 1+z1/8+z2/3+z2^2/8+v2
  y1 <- 1+z1+10*y2+10*y2^2+u1

### fit two stage least square
  iv.fit <- ivreg(y1~z1+y2+I(y2^2)|z1+z2+I(z2^2)) # 2SLS
  summary(iv.fit)

### apply control function method
  cf.fit <- cf(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Control function
  cf.fit$coefficients
  cf.fit$vcov

### apply pretest estimation
  pret.fit <- pretest(y1~z1+y2+I(y2^2),y2~z1+z2+I(z2^2)) # Pretest estimator
  pret.fit$coefficients
  pret.fit$vcov


```
## References
Guo, Z. and D. S. Small (2016), [Control function instrumental variable estimation of nonlinear
causal effect models](https://www.jmlr.org/papers/volume17/14-379/14-379.pdf), The Journal of Machine Learning Research 17(1), 3448–3482.
