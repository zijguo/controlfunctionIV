### cf.R
### Function: Implements the control function method
###           for estimation and inference of nonlinear
###           treatment effects 
###           as developed in Guo and Small (2016)             
### Maintainer: Taehyeon Koo
### E-mail: tk587@scarletmail.rutgers.edu

### cf
### FUNCTION: Point estimate, covariance matrix for treatment effect with  
###           using a single-sample, individual-level data via CF. 
### INPUT: Outcome formula, such as Y ~ X + D + g_2(D) + ... + g_k(D)
###        Treatment formula, such as D ~ X + Z + h_2(Z) + ... + h_k(Z) 
###        where
###        Y, continuous, non-missing, numeric outcome vector (u by 1 vector)
###        D, continuous or discrete, non-missing, numeric treatment vector (n by 1 vector)        
###        Z, continuous or discrete, non-missing, numeric instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete, non-missing, numeric matrix containing p_x 
###           covariates (n by p_x matrix)
###        g_2,...,g_k, functions of treatment variables, eg. D^2
###        h_2,...,h_k, functions of instrument variables, eg. exp(Z) 
### OUTPUT: a list (a) coefficients (scalar numeric value:
###                                  the estimate of the treatment effect)
###                (b) vcov (numeric matrix:
###                          estimated covariance matrix of coefficients)
library(Formula)

cf <- function(outcome.formula,treatment.formula){
  outcome.formula <- Formula(outcome.formula)
  treatment.formula <- Formula(treatment.formula)

  cf.fstg <- lm(treatment.formula)
  e1 <- resid(cf.fstg)
  cf.formula <- as.formula(paste(format(outcome.formula),"+","e1"))
  cf.fit <- lm(cf.formula)

  cf.coef <- coef(cf.fit)
  cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
  cf.vcov <- vcov(cf.fit)
  cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]

  cf.val <- list(
    "coefficients" = cf.coef,
    "vcov" = cf.vcov
  )
  return(cf.val)
}

### Pretest
### FUNCTION: Point estimate, covariance matrix for treatment effect with  
###           using a single-sample, individual-level data via Pretest approach. 
### INPUT: Outcome formula, such as Y ~ X + D + g_2(D) + ... + g_k(D),
###        Treatment formula, such as D ~ X + Z + h_2(Z) + ... + h_k(Z),
###        verbose, which estimator is chosen.
###        where
###        Y, continuous, non-missing, numeric outcome vector (u by 1 vector)
###        D, continuous or discrete, non-missing, numeric treatment vector (n by 1 vector)        
###        Z, continuous or discrete, non-missing, numeric instrument matrix containing p_z 
###           instruments (n by p_z matrix)
###        X, optional continuous or discrete, non-missing, numeric matrix containing p_x 
###           covariates (n by p_x matrix)
###        g_2,...,g_k, functions of treatment variables
###        h_2,...,h_k, functions of instrument variables
###   
### OUTPUT: a list (a) coefficients (scalar numeric value:
###                                  the estimate of the treatment effect)
###                (b) vcov (numeric matrix: estimated covariance matrix of coefficients)                
###                (c) Hausman_statistic (scalar numeric value : 
###                                       test statistic of the validity of cf)                                  
###                (d) p_value (scalar numeric value : 
###                             asymptotic chi square p-value of Hausman statistic)
pretest <- function(outcome.formula,treatment.formula,verbose=FALSE){
  outcome.formula <- Formula(outcome.formula)
  treatment.formula <- Formula(treatment.formula)

  iv.temp<- attr(treatment.formula,"rhs")
  iv.formula <- as.formula(paste(format(outcome.formula),"|",iv.temp))
  tsls.model <- ivreg(iv.formula)
  iv.coef <- coef(tsls.model)
  iv.vcov <- vcov(tsls.model)

  cf.fstg <- lm(treatment.formula)
  e1 <- resid(cf.fstg)
  cf.formula <- as.formula(paste(format(outcome.formula),"+","e1"))
  cf.fit <- lm(cf.formula)

  cf.coef <- coef(cf.fit)
  cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
  cf.vcov <- vcov(cf.fit)
  cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]

  diff <- t(iv.coef-cf.coef)%*%solve(iv.vcov-cf.vcov)%*%(iv.coef-cf.coef)
  prob.larger.than.diff <- pchisq(diff, df = 1, lower.tail = FALSE)
  prob.larger.than.diff
  alpha <- 0.05
  if (prob.larger.than.diff>alpha) {
    if (verbose){
      cat("Level",alpha, "Pretest estimator is Control function estimator.","\n")
      }
    pretest.coef <- cf.coef
    pretest.vcov <- cf.vcov
  } else{
      if (verbose){
        cat("Level",alpha, "Pretest estimator is two stage least square estimator.","\n")
      }
    pretest.coef <- iv.coef
    pretest.vcov <- iv.vcov
  }

  pretest.val <- list(
    "coefficients" = pretest.coef,
    "vcov" = pretest.vcov,
    "Hausman_statistic" = diff,
    "p_value" = prob.larger.than.diff
  )
  return(pretest.val)
}
