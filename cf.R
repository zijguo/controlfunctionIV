#' @title Control-Function
#' @description This function implements the control function method for estimation and inference of nonlinear treatment effects as developed in Guo and Small (2016)
#' @param outcome.formula  formula specification(s) of the regression relationship between outcome, covariates and the endogeneity variables, such as \code{Y ~ X + D + g_2(D) + ... + g_k(D)}.
#' @param treatment.formula formula specification(s) of the regression relationship between endogeneity variables, covariates and the instrumental variables, such as \code{D ~ X + Z + h_2(Z) + ... + h_k(Z)}.
#'
#' @return
#'    \item{\code{coefficients}}{scalar numeric value: the estimate of the treatment effect }
#'    \item{\code{vcov}}{numeric matrix: estimated covariance matrix of coefficients}
#' @export
#'
#'
#' @importFrom stats as.formula lm pchisq quantile resid vcov
#'
#' @examples
#'
#' ### Generate the data
#' library(MASS)
#' library(RobustIV)
#' n <- 10000
#' mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
#' X <- rnorm(n); Z <- rnorm(n)
#' err <- mvrnorm(n,mu=mu,Sigma = V)
#' u1 <- err[,1]; v2 <- err[,2]
#' D <- 1+X/8+Z/3+Z^2/8+v2
#' Y <- 1+X+10*D +10*D^2+u1
#'
#' ### Implement the control function method
#' cf(Y~X+D+I(D^2),D~X+Z+I(Z^2))
#'
#'


cf <- function(outcome.formula,treatment.formula){
  outcome.temp <- gsub(pattern='\\s',replacement="",x=outcome.formula)
  treatment.temp <- gsub(pattern='\\s',replacement="",x=treatment.formula)

  cf.fstg <- lm(treatment.formula)
  e1 <- resid(cf.fstg)
  cf.formula <- as.formula(paste(outcome.temp[2],"~",outcome.temp[3],"+","e1"))
  cf.fit <- lm(cf.formula)

  cf.coef <- coef(cf.fit)
  cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
  cf.vcov <- vcov(cf.fit)
  cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]

  return(list(coefficients = cf.coef, vcov = cf.vcov))
}




#' @title Prestest estimator
#' @description This function implements the pretest approach  for estimation and inference of nonlinear treatment effects as developed in Guo and Small (2016)
#'
#' @param outcome.formula formula specification(s) of the regression relationship between outcome, covariates and the endogeneity variables, such as \code{Y ~ X + D + g_2(D) + ... + g_k(D)}.
#' @param treatment.formula formula specification(s) of the regression relationship between endogeneity variables, covariates and the instrumental variables, such as \code{D ~ X + Z + h_2(Z) + ... + h_k(Z)}.
#' @return
#'    \item{\code{coefficients}}{scalar numeric value: the estimate of the treatment effect}
#'    \item{\code{vcov}}{numeric matrix: estimated covariance matrix of coefficients}
#'    \item{\code{Hausman.stat}}{scalar numeric value : test statistic of the validity of control function}
#'    \item{\code{p.value}}{scalar numeric value : asymptotic chi square p-value of Hausman statistic}
#' @export
#'
#' @importFrom AER ivreg
#'
#' @examples
#'
#' ### Generate the data
#' library(MASS)
#' library(RobustIV)
#' n <- 10000
#' mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
#' X <- rnorm(n); Z <- rnorm(n)
#' err <- mvrnorm(n,mu=mu,Sigma = V)
#' u1 <- err[,1]; v2 <- err[,2]
#' D <- 1+X/8+Z/3+Z^2/8+v2
#' Y <- 1+X+10*D +10*D^2+u1
#'
#' ### Implement the pretest method
#' pretest(Y~X+D+I(D^2),D~X+Z+I(Z^2))
#'
#'
pretest <- function(outcome.formula,treatment.formula){

  outcome.temp <- gsub(pattern='\\s',replacement="",x=outcome.formula)
  treatment.temp <- gsub(pattern='\\s',replacement="",x=treatment.formula)
  iv.formula <- as.formula(paste(outcome.temp[2],"~",outcome.temp[3],"|",treatment.temp[3]))

  tsls.model <- ivreg(iv.formula)
  iv.coef <- coef(tsls.model)
  iv.vcov <- vcov(tsls.model)

  cf.fit <- cf(outcome.formula,treatment.formula)

  cf.coef <- cf.fit$coefficients
  # cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
  cf.vcov <- cf.fit$vcov
  # cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]

  diff <- t(iv.coef-cf.coef)%*%solve(iv.vcov-cf.vcov)%*%(iv.coef-cf.coef)
  prob.larger.than.diff <- pchisq(diff, df = 1, lower.tail = FALSE)
  prob.larger.than.diff
  alpha <- 0.05
  if (prob.larger.than.diff>alpha) {

    cat("Level",alpha, "Pretest estimator is Control function estimator.","\n")
    pretest.coef <- cf.coef
    pretest.vcov <- cf.vcov
  } else{

    cat("Level",alpha, "Pretest estimator is two stage least square estimator.","\n")
    pretest.coef <- iv.coef
    pretest.vcov <- iv.vcov
  }

  pretest.val <- list(
    coefficients = pretest.coef,
    vcov = pretest.vcov,
    Hausman.stat = diff,
    p.value = prob.larger.than.diff
  )
  return(pretest.val)
}

#' @title Two stage least square estimator
#' @description This function implements the two stage least square method for estimation and inference of (non)linear treatment effects. This function just works like a \code{ivreg} in AER.
#'
#' @param outcome.formula formula specification(s) of the regression relationship between outcome, covariates and the endogeneity variables, such as \code{Y ~ X + D + g_2(D) + ... + g_k(D)}.
#' @param treatment.formula formula specification(s) of the regression relationship between endogeneity variables, covariates and the instrumental variables, such as \code{D ~ X + Z + h_2(Z) + ... + h_k(Z)}.
#'
#' @return
#'    \item{\code{coefficients}}{scalar numeric value: the estimate of the treatment effect}
#'    \item{\code{vcov}}{numeric matrix: estimated covariance matrix of coefficients}
#' @export
#' @importFrom AER ivreg
#'
#' @examples
#' #' ### Generate the data
#' library(MASS)
#' library(RobustIV)
#' n <- 10000
#' mu <- rep(0,2); V <- cbind(c(1,0.5),c(0.5,1))
#' X <- rnorm(n); Z <- rnorm(n)
#' err <- mvrnorm(n,mu=mu,Sigma = V)
#' u1 <- err[,1]; v2 <- err[,2]
#' D <- 1+X/8+Z/3+Z^2/8+v2
#' Y <- 1+X+10*D +10*D^2+u1
#'
#' ### Implement the two stage least square method
#' tsls(Y~X+D+I(D^2),D~X+Z+I(Z^2))
#'
tsls <- function(outcome.formula,treatment.formula){

  outcome.temp <- gsub(pattern='\\s',replacement="",x=outcome.formula)
  treatment.temp <- gsub(pattern='\\s',replacement="",x=treatment.formula)
  iv.formula <- as.formula(paste(outcome.temp[2],"~",outcome.temp[3],"|",treatment.temp[3]))

  tsls.model <- ivreg(iv.formula)
  iv.coef <- coef(tsls.model)
  iv.vcov <- vcov(tsls.model)


  return(list(coefficients = iv.coef, vcov = iv.vcov))
}
