#' @title Prestest estimator
#' @description This function implements the pretest estimator by comparing the control function and the TSLS estimators.
#'
#' @param formula A formula describing the model to be fitted.
#' @param alpha The significant level. (default = \code{0.05})
#'
#' @details For example, the formula \code{Y ~ D + I(D^2)+X|Z+I(Z^2)+X} describes the models
#' \eqn{Y = \alpha_0 + D\beta_1 + D^2\beta_2 + X\phi + u}
#' and
#' \eqn{D = \gamma_0 + Z\gamma_1 + Z^2\gamma_2 + X\psi + v}.
#' Here, the outcome is \code{Y}, the endogenous variables is \code{D}, the baseline covariates are \code{X}, and the instrument variables are \code{Z}. The formula environment follows
#' that in the ivreg function in the AER package. The endogenous variable \code{D} must be in the first term of the formula for the outcome model.
#' @return
#'    \code{pretest} returns an object of class "pretest", which is a list containing the following components:
#'    \item{\code{coefficients}}{The estimate of the coefficients in the outcome model.}
#'    \item{\code{vcov}}{The estimated covariance matrix of coefficients.}
#'    \item{\code{Hausman.stat}}{The Hausman test statistic used to test the validity of the extra IV generated by the control function.}
#'    \item{\code{p.value}}{The p-value of the Hausman test.}
#'    \item{\code{cf.check}}{The indicator that the extra IV generated by the control function is valid.}
#' @export
#'
#'
#' @examples
#' data("nonlineardata")
#' Y <- log(nonlineardata[,"insulin"])
#' D <- nonlineardata[,"bmi"]
#' Z <- as.matrix(nonlineardata[,c("Z.1","Z.2","Z.3","Z.4")])
#' X <- as.matrix(nonlineardata[,c("age","sex")])
#' pretest.model <- pretest(Y~D+I(D^2)+X|Z+I(Z^2)+X)
#' summary(pretest.model)
#'
#' @references {
#' Guo, Z. and D. S. Small (2016), Control function instrumental variable estimation of nonlinear causal effect models, \emph{The Journal of Machine Learning Research} 17(1), 3448–3482. \cr
#' }
#'
pretest <- function(formula,alpha = 0.05){


  tsls.model <- AER::ivreg(formula)
  n <- tsls.model$n
  iv.coef <- coef(tsls.model)
  iv.vcov <- vcov(tsls.model)

  cf.model <- cf(formula)

  cf.coef <- cf.model$coefficients
  cf.vcov <- cf.model$vcov

  diff <- t(iv.coef-cf.coef)%*%solve(iv.vcov-cf.vcov)%*%(iv.coef-cf.coef)
  prob.larger.than.diff <- pchisq(diff, df = 1, lower.tail = FALSE)
  prob.larger.than.diff
  if (prob.larger.than.diff>alpha) {
    cf.check = TRUE
    pretest.coef <- cf.coef
    pretest.vcov <- cf.vcov
  } else{

    cf.check = FALSE
    pretest.coef <- iv.coef
    pretest.vcov <- iv.vcov
  }

  pretest.val <- list(
    coefficients = pretest.coef,
    vcov = pretest.vcov,
    Hausman.stat = diff,
    p.value = prob.larger.than.diff,
    alpha = alpha,
    cf.check = cf.check,
    n = n
  )
  class(pretest.val) = 'pretest'
  return(pretest.val)
}
#' Summary of pretest
#'
#' @description Summary function for pretest
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.pretest<- function(object,...){
  pretest <- object
  cat(rep("_", 30), "\n")
  cat("\nHausman Statistic : ",pretest$Hausman.stat,"\n")
  cat("\nP value = ",pretest$p.value,"\n")
  if (pretest$cf.check) {
    cat("\nH0 : Augmented instrumental variables from Control function are valid, is not rejected.","\n")
    cat("\nLevel",pretest$alpha, "Pretest estimator is Control function estimator.","\n")
  } else{
    cat("\nH0 : Augmented instrumental variables from Control function are valid, is rejected.","\n")
    cat("\nLevel",pretest$alpha, "Pretest estimator is Two Stage Least Square estimator.","\n")
  }
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of Pretest Estimators:\n\n")
  coeff <- pretest$coefficients
  std <- sqrt(diag(pretest$vcov))
  t.value <- abs(coeff/std)
  pr.t <- 1-pt(t.value,df = (pretest$n)-1)
  cmat <- cbind(coeff,std,t.value,pr.t)
  colnames(cmat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  printCoefmat(cmat, digits = max(3L, getOption("digits") - 3L))
}
