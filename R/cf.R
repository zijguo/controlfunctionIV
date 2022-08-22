#' @title Control-Function
#' @description Implement the control function method for the inference of nonlinear treatment effects.
#' @param formula A formula describing the model to be fitted.
#' @param d1 The baseline treatment value.
#' @param d2 The target treatment value.
#' @details For example, the formula \code{Y ~ D + I(D^2)+X|Z+I(Z^2)+X} describes the models
#' \eqn{Y = \alpha_0 + D\beta_1 + D^2\beta_2 + X\phi + u}
#' and
#' \eqn{D = \gamma_0 + Z\gamma_1 + Z^2\gamma_2 + X\psi + v}.
#' Here, the outcome is \code{Y}, the endogenous variables is \code{D}, the baseline covariates are \code{X}, and the instrument variables are \code{Z}. The formula environment follows
#' that in the ivreg function in the AER package. The endogenous variable \code{D} must be in the first term of the formula for the outcome model.
#' If either one of \code{d1} or \code{d2} is missing or \code{NULL}, \code{CausalEffect} is calculated assuming that the baseline value \code{d1} is the median of the treatment and the target value \code{d2} is \code{d1+1}.
#' @return
#'    \code{cf} returns an object of class "cf", which is a list containing the following components:
#'    \item{\code{coefficients}}{The estimate of the coefficients in the outcome model.}
#'    \item{\code{vcov}}{The estimated covariance matrix of coefficients.}
#'    \item{\code{CausalEffect}}{The causal effect when the treatment changes from \code{d1} to \code{d2}.}
#'    \item{\code{CausalEffect.sd}}{The standard error of the causal effect estimator.}
#'    \item{\code{CausalEffect.ci}}{The 95\% confidence interval of the causal effect.}
#'
#' @export
#'
#'
#' @importFrom Formula as.Formula
#' @rawNamespace importFrom("stats", "coef", "delete.response", "lm", "model.matrix", "model.response", "pchisq", "printCoefmat", "pt", "qnorm", "quantile", "terms", "update", "vcov")
#' @examples
#' data("nonlineardata")
#' Y <- log(nonlineardata[,"insulin"])
#' D <- nonlineardata[,"bmi"]
#' Z <- as.matrix(nonlineardata[,c("Z.1","Z.2","Z.3","Z.4")])
#' X <- as.matrix(nonlineardata[,c("age","sex")])
#' cf.model <- cf(Y~D+I(D^2)+X|Z+I(Z^2)+X)
#' summary(cf.model)
#'
#' @references {
#' Guo, Z. and D. S. Small (2016), Control function instrumental variable estimation of nonlinear causal effect models, \emph{The Journal of Machine Learning Research} 17(1), 3448–3482. \cr
#' }
#'
#'


cf <- function(formula,d1 = NULL,d2 = NULL){
  if(!inherits(formula,"formula")) {
    stop("method is only for formula objects!")
  }
  # code gratefully lifted from ivreg() (package AER) and ivmodelFormula (package ivmodel).

  mf = match.call()
  m <- match(c("formula"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  formula <- as.Formula(formula)

  stopifnot(length(formula)[1] == 1L, length(formula)[2] %in%
              1:2)
  has_dot <- function(formula) inherits(try(terms(formula),silent = TRUE), "try-error")
  if (has_dot(formula)) {
    f1 <- formula(formula, rhs = 1)
    f2 <- formula(formula, lhs = 0, rhs = 2)
    if (!has_dot(f1) & has_dot(f2)) {
      formula <- as.Formula(f1, update(formula(formula, lhs = 0, rhs = 1), f2))
    }
  }
  mf$formula <- formula
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  Y <- model.response(mf, "numeric");
  n <- length(Y)
  Y = matrix(as.numeric(Y),length(Y),1)
  mt <- terms(formula)
  mtX <- terms(formula, rhs = 1)
  X <- model.matrix(mtX, mf)

  mtZ <- delete.response(terms(formula, rhs = 2))
  Z <- model.matrix(mtZ, mf)

  if("(Intercept)" %in% colnames(X)) {
    intercept=TRUE
    X = X[,!(colnames(X) %in% "(Intercept)"),drop=FALSE]
    Z = Z[,!(colnames(Z) %in% "(Intercept)"),drop=FALSE]
    if(dim(Z)[2] < 1) stop("There aren't any instruments!")
  } else{
    intercept=FALSE
  }

  # Parse X and Z into D, X, and Z

  whichD = !(colnames(X) %in% colnames(Z))
  d = X[,whichD,drop=FALSE]
  if (sum(!whichD) == 0) { # no covariates case
    if (intercept) { # intercept
      first.model <- lm(d~Z)
      e1 <- first.model$residuals[,1] # so the first term of d should be "D"
      second.model <- lm(Y~d+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef)[2:length(cf.coef)] = colnames(d)
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov)[2:length(cf.coef)] = colnames(cf.vcov)[2:length(cf.coef)] = colnames(d)

      if (!is.null(d1)&!is.null(d2)) {
        assign(colnames(d)[1],d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        assign(colnames(d)[1],d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[-1]
      } else {
        # cat("There is no baseline or target treatment value.\n")
        # cat("Calculate the causal effect of increasing of 1 in the median of treatment.\n")
        d1 <- median(d[,1])
        assign(colnames(d)[1],d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        d2 <- d1+1
        assign(colnames(d)[1],d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[-1]
      }

      CausalEffect.sd <- sqrt((d.target-d.base)%*%(cf.vcov[-1,-1])%*%cbind(d.target-d.base))
    } else { # no intercept
      first.model <- lm(d~0+Z)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~0+d+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef) = colnames(d)
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov) = colnames(cf.vcov) = colnames(d)
      if (!is.null(d1)&!is.null(d2)) {
        assign(colnames(d)[1],d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        assign(colnames(d)[1],d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef
      } else {
        # cat("There is no baseline or target treatment value.\n")
        # cat("Calculate the causal effect of increasing of 1 in the median of treatment.\n")
        d1 <- median(d[,1])
        assign(colnames(d)[1],d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        d2 <- d1+1
        assign(colnames(d)[1],d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef
      }
      CausalEffect.sd <- sqrt((d.target-d.base)%*%cf.vcov%*%cbind(d.target-d.base))
    }


  } else { # when there is covariates
    X = X[,!whichD,drop=FALSE]
    whichZ = !(colnames(Z) %in% colnames(X))
    Z = Z[,whichZ,drop=FALSE]
    if (intercept) { # intercept
      first.model <- lm(d~Z+X)
      e1 <- first.model$residuals[,1] # so the first term of d should be "D"
      second.model <- lm(Y~d+X+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef)[2:length(cf.coef)] = c(colnames(d),colnames(X))
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov)[2:length(cf.coef)] = colnames(cf.vcov)[2:length(cf.coef)] = c(colnames(d),colnames(X))
      d.name <- colnames(d)[1]
      if (!is.null(d1)&!is.null(d2)) {
        assign(d.name,d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        assign(d.name,d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[-c(1,seq(length(cf.coef)-ncol(X)+1,length =ncol(X)))]
      } else{
        # cat("There is no baseline or target treatment value.\n")
        # cat("Calculate the causal effect of increasing of 1 in the median of treatment.\n")
        d1 <- median(d[,1])
        assign(d.name,d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        d2 <- d1+1
        assign(d.name,d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[-c(1,seq(length(cf.coef)-ncol(X)+1,length =ncol(X)))]
      }
      CausalEffect.sd <- sqrt((d.target-d.base)%*%(cf.vcov[-c(1,seq(length(cf.coef)-ncol(X)+1,length =ncol(X))),
                                                          -c(1,seq(length(cf.coef)-ncol(X)+1,length =ncol(X)))])%*%cbind(d.target-d.base))
    } else { # no intercept
      first.model <- lm(d~0+Z+X)
      e1 <- first.model$residuals[,1]
      second.model <- lm(Y~0+d+X+e1)
      cf.coef <- second.model$coefficients
      cf.coef <- cf.coef[-which(names(cf.coef)=="e1")]
      names(cf.coef) = c(colnames(d),colnames(X))
      cf.vcov <- vcov(second.model)
      cf.vcov <- cf.vcov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]
      rownames(cf.vcov) = colnames(cf.vcov) = c(colnames(d),colnames(X))
      d.name <- colnames(d)[1]
      if (!is.null(d1)&!is.null(d2)) {
        assign(d.name,d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        assign(d.name,d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[(1:length(cf.coef))[-seq(length(cf.coef)-ncol(X)+1,length =ncol(X))]]
      } else{
        # cat("There is no baseline or target treatment value.\n")
        # cat("Calculate the causal effect of increasing of 1 in the median of treatment.\n")
        d1 <- median(d[,1])
        assign(d.name,d1)
        d.base <- sapply(colnames(d),function(x){eval(parse(text=x))})
        d2 <- d1+1
        assign(d.name,d2)
        d.target <- sapply(colnames(d),function(x){eval(parse(text=x))})
        CausalEffect <- (d.target-d.base)%*%cf.coef[(1:length(cf.coef))[-seq(length(cf.coef)-ncol(X)+1,length =ncol(X))]]
      }
      CausalEffect.sd <- sqrt((d.target-d.base)%*%(cf.vcov[(1:length(cf.coef))[-seq(length(cf.coef)-ncol(X)+1,length =ncol(X))],
                                                     (1:length(cf.coef))[-seq(length(cf.coef)-ncol(X)+1,length =ncol(X))]])%*%cbind(d.target-d.base))

    }

  }
  CausalEffect.ci <- c(CausalEffect-qnorm(0.975)*CausalEffect.sd,CausalEffect+qnorm(0.975)*CausalEffect.sd)
  out <- list(coefficients=cf.coef,vcov =cf.vcov,CausalEffect=CausalEffect,CausalEffect.sd =CausalEffect.sd,CausalEffect.ci =CausalEffect.ci,n=n,d1 = d1,d2 = d2)
  class(out) = 'cf'
  return(out)
}

#' Summary of cf
#'
#' @description Summary function for cf
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.cf<- function(object,...){
  cf <- object
  cat(rep("_", 30), "\n")
  cat("\nCoefficients of Control Function Estimators:\n\n")
  coeff <- cf$coefficients
  std <- sqrt(diag(cf$vcov))
  t.value <- abs(coeff/std)
  pr.t <- 1-pt(t.value,df = (cf$n)-1)
  cmat <- cbind(coeff,std,t.value,pr.t)
  colnames(cmat) <- c("Estimate", "Std.Err", "t value", "Pr(>|t|)")
  printCoefmat(cmat, digits = max(3L, getOption("digits") - 3L))
  cat(rep("_", 30), "\n")
  if (!is.null(cf$d2)&!is.null(cf$d1)) {
    cat("\nThe estimate causal effect when changing treatment from",cf$d1,"to",cf$d2,": ",cf$CausalEffect,"\n" )
    cat("\nStandard error of the estimate of the causal effect:",cf$CausalEffect.sd,"\n")
    cat("\n95% confidence interval for the causal effect: [",cf$CausalEffect.ci[1],",",cf$CausalEffect.ci[2],"]\n",sep = "")
  }
}



