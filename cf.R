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
  cf.vcov <- cf.cov[-which(rownames(cf.vcov)=="e1"),-which(colnames(cf.vcov)=="e1")]

  cf.val <- list(
    "coefficients" = cf.coef,
    "vcov" = cf.vcov,
  )
  return(cf.val)
}

pretest <- function(outcome.formula,treatment.formula){
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
    cat("Level",alpha, "Pretest estimator is Control function estimator.")
    pretest.coef <- cf.coef
    pretest.vcov <- cf.vcov
  } else{
    cat("Level",alpha, "Pretest estimator is two stage least square estimator.")
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
