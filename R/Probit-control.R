#' @title Causal inference in probit outcome models with possibly invalid IVs
#' @description Perform causal inference in the probit outcome model with possibly invalid IVs.
#'
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param invalid If \code{TRUE}, the method is robust to the presence of possibly invalid IVs; If \code{FALSE}, the method assumes all IVs to be valid. (default = \code{TRUE})
#' @param d1 A treatment value for computing CATE(d1,d2|w0).
#' @param d2 A treatment value for computing CATE(d1,d2|w0).
#' @param w0 A vector of the instruments and baseline covariates for computing CATE(d1,d2|w0).
#' @param bs.Niter The bootstrap resampling size for constructing the confidence interval.
#'
#' @return
#'     \code{ProbitControl} returns an object of class "SpotIV", which is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of the model parameter in front of the treatment.}
#'     \item{\code{beta.sdHat}}{The estimated standard error of betaHat.}
#'     \item{\code{cateHat}}{The estimate of CATE(d1,d2|w0).}
#'     \item{\code{cate.sdHat}}{The estimated standard deviation of \code{cateHat}.}
#'     \item{\code{SHat}}{The estimated set of relevant IVs.}
#'     \item{\code{VHat}}{The estimated set of relevant and valid IVs.}
#'     \item{\code{Maj.pass}}{The indicator that the majority rule is satisfied.}
#'
#' @importFrom stats binomial glm median pnorm sd
#' @export
#'
#'
#'
#' @examples
#' data("nonlineardata")
#' Y <- nonlineardata[,"CVD"]
#' D <- nonlineardata[,"bmi"]
#' Z <- as.matrix(nonlineardata[,c("Z.1","Z.2","Z.3","Z.4")])
#' X <- as.matrix(nonlineardata[,c("age","sex")])
#' d1 <- median(D)+1
#' d2 <- median(D)
#' w0 <- c(rep(0,4), 30, 1)
#' Probit.model <- ProbitControl(Y,D,Z,X,invalid = TRUE,d1 =d1, d2 = d2,w0 = w0)
#' summary(Probit.model)
#'
#'
#' @references {
#' Li, S., Guo, Z. (2020), Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables, Preprint \emph{arXiv:2010.09922}.\cr
#' }
#'
#'




ProbitControl<- function(Y, D, Z, X=NULL, intercept=TRUE, invalid=TRUE,
                         d1=NULL, d2=NULL , w0=NULL, bs.Niter=40){
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),is.vector(Y)||(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  if (is.vector(Y)) {
    Y <- cbind(Y)
  }
  Y = as.numeric(Y)
  stopifnot(length(table(Y))==2)
  # Check D
  stopifnot(!missing(D),(is.numeric(D) || is.logical(D)),is.vector(D)||(is.matrix(D) || is.data.frame(D)) && ncol(D) == 1)
  stopifnot(all(!is.na(D)))
  if (is.vector(D)) {
    D <- cbind(D)
  }
  D = as.numeric(D)

  # Check Z
  stopifnot(!missing(Z),(is.numeric(Z) || is.logical(Z)),(is.vector(Z) || is.matrix(Z)))
  stopifnot(all(!is.na(Z)))
  if (is.vector(Z)) {
    Z <- cbind(Z)
  }
  # Check dimesions
  stopifnot(nrow(Y) == nrow(D), nrow(Y) == nrow(Z))

  # Check X, if present
  if(!missing(X)) {
    stopifnot((is.numeric(X) || is.logical(X)),(is.vector(X))||(is.matrix(X) && nrow(X) == nrow(Z)))
    stopifnot(all(!is.na(X)))
  }

  # All the other argument
  stopifnot(is.logical(invalid))
  stopifnot(is.logical(intercept))

  pz<- ncol(Z)
  px<-0
  if(!is.null(X)){
    Z<-cbind(Z,X)
    X <- cbind(X)
    px<-ncol(X)
  }
  stopifnot(length(w0)==ncol(Z))

  if(intercept){ Z<-cbind(Z,1) }
  n <- length(Y)
  Maj.pass=T

  #first stage
  gam.re<-lm(D ~ Z -1)
  gam.hat<-gam.re$coef
  v.hat<-D-Z%*%gam.hat
  sig.v.hat<- mean(v.hat^2)
  gam.cov<-n*vcov(gam.re)[1:pz,1:pz]

  if(invalid){ #
    #reduced form
    Gam.re<-glm(Y~cbind(Z,v.hat)-1, family=binomial(link='probit'))
    Gam.hat<-Gam.re$coef[-length(Gam.re$coef)]
    lam.hat<-Gam.re$coef[length(Gam.re$coef)]
    Gam.cov<-as.matrix(n*vcov(Gam.re)[1:pz,1:pz]+lam.hat^2*gam.cov)
    Cov.gGam<-rbind(cbind(gam.cov,lam.hat^2*gam.cov),
                    cbind(lam.hat^2*Gam.cov, Gam.cov))
    #applying the majority rule
    Select.re<-Majority.test(n=n,ITT_Y=Gam.hat[1:pz], ITT_D=gam.hat[1:pz], Cov.gGam=Cov.gGam)
    SHat<-Select.re$SHat
    VHat <- Select.re$VHat
    if(length(Select.re$VHat)<= length(SHat)/2){
      message('Majority rule fails.','\n')
      Maj.pass=F
    }
    beta.hat<-median(Gam.hat[SHat]/gam.hat[SHat])
    pi.hat<- Gam.hat - gam.hat * beta.hat
    kappa.hat<-lam.hat-beta.hat #coef of v.hat
  }else{ # valid IV method
    SHat=1:pz
    if((!intercept) && px==0){
      coef.re<-glm(Y~D+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz))
      kappa.hat<-coef.re[length(coef.re)]
    }else{
      if(intercept && (px>0)){
        X<-cbind(X,1)
      }else if(intercept && px==0){
        X<-matrix(rep(1,n),ncol=1)
      }
      coef.re<-glm(Y~D+X+v.hat-1)$coef
      beta.hat<- coef.re[1]
      pi.hat<-c(rep(0,pz),coef.re[2:(1+ncol(X))])
      kappa.hat<-coef.re[length(coef.re)]
    }
    VHat = SHat
  }

  cace.hat<-NA; cace.sd<-NA
  if(!is.null(d1) && !is.null(d2) && !is.null(w0)){ #compute cate and its sd
    if(intercept){ w0=c(w0,1)}
    cace.hat = mean(pnorm(as.numeric(d1*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))-
      mean(pnorm(as.numeric(d2*beta.hat+w0%*%pi.hat) + v.hat*kappa.hat))
    #bootstrap for sd
    bs.lst<-list()
    for(i in 1:bs.Niter){
      sample.true <- F
      while(sample.true==F)
        {tryCatch({
            bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),];
            bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, d1=d1,d2=d2,
                            w0=w0, SHat=SHat, invalid=invalid, intercept=intercept);
            sample.true<-T;
            },error=function(e){
            },finally={})
        }
      }
    cace.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[1]))-cace.hat)^2))
    beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
  }else{ #only compute the sd of beta.hat
    bs.lst<-list()
    for(i in 1:bs.Niter){
      sample.true <- F
      while(sample.true==F)
      {tryCatch({
        bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),];
        bs.lst[[i]]<-Probit.boot.fun(data=bootstrap_data, pz=pz, d1=d1,d2=d2,
                                     w0=w0, SHat=SHat, invalid=invalid, intercept=intercept);
        sample.true<-T;
      },error=function(e){
      },finally={})
      }
      beta.sd<-sqrt(mean((unlist(lapply(bs.lst, function(x) x[2]))-beta.hat)^2))
    }
  }
  VHat <- as.numeric(VHat)
  if (!is.null(colnames(Z))) {
    SHat = colnames(Z)[SHat]
    VHat = colnames(Z)[VHat]
  }
  Probit.model <- list(betaHat=beta.hat, beta.sdHat=beta.sd, cateHat=cace.hat, cate.sdHat= cace.sd, SHat=SHat, VHat = VHat, Maj.pass=Maj.pass)
  class(Probit.model) <- "SpotIV"
  return(Probit.model)
}

