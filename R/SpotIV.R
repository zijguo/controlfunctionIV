#' @title SpotIV method for causal inference in semi-parametric outcome model
#' @description Perform causal inference in the semi-parametric outcome model with possibly invalid IVs.
#' @param Y The outcome observation, a vector of length \eqn{n}.
#' @param D The treatment observation, a vector of length \eqn{n}.
#' @param Z The instrument observation of dimension \eqn{n \times p_z}.
#' @param X The covariates observation of dimension \eqn{n \times p_x}.
#' @param intercept Whether the intercept is included. (default = \code{TRUE})
#' @param invalid If TRUE, the method is robust to the presence of possibly invalid IVs; If FALSE, the method assumes all IVs to be valid. (default = \code{TRUE})
#' @param d1 A treatment value for computing CATE(d1,d2|w0).
#' @param d2 A treatment value for computing CATE(d1,d2|w0).
#' @param w0 A vector of the instruments and baseline covariates for computing CATE(d1,d2|w0).
#' @param M.est If \code{TRUE}, \code{M} is estimated based on BIC, otherwise \code{M} is specified by input value of \code{M}. (default = \code{TRUE})
#' @param M The dimension of indices in the outcome model, from 1 to 3. (default = \code{2})
#' @param bs.Niter The bootstrap resampling size for constructing the confidence interval.
#' @param bw A (M+1) by 1 vector bandwidth specification. (default = \code{NULL})

#' @return
#'     \code{SpotIV} returns an object of class "SpotIV", which "SpotIV" is a list containing the following components:
#'     \item{\code{betaHat}}{The estimate of the model parameter in front of the treatment.}
#'     \item{\code{cateHat}}{The estimate of CATE(d1,d2|w0).}
#'     \item{\code{cate.sdHat}}{The estimated standard error of cateHat.}
#'     \item{\code{SHat}}{The set of relevant IVs.}
#'     \item{\code{VHat}}{The set of relevant and valid IVs.}
#'     \item{\code{Maj.pass}}{The indicator that the majority rule is satisfied.}
#' @import dr
#' @import orthoDr
#' @importFrom stats binomial glm median pnorm sd
#' @export
#'

#' @examples
#'\donttest{
#' data("nonlineardata")
#' Y <- nonlineardata[,"CVD"]
#' D <- nonlineardata[,"bmi"]
#' Z <- as.matrix(nonlineardata[,c("Z.1","Z.2","Z.3","Z.4")])
#' X <- as.matrix(nonlineardata[,c("age","sex")])
#' d1 <- median(D)+1
#' d2 <- median(D)
#' w0 <- c(rep(0,4), 30, 1)
#' SpotIV.model <- SpotIV(Y,D,Z,X,invalid = TRUE,d1 =d1, d2 = d2,w0 = w0)
#' summary(SpotIV.model)
#' }
#'
#' @references {
#' Li, S., Guo, Z. (2020), Causal Inference for Nonlinear Outcome Models with Possibly Invalid Instrumental Variables, Preprint \emph{arXiv:2010.09922}.\cr
#' }
#'
#'
#'

SpotIV<- function(Y, D, Z, X=NULL, intercept=TRUE, invalid=TRUE,  d1, d2 , w0,
                  M.est=TRUE, M=2, bs.Niter=40, bw=NULL){
  stopifnot(!missing(Y),(is.numeric(Y) || is.logical(Y)),is.vector(Y)||(is.matrix(Y) || is.data.frame(Y)) && ncol(Y) == 1)
  stopifnot(all(!is.na(Y)))
  if (is.vector(Y)) {
    Y <- cbind(Y)
  }
  Y = as.numeric(Y)

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

  # Check d1,d2, and w0
  stopifnot(!missing(d1),!missing(d2),!missing(w0),is.numeric(d1),is.numeric(d2),is.numeric(w0))

  # All the other argument
  stopifnot(is.logical(intercept))
  stopifnot(is.logical(invalid))
  stopifnot(is.logical(M.est))

  pz<- ncol(Z)
  px<-0
  if(!is.null(X)){
    Z<-cbind(Z,X)
    X <- cbind(X)
    px<-ncol(X)
  }
  stopifnot(length(w0)==ncol(Z))
  n <- length(Y)
  Maj.pass=T

  if (invalid) {
    V <- NULL
  } else {
    V <- 1:pz
  }

  #first-stage regression
  if(intercept){
    gam.re<- lm(D ~ Z)
    gam.hat<-gam.re$coef[-1]
    var.gam <- n * as.matrix(vcov(gam.re)[2:(pz+1),2:(pz+1)])
  }else{
    gam.re<- lm(D ~ Z-1)
    gam.hat<-gam.re$coef
    var.gam <- n * vcov(gam.re)[1:pz,1:pz]
  }

  v.hat <- D-Z%*%gam.hat

  #voting and applying the majority rule
  if(is.null(V)){ # allowing for invalid IVs and do IV selection
    #get reduced-from
    SIR.re <-SIR.est(X.cov=cbind(Z,v.hat), Y, M= M, M.est=M.est)
    M <- ncol(SIR.re$theta.hat)
    Gam.hat<-as.matrix(SIR.re$theta.hat[1:ncol(Z),], ncol=M) #estimate Gamma
    ##voting
    Select.re<-Majority.test(n=n,ITT_Y=Gam.hat[1:pz,1], ITT_D=gam.hat[1:pz],
                         Cov.gGam=rbind(cbind(var.gam, diag(0,pz)),
                                            cbind(diag(0,pz), SIR.re$vGam[1:pz,1:pz])))
    SHat<-Select.re$SHat
    if(length(Select.re$VHat)< length(SHat)/2){
      message('Majority rule fails.','\n')
      Maj.pass=F
    }
    VHat <- as.numeric(Select.re$VHat)
    beta.hat<-sapply(1:M, function(m) median(Gam.hat[SHat,m]/gam.hat[SHat]))
    beta.hat <- matrix(beta.hat,nrow=1,ncol=M)
    pi.hat<- Gam.hat - gam.hat %*% beta.hat
  }else{ ### assume all IVs to be valid
    SIR.re <-SIR.est(X.cov=cbind(Z%*%gam.hat,Z[,-V],v.hat), Y, M= M, M.est=M.est)
    M <- ncol(SIR.re$theta.hat)
    beta.hat<-matrix(SIR.re$theta.hat[1,],nrow=1,ncol=M)
    pi.hat<-matrix(0,nrow=ncol(Z), ncol=M)
    pi.hat[V,]<-0
    pi.hat[-V,]<-SIR.re$theta.hat[2:(ncol(Z)-length(V)+1),]
    SHat <- V
    VHat <- V
  }
  ##estimate cace##
  if(is.null(bw)){
    bw.z<- apply(cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat), 2,
                 function(x) 0.9*n^(-1/6)*min(sd(x), (quantile(x,0.75)-quantile(x,0.25))/1.34))
  }

  asf.dw<-Spot.ASF.est(d1=d1,d2=d2, z0=w0, beta.hat= beta.hat, pi.hat= pi.hat,
                  Y, D, Z, v.hat=v.hat, bw.z=bw.z)
  cace.hat = asf.dw$cace.hat #cace
  bw.z=asf.dw$bw.z
  boot_b<-list()
  for(i in 1: bs.Niter){
    bootstrap_data<-cbind(Y,D,Z)[sample(n,n,replace=T),]
    boot_b[[i]]<-Spot.boot.fun(data=bootstrap_data, M=M, d1=d1,d2=d2, w0=w0, SHat=SHat,
                               bw.z=bw.z, V=V, intercept=intercept, pz=pz)
  }
  cace.sd<-sqrt(mean((unlist(lapply(boot_b, function(x) x[1]))-cace.hat)^2))
  VHat <- as.numeric(VHat)
  if (!is.null(colnames(Z))) {
    SHat = colnames(Z)[SHat]
    VHat = colnames(Z)[VHat]
  }
  SpotIV.model <- list(betaHat = beta.hat, cateHat=cace.hat, cate.sdHat= cace.sd,
                       SHat=SHat, VHat = VHat, Maj.pass=Maj.pass)
  class(SpotIV.model) <- "SpotIV"
  return(SpotIV.model)
}
