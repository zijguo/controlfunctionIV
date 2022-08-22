#' Summary of SpotIV
#'
#' @description Summary function for SpotIV and ProbitControl
#' @keywords internal
#' @return No return value, called for summary.
#' @export
summary.SpotIV <- function(object,...){
  SpotIV <- object
  cat("\nRelevant Instruments:", SpotIV$SHat, "\n");
  cat("\nValid Instruments:", SpotIV$VHat,"\n","\nThus, Majority rule",ifelse(SpotIV$Maj.pass,"holds.","does not hold."), "\n");
  cat(rep("_", 30), "\n")
  cat("\nBetaHat:",SpotIV$betaHat,"\n");
  if (!is.null(SpotIV$beta.sdHat)) {
    cat("\nStandard error of BetaHat:",SpotIV$beta.sdHat,"\n");
    ci.beta <- c(SpotIV$betaHat-qnorm(0.975)*SpotIV$beta.sdHat,SpotIV$betaHat+qnorm(0.975)*SpotIV$beta.sdHat)
    cat("\n95% Confidence Interval for Beta: [", ci.beta[1], ",", ci.beta[2], "]", "\n", sep = '');
  }
  cat("\nCATEHat:",SpotIV$cateHat,"\n");
  cat("\nStandard error of CATEHat:",SpotIV$cate.sdHat,"\n");
  ci <- c(SpotIV$cateHat-qnorm(0.975)*SpotIV$cate.sdHat,SpotIV$cateHat+qnorm(0.975)*SpotIV$cate.sdHat)
  cat("\n95% Confidence Interval for CATE: [", ci[1], ",", ci[2], "]", "\n", sep = '');
}


SIR.est<- function(X.cov,Y, M=2, M.est=TRUE){
  p<- ncol(X.cov)
  n<-length(Y)
  if(M.est){
    nslice=ifelse(length(table(Y))==2,2,8)
    SIR.re<-dr::dr(Y ~ X.cov -1, method='sir', numdir=M, nslices=nslice)
    evalues=SIR.re$evalues
    nobs.slice <- median(SIR.re$slice.info$slice.sizes)
    M <- which.max(
      sapply(1:M, function(m) sum(log(evalues[(m+1):p]+1)-evalues[(m+1):p])*n/2-log(n)*m*(2*p-m+1)/4/nobs.slice))
  }
  if(length(table(Y))==2){
    init.re<- glm(Y~X.cov-1,family=binomial(link='logit'))
  }else{
    init.re<-lm(Y~X.cov-1)
  }
  Gam.init<-init.re$coef
  vGam<-vcov(init.re)*n
  if(M==1){
    theta.hat <- orthoDr::orthoDr_reg(x=X.cov, y = Y, B.initial =as.matrix(Gam.init,ncol=1), ndr=1)$B
  }else{
    theta.hat<-dr::dr(Y ~ X.cov -1, method='sir', numdir=M)$evectors[,1:M] #using a faster computation
  }
  list(theta.hat=theta.hat, vGam=vGam)
}

Majority.test <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01, majority=TRUE) {
  Var.comp.est <- function(Cov.mat, gam.hat, j){
    diag(Cov.mat) + (gam.hat/gam.hat[j])^2 * Cov.mat[j,j] - 2*gam.hat/gam.hat[j] * Cov.mat[j,]
  }
  # Check ITT_Y and ITT_D
  stopifnot(!missing(ITT_Y),!missing(ITT_D),length(ITT_Y) == length(ITT_D))
  stopifnot(all(!is.na(ITT_Y)),all(!is.na(ITT_D)))
  ITT_Y = as.numeric(ITT_Y); ITT_D = as.numeric(ITT_D)
  # Check Sigmas
  stopifnot(!missing(Cov.gGam))

  # Other Input check
  stopifnot(is.numeric(tuning),length(tuning) == 1, tuning >=2)

  # Constants
  pz = length(ITT_Y)

  # First Stage
  Var.gam.hat <- diag(Cov.gGam)[1:pz]
  SHat = (1:pz)[abs(ITT_D) >= (sqrt(Var.gam.hat) * sqrt(tuning*log(pz)/n))]
  if(length(SHat) == 0) {
    warning("First Thresholding Warning: IVs individually weak. TSHT with these IVs will give misleading CIs, SEs, and p-values. Use more robust methods.")
    warning("Defaulting to treating all IVs as strong.")
    SHat= 1:pz
  }
  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  # Second Stage
  # pi.candidate is the estimated value of pi across different candidates
  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat
  for(j in SHat) {
    beta.j = ITT_Y[j] / ITT_D[j]
    pi.j = ITT_Y - ITT_D * beta.j
    #compute three components in eq(33)
    Sig1.j <- Var.comp.est(as.matrix(Cov.gGam[1:pz,1:pz]), ITT_D, j)
    Sig2.j <- Var.comp.est(as.matrix(Cov.gGam[(pz+1):(2*pz),(pz+1):(2*pz)]), ITT_D, j)
    Sig3.j <-  Var.comp.est(as.matrix(Cov.gGam[1:pz,(pz+1):(2*pz)]), ITT_D, j)
    sigmasq.j <- beta.j^2 *Sig1.j +  Sig2.j - 2* beta.j * Sig3.j
    PHat.bool.j <- abs(pi.j) <= sqrt(sigmasq.j) * tuning * sqrt(log(pz)/n)
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  # Voting
  diag(VHats.bool) <- rep(TRUE, nCand)
  VM = rowSums(VHats.bool)
  # cat(VM,'\n')
  VHat = rownames(VHats.bool)[VM > (0.5 * length(SHat))] # Majority winners
  return(list(VHat = VHat,SHat=SHat))
}



Spot.boot.fun<-function(data, M, d1, d2,w0, SHat,
                        bw.z=NULL, V=NULL, intercept, pz){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  if(intercept){
    gam.re<- lm(D ~ Z)
  }else{
    gam.re<- lm(D ~ Z-1)
  }
  gam.bs<-gam.re$coef[1:ncol(Z)]
  v.bs <- D- Z%*%gam.bs
  if(is.null(V)){
    SIR.bs.re <- SIR.est(cbind(Z,v.bs), Y, M=M, M.est=F)
    Gam.bs<-as.matrix(SIR.bs.re$theta.hat[1:ncol(Z),], ncol=M)
    beta.bs<-sapply(1:M, function(m) median(Gam.bs[SHat,m]/gam.bs[SHat]))
    pi.bs <- Gam.bs - gam.bs %*% matrix(beta.bs,nrow=1,ncol=M)
  }else{
    SIR.bs <-SIR.est(X.cov=cbind(Z%*%gam.bs,Z[,-V],v.bs), Y, M= M, M.est=F)
    M <- ncol(SIR.bs$theta.hat)
    beta.bs<-SIR.bs$theta.hat[1,]
    pi.bs<-matrix(0,nrow=ncol(Z), ncol=M)
    pi.bs[V,]<-0
    pi.bs[-V,]<-SIR.bs$theta.hat[2:(ncol(Z)-length(V)+1),]
  }
  asf.dw<-Spot.ASF.est(d1=d1,d2=d2, z0=w0, beta.hat= beta.bs, pi.hat= pi.bs,
                       Y=Y, D=D, Z=Z, v.hat=D-Z%*%gam.bs, bw.z=bw.z)
  asf.dw$cace.hat

}

box.ker<-function(xx, X, h){
  apply(X, 1, function(xi) max(abs(xi-xx)/h)<1)/prod(h)
}

Spot.ASF.est <- function(d1, d2, z0, beta.hat, pi.hat,  Y, D, Z, v.hat, bw.z=NULL){
  beta.hat=as.vector(beta.hat)
  M=length(beta.hat)
  n=length(D)
  v.hat <- as.matrix(v.hat,ncol=1)
  D<-as.matrix(D,ncol=1)
  index <- cbind(D%*%beta.hat+ Z%*%pi.hat,v.hat)

  #cace
  index1.z <- cbind(apply(matrix(d1%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat)
  index2.z <- cbind(apply(matrix(d2%*%beta.hat+z0%*%pi.hat, nrow=1), 2, function(x) rep(x,n)),v.hat)
  #cat(index1.z[1,1],index2.z[1,1],'\n')
  q1<-quantile(index[,1], 0.975)
  q2<-quantile(index[,1], 0.025)

  if(sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0 & sum(index1.z[,1]<= q1 & index1.z[,1]>= q2)>0){
    index1.z<- index1.z[index1.z[,1]<= q1 & index1.z[,1]>= q2,]
    index2.z<- index2.z[index2.z[,1]<= q1 & index2.z[,1]>= q2,]
  }
  asf.dw1<- NA
  asf.dw2<-NA
  bw.z<-bw.z/1.5
  while(is.na(asf.dw1) | is.na(asf.dw2)){
    bw.z<-bw.z*1.5
    if(M==1){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] ])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2]])), na.rm=T)
    }else if(M==2){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3]])), na.rm=T)
    }else if(M==3){
      asf.dw1<-mean(apply(index1.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
      asf.dw2<-mean(apply(index2.z,1, function(xx) mean(Y[abs(xx[1]-index[,1])<=bw.z[1] & abs(xx[2]-index[,2])<=bw.z[2] & abs(xx[3]-index[,3])<=bw.z[3] & abs(xx[4]-index[,4])<=bw.z[4]])), na.rm=T)
    }
  }


  cace.hat <- asf.dw1-asf.dw2


  list(cace.hat = cace.hat, ace.hat=0, bw=bw.z, bw.z=bw.z)
}

