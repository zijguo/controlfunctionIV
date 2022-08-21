Probit.boot.fun<-function(data, pz,d1=NULL, d2=NULL,w0=NULL, SHat, invalid=invalid, intercept=intercept){
  Y<-data[,1]
  D<- data[,2]
  Z<-data[,-c(1,2)]
  gam.bs<-lm(D~Z-1)$coef
  v.bs <- D- Z%*%gam.bs
  if(invalid){
    Gam.bs.re<-glm(Y~cbind(Z,v.bs)-1, family=binomial(link='probit'))
    Gam.bs<-Gam.bs.re$coef[-length(Gam.bs.re$coef)]
    lam.bs<-Gam.bs.re$coef[length(Gam.bs.re$coef)]
    beta.bs<-median(Gam.bs[SHat]/gam.bs[SHat])
    pi.bs <- Gam.bs - gam.bs *beta.bs
    kappa.bs<-lam.bs-beta.bs
  }else{
    if(pz==ncol(Z)){
      coef.bs<-glm(Y~D+v.bs-1)$coef
      beta.bs<- coef.bs[1]
      pi.bs<-c(rep(0,pz))
      kappa.bs<-coef.bs[length(coef.bs)]
    }else{
      X<- as.matrix(Z[,-(1:pz)])# the last column is the intercept
      coef.bs<-glm(Y~D+X+v.bs-1)$coef
      beta.bs<- coef.bs[1]
      pi.bs<-c(rep(0,pz),coef.bs[2:(1+ncol(X))])
      kappa.bs<-coef.bs[length(coef.bs)]
    }

  }

  cace.bs<-NA
  if(!is.null(d1) && !is.null(d2) && !is.null(w0)){
    cace.bs = mean(pnorm(as.numeric(d1*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))-
      mean(pnorm(as.numeric(d2*beta.bs+w0%*%pi.bs) + v.bs*kappa.bs))
  }
  c(cace.bs, beta.bs)
}

Majority.test <- function(n, ITT_Y,ITT_D, Cov.gGam, tuning = 2.01) {
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

