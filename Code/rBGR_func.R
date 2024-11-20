library(extraDistr)
library(truncnorm)
Rcpp::sourceCpp("./code/expTail/rBGR_RcppGibbsFunc.cpp")

reScl.xi.eta<-function(param_list_){
  l.eta_<-param_list_$l.eta_
  l.xi_<-param_list_$l.xi_
  q<-ncol(l.eta_)
  for(k in 1:q){
    l.xibar=apply(matrix(abs(l.xi_[,k]),nrow=nrow(l.xi_)),1,mean)
    l.eta_[,k] = l.eta_[,k] * l.xibar
    l.xi_[,k] = l.xi_[,k] / l.xibar
  }
  return(list(l.xi_ = l.xi_, l.eta_ = l.eta_))
}

update_m<-function(l.xi_){
  l.m_=matrix(2 * rbinom(prod(dim(l.xi_)),1, prob= 1/(1+exp(-2*l.xi_)) ) -1,
              nrow=nrow(l.xi_))
  return(list(l.m_=l.m_))
}

update_nu<-function(a_nu, b_nu, param_list_){
  l.eta_ = param_list_$l.eta_
  l.Z_ = param_list_$l.Z_
  p<-dim(l.eta_)[1]
  q<-dim(l.eta_)[2]
  l.nu_<-matrix(1/rgamma(p*q,shape=a_nu+0.5, rate=b_nu+(l.eta_^2/2/l.Z_)),nrow=p,ncol=q)
  return(list(l.nu_ = l.nu_))
}

update_Z<-function(v0_, param_list_){
  rho_= param_list_$rho_
  l.eta_ = param_list_$l.eta_
  l.nu_ = param_list_$l.nu_
  p<-dim(l.eta_)[1]
  q<-dim(l.eta_)[2]
  l.odds<-rho_*sqrt(v0_)*exp((1-v0_)*l.eta_^2/(2*l.nu_*v0_))/(1-rho_)
  l.odds<-ifelse(l.odds==Inf,10^{5},l.odds)
  l.Z_<-matrix(ifelse(rbinom(p*q,1,l.odds/(1+l.odds)),1,v0_),nrow=p)
  return(list(l.Z_ = l.Z_))
}

update_rho<-function(a_rho, b_rho, l.Z_){
  res.rho <- rbeta(1,a_rho + sum(l.Z_==1),
                   b_rho + sum(l.Z_!=1))
  return(res.rho)
}

update_Drho<-function(ad_rho, bd_rho, Dj_){
  res_rho <- rbeta(1,ad_rho + sum(Dj_==1),bd_rho + sum(Dj_!=1))
  return(res_rho)
}

update_D<-function(dfj_, meanj_, rhoj_, Dj_, hyperVarj_, tail_ = c("Exp","Poly")){
  n<-length(dfj_)
  # df_rescl<-(df-mean(df))/sd(df)
  # set.seed(ite)
  D_old = Dj_
  acpt<-D_new<-rep(NA,n)
  # acpt<-rep(NA,n)
  D_newNormIdx<- rbinom(n,1,rhoj_)
  
  for(i in 1:n){
    # propose from prior
    if(D_newNormIdx[i]==1){
      Di_prop <- 1
    }else{
      if(tail_=="Exp"){
        Di_prop <- rexp(1, rate = hyperVarj_/2)
      }
      if(tail_=="Poly"){
        Di_prop <- extraDistr::rinvchisq(1, nu = 4, tau = hyperVarj_)
      }
    }
    #llh
    llh_diff <- dnorm(dfj_[i], meanj_[i], sd=sqrt(Di_prop),log=T) - dnorm(dfj_[i], mean = meanj_[i], sd=sqrt(D_old[i]),log=T)
    MH.u<-runif(1)
    if(llh_diff>log(MH.u)){
      # acpt[i]=1
      D_new[i]<-Di_prop
    }else{
      # acpt[i]=0
      D_new[i]<-D_old[i]
    }
  }
  return(D_new)
}

update_rndScly<-function(y_, mean_, dy_, rho_, hyperVar_, tail=c("Poly","Exp")){
  n<-length(y_)
  dy_new<-update_D(dfj_ = y_, meanj_ = mean_, rhoj_ = rho_, Dj_ = dy_, hyperVarj_ = hyperVar_, tail_ = tail)
  rho_new <- update_Drho(ad_rho=1, bd_rho = 1, Dj_ = dy_new)
  if(tail=="Exp"){
    hyperVar_new<- rgamma(1, shape = sum(dy_new!=1)+1, rate= sum(dy_new[dy_new!=1])/2+1.78)
  }else{
    hyperVar_new<- rgamma(1, shape= sum(dy_new!=1)*4/2, rate=(4/2)*sum(1/dy_new[dy_new!=1]))
  }
  return(list(d=dy_new, hyperVar = hyperVar_new, rho = rho_new))
}

update_rndSclx<-function(x_, DMat, rhoVec, hyperVarVec, tailVec){
  p<-ncol(x_)
  DMat_new<-matrix(NA,nr=nrow(DMat),ncol=p)
  rhoVec_new<-hyperVarVec_new<-rep(NA,p)
  for(j in 1:p){
    dj_new<-update_D(dfj_ = x_[,j], rhoj_ = rhoVec[j], Dj_ = DMat[,j], meanj_ = rep(0, nrow(DMat)),
                     hyperVarj_  =hyperVarVec[j], tail_ = tailVec[j])
    rho_j <- update_Drho(ad_rho=1, bd_rho = 1, Dj_ = dj_new)
    if(tailVec[j]=="Exp"){
      hyperVarj<- rgamma(1, shape = sum(dj_new!=1)+1, rate= sum(dj_new[dj_new!=1])/2+1.78)
    }
    if(tailVec[j]=="Poly"){
      hyperVarj<- rgamma(1, shape= sum(dj_new!=1)*4/2, rate=(4/2)*sum(1/dj_new[dj_new!=1]))
    }
    hyperVarVec_new[j]<-hyperVarj
    rhoVec_new[j]<-rho_j
    DMat_new[,j]<-dj_new
  }
  return(list(D=DMat_new, rho=rhoVec_new, hyperVar = hyperVarVec_new))
}

update_t_MH<-function(y_, x_, t_, theta_, d0_, tMax, tsd=1){
  n<-nrow(x_)
  t_old = t_
  betaOld<- theta_ * matrix(as.numeric(abs(theta_)>t_old),nrow=n)
  regressorOld <- apply(x_ * betaOld,1,sum)
  
  t_prop<-runif(1,max=tMax)
  betaProp<- theta_ * matrix(as.numeric(abs(theta_)>t_prop),nrow=n)
  regressorProp <- apply(x_ * betaProp,1,sum)
  
  llh_diff <- sum(dnorm(y_, mean=regressorProp, sd=sqrt(d0_),log=T)) - 
    sum(dnorm(y_, mean=regressorOld, sd=sqrt(d0_),log=T))
  MH.u<-runif(1)
  if(llh_diff>log(MH.u)){
      # acpt[i]=1
    t.res<-t_prop
  }else{
      # acpt[i]=0
    t.res<-t_old
  }
  return(t.res)
}

rBGR_mcmc_Int<-function(y, x, u, iteNum, burnIn, seed_, priorB, sigProp,  tail){
  n<-dim(x)[1]
  p_rest<-dim(x)[2]
  
  u.int<-cbind(1,u)
  q.int<-dim(u.int)[2]
  a_nu = 2; b_nu = 2;
  a_rho = 1;b_rho=1;
  ad_rho =1; bd_rho=1;
  a_s2=1; b_s2=1
  v0_=10^{-5}
  set.seed(seed_)
  
  param_list<-list()
  
  param_list$l.xi_<-matrix(rnorm(p_rest*q.int), nrow=p_rest, ncol=q.int)
  param_list$l.m_<-matrix(1, nrow=p_rest, ncol=q.int)
  param_list$l.eta_<-matrix(rnorm(p_rest*q.int),nrow=p_rest,ncol=q.int)
  param_list$l.nu_<-matrix(1,nrow=p_rest,ncol=q.int)
  param_list$l.Z_<-matrix(1,nrow=p_rest,ncol=q.int)
  
  # param_list$mu.xi_<-matrix(rnorm(p_rest), nrow=p_rest,ncol=1)
  # param_list$mu.m_<-matrix(1, nrow=p_rest,ncol=1)
  # param_list$mu.eta_<-matrix(rnorm(p_rest),nrow=p_rest,ncol=1)
  # param_list$mu.nu_<-matrix(1,nrow=p_rest,ncol=1)
  # param_list$mu.Z_<-matrix(1,nrow=p_rest,ncol=1)
  
  param_list$t_<-runif(1,0,1)
  param_list$rho_<-0.5
  param_list$priorB_<-priorB
  
  param_list$D_<-matrix(1,nrow=n,ncol=p_rest+1)
  param_list$D_rho <- runif(p_rest+1)
  param_list$hyperVar=rep(1, p_rest+1)
  
  u.int<-cbind(1,u)
  l.xi.Arr<-l.eta.Arr<-array(NA,dim=c(iteNum-burnIn,p_rest,q.int))
  mu.Z.Mat<-mu.xi.Mat<-mu.eta.Mat<-matrix(NA, nrow=iteNum-burnIn, ncol=p_rest)
  l.Z.Arr<-array(NA,dim=c(iteNum-burnIn,p_rest,q.int))
  tVec<-rep(NA,iteNum-burnIn)
  # D.Arr<-array(NA,dim=c(iteNum-burnIn,n,p_rest+1))
  llh<-rep(NA,iteNum)
  t.all<-rep(NA,iteNum)
  
  set.seed(seed_)
  for(ite in 1:iteNum){
    if(ite%%1000==0){
      print(ite)
    }
    l.eta.res<-Gibbs_l_etaC(y_ = y/sqrt(param_list$D_[,(p_rest+1)]), x_ = x/sqrt(param_list$D_[,-(p_rest+1)]), u_ = u.int, 
                            d0_ = rep(1,n), t_ = param_list$t_,
                            l_eta_ = param_list$l.eta_, 
                            l_xi_ = param_list$l.xi_, 
                            l_Z_ = param_list$l.Z_, 
                            l_nu_ = param_list$l.nu_)
    param_list$l.eta_<-l.eta.res
    
    l.xi.res<-Gibbs_l_xiC(y_ = y/sqrt(param_list$D_[,(p_rest+1)]), x_ = x/sqrt(param_list$D_[,-(p_rest+1)]), u_ = u.int, 
                          d0_ = rep(1,n), t_ = param_list$t_,
                          l_eta_ = param_list$l.eta_, 
                          l_xi_ = param_list$l.xi_, 
                          l_Z_ = param_list$l.Z_, 
                          l_nu_ = param_list$l.nu_, 
                          l_m_ = param_list$l.m_)
    param_list$l.xi_<-l.xi.res
    
    xi.eta.res<-reScl.xi.eta(param_list_ = param_list)
    param_list$l.xi_<-xi.eta.res$l.xi_
    param_list$l.eta_<-xi.eta.res$l.eta_
      
    m.res<-update_m(l.xi_ = param_list$l.xi_)
    param_list$l.m_<-m.res$l.m_
    
    Z.res<-update_Z(v0_=v0_, param_list_ = param_list)
    param_list$l.Z_<-Z.res$l.Z_
    
    nu.res<-update_nu(a_nu=a_nu, b_nu=b_nu, param_list_ = param_list)
    param_list$l.nu_<-nu.res$l.nu_
    
    rho.res<-update_rho(a_rho=a_rho, b_rho=b_rho, l.Z_ = param_list$l.Z_)
    param_list$rho_<-rho.res
    
    alpha_ <- param_list$l.xi_ * param_list$l.eta_
    theta_<- u.int %*% t(alpha_)
    t.res<-update_t_MH(y_ = y/sqrt(param_list$D_[,(p_rest+1)]), x_ = x/sqrt(param_list$D_[,-(p_rest+1)]), 
                       t_ = param_list$t_, theta_ = theta_, d0_ = rep(1,n), tMax =1 , tsd=1)
    param_list$t_<-t.res
    
    beta_<- theta_ * matrix(as.numeric(abs(theta_)>param_list$t_),nrow=n)
    regressor_ <- apply(x/sqrt(param_list$D_[,-(p_rest+1)]) * beta_,1,sum)
    
    d_y.ls<-update_rndScly(y_ = y, mean_ = rep(0,n), rho_ = param_list$D_rho[p_rest+1], 
                           dy_ = param_list$D_[,p_rest+1], hyperVar_ = param_list$hyperVar[p_rest+1],
                           tail = tail[p_rest+1])
    param_list$D_[,p_rest+1]<-d_y.ls$d
    param_list$D_rho[p_rest+1]<-d_y.ls$rho
    param_list$hyperVar[p_rest+1]<-d_y.ls$hyperVar
    
    d_x.ls<-update_rndSclx(x_ = x, DMat = param_list$D_[,1:p_rest],
                           rhoVec = param_list$D_rho[1:p_rest], 
                           hyperVarVec = param_list$hyperVar[1:p_rest],tailVec = tail[1:p_rest])
    param_list$D_[,1:p_rest]<-d_x.ls$D
    param_list$D_rho[1:p_rest]<-d_x.ls$rho
    param_list$hyperVar[1:p_rest]<-d_x.ls$hyperVar
    
    regressor_ <- apply(x/sqrt(param_list$D_[,1:p_rest]) * beta_,1,sum)
    llh[ite]<- sum(dnorm(y/sqrt(param_list$D_[,(p_rest+1)]), mean=regressor_, log=T))
    t.all[ite]<- param_list$t_
    if(ite >burnIn){
      l.xi.Arr[(ite - burnIn),,]=param_list$l.xi_
      l.eta.Arr[(ite - burnIn),,]=param_list$l.eta_
      l.Z.Arr[(ite - burnIn),,]=param_list$l.Z_
      # D.Arr[(ite - burnIn),,] = param_list$D_
      tVec[ite - burnIn] = param_list$t_
    }
  }
  return(list(llh=llh, t=tVec, t.all=t.all,
              l.eta=l.eta.Arr, l.xi=l.xi.Arr, l.Z=l.Z.Arr))
}
