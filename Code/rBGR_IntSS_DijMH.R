library(extraDistr)
library(truncnorm)
Rcpp::sourceCpp("./rBGR_RcppGibbsFunc.cpp")

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

update_d0<-function(y_, x_, u_, param_list_){
  n<-length(y_)
  dgf0_=param_list_$dgf0_
  xi_<-param_list_$l.xi_
  eta_<-param_list_$l.eta_
  alpha_ <- xi_ * eta_
  t_=param_list_$t_
  sigma2_=param_list_$sigma2_
  d0.res<-rep(NA,n)
  theta_<- u_ %*% t(alpha_)
  beta_ <-theta_ * matrix(as.numeric(abs(theta_)>t_),nrow=n)
  mean_<-apply(x_ * beta_,1,sum)
  for(i in 1:n){
    d0.res[i]<-extraDistr::rinvchisq(1,dgf0_+1, (dgf0_ * sigma2_+ (y_[i]-mean_[i])^2)/(dgf0_+1))
  }
  return(d0.res)
}

update_sigma2<-function(dgf0_,d0_){
  n<-length(d0_)
  return(rgamma(1,shape = n*dgf0_/2, rate=dgf0_*sum(1/d0_)/2))
}

# update_Dj_nodeWise<-function(xj_, dgfj_){
#   xj.vec<-as.numeric(as.matrix(xj_))
#   n<-length(xj.vec)
#   xj_rescl<-(xj.vec-mean(xj.vec))/sd(xj.vec)
#   Dj_<-rep(NA,n)
#   for(i in 1:n){
#     Dj_[i] <- extraDistr::rinvchisq(1,dgfj_+1,(dgfj_+xj_rescl[i]^2)/(dgfj_+1))
#   }
#   return(Dj_)
# }

update_Dj_MH<-function(xj_, u_, j_, param_list_){
  xj.vec<-as.numeric(as.matrix(xj_))
  n<-length(xj.vec)
  xj_rescl<-(xj.vec-mean(xj.vec))/sd(xj.vec)
  
  Dj_old = param_list_$Dj_[,j_]
  prop_dgfj_<-param_list_$prop_dgfj_[j_]
  prior_dgfj_<-param_list_$prior_dgfj_[j_]
  Dj_rho_<- param_list_$Dj_rho_[j_]
  
  Dj_<-rep(NA,n)
  Dij_normIdx<- rbinom(n,1,Dj_rho_)
  
  xi_<-param_list_$l.xi_
  eta_<-param_list_$l.eta_
  alpha_ <- xi_ * eta_
  theta_<- u_ %*% t(alpha_)
  beta_ <-theta_ * matrix(as.numeric(abs(theta_)>param_list_$t_),nrow=n)
  for(i in 1:n){
    Dij_old<-Dj_old[i]
    if(Dij_normIdx[i]==1){
      Dij_prop <- 1
    }else{
      Dij_prop <- extraDistr::rinvchisq(1, prop_dgfj_,1)
    }
    # prior 
    prior_diff <-extraDistr::dinvchisq(Dij_prop, prior_dgfj_,1) - extraDistr::dinvchisq(Dij_old, prior_dgfj_,1)
    
    # propProb
    old2Prop <- ifelse(Dij_old==1, Dj_rho_, (1-Dj_rho_)*extraDistr::dinvchisq(Dij_old,prop_dgfj_,1))
    prop2Old <- ifelse(Dij_prop==1, Dj_rho_, (1-Dj_rho_)*extraDistr::dinvchisq(Dij_prop,prop_dgfj_,1))
    prop_diff <- old2Prop - prop2Old
    
    #llh
    llh_diff <- beta_[i]* xj_rescl[i] * (1/sqrt(Dij_prop)-1/sqrt(Dij_old))
    
    lr <- prior_diff + llh_diff + prop_diff
    
    MH.u<-runif(1)
    if(lr>MH.u){
      Dj_[i]<-Dij_prop
    }else{
      Dj_[i]<-Dij_old
    }
  }
  return(Dj_)
}

update_Dj_rho<-function(ad_rho, bd_rho, Dj_){
  res.Dj.rho <- rbeta(1,ad_rho + sum(Dj_==1),
                   bd_rho + sum(Dj_!=1))
  return(res.Dj.rho)
}

update_Gibbs_t<-function(y_, x_, u_, d0_, l_eta_, l_xi_, priorB_, sigProp_){
  n<-length(y_)
  p<-ncol(x_)
  l.alpha_ <- l_xi_ * l_eta_
  theta_<- u_ %*% t(l.alpha_)
  g_k<-matrix(0,nr=nrow(theta_), nc=ncol(theta_))
  for(i in 1:n){
    PiVec<- (theta_[i,]^2 * x_[i,]^2) - (2*theta_[i,]*x_[i,]*y_[i])
    theta_I<-theta_[i,]
    QijVec<-2*theta_[i,]*x_[i,]
    rmnOrdSet<-ordIdx<-order(abs(theta_I))
    for(jPrime in ordIdx){
      rmnOrdSet<-setdiff(rmnOrdSet,jPrime)
      # QijVec[jPrime]<- QijVec[jPrime] * prod(theta_[i,rmnOrdSet]) * prod(x_[i,rmnOrdSet])
      QijVec[jPrime]<- QijVec[jPrime] * sum(theta_[i,rmnOrdSet] * x_[i,rmnOrdSet])
    }
    g_k[i,]<-g_k[i,] -(PiVec + QijVec)/(2*d0_[i])
  }
  
  absTheta<-as.numeric(abs(theta_))
  ordAbsTheta<-order(absTheta)
  g_kVec <- as.numeric(g_k)
  g_kVec <- g_kVec[ordAbsTheta]
  cndtThre <- absTheta[ordAbsTheta]
  
  # sigLBCutOff<-quantile(absTheta, prob= 1 - sigProp_)
  # sigLBIdx <- cndtThre > sigLBCutOff
  # cndtThre<-cndtThre[sigLBIdx]
  # g_kVec <- g_kVec[sigLBIdx]
  
  sigLB_Prop<-1-sigProp_
  sigLB_cutOff <- quantile(cndtThre, sigLB_Prop)
  sigLB_cutOffIdx <- cndtThre > sigLB_cutOff
  cndtThre <- cndtThre[sigLB_cutOffIdx]
  g_kVec <- g_kVec[sigLB_cutOffIdx]
  
  if(priorB_ < max(absTheta)){
    priorB_cutOffIdx<- cndtThre<priorB_
    cndtThre<-cndtThre[priorB_cutOffIdx]
    g_kVec <- g_kVec[priorB_cutOffIdx]
  }
  cndtThre<-c(cndtThre,priorB_)
  g_kVec <- c(g_kVec,0)
  # g_kVec <- g_kVec - log(priorB_)
  
  F.vec<-rev(cumsum(rev(g_kVec)))
  logProb_ <- F.vec - corrcoverage::logsum(F.vec)
  # ggplot(data.frame(x=cndtThre,y=exp(logProb_)),aes(x=x,y=y))+geom_line()
  selIdx<- sample(1:length(cndtThre), 1, replace = T, prob = exp(logProb_))
  if(selIdx == 1){
    t.res<-runif(1,min=0,max=cndtThre[1])
  }else{
    t.res<-runif(1,min=cndtThre[selIdx-1],max=cndtThre[selIdx])
  }
  return(t.res)
}

rBGR_mcmc_Int<-function(y, x, u, N, burnin, seed_, dgf0, priorB, sigProp, prior_dgfj, prop_dgfj){
  n<-dim(x)[1]
  p.ord<-dim(x)[2]
  
  u.int<-cbind(1,u)
  q.int<-dim(u.int)[2]
  a_nu = 2; b_nu = 2;
  a_rho = 1;b_rho=1;
  ad_rho =1; bd_rho=1;
  a_s2=1; b_s2=1
  v0_=10^{-5}
  set.seed(seed_)
  
  param_list<-list()
  
  param_list$l.xi_<-matrix(rnorm(p.ord*q.int), nrow=p.ord, ncol=q.int)
  param_list$l.m_<-matrix(1, nrow=p.ord, ncol=q.int)
  param_list$l.eta_<-matrix(rnorm(p.ord*q.int),nrow=p.ord,ncol=q.int)
  param_list$l.nu_<-matrix(1,nrow=p.ord,ncol=q.int)
  param_list$l.Z_<-matrix(1,nrow=p.ord,ncol=q.int)
  
  # param_list$mu.xi_<-matrix(rnorm(p.ord), nrow=p.ord,ncol=1)
  # param_list$mu.m_<-matrix(1, nrow=p.ord,ncol=1)
  # param_list$mu.eta_<-matrix(rnorm(p.ord),nrow=p.ord,ncol=1)
  # param_list$mu.nu_<-matrix(1,nrow=p.ord,ncol=1)
  # param_list$mu.Z_<-matrix(1,nrow=p.ord,ncol=1)
  
  param_list$t_<-runif(1,0,1)
  param_list$rho_<-0.5
  param_list$sigma2_=1
  
  param_list$priorB_<-priorB
  param_list$d0_<-rep(1,n)
  param_list$Dj_<-matrix(1,nrow=n,ncol=p.ord)
  
  param_list$dgf0_<-dgf0
  param_list$prior_dgfj_<-rep(prior_dgfj,p.ord)
  param_list$prop_dgfj_<-rep(prop_dgfj,p.ord)
  param_list$Dj_rho_ <- rep(0.5,p)
  
  u.int<-cbind(1,u)
  
  l.xi.Arr<-l.eta.Arr<-array(NA,dim=c(N-burnin,p.ord,q.int))
  mu.Z.Mat<-mu.xi.Mat<-mu.eta.Mat<-matrix(NA, nrow=N-burnin, ncol=p.ord)
  l.Z.Arr<-array(NA,dim=c(N-burnin,p.ord,q.int))
  tVec<-rep(NA,N-burnin)
  Dj.Arr<-array(NA,dim=c(N-burnin,n,p.ord))
  llh<-rep(NA,N)
  t.all<-rep(NA,n)
  
  set.seed(seed_)
  for(ite in 1:N){
    if(ite%%100==0){
      print(ite)
    }
    l.eta.res<-Gibbs_l_etaC(y_ = y, x_ = x/sqrt(param_list$Dj_), u_ = u.int, 
                            d0_ = param_list$d0_, t_ = param_list$t_,
                            l_eta_ = param_list$l.eta_, 
                            l_xi_ = param_list$l.xi_, 
                            l_Z_ = param_list$l.Z_, 
                            l_nu_ = param_list$l.nu_)
    param_list$l.eta_<-l.eta.res
    
    l.xi.res<-Gibbs_l_xiC(y_ = y, x_ = x/sqrt(param_list$Dj_), u_ = u.int, 
                          d0_ = param_list$d0_, t_ = param_list$t_,
                          l_eta_ = param_list$l.eta_, 
                          l_xi_ = param_list$l.xi_, 
                          l_Z_ = param_list$l.Z_, 
                          l_nu_ = param_list$l.nu_, 
                          l_m_ = param_list$l.m_)
    param_list$l.xi_<-l.xi.res
    
    xi.eta.res<-reScl.xi.eta(param_list_ = param_list)
    param_list$l.xi_<-xi.eta.res$l.xi_
    param_list$l.eta_<-xi.eta.res$l.eta_
    
    t.res<-update_Gibbs_t(y_ = y, x_ = x/sqrt(param_list$Dj_), u_ = u.int, 
                          d0_ = param_list$d0_,
                          l_eta_ = param_list$l.eta_, 
                          l_xi_ = param_list$l.xi_, 
                          priorB_ = param_list$priorB, sigProp_=sigProp)
    param_list$t_<-t.res
      
    m.res<-update_m(l.xi_ = param_list$l.xi_)
    param_list$l.m_<-m.res$l.m_
    
    Z.res<-update_Z(v0_=v0_, param_list_ = param_list)
    param_list$l.Z_<-Z.res$l.Z_
    
    nu.res<-update_nu(a_nu=a_nu, b_nu=b_nu, param_list_ = param_list)
    param_list$l.nu_<-nu.res$l.nu_
    
    rho.res<-update_rho(a_rho=a_rho, b_rho=b_rho, l.Z_ = param_list$l.Z_)
    param_list$rho_<-rho.res
    
    # sigma2.res<-update_sigma2(y_ = y, x_ = x, u_ = u, 
    #                           param_list_ = param_list, 
    #                           a_s2 = a_s2, b_s2 = b_s2)
    # param_list$sigma2_<-sigma2.res
    
    d0.res<-update_d0(y_=y, x_=x/sqrt(param_list$Dj_), u_ = u.int, param_list_ = param_list)
    param_list$d0_<-d0.res
  
    sigma2.res<-update_sigma2(dgf0_=param_list$dgf0_, d0_ = param_list$d0_)
    param_list$sigma2_<-sigma2.res
    
    
    for(j in 1:p.ord){
      # dj.res<-update_Dj_nodeWise(xj_ = x[,j],dgfj_ = param_list$dgf0_)
      dj.res <- update_Dj_MH(xj_ = x[,j], u_ = u.int, j_ = j, param_list_ = param_list)
      param_list$Dj_[,j]<-dj.res
      
      Dj.rho.res <- update_Dj_rho(ad_rho, bd_rho, Dj_ = dj.res)
      param_list$Dj_rho_[j]<-Dj.rho.res
    }
    
    
    xi_<-param_list$l.xi_
    eta_<-param_list$l.eta_
    alpha_ <- xi_ * eta_
    theta_<- u.int %*% t(alpha_)
    # mu.alpha_ <- param_list$mu.xi_ * param_list$mu.eta_
    # theta_<- mu.alpha_ * rep(1,n) + u %*% t(l.alpha_)
    # theta_<- u %*% t(l.alpha_)
    beta_<- theta_ * matrix(as.numeric(abs(theta_)>param_list$t_),nrow=n)
    regressor_ <- apply(x/sqrt(param_list$Dj_) * beta_,1,sum)
    llh[ite]<- sum(dnorm(y, mean=regressor_, sd=sqrt(param_list$d0_),log=T))
    t.all[ite]<- param_list$t_
    if(ite >burnin){
      l.xi.Arr[(ite - burnin),,]=param_list$l.xi_
      l.eta.Arr[(ite - burnin),,]=param_list$l.eta_
      l.Z.Arr[(ite - burnin),,]=param_list$l.Z_
      Dj.Arr[(ite - burnin),,] = param_list$Dj_
      tVec[ite - burnin] = param_list$t_
    }
  }
  return(list(llh=llh, Dj = Dj.Arr, t=tVec, t.all=t.all,
              l.eta=l.eta.Arr, l.xi=l.xi.Arr, l.Z=l.Z.Arr))
}
