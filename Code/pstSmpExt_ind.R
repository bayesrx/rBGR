library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(igraph)
library(ggraph)
library(here)


# cancer.type<-"OV"
cancer.type<-"LUAD"
Y.df<-as.matrix(readRDS(here("..","TCGA_Data","Proteomics_Data_Y",paste0("Y_",cancer.type,".rds"))))
X.df<-as.matrix(readRDS(here("..","TCGA_Data","mRNA_Immune_X",paste0("X_",cancer.type,".rds"))))

Y.mean<-apply(Y.df,2,mean)
Y.mean.mat<-matrix(rep(Y.mean,nrow(Y.df)),byrow=T, nrow=nrow(Y.df))
Y.df.cen<-Y.df - Y.mean.mat
Y.df<-Y.df.cen
rm(Y.df.cen)

X.df<-log(X.df)
X.sd<-apply(X.df,2,sd)
X.sd.mat<-matrix(rep(X.sd,nrow(X.df)),nrow=nrow(X.df),byrow=T)
X.mean<-apply(X.df,2,mean)
X.mean.mat<-matrix(rep(X.mean,nrow(X.df)),nrow=nrow(X.df),byrow=T)
X.df.reScl<-(X.df-X.mean.mat)/X.sd.mat
X.df<-X.df.reScl[,c("T cells","Monocytes","Neutrophils")]
rm(X.df.reScl)
X.int<-cbind(1,X.df)
rm(X.df)
immune.vec<-c("Intercept","T cells","Monocytes","Neutrophils")


n<-nrow(Y.df)
p<-ncol(Y.df)
q<-ncol(X.int)

f.path<-here("..","MCMC_Res",cancer.type)
f.ls<-list.files(f.path)

#### Calculate patient-specific (asymmetric) network ####
f.name<-paste0("IntSS_yOrder",1,".RDS")
fit.mcmc<-readRDS(here(f.path,f.name))
iteNum<-length(fit.mcmc$t)

alpha.mcmc<-array(0,dim=c(iteNum,p,p,q))
t.mcmc<-array(0,dim=c(iteNum,p))
print(iteNum)
for(j in 1:p){
  print(j)
  f.name<-paste0("IntSS_yOrder",j,".RDS")
  if(!(f.name %in% f.ls)){
    print(paste("yOrder:", j))
  }else{
    fit.mcmc<-readRDS(here(f.path,f.name))
  }
  rng<-1:iteNum
  alpha.j<-fit.mcmc$l.eta[rng,,]*fit.mcmc$l.xi[rng,,]
  alpha.mcmc[,j,-j,]<-alpha.j
  t.mcmc[,j]<-fit.mcmc$t[rng]
}

print("Coefficient Symmetrization")
selIdx<-readRDS(here(f.path,"sym_intSS_PopRes.RDS"))$sel # selIdx=1: j->k; selIdx=2: k->j
for(l in 1:iteNum){
  if(l%%100==0)print(l)
  alpha_<-alpha.mcmc[l,,,]
  for(h in 1:q){
    alpha.h<-alpha_[,,h]
    selIdx.h<-selIdx[selIdx$h==h,]
    for(rIdx in 1:nrow(selIdx.h)){
      j<-selIdx.h[rIdx,1]
      k<-selIdx.h[rIdx,2]
      pos<-selIdx.h[rIdx,4]
      # Z_<-selIdx.h[rIdx,5]
      if(pos==1){
        # alpha.h[j,k]<-alpha.h[j,k]*Z_
        alpha.h[k,j]<-alpha.h[j,k]
      }else{
        # alpha.h[k,j]<-alpha.h[k,j]*Z_
        alpha.h[j,k]<-alpha.h[k,j]
      }
    }
    alpha_[,,h]<-alpha.h
  }
  alpha.mcmc[l,,,]<-alpha_
}
saveRDS(list(alpha=alpha.mcmc,t=t.mcmc),
        here(f.path,"IteWise_intSS_SymParam_noZ.RDS"))


print("Calculate nonZero Posterior Probability")
ePos<-eNeg<-eZero<-array(0, dim=c(p,n,p))
for(j in 1:p){
  print(j)
  negTmp<-zeroTmp<-posTmp<-array(0,dim=c(n,p-1))
  for(l in 1:iteNum){
    alpha.k<-alpha.mcmc[l,j,-j,]
    theta.j.X<-X.int %*% t(alpha.k)
    beta.j.X<-theta.j.X * as.numeric(abs(theta.j.X)>t.mcmc[l,j])
    posTmp<-posTmp+as.numeric(beta.j.X>0)
    negTmp<-negTmp+as.numeric(beta.j.X<0)
    zeroTmp<-zeroTmp+as.numeric(beta.j.X==0)
  }
  posTmp<-posTmp/iteNum
  negTmp<-negTmp/iteNum
  zeroTmp<-zeroTmp/iteNum
  ePos[j,,-j]<-posTmp
  eNeg[j,,-j]<-negTmp
  eZero[j,,-j]<-zeroTmp
}
saveRDS(list(Pos=ePos,Neg=eNeg,Zero=eZero),
        here(f.path,"indEdge_intSS_aSymArr_ePP.RDS"))

PrBeta<-matrix(NA,nrow=n*(p-1)*p/2,ncol=9)
colnames(PrBeta)<-c("j","k","i","PrNeg_jk","PrNeg_kj","PrZero_jk","PrZero_kj","PrPos_jk","PrPos_kj")
rIdx<-0
for(i in 1:n){
  net.i.pos<-ePos[,i,]
  net.i.neg<-eNeg[,i,]
  net.i.zero<-eZero[,i,]
  for(j in 1:(p-1)){
    for(k in (j+1):p){
      rIdx<-rIdx+1
      PrBeta[rIdx,]<-c(j,k,i, net.i.neg[j,k], net.i.neg[k,j],
                       net.i.zero[j,k], net.i.zero[k,j],
                       net.i.pos[j,k], net.i.pos[k,j])
    }
  }
}
saveRDS(PrBeta,here(f.path,"indEdge_intSS_aSym_ePP.RDS"))

PrBeta_nonZero<-1-PrBeta[,c("PrZero_jk","PrZero_kj")]
selIdx.edge<-unlist(apply(PrBeta_nonZero,1,which.max))
undir.PrBeta<-PrBeta_nonZero[cbind(1:nrow(PrBeta_nonZero),selIdx.edge)]
undir.edge<-data.frame(PrBeta[undir.PrBeta>0.5,])
saveRDS(undir.edge, here(f.path,"indEdge_intSS_sym_ePP.RDS"))

