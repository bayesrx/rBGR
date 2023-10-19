library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
library(igraph)
library(ggraph)
library(here)

fit.mcmc.resExtract<-function(mcmc.ls,rng){
  l.eta.post<-mcmc.ls$l.eta[rng,,]
  l.xi.post<-mcmc.ls$l.xi[rng,,]
  Z.post<-mcmc.ls$l.Z[rng,,]
  alpha_<-apply(l.eta.post*l.xi.post,c(2,3),mean)
  Z_<-apply(Z.post,c(2,3),mean)
  return(list(alpha=alpha_,Z=Z_))
}

cancer.type<-"OV"
# cancer.type<-"LUAD"
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
X.df<-cbind(1,X.df)

n<-nrow(Y.df)
p<-ncol(Y.df)
q<-ncol(X.df)

#### population level ####
f.path<-here("..","MCMC_Res",cancer.type)
rBGR.ls<-list.files(f.path)
res.ls<-list()

alpha.arr<-Z.arr<-array(0,dim=c(p,p,q))
for(j in 1:p){
  print(j)
  f.name<-paste0("IntSS_yOrder",j,".RDS")
  
  if(!(f.name %in% rBGR.ls)){
    print(paste(f.name))
  }else{
    fit.mcmc<-readRDS(here(f.path,f.name))
  }
  rng<-1:length(fit.mcmc$t)
  tmp.ls<-fit.mcmc.resExtract(mcmc.ls = fit.mcmc, rng = rng)
  Z.rBGR.j<- tmp.ls$Z
  alpha.rBGR.j<- tmp.ls$alpha
  for(h in 1:q){
    Z.arr[j,-j,h]<-Z.rBGR.j[,h]
    alpha.arr[j,-j,h]<-alpha.rBGR.j[,h]
  }
}
saveRDS(list(Z=Z.arr,alpha=alpha.arr),here(f.path,"asym_intSS_PopSum.RDS"))

alpha.min.df<-matrix(NA,ncol=5,nrow=0)
alpha.arr.sym<-Z.arr.sym<-array(0,dim=c(p,p,q))
selIdx.df<-matrix(NA,ncol=5,nrow=q*p*(p-1)/2)
rIdx<-0
for(j in 1:(p-1)){
  for(k in (j+1):p){
    for(h in 1:q){
      rIdx<-rIdx+1
      can.Z<-c(Z.arr[j,k,h],Z.arr[k,j,h])
      can.alpha<-c(alpha.arr[j,k,h],alpha.arr[k,j,h])
      sel.idx<-which.min(can.Z)
      
      alpha.arr.sym[j,k,h]<-can.alpha[sel.idx]
      Z.arr.sym[j,k,h]<-can.Z[sel.idx]
      Z.min<-can.Z[sel.idx]
      selIdx.df[rIdx,]<-c(j,k,h,sel.idx,as.numeric(Z.min>0.5))
      
      # plt.df.min<-rbind(plt.df.min,c(j,k,as.numeric(Z.min>0.5),h))
      if(Z.min>0.5){
        # if(sign(alpha.arr[j,k,h]*alpha.arr[k,j,h])!=1){
        #   print(paste0("Inconsistent sign for j:",j,"; k:",k,"; h:", h))
        # }
        tmp<-c(j,k,h,can.alpha[sel.idx],can.Z[sel.idx])
        alpha.min.df<-rbind(alpha.min.df,tmp)
      }
    }
  }
}
colnames(selIdx.df)<-c("j","k","h","selIdx","Z")
colnames(alpha.min.df)<-c("j","k","h","alpha","Z")
selIdx.df<-data.frame(selIdx.df)
saveRDS(list(alpha=alpha.arr.sym,Z=Z.arr.sym,
             sel=selIdx.df,coef=alpha.min.df),here(f.path,"sym_intSS_PopRes.RDS"))
