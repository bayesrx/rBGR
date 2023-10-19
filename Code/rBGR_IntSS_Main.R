library(here)
source(here("rBGR_IntSS_DijMH.R"))

cancer.type<-"OV"
args <- commandArgs(TRUE)
#args<-c("1")
print(args)
y.ord<-as.numeric(args[1])

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

n<-nrow(Y.df)
p<-ncol(Y.df)
q<-ncol(X.df)

nu.true<-3


iteNum=20000
# iteNum=200
system.time({
  print(iteNum)
  fit.mcmc<-rBGR_mcmc_Int(y=Y.df[,y.ord], x=Y.df[,-y.ord], u=X.df,
                          N=iteNum, burnin=(iteNum * 0.9), 
                          seed_=y.ord *5679 , priorB = 1, sigProp = 0.2,
                          dgf0 = nu.true, prior_dgfj = nu.true, prop_dgfj = nu.true)
})
saveRDS(fit.mcmc, here("..","MCMC_Res",paste0(cancer.type,"_",y.ord,".RDS")))





